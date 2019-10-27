#!/usr/local/bin/python2.7
# encoding: utf-8
"""
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe Profil

@author:     Michel Le Berre, Pierre Puiseux
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 31 Jan 2013
"""
from spleen.utilitaires.utilitairesdivers import (className, debug, rdebug, dist, rcercle,diff, debug, XY, rdebug, dist2)
import numpy as np
import cPickle, os
from parametresprofil import ProfsParamNew, ProfsParam, ProfsParam1
from naca import naca4,naca5
from numpy import (zeros,abs, arctan2,asarray, hstack, vstack, linspace,)
from scipy.optimize import newton
from matplotlib import pyplot as plt
from matplotlib.widgets import CheckButtons
plt.rcParams["figure.figsize"] = (20,10)
from spleen.preferences import ProfilPrefs
from spleen.splineabstraite import arrange
from spleen.splinecomposee import NSplineComposee
from spleen.utilitaires.utilitairesprofil import computeCordeAndNBA
from plotly.offline import plot
import plotly.graph_objs as go

class Profil(NSplineComposee):
    prefs = ProfilPrefs
    class Default(dict) :
        u"""Un dictionnaire avec les valeurs par défaut.
        1. Il initialise les clés intermédiaires nécessaires à la
            construction : (qui commencent par '_', dans le dict.__init__)
        2. la property Default.dump fournit les attributs par defaut exposés
            dans l'API, i.e. après construction Profil(), le Profil expose
            les attributs K de dump.
        >>> P = Profil()
        >>> print P.Default().dump
        ... {'name'   : 'Profil', 'ruptures': [0, 1, -1],
        ...  'cpoints': array([[1., 0.],[0., 0.],[1., 0.]]),
        ...  'methode': [('cubic', ((2, 0.0, 0.0), (1, 0.0, -1.0))),
        ...              ('cubic', ((1, 0.0, -1.0), (2, 0.0, 0.0)))],
        ...  'role'   : 'Profil', 'mode': ['courbure', 'courbure'],
        ...  'nbpd'   : [1000, 1000], 'nbpe': [50, 50]}
        """
        def __init__(self) :
#             prefs = ProfilPrefs
            self.cpoints = asarray([[1.0, 0.0],[0.0, 0.0],[1.0, 0.0],])
                        #les points de séparation des splines
            self.ruptures = [0,1,-1]
                        #les methodes pour chaque brin de spline
            self.methode = [('cubic', ((2, 0.0,  0.0), (1, 0.0, -1.0))),
                        ('cubic', ((1, 0.0, -1.0), (2, 0.0,  0.0)))]
            self.mode    = ['courbure','courbure']
            self.nbpe    = [50, 50]
            self.nbpd    = [1000, 1000]
            self.name      = 'Profil'
            self.role      = 'Profil'
            self.nbspline  = 0,# nb appels à computeSpline
            self.nbdpoint  = 0,# nb calculs dpoints
            self.nbech     = 0,# nb echantillonnages (epoints)
            dict.__init__(self,
                   _cpoints  = self.cpoints,
                   _ruptures = self.ruptures,
                   _methodes = self.methode,
                   _nbpds    = self.nbpd,
                   _modes    = self.mode,
                   _nbpes    = self.nbpe,
                   name      = self.name,
                   role      = self.role,
                   nbspline  = self.nbspline,
                   nbdpoint  = self.nbdpoint,
                   nbech     = self.nbech,
                   )
        @property
        def dump(self):
            d = {}#est comme self, sans les _ : '_xxx'=>'xxx'
            K = ('cpoints','ruptures','methode','nbpe','mode','nbpd',
                        'name','role')
            for key in K :
                d[key] = getattr(self,key)
            return d

#             defautdump = dict(cpoints=cpoints)
    u"""
    Un profil est une NSplineComposee comportant deux NSplineSimple, l'une pour
    l'extrados, l'autre pour l'intrados.
    Avec un point de jonction (rupture) au BA.
    Les deux splines se joignent au BA avec une tangente verticale.

    Comme pour toute spline, il y a trois sortes de points :
    - cpoints = points de contrôle avec
            nba = numéro du BA dans les points de contrôle

    - epoints = points échantillonnés. Ce sont ceux qui servent à la coupe.
        iba = numéro du BA dans les points échantillonnés
        iouvext, iouvint = num. des pts ouverture dans les points échantillonnés
        profparam = concerne les points échantillonnés
        nbpe = [nbpeext, nbpeint] nb de pts d'échantillonnage sur chaque spline

    - dpoints = points de discrétisation, uniquement pour visualisation.
    """

    def __init__(self, **dump):
        u"""
        """
        self.conforme = None
        self.nb_echantillonnages = 0
        self._profparam = None
        try :
            naca = dump.pop('naca')
            if len(naca) == 2 : #liste non vide
                serie = naca[0].lower().strip('naca')
                if len(serie) == 4 : points = naca4(serie, naca[1])
                if len(serie) == 5 : points = naca5(serie, naca[1])
                dump['cpoints'] = points
                dump['name'] = 'naca'+serie
#                 corde, nba = computeCordeAndNBA(points)
#                 dump['ruptures'] = [0, nba, -1]
        except KeyError as msg :
            pass
#         for key in ('cpoints','points') :#P = Profil(cpoints=[[x,y],...])
#             if dump.has_key(key) and not dump.has_key('ruptures'):
#                 corde, nba = computeCordeAndNBA(dump[key])
#                 dump['ruptures'] = [0, nba, -1]

        #Appelle setDefaultValues() puis load() puis _update()
        super(Profil, self).__init__(**dump)
#         il faut calculer nba (i.e. ruptures) AVANT le _update
#         self._update()
        if len(self.splines) == 1:
            '''c'est un profil brut, sans BA, une seule spline, il faut tout
               refaire c'est le cas de tous les projets antérieurs à 14 dec 2017
            '''
            corde, nba = computeCordeAndNBA(self.cpoints)
            self.precision = [1000]
            self.mode = ['courbure']
            Re, Ri = corde/2, corde/4#C'est ce qui marche le mieux
            self.split(nba,
                       ('cubic',((2, 0, 0), (1, 0, -abs(Re)))),#extrados
                       ('cubic',((1, 0, -abs(Ri)), (2, 0, 0)))#intrados
                       )
            '''L'ouverture n'est pas précisée dans le dump
            (qui est chargé par super(Profil, self).__init__(**dump)),
            elle est donc définie par profparam.???
            Il devrait y avoir les valeurs par defaut.... Je les met ici'''
            iouvext, iouvint = self.profparam.iouverture
            self.pouverture = self.prefs.pouverture
#             try :
#                 ouvint,  ouvext  = self.cpoints[iouvint], self.cpoints[iouvext]
#                 pouvint=100*ouvint[0]/self.corde
#                 pouvext=100*ouvext[0]/self.corde
#                 self.pouverture = (-pouvext, -pouvint)
#             except IndexError as msg :
#                 pass
# #                 debug(msg)

            '''Les noms'''
            self.splines[0].role = 'Extrados'
            self.splines[0].name = 'Extrados-'+self.name
            self.splines[1].role = 'Intrados'
            self.splines[1].name = 'Intrados-'+self.name

#     def setDefaultValues(self):
#         """
#         Valeurs par defaut:
#         Pour constructeur vide
#         """
#         for key, value in self.default.iteritems() :
#             setattr(self, key, value)
#         self._cpoints    = asarray([[1.0, 0.0],[0.0, 0.0],[1.0, 0.0],])
#         self._ruptures   = [0,1,2]#les points de séparation des splines
#         #les methodes pour chaque brin de spline
#         self._methodes   = [('cubic', ((2, 0.0,  0.0), (1, 0.0, -1.0))),
#                             ('cubic', ((1, 0.0, -1.0), (2, 0.0,  0.0)))]
#         self._nbpds      = [1000, 1000]
#         self._modes      = ['courbure','courbure']#polyligne
#         self._nbpes      = [50, 50]
#         self.name        = className(self)#+'(%d)'%(id(self))
#         self.role        = className(self)
#         self.nbspline    = 0# nb appels à computeSpline
#         self.nbdpoint    = 0# nb calculs dpoints
#         self.nbech       = 0# nb echantillonnages (epoints)


    def _getT(self, x, t0=None, nbit=False, dt=[-0.1, 1.1], tol=1.0e-10, maxiter=50):
        u"""
        :return t: float, la valeur du paramètre t correspondant à l'abscisse
            |x|/100. Si t est non trouvé, ou si x==0.0 retourne np.nan
            Plus précisement la recherche se fait dans

                - S = extrados si x>0
                - S = intrados si x<0

            le t retourné est l'UNIQUE t tel que |x|/100 == S.sx(t).
            On suppose que l'intrados et l'extrados sont des FONCTIONS de x,
            c'est à dire qu'à un x donné dans [0,100], ne correspond qu'un seul
            point de l'[in,ex]trados.
            Si (par exemple à l'intrados) à x correspondent DEUX points,
            alors _getT(x) retourne est le premier point trouvé.
            Il faut donc utiliser t0 pour trouver le second, mais ça n'est pas
            prévu pour...
            Non testé dans ce dernier cas.
        :param x: float, 0<=x<=100. la position du point (dont on cherche
            l'abcisse curviligne t), en % de corde.
            La recherche se fait sur l'intrados si x<0 et sur l'extrados si x>0.
            si x==0.0, retourne nan
        :param t0 : float, la recherche se fait par itérations de Newton, elle
            démarre à t=t0
        :param nbit : bool
            - si nbit est True, retourne le couple (t, nb iterations)
            - si nbit est False, retourne t.
        :param dt : [tmin, tmax], est l'intervalle de recherche.
        :param tol: float, paramètre passé à np.newton(),
            (=The allowable error of the zero value.)
        :param maxiter: int, paramètre passé à np.newton(),
            (=Maximum number of iterations.)

        """
        if x == 0.0 : return (np.nan, 0) if nbit else np.nan
        absx = abs(x)/100.0
        #abscisse Px du point de la corde situé à absx% de corde
        Ax, Fx = self[self.nba][0], self[0][0]
        Px = (Ax + absx*(Fx-Ax))
        if absx > 1.0 :
            raise ValueError("∣x∣=%g devrait etre dans [0,100]. C'est une abscisse en %% de corde"%abs(x))
        S = self.splines[0] if x>0 else self.splines[1]
        if t0 is None :
            if x>0 : t0 = 1.0 - absx
            elif x<0 : t0 = absx
        k = 0
        while 1 :
            t = newton(lambda t: S.sx(t)-Px, t0, lambda t:S.sx(t,1),
                       tol=tol, maxiter=50, fprime2=lambda t:S.sx(t,2))
            k += 1
            return (t,k) if nbit else t
            # if not dt[0] <= t <= dt[1] : #dépasse les bornes
            #     return (np.nan, k) if nbit else np.nan
            # t0 += 0.1
    @property
    def nbpe(self):
        u"""
        :return: le nb de points d'échantillonnage.
        Dès que self.profparam existe, c'est lui qui décide de nbpe."""
        try :
            self.splines[0].nbpe = 1+self.iba
            self.splines[1].nbpe = self.profparam.nptprof-self.iba
            self._nbpe = [1+self.iba, self.profparam.nptprof-self.iba]
            return self._nbpe
        except AttributeError:#profparam n'existe pas encore
            return [s.nbpe for s in self.splines]
        except IndexError:#la spline[1] n'existe pas encore
            return [s.nbpe for s in self.splines]

    @nbpe.setter
    def nbpe(self, nbpe):
        u"""Inutilisé, sert seulement à éviter un AttributeError"""
        self._nbpe = nbpe

    def arrange(self, dump):
        u"""
        dump est un dictionnaire
        - comme celui retourné par self.toDump()
        - ou bien un dump avec 'cpoints' (quasiment) seul
        - ou bien un dump de SplineSimple  => conversion en Profil
        - ou bien un dump de SplineComposee => non traité
        On arrange le dump pour qu'il puisse convertir en Profil à peu
        près n'importe quoi
        """
        arrange(dump)#des cles mal nommées
        key='cpoints' #P = Profil(cpoints=[[x,y],...])
        if dump.has_key(key) and not dump.has_key('ruptures'):
            _, nba = computeCordeAndNBA(dump[key])
            dump['ruptures'] = [0, nba, -1]

        keys = dump.keys()
        key = 'classname'
        if key in keys :
            clsname = dump[key]
            if 'Profil' in clsname :
                return
            elif clsname == 'NSplineSimple' :
                debug(u"Conversion %s => %s"%(clsname, className(self)))
                cpoints = dump.pop('cpoints')#a l'abri
                name = dump.pop('name')
                dump.update(self.Default().dump)
                _, nba = computeCordeAndNBA(cpoints)
                dump.update(ruptures=[0, nba, len(cpoints)-1], cpoints=cpoints,
                            name=name, classname=className(self))
                return
        elif dump.has_key('cpoints') :
            #pas de key 'classname', appel en provenance de open('xxx.gnu')
            #ou constructeur vide ou autre
            cpoints = dump.pop('cpoints')#a l'abri
            try : name = dump.pop('name')
            except KeyError : name=className(self)

            dump.update(self.Default().dump)
            _, nba = computeCordeAndNBA(cpoints)
            dump.update(ruptures=[0, nba, len(cpoints)-1], cpoints=cpoints,
                        name=name, classname=className(self))

    def load(self, dump):
        u"""
        dump est un dictionnaire
        - comme celui retourné par self.toDump()
        - ou bien un dump avec 'cpoints' (quasiment) seul
        - ou bien un dump de SplineSimple  => conversion en Profil
        - ou bien un dump de SplineComposee => non traité
        """
        self.arrange(dump)
#         rdebug(nbpe=self.nbpe)
        super(Profil, self).load(dump)
#         rdebug(self_nbpe=self.nbpe, list_nbpes=[s.nbpe for s in self.splines])

        self.splines[0].name = self.name+'#Extrados'
        self.splines[0].role = 'extrados'

        self.splines[1].name = self.name+'#Intrados'
        self.splines[1].role = 'intrados'
        # dump = nom de fichier
#         if isinstance(dump,(str, unicode)) and os.path.isfile(dump):
#             f = open(dump,'r')
#             d = cPickle.load(f)
#             self.load(d)
        try :
            #Constructeur recopie :si dump possede la cle 'profil'
            #On le dumpe => pdump
            #puis les autres cles de dump remplacent celles de pdump, le cas echeant
            profil = dump.pop('profil')
            if not isinstance(profil, Profil) :
                raise TypeError("parametre 'profil' : une instance de Profil est attendue au lieu de : %s"%profil.__class__.__name__)
            pdump = profil.toDump()
            for key, value in dump.iteritems():#les autres cles de dump remplacent celles de pdump
                pdump[key] = value
#             dump = pdump
            self.load(pdump)#ne pas oublier, pour appeler super(Profil, self).load(dump)
            #ATTENTION, ne pas oublier le return, sinon il repasse deux fois dans le try suivant,
            #d'abord avec la cle 'profparam', puis sans la cle 'profparam' => valeurs par defaut
            return
        except KeyError :
            pass

        try :
            dparams = dump.pop('profparam')
#             debug('profparam', dparams)
            #on ne sait pas si c'est
            #un dictionnaire (projets recents > 25.02.2017) ou
            #un objet de type ProfsParam (projets tres anciens < xx.01.2017) ou
            #un objet de type ProfsParamNew (xx.01.2017 <projet< 25.02.2017 ou
            #une chaîne de caractères '[nptprof=86, iouverture=(46, 47), iba=37]'
            #une liste (certains autres projets?) ou
            self.profparam = ProfsParam1(**dparams)
        except TypeError as msg:
#             debug('TypeError', msg)
            self.profparam = dparams.castToProfsParam1()
        except KeyError as msg :#pas de profparam dans cette liste d'arguments on prend ceux par defaut
#             debug('KeyError', msg)
            self.profparam = ProfsParam1(*self.prefs.profparam1)
#         debug(name=self.name, profparams=self.profparam)
        try :
            self.pouverture = dump.pop('pouverture')
        except KeyError :
            self.pouverture = ProfilPrefs.pouverture
#         rdebug(self.nbpe)

    def dump(self,filename,format_='new'):
        cPickle.dump(self.toDump(format_),open(filename,'wb'))

    def toDump(self, format_='new'):
        u'''
        :return: lors d'un dump d'une base de profils, le profil self est sauvée sous cette forme.
        :rtype: dict.
        :param format_: str, 'new' ou bien 'old'

            - si format=='old', ce sont les points échantillonés qui sont sauvegardés,
            en plus des profparam, le nom, la position de l'ouverture en % de corde
            et le rôle
            - si format=='new', c'est la spline complète qui est sauvegardée,
            voir NSplineComposee.toDump().

        '''
        if format_=='new' :
            dump = super(Profil, self).toDump()
        elif format_=='old':
            #On ne garde que les points échantillonnés
            dump = {
                'points' : self.epoints,
                'profparam' : self.profparam,#.toDump(),
                'name' : self.name,
                'role' : self.role,
                    }
        try : dump['pouverture'] = self.pouverture#ouverture en % de corde
        except AttributeError : pass
        dump['profparam'] = self.profparam.toDump()
        return dump

#     def __repr__(self, format_='spl'):
#         if format_== 'spl' :
#             S = repr(self.toDump())
#         elif format_ in ('gnu', 'txt') :
#             S='#\n### Profil %s\n#'%self.name
#             for x,y in self.epoints :
#                 S += '\n'+str(x)+', '+str(y)

#     def toDumpOld(self,complete=False):
#         '''Lors d'un dump oldschool d'une base de profils, le profil self est sauvée sous cette forme'''
#         return todump

    def hardScale(self, scale, centre=None):
        u'''
        modification in situ de self, par application d'une homothétie
        de centre 'centre' et de rapport 'scale'.
        Les points de contrôle eux même sont modifiés.
        :param scale: float ou (float, float)

            - si scale est un float, le rapport d'homothétie est (hx, hy)=(scale, scale)
            - si scale est un (float, float), le rapport d'homothétie est (hx, hy)=(scale[0], scale[1])

            chaque point de contrôle P=(x,y) est transformé en P'=(x',y') avec
                x' = centre[0] + hx*(x-centre[0]), et
                y' = centre[1] + hy*(y-centre[1])
        '''
#         rdebug(self.profparam)
#         exit()
        if isinstance(scale, (int,float)):
            scale = (scale,scale)
#         rdebug(scale)
        super(Profil, self).hardScale(scale, centre)
#         """Mise à jour ouverture"""
#         prof=Profil(profil=self)
#         prof.hardScale(scale)
#         prof.name='%.2gx'%scale
        return

    def scaled(self, scale, centre=None):
        '''
        Retourne une COPIE de self, mise à l'échelle scale.
        Pas de modif de self
        '''
#         rdebug(self.profparam)
#         exit()
        if isinstance(scale, (int,float)):
            scale = (scale,scale)
#         rdebug(scale)
        prof = Profil(profil=self)
        prof.hardScale(scale, centre)
        prof.name = self.name+'%.2gx,%.2gx'%scale
        return prof

    @property
    def pouverture(self):
        u"""Position des ouvertures en % de corde"""
        return self._pouverture

    @pouverture.setter
    def pouverture(self, pouverture):
        u"""Position des ouvertures en % de corde
        pouverture = (xe, xi)
        xe et xi sont les positions de l'ouverture en % de corde avec la convention :
        le point (x=xe ou x=xi) est
        - sur l'extrados si x>0
        - sur l'intrados si x<0
        """
        if pouverture[0] < pouverture[1] :
            msg = u"%s : attention, recouvrement extrados/intrados : %s"%(self.name, str(pouverture))
            debug(msg)
        if pouverture[0]*pouverture[1] < 0 :
            msg = u"%s : attention, les points ouverture ne sont pas sur le meme trados : %s"%(self.name, str(pouverture))
            debug(msg)
        self._pouverture = pouverture
#         stack()
        if hasattr(self, '_epoints') :
            #pour que self soit rééchantillonné au besoin
            del self._epoints
    @property
    def recouvrement(self):
        u"""L'intrados tissu recouvre l'extrados tissu"""
        pam, pav = self.pouverture
        return pam < pav
    @property
    def touverture(self):
        u"""
        Retourne (kam, tam), (kav,tav) où
        - kam est le numero de la spline contenant le point amont de l'ouverture (0=extrados, 1=intrados)
        - tam est le paramètre (abscisse curviligne) du point amont de l'ouverture sur sa spline.
        - kam et tav : idem pour le point aval.
        >>> (kam, tam), (kav,tav) = self.touverture
        >>> xam, yam = self[kam](tam)#les coordonnées du point amont
        >>> xav, yav = self[kav](tav)#les coordonnées du point aval
        """
        pam, pav = self.pouverture[0], self.pouverture[1]#ouverture amont/aval en % de corde

        if pam < 0.0 : #point amont à l'intrados
            kam, tam = 1, self._getT(pam)
            if pav < 0.0 : #Le cas usuel : tout à l'intrados
#                 debug(loc_ouv='intrados,intrados')
                kav, tav = 1, self._getT(pav)
            elif pav == 0.0 :#point aval au BA self._getT(0) retourn nan Recouvrement
                kav, tav = 1, 0.0
            else : # pav > 0.0 : extrados => recouvrement
                kav, tav = 0, self._getT(pav)
        elif pam == 0.0 : #point amont au BA, self._getT(0) retourn nan
            if pav < 0.0 :
                kav, tav = 1, self._getT(pav)
                kam, tam = 1, 0.0
            elif pav == 0.0 :
                kav, tav = 1, 0.0
                kam, tam = 1, 0.0
            else : # pav > 0.0 : tout extrados
                kav, tav = 0, self._getT(pav)
                kam, tam = 0, 1.0
        else : #pam > 0.0 : #point amont à l'extrados
            kam, tam = 0, self._getT(pam)
            if pav < 0.0 : #point aval intrados
                kav, tav = 1, self._getT(pav)
            elif pav == 0.0 :#self._getT(0) retourn nan
                kav, tav = 1, 0.0
            else : # pav > 0.0 : extrados => recouvrement
                kav, tav = 0, self._getT(pav)
        return (kam, tam), (kav, tav)
#         return (0 if pam>0 else 1, self._getT(pam)), (0 if pav>0 else 1,self._getT(pav))#l'abscisse curv de pam/pav dans la spline S
    @property
    def ouverture(self):
        u"""Retourne les coord des deux points l'ouverture"""
        (kam, tam), (kav,tav) = self.touverture
        return self.splines[kam](tam), self.splines[kav](tav)

    @property
    def douverture(self, nbp=100):
        u"""Retourne un tableau numpy de nbp points discrétisés de l'ouverture"""
        touv = (kam, tam), (kav,tav) = self.touverture
#         debug(pouv=self.pouverture, touv=touv)
#         if kam!=kav :
        if kam == kav : #deux points sur le meme trados
            return self.splines[kam].echantillonner(nbp, ta=tam, tb=tav, mode='linear', onlyT=False)
        else :
            rdebug('TODO. Pour le moment les deux points de l\'ouverture doivent etre sur le meme trados' )
            return zeros((nbp,2))

    def echantillonner(self):
        u"""
        Echantillonnage de tout le profil, compte tenu de l'ouverture,
        selon les prescriptions de self.profparam.
        ********************************************************************
        *******On suppose que l'ouverture est entièrement à l'intrados******
        ********************************************************************
        Le profil est constitué de 2 splines.
        - la spline extrados, du BF au BA (vrai BA, le point le plus éloigné du BF)
        - la spline intrados, du vrai BA au BF
        les deux splines sont echantillonnées en mode 'cos', i.e. beaucoup de
        points au début et à la fin (BA et BF) et points espacés au milieu.
        Les epoints sont ensuite modifiés pour satisfaire les prescriptions
        de profparam
        Les points d'échantillonnage sont constitués de 4 morceaux :
        - extrados : du BF au vrai BA => Text, eext
        - retour : du BA au 1er point de l'ouverture => Tret
        - ouverture : l'ouverture = Touv
        - intrados : du dernier point ouverture au BF => Tfin
        """
#         epoints = super(Profil, self).echantillonner(nbp, mode)
#         debug(u'%s : <%s> %d-eme echantillonnage'%(self.name, className(self), self.nb_echantillonnages))

        mode=['cos','courbure']
        self.nb_echantillonnages += 1
        _, ri = self.rba
        if abs(ri) < self.corde/100 :
            raise RuntimeWarning('Rayon de bord d\'attaque quasi nul (%.2g) à l\'intrados, => problemes d\'echantillonnage'%ri)

        """extrados : echantillonnage standard"""

        pp = self.profparam
        Se = self.splines[0]#extrados
        Text = Se.echantillonner(nbp=1+pp.iba, ta=0, tb=1, mode=mode[0])
        eext = Se(Text)#points echantillon extrados
        #pour debug
#         eint = np.zeros((pp.nptint-pp.nptret,2))
#         self._epoints = np.vstack((eext,eint))
#         return self._epoints
        #fin pour debug

        u"""Intrados : les ouvertures doivent coincider avec des points echantillon fixés par profparam"""
        pam, pav = self.pouverture[0], self.pouverture[1]#ouverture amont/aval en % de corde
        if pam>0 or pav>0 :
            raise NotImplementedError(u"Ouverture sur l'extrados : pouverture=%s"%str(self.pouverture))

        Si = self.splines[1]#pam<0, et Si = intrados
        (kam, tam), (kav, tav) = self.touverture#l'abscisse curv de pam/pav dans la spline Si
        # Il faut nbpret points dans [BA, am], y compris BA et ouvamont
        if kam==0 or kav==0 :#kam, kav=numéro de la spline
            raise NotImplementedError('Les deux points d\'ouverture doivent être sur l\'intrados')

        nbpret = pp.nptret
#         debug(nbpret=nbpret)
        nbpouv = pp.nptouv+2#les points ouverture ne sont pas compris dans pp.nptouv
        nbpfin = pp.nptprof - pp.iouvint
#         debug(touv=(tam,tav), nbpret=nbpret, nbpouv=nbpouv, nbpfin=nbpfin)
#         debug('appel echantillonner retour', nbpret=nbpret, ta=0, tb=tam, mode=mode[1])
        try :
            #echantillonnage retour BA=>am
            Tret = Si.echantillonner(nbp=nbpret, ta=0, tb=tam, mode='courbure')
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
        except ValueError :#as msg :
            rdebug('nbpret=%d<0 ?'%nbpret)
            raise
#         debug(eret=eret, len_eret=len(eret))
#         debug('appel echantillonner ouverture', nbpouv=nbpouv, ta=tam, tb=tav, mode=mode[1])
        try :
            Touv = Si.echantillonner(nbp=nbpouv, ta=tam, tb=tav, mode='courbure')#echantillonnage ouverture
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
#         rdebug(eouv=eouv, len_eouv=len(eouv))
        #Tfin = les t de l'intrados tissu
        Tfin = Si.echantillonner(nbp=nbpfin, ta=tav, tb=1, mode='cos')#echantillonnage ouverture=>BF
#         rdebug('===> Nb points', ret=len(eret), ouv=len(eouv),efin=len(efin))
#         debug(Tret=Tret, Touv=Touv,Tfin=Tfin)
        # les trois absc. curv. Tret[0], Touv[0] et Touv[-1] désignent des points doubles
        # - Tret[0]  est l'absc. curv. du BA ; S(Tret[0]) il est deja dans eext
        # - Touv[0]  est l'absc. curv. de l'extrémité amont de l'ouverture, elle coincide avec Tret[-1]
        # - Touv[-1] est l'absc. curv. de l'extrémité avale de l'ouverture, elle coincide avec Tfin[0]
        # Tfin contient les absc. curv. des points intrados tissu
        # Tint = les t de l'intrados théorique
        Tint = hstack([Tret[1:],Touv[1:-1],Tfin])
        self._techint = Tint
        self._techext = Text
        self._tech = [Text, Tint]
#         debug(T_intrados=Tint)
        eint = Si(Tint)#les points intrados théorique echantillonnés
        self._epoints = vstack((eext,eint))
#         debug(nb_points_echantillonnage=len(self._epoints))
        return self._epoints


    def echantillonnerOld(self):
        u"""
        Echantillonnage de tout le profil, par tranche, compte tenu de l'ouverture,
        selon les prescriptions de self.profparam.
        ********************************************************************
        *******On suppose que l'ouverture est entièrement à l'intrados******
        ********************************************************************
        Le profil est constitué de 2 splines.
        - la spline extrados, du BF au BA (vrai BA, le point le plus éloigné du BF)
        - la spline intrados, du vrai BA au BF
        Les points d'échantillonnage sont constitués de 4 morceaux :
        - extrados : du BF au vrai BA => Text, eext
        - retour : du BA au 1er point de l'ouverture => Tret
        - ouverture : l'ouverture = Touv
        - intrados : du dernier point ouverture au BF => Tfin
        """
#         epoints = super(Profil, self).echantillonner(nbp, mode)
#         debug(u'%s : <%s> %d-eme echantillonnage'%(self.name, className(self), self.nb_echantillonnages))

        mode=['courbure','courbure']
        self.nb_echantillonnages += 1
        _, ri = self.rba
        if abs(ri) < self.corde/100 :
            raise RuntimeWarning('Rayon de bord d\'attaque quasi nul (%.2g) à l\'intrados, => problemes d\'echantillonnage'%ri)

        """extrados : echantillonnage standard"""

        pp = self.profparam
        Se = self.splines[0]#extrados
        Text = Se.echantillonner(nbp=1+pp.iba, ta=0, tb=1, mode=mode[0])
        eext = Se(Text)#points echantillon extrados
        #pour debug
#         eint = np.zeros((pp.nptint-pp.nptret,2))
#         self._epoints = np.vstack((eext,eint))
#         return self._epoints
        #fin pour debug

        u"""Intrados : les ouvertures doivent coincider avec des points echantillon fixés par profparam"""
        pam, pav = self.pouverture[0], self.pouverture[1]#ouverture amont/aval en % de corde
        if pam>0 or pav>0 :
            raise NotImplementedError(u"Ouverture sur l'extrados : pouverture=%s"%str(self.pouverture))

        Si = self.splines[1]#pam<0, et Si = intrados
        (kam, tam), (kav, tav) = self.touverture#l'abscisse curv de pam/pav dans la spline Si
        # Il faut nbpret points dans [BA, am], y compris BA et ouvamont
        if kam==0 or kav==0 :#kam, kav=numéro de la spline
            raise NotImplementedError('Les deux points d\'ouverture doivent être sur l\'intrados')

        nbpret = pp.nptret
#         debug(nbpret=nbpret)
        nbpouv = pp.nptouv+2#les points ouverture ne sont pas compris dans pp.nptouv
        nbpfin = pp.nptprof - pp.iouvint
#         debug(touv=(tam,tav), nbpret=nbpret, nbpouv=nbpouv, nbpfin=nbpfin)
#         debug('appel echantillonner retour', nbpret=nbpret, ta=0, tb=tam, mode=mode[1])
        try :
            #echantillonnage retour BA=>am
            Tret = Si.echantillonner(nbp=nbpret, ta=0, tb=tam, mode=mode[1])
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
        except ValueError :#as msg :
            rdebug('nbpret=%d<0 ?'%nbpret)
            raise
#         debug(eret=eret, len_eret=len(eret))
#         debug('appel echantillonner ouverture', nbpouv=nbpouv, ta=tam, tb=tav, mode=mode[1])
        try :
            Touv = Si.echantillonner(nbp=nbpouv, ta=tam, tb=tav, mode=mode[1])#echantillonnage ouverture
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
#         rdebug(eouv=eouv, len_eouv=len(eouv))
        #Tfin = les t de l'intrados tissu
        Tfin = Si.echantillonner(nbp=nbpfin, ta=tav, tb=1, mode=mode[1])#echantillonnage ouverture=>BF
#         rdebug('===> Nb points', ret=len(eret), ouv=len(eouv),efin=len(efin))
#         debug(Tret=Tret, Touv=Touv,Tfin=Tfin)
        # les trois absc. curv. Tret[0], Touv[0] et Touv[-1] désignent des points doubles
        # - Tret[0]  est l'absc. curv. du BA ; S(Tret[0]) il est deja dans eext
        # - Touv[0]  est l'absc. curv. de l'extrémité amont de l'ouverture, elle coincide avec Tret[-1]
        # - Touv[-1] est l'absc. curv. de l'extrémité avale de l'ouverture, elle coincide avec Tfin[0]
        # Tfin contient les absc. curv. des points intrados tissu
        # Tint = les t de l'intrados théorique
        Tint = hstack([Tret[1:],Touv[1:-1],Tfin])
        self._techint = Tint
        self._techext = Text
        self._tech = [Text, Tint]
#         debug(T_intrados=Tint)
        eint = Si(Tint)#les points intrados théorique echantillonnés
        self._epoints = vstack((eext,eint))
#         debug(nb_points_echantillonnage=len(self._epoints))
        return self._epoints
    @property
    def techext(self):
        return self.tech[0]
    @property
    def techint(self):
        return self.tech[1]

    @property
    def info(self):
        infos = super(Profil,self).info
        infos1 = [infos[k] for k in (-3,-2,-1)]
        infos = [infos[k] for k in (0,1,2,3,4,6,7,8,9)]
        infos.extend(infos1)
        i = u'epaisseur'
        try :
            i1 = u"%g"%(self.height)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)
        i = u'longueur'
        try :
            i1 = u"%g"%(self.longueur('r'))
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)

#         ni = ()
#         infos = infos + [u"    *Extrados :\n     ----------"] + self.splines[0].info
#         infos = infos + [u""]
#         infos = infos + [u"    *Intrados :\n     ----------"] + self.splines[1].info
#         infos = infos + [u""]
#         infos = infos + [u"*divers (%s) :"%self.name+u"\n     -------"]
        infos = infos + [u'%25s = '%u"pouverture"+u"%s"%(str(self.pouverture))]
        infos = infos + [u'%25s = '%u"corde"+u"%s"%(str(self.corde))]
        infos = infos + [u'%25s = '%u"profparam"+u"%s"%(str(self.profparam))]
        return infos
#
#         infos = infos + [u"fin <SplineProfil>"]

    @property
    def profparam(self) :
        return self._profparam
#
    @profparam.setter
    def profparam(self,profparam):
        u'''
        Les paramètres du profil :
        - ProfsParam1 (nptprof, iouvext, iouvint, iba,nbpbf,pourcentbf)
        - ou ProfsParamNew(nptprof, iouvext, iouvint, iba) => old fashion.
        - ou ProfsParam(nptext, nptint, iba) => very old fashion.
        Utiles uniquement pour l'échantillonnage.
        Quand on affecte le profparam, on doit répercuter ces modifs sur les
        splines intrados et extrados'''
        if not isinstance(profparam, (ProfsParam1, ProfsParamNew, ProfsParam)) :
            raise(TypeError,u"Une instance de ProfsParam1 est attendue,impossible d'affecter un %s : %s."%(profparam.__class__.__name__,str(profparam)))
        if hasattr(self, '_profparam')\
                and self._profparam is not None\
                and not profparam == self._profparam :
            try : del self._epoints
            except : pass

        pp = profparam.castToProfsParam1()
        self._profparam = pp
        self.splines[0].nbpe = 1+pp.iba
        self.splines[1].nbpe = pp.nptprof-pp.iba
        self._nbpe = [1+pp.iba, pp.nptprof-pp.iba]

#         s = self.splines[0]
#         s.
    def appendPoint(self, p):
        raise RuntimeError("Impossible d'ajouter un point a la fin d'un profil.")

#     def insertPoint(self,pos,k=None):
#         i = super(Profil,self).insertPoint(pos,k)
# #         iouvext, iouvint = self.profparam.iouverture
# #         if 0 <=i<= iouvext:
# # #             self.profparam.iouvext += 1
# #             self.profparam.iouverture = iouvext+1, iouvint+1
# #         elif self.profparam.iouvint <= i :
# #             pass
# #         self._update()
#         return i

#     def removePoint(self, pnum):
#         u"""
#         - Suppression du point
#         - mise à jour des profparam => NOOOON
#         """
#         point = super(Profil,self).removePoint(pnum)
# #         iouvext, iouvint = self.profparam.iouverture
# #         if 0 <= pnum <= iouvext:
# #             self.profparam.iouverture = iouvext-1, iouvint-1
# # #             self.profparam.iouvext-=1
# #         self._update() #update est deja appelé par super(Profil,self).removePoint(pnum)
#         return point

    def extrados(self, genre='tissu'):
        if genre=='spline' :
            return self.splines[0]
        if genre.lower()=='tissu':#délimité par l'ouverture
#            print whoami(self), 'self.ouverture[0], self.nba', self.ouverture[0], self.nba
            nfin=self.iouverture[0]
        else :#délimité par le Bord d'attaque
            nfin=self.nba
        return self.cpoints[:1+nfin,:]

    def intrados(self,genre='tissu'):
        if genre=='spline' :
            return self.splines[1]
        if genre.lower()=='tissu':#délimité par l'ouverture
            ndeb=self.iouverture[1]
        else :#délimité par le Bord d'attaque
            ndeb=self.nba
#        points = self.points
        return self.cpoints[ndeb:,:]

#     def copy(self):
#         return Profil(**self.toDump())

    def profilDilate(self, coefext, coefint):
        '''
        Calcule et retourne un profil dilaté.
        '''
        self.dilatext = 1+coefext/100.0
        self.dilatint = 1+coefint/100.0
        dilate = self.copy()
        dilate.splines[0].hardScale((1.0, self.dilatext), centre=asarray([0,0]))
        dilate.splines[1].hardScale((1.0, self.dilatint), centre=asarray([0,0]))
        return dilate
        # XY = self.points.copy()
        # XY[0:self.iba,1]*=self.dilatext
        # XY[self.iba:,1]*=self.dilatint
        # return XY


    def profilDilateOld(self,coefext,coefint):
        '''
        Calcul d'un profil dilaté
        '''
        self.dilatext=1+coefext/100.0
        self.dilatint=1+coefint/100.0
        XY = self.points.copy()
        XY[0:self.iba,1]*=self.dilatext
        XY[self.iba:,1]*=self.dilatint
        return XY

#     def lisserNew(self, profparam, tension=None,degre=(3,3)):
#         """faire quelque chose comme
#         spline = SplineProfil(self.points)
#         spline.echantillonner(profparam)
#         self.profilLisse =spline.epoints"""

    @property
    def profilLisse(self) :
        return self.dpoints

#     def validate(self):
#         u'''inutile        '''
#         return self.update()

#     def update(self):
#         return super(Profil,self)._update()

    def normalise(self, scale=1.0):
#        alert(self)
        facteur = scale/self.corde
        cba, cbf = asarray(self[self.nba]), asarray(self[0])
        u = cbf-cba
        alfa = arctan2(u[1], u[0])
        if abs(alfa) > 1.0e-5 :
#             debug( '##### u, alfa', u, 180*alfa/np.pi)
#             debug(alfa=alfa, cba=cba)
            self.hardRotate(-alfa, cba, 'radians')
        if facteur != 1.0 :#and facteur != (1.0, 1.0):
            self.hardScale((facteur, facteur),cba)
        self.translate(-cba)
        self.cpoints[0] = (1,0)#Ca fait un update...
        self.cpoints[self.nba] = (0,0)#Ca fait un update...
#         self._update()
        return cba, cbf

    def normalised(self, scale=1.0):
        u"""Calcule et retourne le QPolygonF normalisé"""
        prof = Profil(profil=self)
        prof.normalise(scale)
        return prof

    def ligneMoyenne(self):
        u"""Calcule et retourne la ligne moyenne"""
        T = linspace(0,1,100)
        return 0.5*(self.splines[0](T) + self.splines[1](T)[::-1])

    def verification(self):
        u'''
        Vérifie que la nervure est bien dans l'ordre
            BF>Extrados>BA>intrados,
            avec les yextrados >= 0  et
            yintrados <= 0
        '''
        def fmt(question, reponse, negatif=False):
            if negatif :
                alert = u'  ' if not reponse else u'=>'
            else :
                alert = u'  ' if reponse else u'=>'
            return u'%s %30s %s'%(alert, question, reponse)
        points=self.cpoints
        p = Profil(cpoints=points)
#         debug(p)
        p.normalise()
#         debug( self.name, p.points)
#         debug( p.points)
        extr=p.extrados(genre='theorique')
        intr=p.intrados(genre='theorique')
        reponse=[u"\n****\nVerification de %s (apres normalisation): "%self.name,(39+len(self.name))*'=']
        #BA-BF
        oks=[]
        ok=dist(p[0],p[-1])<=1.0e-8
        oks.append(ok)
        reponse.append(fmt(u"Profil ferm'e ?", ok))
#         reponse.append(s+u"Profil ferm'e : = %s"%ok)
        ok=dist(p[0],(1.0,0.0))<=1.0e-8
        oks.append(ok)
        reponse.append(fmt(u"BF =", str(p[0])))
        ok=dist(p[p.nba],(0.0,0.0))<=1.0e-8
        oks.append(ok)
        reponse.append(fmt(u"BA =", str(p[p.nba])))

        # Position relative intrados, extrados
        minex=np.argmin(extr[:,1])
        maxex=np.argmax(extr[:,1])
        minin=np.argmin(intr[:,1])
        maxin=np.argmax(intr[:,1])
        yminex=np.amin(extr[:,1])
        ymaxex=np.amax(extr[:,1])
        yminin=np.amin(intr[:,1])
        ymaxin=np.amax(intr[:,1])
        rep = fmt(u"indices extrados (min, max) =", u"(%-3d,%-3d)"%(minex,maxex))
        reponse.append(rep)
        rep = fmt(u"y extrados : (ymin, ymax) =", u"(%-5.3g,%-5.3g)"%(yminex,ymaxex))
        reponse.append(rep)
        rep = fmt(u"indices intrados (min, max) =", u"(%-3d,%-3d)"%(minin,maxin))
        reponse.append(rep)
        rep = fmt(u"y intrados : (ymin, ymax) =", u"(%-5.3g,%-5.3g)"%(yminin,ymaxin))
        reponse.append(rep)
#         reponse.append(u'    Extrados :  (min, max)=(%d, %d) ; (ymin, ymax)=(%.2g, %.2g)'%(minex,maxex,yminex,ymaxex))
#         reponse.append(u'    Intrados :  (min, max)=(%d, %d) ; (ymin, ymax)=(%.2g, %.2g)'%(minin,maxin,yminin,ymaxin))

        # Yi intrados negatif,
        ok=(ymaxin<=1.0e-8)
        oks.append(ok)
        reponse.append(fmt(u'Y intrados <=0 ?', ok))

        # Yi extrados positif
        ok=(yminex>=-1.0e-8)
        oks.append(ok)
        reponse.append(fmt(u'Y extrados >=0 ?',ok))

        # Extrados plus épais que intrados
        ok=(abs(ymaxex)>=abs(yminin))
        oks.append(ok)
        reponse.append(fmt(u'Extr. plus epais intr. ?' ,ok))
#         reponse.append(s+u'Extrados plus epais que intrados ? %s'%ok)

        # Les xi décroissent pour l'extrados
        dEx=diff(extr)
        dxex=dEx[:,0]#extr[1:,0] - extr[:-1,0]
        # dyex=dEx[:,1]#extr[1:,0] - extr[:-1,0]
        ok=np.all(dxex<0)
        oks.append(ok)
        reponse.append(fmt(u"Pli a l'extrados ?", not ok, negatif=True))
        # Les xi croissent pour l'intrados
        dIn=diff(intr)
        dxin=dIn[:,0]#intr[1:,0] - intr[:-1,0]
        # dyin=dIn[:,1]#intr[1:,0] - intr[:-1,0]
        ok=np.all(dxin>0)
        oks.append(ok)
        reponse.append(fmt(u"Pli a l'intrados ?", not ok, negatif=True))
        # Convexité de l'extrados => pas significatif.
        if 0 :
            d2Ex=dEx[1:]-dEx[:-1]
            ok=np.all(d2Ex[:,1]>=0)
            if not ok :
                debug('dY/dX   extrados\n',-dEx[:,1]/dEx[:,0])
                debug('d2Y/dX2 extrados\n',d2Ex[:,1]/(dEx[1:,0])*dEx[:-1,0])
            s='' if ok else u'*** '
            reponse.append(s+u"Extrados concave : %s"%(ok))
            oks.append(ok)
            # Petit creux en BF à l'intrados ?
            d2In=dIn[1:]-dIn[:-1]
            ok=np.all(d2In[:,1]>=0)
            s='' if ok else u'*** '
            if not ok :
                debug('dY/dX   Intrados\n',dIn[:,1]/dIn[:,0])
                debug('d2Y/dX2 Intrados\n',d2In[:,1]/(dIn[1:,0])*dIn[:-1,0])
            reponse.append(s+u"Intrados convexe : %s"%(ok))
            oks.append(ok)

        oks=asarray(oks)
        if np.all(oks) :
            self.conforme=True
#             reponse.append(str(oks))
            reponse.append(u'Tout semble OK\n****')
        else :
            self.conforme=False
#             reponse.append(str(oks))
            reponse.append(
                u"""
        Profil douteux.
        - Verifiez qu'il est dans l'ordre : BF>extrados>BA>intrados.
        - Verifiez qu'il n'a pas de pli.
        - Verifiez que l'extrados est plus epais de l'intrados.
        - Verifiez qu'il est ferme (dernier point=premier point)
                """)
        return reponse

    @property
    def iouverture(self):
        u"""Les deux numéros de points qui constituent l'ouverture"""
        return self.profparam.iouverture
    @iouverture.setter
    def iouverture(self,tranche):
        self.profparam.iouverture = (tranche[0], tranche[1])

    @property
    def rba(self):
        u"""
        Retourne les rayons de bord d'attaque extrados et intrados,
        i.e. les dérivées de y(t) au BA.
        Dans ce qui suit, l'attribut methode est la méthode de la
            spline 0 (extrados) ou 1 (intrados) avec :
        - pour l'extrados la methode est self.splines[0].methode =
            ('cubic',((k=2, derivéee k-ieme de xBF,yBF), (k=1, dérivée de xBA,yBA))
        - pour l'intrados la methode est self.splines[1].methode =
            ('cubic',((k=1, xBA',yBA'), (k=2, xBF",yBF"))
        """
        re, ri = self.splines[0].methode[1][1][2], self.splines[1].methode[1][0][2]
#         rdebug('Rayon bord d\'attaque : rext, rint = %.2f, %.2f'%(re,ri))
#         if re != ri :pass
        return re, ri
    @rba.setter
    def rba(self, rba):
        re, ri = rba
        sext, (bf, ba) = self.splines[0].methode
        ba = (1,0,re)
        self.splines[0].methode = (sext,(bf,ba))#appelle _update
        sint, (ba, bf) = self.splines[1].methode
        ba = (1,0,ri)
        self.splines[1].methode = (sint,(ba,bf))#appelle _update
        self._update()

    @property
    def erba(self):
        u"""Retourne le rayon de BA calculé avec les points echantillonnés"""
        P = self.epoints
        corde = self.corde
        nba = self.profparam.iba
        Re = rcercle(P[nba-2],P[nba-1],P[nba]) if nba>=2 else corde/2
        R0 = rcercle(P[nba-1],P[nba],  P[nba+1]) if nba>=1 and len(P)>nba+1 else corde/3
        Ri = rcercle(P[nba],  P[nba+1],P[nba+2]) if len(P)>nba+2 else corde/4
#         debug(Re=Re, Ri=Ri, R0=R0, corde=corde)
        return (Re,R0,Ri)

    @property
    def drba(self):
        u"""Retourne le rayon de BA calculé avec les points discrétisés"""
        P = self.dpoints
        corde = self.corde
        nba = len(self.splines[0].dpoints)-1
        Re = rcercle(P[nba-2],P[nba-1],P[nba]) if nba>=2 else corde/2
        R0 = rcercle(P[nba-1],P[nba],  P[nba+1]) if nba>=1 and len(P)>nba+1 else corde/3
        Ri = rcercle(P[nba],  P[nba+1],P[nba+2]) if len(P)>nba+2 else corde/4
#         debug(Re=Re, Ri=Ri, R0=R0, corde=corde)
        return (Re,R0,Ri)

    @property
    def nba(self):
        u"""
        Retourne le numéro du dernier point de contrôle de l'extrados.
        C'est le BA. Il a pour coordonnées (0,0) pour un ProfilNormalise.
        Par construction, c'est le point de controle le plus éloigné de BF=self.cpoints[0]
        Il ne peut changer que s'il y a insertion ou suppression de point à l'extrados
        """
        if len(self)<=1 : return 0
        return len(self.splines[0])-1

    @property
    def corde(self):
        if len(self)>0 :
            return dist(self[0],self[self.nba])
        else :
            return None

    @property
    def points(self):
        """Retourne les points echantillinnés.
        Pour les autres parties de Axile, compatibilité ascendante"""
        try :
            return self.epoints
        except TypeError :#on peut pas echantillonner
#             rdebug('impossible d\'echantillonner : nb points de controle=%d'%len(self))
            return np.zeros((0,2))
#     @property
#     def qcpolygon(self):
#         return qpolygonFrom(self.cpoints)
# #
    # @property
    # def qepolygon(self):
    #     return qpolygonFrom(self.epoints)
#
    @property
    def iba(self):
        return self.profparam.iba
    @iba.setter
    def iba(self,index):
        u"Bord d'attaque imposé"
        self.profparam.iba=index

    def plot0(self, figure=None, aspect={}, titre=None, more=[], texts=[],
             show=True, numbers=['3p'], filename='temp-plot.html'):
        """
        :param figure: une instance de matplotlib.figure.Figure
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
        :param texts : list, [(x,y, 'text', style=, fontsize=...), ...] texts en x,y
        :param numbers: list ou set ou tuple ['3c','12e','100d'] numéroter les
            points controle, echantillon, discretisation
            - '3c' pour numéroter les points de Contrôle par pas de 3,
            - '100d' les points de Discrétisation par pas de 100,
            - '24e' les points d'Echantillonnage par pas de 24.
        :param show: bool, True pour affichage de TOUS les plots en attente, False sinon
        :param aspect: dict, cf defaultaspect ci-dessous.
        :return figure: la figure passée en argument ou nouvelle figure si None

        """
    #     enable_mpl_offline()
        defaultaspect = {
                    'ce':'ro',#control extr:        red
                    'de':'r-',#discretisation extr : blue
                    'ee':'k.',#echantillon extr :    green
                    'ci':'bo',#control intr :        red
                    'di':'b-',#discretisation intr : blue
                    'ei':'k.',#echantillon intr :    green
                    'r' :'k*',#rupture:         black
                    'o' :'r*',#ouverture:         bleu
                    }
        defaultaspect.update(aspect)
        C = self.cpoints
        D = self.dpoints
        E = self.epoints
        CX, CY = XY(C)
        DX, DY = XY(D)
        EX, EY = XY(E)
        xouvext, youvext = self.epoints[self.profparam.iouvext]
        xouvint, youvint = self.epoints[self.profparam.iouvint]
        goOuv = go.Scatter({'x':(xouvext,xouvint), 'y':(youvext,youvint)},
                           mode='markers',
                           marker={'color':'green','size':12,'symbol':"diamond"},
                           name=u'Ouverture',
                           hoverinfo='x+y')
        goC = go.Scatter({'x':CX, 'y':CY, 'text':range(len(C))},
                         mode='markers',
                         marker={'color':'red','size':12},
                         name=u'Points de contrôle',
                         hoverinfo='text')
        goE = go.Scatter({'x':EX, 'y':EY, 'text':range(len(E))},
                         mode='markers',
                         marker={'color':'blue','size':8,'symbol':'star'},
                         name=u'Points échantillonnés',
                         hoverinfo='text')
        goD = go.Line({'x':DX, 'y':DY},
                      mode='lines',
                      line={'color':'green','width':1},
                      name=u'Spline',
                      hoverinfo='none')
        data = [goC,goD,goE,goOuv]
        #Pour aspect ratio
        layout = go.Layout(width = 2200, height = 1200,
                           title = titre,
                           xaxis = dict(nticks = 20,domain = [0, 1],
                                        title = 'x'),
                           yaxis = dict(scaleanchor = "x", nticks = 10,
                                        domain = [0, 1],),
                           legend = dict(x=0.9,y=0.9),
                           hovermode = 'closest')

        figure = go.Figure(data=data, layout=layout)
        plot(figure, show_link=True,
             link_text=u'Exporter vers plot.ly (pour édition)', filename=filename)
    #     fig['layout'].update(scene=dict(aspectmode="data"))
        return

    def plot(self, figure=None, aspect={}, titre=None, more=[], texts=[],
             show=True, buttons=True, numbers=['3p']):
        """
        :param figure: une instance de matplotlib.figure.Figure
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
        :param texts : list, [(x,y, 'text', style=, fontsize=...), ...] texts en x,y
        :param numbers: list ou set ou tuple ['3c','12e','100d'] numéroter les
            points controle, echantillon, discretisation
            - '3c' pour numéroter les points de Contrôle par pas de 3,
            - '100d' les points de Discrétisation par pas de 100,
            - '24e' les points d'Echantillonnage par pas de 24.
        :param show: bool, True pour affichage de TOUS les plots en attente, False sinon
        :param aspect: dict, cf defaultaspect ci-dessous.
        :return figure: la figure passée en argument ou nouvelle figure si None

        """
        defaultaspect = {
                    'ce':'ro',#control extr:        red
                    'de':'r-',#discretisation extr : blue
                    'ee':'k.',#echantillon extr :    green
                    'ci':'bo',#control intr :        red
                    'di':'b-',#discretisation intr : blue
                    'ei':'k.',#echantillon intr :    green
                    'r' :'k*',#rupture:         black
                    'o' :'r*',#ouverture:         bleu
                    }
        defaultaspect.update(aspect)
        if figure is None :#pas de figure => on la crée
            figure = plt.figure('plot(self)')

        axes = figure.get_axes()
        if axes : ax = axes[0]
        else : ax = figure.subplots()

        if titre : ax.set_title(titre)
        ax.axis('equal')
        C = self.cpoints
        D = self.dpoints
        E = self.epoints

        for chars in numbers :
            if len(chars)>1 :
                step = eval(chars[:-1])
            else :
                step = 1#numerotation de un point sur step
            char = chars[-1]
            if char=='c' :
                for k, p in enumerate(C[::step]) :
                    ax.text(p[0], p[1], '%d'%(step*k))
            elif char=='e' :
                for k, p in enumerate(E[::step]) :
                    ax.text(p[0], p[1], '%d'%(step*k))

        for txt in texts :
            ax.text(*txt)

        Cex = self.splines[0].cpoints
        Dex = self.splines[0].dpoints
        Eex = self.splines[0](self.techext)
        Cin = self.splines[1].cpoints
        Din = self.splines[1].dpoints
        Ein = self.splines[1](self.techint)
        R = asarray(self[self.ruptures])#des points
        xouvext, youvext = self.epoints[self.profparam.iouvext]
        xouvint, youvint = self.epoints[self.profparam.iouvint]
        fmtr  = defaultaspect['r']
        fmto  = defaultaspect['o']
        fmtce = defaultaspect['ce']
        fmtde = defaultaspect['de']
        fmtee = defaultaspect['ee']
        fmtci = defaultaspect['ci']
        fmtdi = defaultaspect['di']
        fmtei = defaultaspect['ei']
        if fmtce : econtrol,     = ax.plot(Cex[:,0], Cex[:,1], fmtce, lw=1,
                                            label=u'Contrôle extrados')
        if fmtde : espline,      = ax.plot(Dex[:,0], Dex[:,1], fmtde, lw=1,
                                            label=u'Spline extrados')
        if fmtee : eechantillon, = ax.plot(Eex[:,0], Eex[:,1], fmtee, lw=1,
                                            label=u'Échantillon extrados')
        if fmtci : icontrol,     = ax.plot(Cin[:,0], Cin[:,1], fmtci, lw=1,
                                            label=u'Contrôle intrados')
        if fmtdi : ispline,      = ax.plot(Din[:,0], Din[:,1], fmtdi, lw=1,
                                            label=u'Spline intrados')
        if fmtei : iechantillon, = ax.plot(Ein[:,0], Ein[:,1], fmtei, lw=1,
                                            label=u'Échantillon intrados')
        if fmtr  : BABF,         = ax.plot(R[:,0], R[:,1], fmtr, lw=1,
                                            markersize=12, label=u'BA-BF')
        if fmto  : ouverture,    = ax.plot([xouvext, xouvint],
                                           [youvext, youvint],
                                           fmto, markersize=12,
                                           label=u'Ouverture')
        for x, y, color, name in more:
            temp, = ax.plot(x, y, color,label=name)
#             if not name : continue
#             if buttons :
#                 draws.append(temp)
#                 labels.append(name)
#                 values.append(True)


        figure.legend()
        figure.subplots_adjust(left=0.2)
        if buttons :
            rax = figure.add_axes([0.05, 0.4, 0.1, 0.15])#
            labels = (u'Extrados spline', u'Extrados control', u'Extrados échantillon',
                      u'Intrados spline', u'Intrados control', u'Intrados échantillon',
                      u'BA-BF', u'Ouverture')
            draws  = (espline, econtrol, eechantillon,
                      ispline, icontrol, iechantillon,
                      BABF,   ouverture)
            values = (True, False, False,
                      True, False, False,
                      False, False)

            for draw, value in zip(draws,values):
                try : draw.set_visible(value)
                except AttributeError : pass
            check = CheckButtons(rax, labels, values)
            def func(label):
                try :
                    draw = draws[labels.index(label)]
                    draw.set_visible(not draw.get_visible())
                except Exception as msg:
                    rdebug(u"label=%s marche po : %s"%(label,str(msg)))
                    pass
                plt.draw()
            check.on_clicked(func)
        if show : plt.show()
        return figure

    def update(self):
        u"""
        Mise à jour de nba et profparam.
        Doit être appelé à chaque modification (suppression, insertion, deplacement) d'un point du profil
        - suppression, insertion de point : on reconstruit profparam entier=> ?????.
        - déplacement de point : seul nba peut changer.
        On appelle le _update de NSplineComposee (i.e. celui de NSplineAbstraite),
        """
        try : del self._techint
        except AttributeError : pass
        try : del self._techext
        except AttributeError : pass

        return super(Profil,self)._update()
    _update = update

if __name__=="__main__":
    from spleen.tests.testsprofil import testMain
    from spleen import spleenconfig
    spleenconfig.TEST_MODE = False
    testMain()
