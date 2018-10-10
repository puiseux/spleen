#!/usr/local/bin/python2.7
# encoding: utf-8
'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe Profil

@author:     Michel Le Berre, Pierre Puiseux
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 31 Jan 2013
'''
import sys,os,math,cPickle
from path import Path
from pprint import pprint
from polyline import NPolyLine
from config import DATA_DIR,VALIDATION_DIR,RUNS_DIR
# from inout.lecteurdata import LecteurData
# from inout.writerpts import writeProfil
from utilitaires import (stack, debug, rdebug, dist2, dist, rcercle)
import numpy as np
# from gui.graphicsbase.graphicscommon import p2s, p2t
from paramgeneraux import ProfsParamNew, ProfsParam
from naca import naca4,naca5
# from model.basicobjects.splinesimple import NSplineSimple
from splinecomposee import NSplineComposee
from lecteurs import pointsFrom
from utilitairesprofil import computeCordeAndNBA
from pprint import pprint
from numpy import (zeros,)
from numpy import asarray as array
from scipy.optimize import newton
from preferences import ProfilPrefs
class Profil(NSplineComposee):
    prefs = ProfilPrefs

    u"""
    Un profil est une NSplineComposee comportant deux NSplineSimple, l'une pour l'extrados,
    l'autre pour l'intrados. Avec un point de jonction (rupture) au BA.
    Les deux splines se joignent au BA avec une tangente verticale.

    Il y a trois sortes de points :
    - cpoints = points de contrôle avec
            nba = ciba = numéro du BA dans les points de contrôle

    - epoints = points échantillonnés. Ce sont ceux qui servent à la coupe.
            iba = eiba = numéro du BA dans les points échantillonnés
            iouvext, iouvint = numéros des points ouverture dans les points échantillonnés
            profparam = concerne les points échantillonnés
            nbpe = [nbpeext, nbpeint] nb de points d'échantillonnage sur chaque spline

    - dpoints = points de discrétisation, uniquement pour visualisation
    """

    def __init__(self, **dump):
        u"""
        """
#         debug(dump_keys=dump.keys())
#         pprint(dump)
        self.conforme = None
        self.nb_echantillonnages = 0
#         self.nb_normalisations = 0
        try :
            naca = dump.pop('naca')
            if len(naca) == 2 : #liste non vide
                serie = naca[0].lower().strip('naca')
                if len(serie) == 4 : points = naca4(serie, naca[1])
                if len(serie) == 5 : points = naca5(serie, naca[1])
                dump['points'] = points
                dump['name'] = 'naca'+serie
        except KeyError :
    #             debug( msg)
            pass
        #Appelle setDefaultValues() puis load() puis _update()
        super(Profil, self).__init__(**dump)
#         il faut calculer nba (i.e. ruptures) AVANT le _update
#         self._update()
        '''c'est un profil brut, sans BA, une seule spline, il faut tout refaire
           c'est le cas de tous les projets antérieurs à 14 dec 2017'''
        if len(self.splines) == 1:
            corde, nba = computeCordeAndNBA(self.cpoints)
            self.precision = [1000]
            self.mode = ['courbure']
#             P = self.cpoints
#             #Rayons de BA [ex,in]trados marche pas tres bien
#             Re = rcercle(P[nba-2],P[nba-1],P[nba]) if nba>=2 else corde/2
#             R0 = rcercle(P[nba-1],P[nba],  P[nba+1]) if nba>=1 and len(self)>nba+1 else corde/3
#             Ri = rcercle(P[nba],  P[nba+1],P[nba+2]) if len(self)>nba+2 else corde/4
# #             debug(Re=Re, Ri=Ri, R0=R0, corde=corde)
            Re, Ri = corde/2, corde/4#C'est ce qui marche le mieux
            self.split(nba,
                       ('cubic',((2, 0, 0), (1, 0, -abs(Re)))),#extrados
                       ('cubic',((1, 0, -abs(Ri)), (2, 0, 0)))#intrados
                       )
#             debug(rayonBA=self.rba)
            '''L'ouverture n'est pas précisée dans le dump
            (qui est chargé par super(Profil, self).__init__(**dump)),
            elle est donc définie par profparam'''
            iouvext, iouvint = self.profparam.iouverture
            try :
                ouvint,  ouvext  = self.cpoints[iouvint], self.cpoints[iouvext]
#                 debug(corde=self.corde, ouvint=ouvint, ouvext=ouvext)
                pouvint=100*ouvint[0]/self.corde
                pouvext=100*ouvext[0]/self.corde
#                 debug(pouvext=pouvext, pouvint=pouvint)
                self.pouverture = (-pouvext, -pouvint)
#                 debug()
            except IndexError:# as msg :
                pass
#                 debug(msg)

            '''Les noms'''
            self.splines[0].role = 'Extrados'
            self.splines[0].name = 'Extrados-'+self.name
            self.splines[1].role = 'Intrados'
            self.splines[1].name = 'Intrados-'+self.name
#     @property
#     def points(self):
#         return self.epoints
    def _getT(self, x, t0=None, nbit=False, dt=[-0.1, 1.1], tol=1.0e-10, maxiter=50):
        u"""
        :return: la valeur du paramètre t correspondant à l'abscisse |x|/100.
            Si t est non trouvé, retourne np.nan
            Plus précisement la recherche se fait dans

                - S = extrados si x>0
                - S = intrados si x<0

            le t retourné est l'UNIQUE t tel que |x|/100 == S.sx(t).
            On suppose que l'intrados et l'extrados sont des FONCTIONS de x,
            c'est à dire qu'à un x donné dans [0,100], ne correspond qu'un seul point de l'[in,ex]trados.
            Si (par exemple à l'intrados) à x correspondent DEUX points,
            alors _getT(x) retourne est le premier point trouvé.
            Il faut donc utiliser t0 pour trouver le second, mais ça n'est pas prévu pour...
            Non testé dans ce dernier cas.
        :param x: float, la position du point (dont on cherche l'abcisse curviligne t),
            en % de corde.
            La recherche se fait sur l'intrados si x<0 et sur l'extrados si x>0.
        :type x: float, 0<=x<=100.
        :param t0 : la recherche se fait par itérations de Newton, elle démarre à t=t0
        :type t0 : float.
        :param nbit :

            - si nbit est True, retourne le couple (t, nb iterations)
            - si nbit est False, retourne t.

        :type nbit : bool.
        :param dt : est l'intervalle de recherche.
        :param tol: float, paramètre passé à np.newton(), (=The allowable error of the zero value.)
        :param maxiter: int, paramètre passé à np.newton(), (=Maximum number of iterations.)

        """
        if x == 0.0 : return (np.nan, 0) if nbit else np.nan
        absx = abs(x)/100.0
        #abscisse Px du point de la corde situé à absx% de corde
        Ax, Fx = self[self.nba][0], self[0][0]
        Px = (Ax + absx*(Fx-Ax))
        if absx > 1.0 :
            raise ValueError("∣x∣=%g devrait etre dans [0,100]. C'est une abscisse en %% de corde"%abs(x))
        if x>=0 : S = self.splines[0]
        else    : S = self.splines[1]
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
            return [1+self.iba, self.profparam.nptprof-self.iba]
        except AttributeError:
            return self._nbpe

    @nbpe.setter
    def nbpe(self, nbpe):
        u"""Inutilisé, sert seulement à éviter un AttributeError"""
        self._nbpe = nbpe

    def load(self, dump):
        u"""dump est un dictionnaire comme celui retourné par self.toDump()"""
#         debug('entree load',dump_keys=dump.keys())
        super(Profil, self).load(dump)
#         debug('entree load apres super.load',dump_keys=dump.keys())
        if isinstance(dump,(Path,)) :
            f = open(dump,'r')
            d = cPickle.load(f)
            self.load(d)
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
            self.profparam = ProfsParamNew(**dparams)
        except TypeError :#as msg:
#             debug('TypeError', msg)
            self.profparam = dparams.castToProfsParamNew()
        except KeyError :#as msg :#pas de profparam dans cette liste d'arguments on prend ceux par defaut
#             debug('KeyError', msg)
            self.profparam = ProfsParamNew(*self.prefs.profparamnew)
#         debug(name=self.name, profparams=self.profparam)
        try :
            self.pouverture = dump.pop('pouverture')
        except KeyError :
#             debug('KeyError', msg)
            self.pouverture = ProfilPrefs.pouverture

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
            dump = {'points' : self.epoints,
                    'profparam' : self.profparam,#.toDump(),
                    'name' : self.name,
                    'role' : self.role,
                    }
        dump['pouverture'] = self.pouverture#ouverture en % de corde
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
        prof=Profil(profil=self)
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
            rdebug(u"Attention, recouvrement extrados/intrados : %s"%(str(pouverture)))
        if pouverture[0]*pouverture[1] < 0 :
            rdebug(u"Attention, les points ouverture ne sont pas sur le meme trados : %s"%(str(pouverture)))
        self._pouverture = pouverture
#         stack()
        if hasattr(self, '_epoints') :
            del self._epoints#pour que les points d'échantillonnage contiennent les points d'ouverture
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
        (kam, tam), (kav,tav) = self.touverture
#         debug(pouv=self.pouverture, touv=touv)
#         if kam!=kav :
        if kam == kav : #deux points sur le meme trados
            return self.splines[kam].echantillonner(nbp, ta=tam, tb=tav, mode='lineaire', onlyT=False)
        else :
            rdebug('TODO. Pour le moment les deux points de l\'ouverture doivent etre sur le meme trados' )
            return zeros((nbp,2))

    def echantillonner(self):
        u"""Echantillonnage de tout le profil, par tranche, compte tenu de l'ouverture,
        selon les prescriptions de self.profparam.
        ********************************************************************
        *******On suppose que l'ouverture est entièrement à l'intrados******
        ********************************************************************
        Le profil est constitué de 2 splines.
        - la spline extrados, du BF au BA (vrai BA, le point le plus éloigné du BF)
        - la spline intrados, du vrai BA au BF
        Les points d'échantillonnage sont constitués de 4 morceaux :
        - extrados : du BF au vrai BA
        - retour : du BA au 1er point de l'ouverture
        - ouverture : l'ouverture
        - intrados : du dernier point ouverture au BF
        """
#         epoints = super(Profil, self).echantillonner(nbp, mode)
        mode=['courbure','courbure']
        self.nb_echantillonnages += 1
        _, ri = self.rba
        if abs(ri) < self.corde/100 :
            raise RuntimeWarning('Rayon de bord d\'attaque quasi nul (%.2g) à l\'intrados, => problemes d\'echantillonnage'%ri)

        """extrados : echantillonnage standard"""

#         if profparam is None :
        pp = self.profparam
#         rdebug(pp=pp)
#         elif isinstance(profparam, ProfsParamNew) :
#             pp = self.profparam = profparam
#         else :
#             raise TypeError("J'attend un ProfsParamNew au lieu de %s"%type(profparam.__class_.__name__))
#        rdebug('name=%s'%self.name,'nb_echantillonages=%d'%self.nb_,'profparam=%s'%pp)
        Se = self.splines[0]
#         pp = self.profparam
        Text = Se.echantillonner(nbp=1+pp.iba, ta=0, tb=1, mode=mode[0], onlyT=True)
#         debug(Text)
        self.techext = Text
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
#         S = self.splines[0] if pam>=0 else self.splines[1]#Normalement pam<0, et S = intrados

        Si = self.splines[1]#pam<0, et Si = intrados
        (kam, tam), (kav, tav) = self.touverture#l'abscisse curv de pam/pav dans la spline Si
#         am, av = S(tam), S(tav)#Le point ouverture amont/aval
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
            Tret = Si.echantillonner(nbp=nbpret, ta=0, tb=tam, mode=mode[1], onlyT=True)
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
        except ValueError :#as msg :
            rdebug('nbpret=%d<0 ?'%nbpret)
            raise
#         debug(eret=eret, len_eret=len(eret))
#         debug('appel echantillonner ouverture', nbpouv=nbpouv, ta=tam, tb=tav, mode=mode[1])
        try :
            Touv = Si.echantillonner(nbp=nbpouv, ta=tam, tb=tav, mode=mode[1], onlyT=True)#echantillonnage ouverture
        except RuntimeWarning :
            rdebug('Verifier que le Rayon de bord d\'attaque est non nul')
            raise
#         rdebug(eouv=eouv, len_eouv=len(eouv))
        #Tfin = les t de l'intrados tissu
        Tfin = Si.echantillonner(nbp=nbpfin, ta=tav, tb=1, mode=mode[1], onlyT=True)#echantillonnage ouverture=>BF
#         rdebug('===> Nb points', ret=len(eret), ouv=len(eouv),efin=len(efin))
#         debug(Tret=Tret, Touv=Touv,Tfin=Tfin)
        # les trois absc. curv. Tret[0], Touv[0] et Touv[-1] désignent des points doubles
        # - Tret[0]  est l'absc. curv. du BA ; S(Tret[0]) il est deja dans eext
        # - Touv[0]  est l'absc. curv. de l'extrémité amont de l'ouverture, elle coincide avec Tret[-1]
        # - Touv[-1] est l'absc. curv. de l'extrémité avale de l'ouverture, elle coincide avec Tfin[0]
        # Tfin contient les absc. curv. des points intrados tissu
        # Tint = les t de l'intrados théorique
        Tint = np.hstack([Tret[1:],Touv[1:-1],Tfin])
        self.techint = Tint
#         debug(T_intrados=Tint)
        eint = Si(Tint)#les points intrados théorique echantillonnés
        self._epoints = np.vstack((eext,eint))
#         debug(nb_points_echantillonnage=len(self._epoints))
        return self._epoints
    @property
    def info(self):
        infos = super(Profil,self).info
        infos = [infos[k] for k in (0,1,9,10,11)]
        infos = infos + [u"*Extrados :\n     ----------------"] + self.splines[0].info
        infos = infos + [u""]
        infos = infos + [u"*Intrados :\n     ----------------"] + self.splines[1].info
        infos = infos + [u""]
        infos = infos + [u"*divers (%s) :"%self.name+u"\n     -------"]
        infos = infos + [u'%20s = '%u"pouverture"+u"%s"%(str(self.pouverture))]
        infos = infos + [u'%20s = '%u"corde"+u"%s"%(str(self.corde))]
        infos = infos + [u'%20s = '%u"profparam"+u"%s"%(str(self.profparam))]
        return infos
#
#         infos = infos + [u"fin <SplineProfil>"]

    @property
    def profparam(self) :
        return self._profparam
#
    @profparam.setter
    def profparam(self,profparam):
        u'''Les paramètres du profil :
        - ProfsParamNew (nptprof, iouvext, iouvint, iba)
        - ou ProfsParam(nptext, nptint, iba) => old fashion.
        Utiles uniquement pour l'échantillonnage'''
        if not isinstance(profparam, (ProfsParamNew, ProfsParam)) :
            raise(TypeError,u"Une instance de ProfsParamNew est attendue,impossible d'affecter un %s : %s."%(profparam.__class__.__name__,str(profparam)))
        if hasattr(self, '_profparam') and not profparam == self._profparam :
            try : del self._epoints
            except : pass
#         debug(profparam)
#         try : debug(avant=self._profparam)
#         except : debug(avant='pas de profparam')
        self._profparam=profparam.castToProfsParamNew()
#         debug(apres=self._profparam)
#         self.emit(SIGNAL('modifProfparam()'))
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
        dilate.splines[0].hardScale((1.0, self.dilatext), centre=array([0,0]))
        dilate.splines[1].hardScale((1.0, self.dilatint), centre=array([0,0]))
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
        cba, cbf = np.array(self[self.nba]), np.array(self[0])
        u = cbf-cba
        alfa = np.arctan2(u[1], u[0])
        if np.abs(alfa) > 1.0e-5 :
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

    def verification(self):
        '''Verifions que la nervure est bien dans l'ordre BF>Extrados>BA>intrados'''

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
        s='' if ok else '**** '
        reponse.append(s+u"Profil ferm'e : = %s"%ok)
        ok=dist(p[0],(1.0,0.0))<=1.0e-8
        oks.append(ok)
        s='' if ok else '**** '
        reponse.append(s+u"BF : = %s"%str(p[0]))
        ok=dist(p[p.nba],(0.0,0.0))<=1.0e-8
        oks.append(ok)
        s='' if ok else '****'
        reponse.append(s+u"BA : = %s"%str(p[p.nba]))

        # Position relative intrados, extrados
        minex=np.argmin(extr[:,1])
        maxex=np.argmax(extr[:,1])
        minin=np.argmin(intr[:,1])
        maxin=np.argmax(intr[:,1])
        yminex=np.amin(extr[:,1])
        ymaxex=np.amax(extr[:,1])
        yminin=np.amin(intr[:,1])
        ymaxin=np.amax(intr[:,1])
        reponse.append(u'Extrados :  (min, max)=(%d, %d) ; (ymin, ymax)=(%.2g, %.2g)'%(minex,maxex,yminex,ymaxex))
        reponse.append(u'Intrados :  (min, max)=(%d, %d) ; (ymin, ymax)=(%.2g, %.2g)'%(minin,maxin,yminin,ymaxin))

        # Yi intrados negatif,
        ok=(ymaxin<=1.0e-8)
        s='' if ok else '****'
        oks.append(ok)
        reponse.append(s+u'Y intrados <=0 : %s'%ok)

        # Yi extrados positif
        ok=(yminex>=-1.0e-8)
        s='' if ok else '****'
        oks.append(ok)
        reponse.append(s+u'Y extrados >=0 : %s'%ok)

        # Extrados plus épais que intrados
        ok=(abs(ymaxex)>=abs(yminin))
        s='' if ok else '****'
        oks.append(ok)
        reponse.append(s+u'Extrados plus epais que intrados : %s'%ok)

        # Les xi décroissent pour l'extrados
        dEx=extr[1:]-extr[:-1]
        dxex=dEx[:,0]#extr[1:,0] - extr[:-1,0]
        # dyex=dEx[:,1]#extr[1:,0] - extr[:-1,0]
        ok=np.all(dxex<0)
        s='' if ok else '****'
        oks.append(ok)
        reponse.append(s+u"les Xi decroissants a l'extrados : %s"%(ok))
        # Les xi croissent pour l'intrados
        dIn=intr[1:]-intr[:-1]
        dxin=dIn[:,0]#intr[1:,0] - intr[:-1,0]
        # dyin=dIn[:,1]#intr[1:,0] - intr[:-1,0]
        ok=np.all(dxin>0)
        s='' if ok else '****'
        oks.append(ok)
        reponse.append(s+u"les Xi croissants a l'intrados : %s"%(ok))

        # Convexité de l'extrados => pas significatif.
        if 0 :
            d2Ex=dEx[1:]-dEx[:-1]
            ok=np.all(d2Ex[:,1]>=0)
            if not ok :
                debug('dY/dX   extrados\n',-dEx[:,1]/dEx[:,0])
                debug('d2Y/dX2 extrados\n',d2Ex[:,1]/(dEx[1:,0])*dEx[:-1,0])
            s='' if ok else '****'
            reponse.append(s+u"Extrados concave : %s"%(ok))
            oks.append(ok)
            # Petit creux en BF à l'intrados ?
            d2In=dIn[1:]-dIn[:-1]
            ok=np.all(d2In[:,1]>=0)
            s='' if ok else '****'
            if not ok :
                debug('dY/dX   Intrados\n',dIn[:,1]/dIn[:,0])
                debug('d2Y/dX2 Intrados\n',d2In[:,1]/(dIn[1:,0])*dIn[:-1,0])
            reponse.append(s+u"Intrados convexe : %s"%(ok))
            oks.append(ok)

        oks=np.array(oks)
        if np.all(oks) :
            self.conforme=True
            reponse.append(str(oks))
            reponse.append('Tout semble OK\n****')
        else :
            self.conforme=False
            reponse.append(str(oks))
            reponse.append(u"Profil douteux. Verifiez qu'il est dans l'ordre : BF>extrados>BA>intrados")
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

    def _isLissable(self):
        ''' Dans certains cas pathologique, il n'y a pas assez de points pour lisser et pour définir l'ouverture'''
        return len(self)>=5 and len(self)-1>=self.nba>=3# \
                    # and len(self.extrados()) >=3 \
                    # and len(self.intrados()) >=3
    lissable=property(_isLissable)

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

    def plot0(self, plt, control=True, nbpd=None, nbpe=None, titre=None, show=True):
#         from matplotlib import pyplot as plt
#         from matplotlib.widgets import CheckButtons
#         plt.figure(numfig)
#         debug(epoints=self.epoints)
        if nbpd is not None : self.precision = nbpd
        else : nbpd = self.precision
        if nbpe is not None : self.nbpe = nbpe
        else : nbpe = self.nbpe
        self.echantillonner()#nbpe)
        D = self.dpoints
        C = self.cpoints
        E = self.epoints
        _, ax = plt.subplots()
        if titre is None : titre = self.name+str(self.rba)
        plt.title(titre)
        ax.plot(D[:,0], D[:,1], 'b-', label=u'spline %s'%str(self.precision), lw=1)
        ax.plot(E[:,0], E[:,1], 'bx', label=u'échantillon %s'%str(self.nbpe), lw=1)
        ax.plot(C[:,0], C[:,1], 'ro', label=u'points contrôle %s'%str(len(self)), lw=1)
        xBA, yBA = self.epoints[self.iba]
        xouvext, youvext = self.epoints[self.profparam.iouvext]
        xouvint, youvint = self.epoints[self.profparam.iouvint]
        ax.plot([xouvext, xouvint],[youvext, youvint],'yo', label=u'ouverture:%d, %d'%(self.profparam.iouvext, self.profparam.iouvint))
        ax.plot([xBA],[yBA],'go', label=u'BA:%d'%self.iba)
#         ax.plot([xBA],[yBA],'gO', label=u'BA:%d'%self.iba)
        plt.legend()
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')
        if show : plt.show()
#         return plt

    def plot(self, plt, control=True, nbpd=None, nbpe=None, titre=None, show=True):
        from matplotlib import pyplot as plt
        from matplotlib.widgets import CheckButtons
        _, ax = plt.subplots()
        if titre is None : titre = self.name#+str(self.rba)
        plt.title(titre)
        ED = self.splines[0].dpoints
        EC = self.splines[0].cpoints
        espline,  = ax.plot(ED[:,0], ED[:,1], 'b-.', lw=1)
        ID = self.splines[1].dpoints
        IC = self.splines[1].cpoints
        ispline,  = ax.plot(ID[:,0], ID[:,1], 'b-.', lw=1)
#         iechant,  = ax.plot(IC[:,0], IC[:,1], 'ro', lw=1)
#         debug(epoints=self.epoints.shape)
        AD = self.dpoints
#         debug(AD)
        if AD is None :
            if show :
                plt.title('Pas de profil')
                plt.show()
                return
            else :
                return
        all, = ax.plot(AD[:,0], AD[:,1], 'g--', lw=1)
        try :
            self.echantillonner()
        except ValueError as msg:
            debug(self, msg=msg)
        AE = self.epoints
        ech, = ax.plot(AE[:,0], AE[:,1], 'bo', lw=1)

        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        econtrol, = ax.plot(EC[:,0], EC[:,1], 'ro', lw=1)
        icontrol, = ax.plot(IC[:,0], IC[:,1], 'ro', lw=1)
        xouvext, youvext = self.epoints[self.profparam.iouvext]
        xouvint, youvint = self.epoints[self.profparam.iouvint]
        xBA, yBA = self.epoints[self.iba]
        ouverture, = ax.plot([xouvext, xouvint],[youvext, youvint],'yo', label=u'ouverture:%d, %d'%(self.profparam.iouvext, self.profparam.iouvint))
        BA, = ax.plot([xBA],[yBA],'go', label=u'BA:%d'%self.iba)
        rax = plt.axes([0.05, 0.4, 0.1, 0.5])
        labels   = ('ext', 'ext control','int', 'int control', 'All', 'Ech', 'BA', 'ouverture')
        entities = (espline,econtrol,ispline,icontrol,all,  ech,  BA,   ouverture)
        valinit  = (False,  False,   False,  False,   True, True, True, True)

        for entity, value in zip(entities,valinit):
            entity.set_visible(value)
        check = CheckButtons(rax,labels,valinit)

        def func(label):
            if label == 'ext':
                espline.set_visible(not espline.get_visible())
            elif label == 'ext control':
                econtrol.set_visible(not econtrol.get_visible())
            elif label == 'int':
                ispline.set_visible(not ispline.get_visible())
            elif label == 'int control':
                icontrol.set_visible(not icontrol.get_visible())
            elif label == 'All':
                all.set_visible(not all.get_visible())
            elif label == 'Ech':
                ech.set_visible(not ech.get_visible())
            elif label == 'BA':
                BA.set_visible(not BA.get_visible())
            elif label == 'ouverture':
                ouverture.set_visible(not ouverture.get_visible())
            plt.draw()
        check.on_clicked(func)
#         return plt
        if show : plt.show()
def testProfil(filename):
#     p = Profil(points=None, parent=None, naca=['2415', 50], name=None)
#     print p.verification()
#     print p
    from matplotlib import pyplot as plt
#     print " => Constructeur presque vide, puis append(5,10):\n    +++++++++++++++++++++++++++++++++++++"
    p = Profil()
#     debug( p)
#     try : p.appendPoint((5,10))
#     except RuntimeError as msg : rdebug(msg)
#     numfig = -1
#     print " => Constructeur filename :\n    +++++++++++++++++++++"
    debug(" => constructeur np.ndarray  :\n    +++++++++++++++++++++++")
    p = Profil(points=pointsFrom(filename))
    p.normalise()
    debug(p)
    for msg in p.verification() :
        print msg
#     p.plot(plt,titre='Profil(points=pointsFrom(filename))')

    debug(" => test toDump/load  :\n    +++++++++++++++++++++++")
    dump = p.toDump()
    debug(p)
#     pprint(dump)
    p = Profil(**dump)
#     p.verification()
#     p.plot(plt, titre='p = Profil(**p.toDump())')

#     p[46] = (3.3, -1.9)
#     p.plot(plt, titre='p[46] = (3.3, -1.9)')

#     p = Profil(points=points)
#     debug( p)
#     print " => constructeur QPolygon  :\n    +++++++++++++++++++++"
#     p = Profil(points=pointsFrom(p.qpolygon))
#     print p
    debug(" => constructeur recopie  :\n    +++++++++++++++++++++")
#     rdebug(pp=p.profparam)
#     try :
#         p.iouverture = p.iba-10,p.iba+10
#         p.plot(plt, titre='PROBLEM:p.iouverture = p.iba-10,p.iba+10')
#     except UnboundLocalError as msg:
#         rdebug('nombre de points de retour nptret=%d<0 ? '%p.profparam.nptret,msg=msg)
    debug(pp=p.profparam)
    p.iouverture = p.iba+5,p.iba+10
    p.pouverture = -1.0,-5.0
    q = Profil(profil=p)
#     debug(q)
#     rdebug(pp=p.profparam)
    debug(qp=q.profparam)
#     exit()
#     return
    debug('scaled')
    debug(p.scaled(2))
    #print "constructeur QPolygon :",p
    #print '########## Profil divers (BA, absCurv, extrados, ligne moyenne...) ##########'
#     p.removePoint(0)
#     curv = absCurv(p)
#     dcurv = curv[1:]-curv[:-1]
#     print dcurv
#     print 'p.absCurv()',absCurv(p)
    debug('########## Profil geometrie ##########')
#     rdebug(p.profparam)
#     exit()
    p.iouverture = p.iba+3,p.iba+10#sinon il sait pas echantillonner
    debug(p)
    p.plot(plt, titre='p.iouverture = p.iba+3,p.iba+10')

    p.hardScale((2,2))
    debug(p)
    pprint(p.toDump())
    p.plot(plt, titre='p.hardScale((2,2))')
    p.hardRotate(30,centre=(2,2))
    p.plot(plt, titre='p.hardRotate(30,centre=(2,2))')
    p.translate((100,100))
    p.plot(plt, titre='p.translate((100,100)')
#     print p.verification()
    debug( p)
    p.normalise()
    debug(p)
    p.plot(plt, titre='normalise()')
    p[1] = (0.6,0.12)
    p.plot(plt, titre='p[1] = (0.6,0.12)')
    p.iouverture = p.iba+2,p.iba+5
    p.plot(plt, titre='ouverture = %d, %d'%(p.iba+2,p.iba+5))
    p.insertPoint((0.23,0.16))
    p.plot(plt, titre='insertPoint((0.23,0.16))')
    p.removePoint(1)
    p.plot(plt, titre='p.removePoint(1)')
    plt.show()
#     return

    p.iouverture=p.nba+2,p.nba+3
    debug( p)
#     exit()
    debug(p.points)
    debug( p)
#     filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
    p=Profil(points=pointsFrom(filename))
    mp=Profil(profil=p)
    centre=array((7,0))
    mp.hardRotate(30,centre)
    mp.hardScale((0.5,0.5),centre)
    p.dump(Path(RUNS_DIR,'profildump.pkl'))
    f=open(Path(RUNS_DIR,'profildump.pkl'),'r')
    d=cPickle.load(f)
    pprint(d)
    p.load(d)
    debug(p)
    p.normalise()
    debug(p)
    rdebug('################## Fin testprofil ##################')
#     exit()

def testElaguer(filename):
    from matplotlib import pyplot as plt
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba)
    debug('Extrados : Se\'(1.0)=%s'%(p.splines[0](1.0,1)))
    debug('Intrados : Si\'(0.0)=%s'%(p.splines[1](0.0,1)))
    debug('sinuosite=',p.splines[1].integraleCourbure(0.01,1,1000))
#     return
    p.elaguer(eps=1, replace=True)
    debug(rba=p.rba)
    debug('sinuosite=',p.splines[1].integraleCourbure(0.01,1,1000))
    p.plot(plt, nbpd=[1000,1000],titre='elagage')
    plt.show()

def testDivers(filename):
    if 'spl' in filename.ext :
        pass
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba, corde=p.corde)
    p.normalise()
    debug(rba=p.rba, corde=p.corde)

def testSaveAndOpen(filename):
    from matplotlib import pyplot as plt
    if '.spl' in filename :
        S = Profil(points=pointsFrom(filename))
#         with open(filename,'r') as f :
#             lines = f.readlines()
#         for line in lines :
#             dump = eval(line)
#             S = Profil(**dump)
        print S
    else :
        S = Profil(points=pointsFrom(filename),
#                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                      mode=['courbure'],
                      precision=[3000])
    S.normalise()
    fname = Path(filename.dirname, filename.namebase+'Test.spl')
    debug(fname)
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = Profil(**dump)
    S1.elaguer(1, True)
    S1.rba = (-1,-1)
    S1.plot(plt,titre='rba=%s'%str(S1.rba))
    plt.show()
    dump  = S.toDump()
    dump1 = S1.toDump()
    debug('dump==dump1 ?', dump==dump1)
    pprint(dump1)
    debug('S.rba   = %s\n'%str(S.rba)+\
          '    S1.rba  = %s\n'%str(S1.rba)+\
          '    S.erba  = %s\n'%str(S.erba)+\
          '    S1.erba = %s\n'%str(S1.erba)+\
          '    S.drba  = %s\n'%str(S.drba)+\
          '    S1.drba = %s'%str(S1.drba))


def testEchantillonner(filename):
    from matplotlib import pyplot as plt
    P = Profil(points=pointsFrom(filename),
#                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
#                       mode=['linear'],
                      precision=[3000])
    P.normalise()
    P.elaguer(2, True)
    P.pouverture = pam, pav = -10, -40#en % de corde
    debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)
    debug(ouverture=P.ouverture)

    debug('P._getT(%.f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
#     P.echantillonner()
    P.plot(plt, titre='echantillonage : pouv=%s, touv=%s'%(str(P.pouverture),str(P.touverture)))
    P.hardScale(2, centre=array([0,0]))
    P.hardRotate(180, centre=(0,0))
#     debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)#,P_t=sint(t))
    debug(ouverture=P.ouverture)
    debug('P._getT(%.1f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
    P.plot(plt,titre='echantillonage : rotation 180')#%(str(P.pouverture),str(P.touverture)))

def testOuverture(filename):
    from matplotlib import pyplot as plt
    P = Profil(points=pointsFrom(filename), name='Ouverture par defaut')
    P.plot0(plt, titre='testOuverture:%s'%P.name)
    pouv = (-20, -50)
    P.pouverture = pouv
    P.plot0(plt, titre='ouverture=%s'%(str(pouv)))

if __name__=="__main__":
    ##########
    sys.path.append(Path('..').abspath)
    p = Profil()#constructeur vide
    debug('Constructeur Vide', p)
    filename = Path(VALIDATION_DIR,'1.gnu')
    filename = Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename = Path(VALIDATION_DIR,'faial2.dxf')
    filename = Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    filename = Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename = Path(VALIDATION_DIR,'simple','profil.gnu')
    filename = Path(RUNS_DIR,'profils','D2v5p2.pts')
    filename = Path(VALIDATION_DIR,'reference.pts')
    filename = Path(VALIDATION_DIR,'P0.spl')
    if 0 : testOuverture(filename)
#     exit()
    if 0 : testProfil(filename)
    if 0 : testEchantillonner(filename)
    if 0 : testElaguer(filename)
    if 0 : testSaveAndOpen(filename)
    if 1 : testDivers(filename)
