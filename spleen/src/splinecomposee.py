#!/usr/local/bin/python2.7
# encoding: utf-8
u'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSplineComposee
Description :
@author:      puiseux
@copyright:   2016 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
import math
from splineabstraite import NSplineAbstract, arrange
from matplotlib.widgets import CheckButtons
from preferences import SplinePrefs
from numpy import (asarray, linspace, log, vstack, zeros, ndarray, empty, nan,
                   abs, where, isnan, isfinite)
from collections import Iterable
from splinesimple import NSplineSimple
from utilitaires import (dist2, dist, segmentPlusProche,debug, rdebug, className)
from matplotlib import pyplot as plt
class NSplineComposee(NSplineAbstract):
    u"""
    Une spline composée de plusieurs splines simples non périodiques (des brins),
    pour points anguleux, ou profils
    Les paramètres doivent être :
    - cpoints est l'ensemble des points de contrôle
    - rupture est une liste des numéros de points anguleux, dans cpoints,
        commençant par kr0=0 et finissant par -1 ou len(cpoints)-1
        une SplineComposee comporte ns=len(ruptures)-1 splines simples (non périodiques)
    - methode est la liste des méthodes des ns splines
    - nbpd la liste des ns nombre de points pour visualiser les ns splines
    - mode est la liste des ns mode de répartition des points d'échantillonnage
    - nbe la liste des ns nombres de points d'échantillonnage
    - Plus les autres paramètres de SplineAbstract.load() (name, role, ...)
    - """

    prefs = SplinePrefs
    class Default(dict) :
        u"""Un dictionnaire avec les valeurs par défaut"""
        def __init__(self) :
#             prefs = SplinePrefs
#             dict.__init__(self,#spline composee par defaut
#                        _cpoints  = zeros((0,2)),
#                        _ruptures = [0],#les points de séparation des splines
#                        #les methodes pour chaque brin de spline
#                        _methodes = [],
#                        _nbpds    = [],
#                        _modes    = [],#polyligne
#                        _nbpes    = [],
#                        name      = 'NSplineComposee',#+'(%d)'%(id(self))
#                        role      = 'NSplineComposee',
#                        nbspline  = 0,# nb appels à computeSpline
#                        nbdpoint  = 0,# nb calculs dpoints
#                        nbech     = 0,# nb echantillonnages (epoints)
#                        )
#
            self.cpoints = zeros((0,2))
                        #les points de séparation des splines
            self.ruptures = [0]
                        #les methodes pour chaque brin de spline
            self.methode = []
            self.mode    = []
            self.nbpe    = []
            self.nbpd    = []
            self.name      = 'NSplineComposee'
            self.role      = 'NSplineComposee'
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

    def __init__(self, **dump):
        u"""
        :param cpoints: doit contenir tous les points de toutes les splines-composantes
        """
        self.splines = []
        #appel defaultet load(dump)
        super(NSplineComposee, self).__init__(**dump)
        #suppression des intermediaires à la construction
        self._update()

    def brinUnique(self, S0):
        u"""S0 de type NSplineSimple => brin unique de self"""
        self.splines   = [S0]
        self._ruptures = [0, -1]
        self._cpoints  = S0._cpoints#.copy()?
        self._methodes = [S0._methode]
        self._modes    = [S0._mode]
        self._nbpes    = [S0._nbpe]
        self._nbpds    = [S0._nbpd]
        self.name      = S0.name#+'-#0'
        self.gparent   = S0.gparent
        self.role      = S0.role
        self._update()

    def load(self, dump):
        u"""
        """
        #on vire tout
        if self.splines :
            while self.splines :
                self._update()#ne pas bouger, a faire avant
                spline = self.splines.pop()
                del spline
                try : del self._cpoints
                except AttributeError : pass
            # on doit réinitialiser aux valeurs par défaut
            #qui ont été supprimées par le _update ci-dessus
            for key, value in self.Default().iteritems() :
                setattr(self, key, value)

        arrange(dump)#renommage de cles mal nommées (precision, classename)
        keys = dump.keys()
        key = 'classname'
        if 'classname' in keys :
            clsname = dump['classname']
            if clsname == 'NSplineSimple' :
                #dump est une spline simple,
                #on la charge comme spline self.splines[0]
                #self est à une seule composante
                debug(u"Conversion %s => %s"%(clsname, className(self)))

                S0 = NSplineSimple(**dump)
                return self.brinUnique(S0)
# #             else :
# #                 msg = u'Pas de conversion %s => %s'%(clsname, className(self))
# #                 raise NotImplementedError, msg
#         elif not hasattr(self, '_ruptures') :
#             #on remet les valeurs par defaut qui on disparu dans le _update plus haut
#             dump1 = {}
#             dump1.update(dump)#On conserve une copy de dump
#             default  = self.Default()
#             for key, value in default.iteritems() :
#                 setattr(self, key, value)
#             dump.update(default.dump)
#             dump.update(dump1)

#       else : #ici on n'a pas trouve 'classname'
        if 'name'      in keys : self.name       = dump.pop('name')
        if 'gparent'   in keys : self.gparent    = dump.pop('gparent')
        if 'role'      in keys : self.role       = dump.pop('role')
        if 'ruptures'  in keys : self._ruptures  = dump.pop('ruptures')
        if 'methode'   in keys : self._methodes  = dump.pop('methode')#type spline
        if 'mode'      in keys : self._modes     = dump.pop('mode')#mode echantillonnage
        if 'nbpe'      in keys : self._nbpes     = dump.pop('nbpe')#nbp echantillonnage
        if 'nbpd'      in keys : self._nbpds     = dump.pop('nbpd')#nbp discretisation
        if 'precision' in keys : self._nbpds     = dump.pop('precision')#idem
        ns = len(self._ruptures)-1#nb splines
        methodes   = self._methodes
        modes      = self._modes
        nbpes      = self._nbpes
        nbpds      = self._nbpds
#         debug((len(methodes), len(nbpds), len(modes), len(nbpes),), ns=ns, )
        if (len(methodes), len(nbpds),len(modes), len(nbpes),) != 4*(ns,) :
            raise ValueError, u"""probleme de dimensions, les parametres 'methode', 'mode',
            'nbpd', 'nbpe' doivent etre des tuple ou des listes a %d elements
            'ruptures' doit etre une liste ou un tuple a %d elements.
            ruptures=%s, methodes=%s, modes=%s, nbpes=%s, nbpds=%s
            """%(ns, 1+ns, self._ruptures, methodes, modes, nbpes, nbpds)
        u"""
        maintenant il y a une methode, mode, ... par brin
        On charge les points de contrôle"""
#
        try : cpoints = self._cpoints
        except AttributeError : cpoints = zeros((1,2))
        for key in ('cpoints', 'points') :#la cle 'cpoints' est prioritaire, en cas
            if dump.has_key(key):
                cpoints = dump.pop(key)
                if isinstance(cpoints, (list, tuple, ndarray)) :#types Python de base
                    #obligé de considerer ce cas, car pour le load d'un pkl cpoints est une liste
                    cpoints = asarray(cpoints)
                if not isinstance(cpoints, ndarray) :
                    msg = u'Les points de controle doivent etre de type ndarray((*,2)), et non pas %s.'%cpoints.__class__.__name__
                    msg = msg + u' \nUtilisez la fonction points=pointsFrom(points) pour transformer a peut pres n\'importe quoi en ndarray((n,2)).'
                    msg = msg + u' \nLes points de controle ne doivent pas avoir de points doubles consecutifs.'
                    msg = msg + u' \nUtilisez points = eliminerPointsDoublesConsecutifs(points) en cas de doute.'
                    raise TypeError, msg
                break#la cle 'cpoints' est prioritaire,
        if len(cpoints)>1 and len(self._ruptures) == 1 :
            #c'est une initialisation avec cpoints, et rupture, ...
            # avec les valeurs de Default. => on met une NSplineSimple unique
            S0 = NSplineSimple(cpoints=cpoints, name=self.name)
            self.brinUnique(S0)
            return

        ok = self._ruptures[0] == 0 and self._ruptures[-1] in (0,-1,len(cpoints)-1)
        for r in self._ruptures[1:-1] :
            ok = ok and 0 < r < len(cpoints)-1
        if not ok :
            msg = [
                    u'Probleme avec les points de rupture %s:'%self._ruptures,
                    u'Le premier devrait etre 0, le dernier devrait etre 0 ou %d ou -1'%(len(cpoints)-1),
                    u'Les autres devraient etre dans l\'ordre strictement croissant, et dans l\'intervalle ]0,%d[ '%(len(cpoints)-1),
                    u"len(cpoints)=%d, ruptures=%s, methodes=%s, modes=%s, nbpes=%s, nbpds=%s"%(len(cpoints), self._ruptures, methodes, modes, nbpes, nbpds)
                   ]
            raise ValueError('\n'.join(msg))
        self._ruptures[-1] = len(cpoints)-1
        for k in range(len(methodes)) :
            deb, fin = self.ruptures[k], self.ruptures[k+1]
            self.splines.append(NSplineSimple(
                                            #Attention a bien mettre une copy de cpoints[deb:1+fin]
                                            #sinon, en changeant (hardScale p.ex.) la spline 1, la spline 2 change aussi
                                            #car elles on un point commun. Ce point serait modifié deux fois
                                            #C'est un bug très difficile à trouver !!!
                                            cpoints   = cpoints[deb:1+fin].copy(),
                                            methode   = methodes[k],
                                            nbpe      = nbpes[k],
                                            mode      = modes[k],
                                            nbpd      = nbpds[k],
                                            name      = self.name+'#%d'%k,
                                            role      = 'spline#%d'%k))
#         try : del self._ruptures
#         except AttributeError : pass

    def elaguer(self, eps=0.5, replace=False):
        A = []
        for spline in self.splines :
            a = spline.elaguer(eps, replace)
            A.append(a)
        self._update()
        try : del self._ruptures
        except AttributeError : pass
        return A

    def echantillonner(self, nbp=None, mode=None):
        if mode is None : mode = self.mode
        if nbp is None : nbp=self.nbpe
        epoints = self.splines[0].epoints
        for n, m, spline in zip(nbp[1:], mode[1:], self.splines[1:]) :
            spline.nbp, spline.mode = n, m
            #on ne met pas le premier epoint sinon il serait en double
            #avec celui de la spline précédente.
            epoints=vstack([epoints,spline.epoints[1:]])
        self._epoints = epoints
        self._tech = [s._tech for s in self.splines]
        return epoints

    @property
    def epoints(self):
        if hasattr(self, '_epoints') :
            return self._epoints
        else :
            self.echantillonner()
            return self._epoints

    @property
    def tech(self):
        u"""Les temps T de l'echantillonnage"""
        if not hasattr(self, '_tech') :
            self.echantillonner()
        return self._tech

    def __pointsFrom(self, listearrays):
        u"""
        C'est en réalité une fonction plutôt qu'une méthode (self n'est pas utilisé)
        mais à utiliser uniquement pour NSplineComposee.
        :param listearrays: une liste de ndarray de shape (*,2)
        Recompose un ndarray((n,2)) de points à partir des [c,e,d]points des brins de la spline.
        Les points de part et d'autre de chaque rupture ne sont présent qu'en un seul exemplaire
        dans le tableau recomposé."""
        if not isinstance(listearrays, (list, tuple)) :
            raise TypeError('Une liste de ndarray((n,2)) est attendu au lieu de %s'%listearrays.__class_.__name__)
        if len(listearrays)==0 :
            return zeros((1,2))
        else :
            points = listearrays[0]
        for a in listearrays[1:]:
            if a is not None and len(a) > 0 :
                points = vstack((points, a[1:]))
        return points

    def __listFrom(self, points):
        u'''
        C'est bien une méthode de NSplineComposee, car self est utilisé.
        :param points: un ndarray de shape (*,2)
        Décompose [c,d,e]points en liste de ndarray((*,2)) pour les splines composantes.
        Les points de ruptures sont en double : dans chaque brin de part et d'autre de la rupture.
        Cette fonction est une sorte de réciproque à droite de __pointsFrom()
        >>> points == self.__pointFrom(self.__listFrom(self.cpoints)) => vrai
        '''
        lp = []
        ruptures = self.ruptures
        for k in range(len(self.methode)) :
            deb, fin = ruptures[k], ruptures[k+1]
            lp.append(points[deb:1+fin])
        return lp

    def isClosed(self, eps=0.0):
        try : return dist2(self[0], self[-1]) <= eps*eps
        except IndexError : return None

    def close(self):
        """
        Fermeture de self._cpoints.
        Appelé automatiquement si self.methode[0] = 'periodic'
        :return : int,
            # 0 si rien touché,
            # 1 si ajustement mineur (<= self.eps) de self[-1]
            # 2 si ajout de point self[-1]=self[0]
        """
        if self.isClosed(0) :
            msg = u"    cpoints deja ferme"
        elif dist2(self[0],self[-1]) <= self.eps:
            #presque fermé, on n'ajoute pas de point
            p1, p2 = self[0], self[-1]
            self._cpoints[-1] = self._cpoints[0]
            msg = u"""    les points self[0]=%s et self[-1]=%s sont presque identiques.
            Ils sont à une distance d=%.3g
             => self[-1] est (legerement) modifie.
            """%(p1,p2, dist(p1,p2))
        else :
            #On ajoute un point
            p1, p2 = self[0], self[-1]
            self._cpoint.append[self[0]]
            msg = u"""    Les points self[0]=%s et self[-1]=%s sont distincts.
            Ils sont à une distance d=%.3g
             => le point self[0] a ete rajoute en fin de spline.
            """%(p1, p2, dist(p1,p2))
        debug(msg)

    def __len__(self):
        return len(self.cpoints)

    @property
    def cpoints(self):
        u"""Les points de contrôle, sous forme ndarray((n,2)),
        Les points de contrôle sont dans les splines composantes.
        Le self._cpoints et self._qcpolygon de NSplineAbstract sont inutiles.
        On n'a pas besoin de les créer"""
        return self.__pointsFrom([s.cpoints for s in self.splines])

    @property
    def dpoints(self):
        u"""
        Les points discretises, recalculés seulement en cas d'update
        On va les chercher dans les splines composantes.
        """
#         self._dpoints =
        return self.__pointsFrom([s.dpoints for s in self.splines])

    @property
    def nbpd(self):
        return [s.nbpd for s in self.splines]

    @nbpd.setter
    def nbpd(self, prec):
        self._nbpds = prec
        for k, spline in enumerate(self.splines) :
            spline.nbpd = prec[k]
    @property
    def nbpe(self):
        return [s.nbpe for s in self.splines]
    @nbpe.setter
    def nbpe(self, nbpe):
        self._nbpes = nbpe
        for k, spline in enumerate(self.splines) :
            spline.nbpe = nbpe[k]

    def __str__(self):
        return u'\n'.join(self.info)

    @property
    def info(self):
        infos=[
                u"<%s>"%className(self),
                u'%25s = '%u'name'+u'%s'%self.name,
                u'%25s = '%u'role'+u'%s'%self.role,
                u'%25s = '%u'nb pts controle'+u"%d"%len(self.cpoints),
                u'%25s = '%u'closed'+u'%s'%self.isClosed(),
                u'%25s = '%u'nb_splines, n_bech'+u"%s, %s"%(self.nbspline, self.nbech),
                u'%25s = '%u'methode'+u"%s"%str(self.methode),
                u'%25s = '%u'nb pts discretisation'+u"%s"%str(self.nbpd),
                u'%25s = '%u'mode echantillonage'+u"%s"%self.mode,
                u'%25s = '%u'nb pts echantillon'+u"%s"%self.nbpe,
#                 u'%25s = '%u'nb updates'+u"%s"%self.nbupdate,
                ]
        i = u'largeur, hauteur'
        try :
            i1 = u"%g, %g"%(self.width, self.height)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)
#
        i = u'position cg'
        try :
            i1 = u"%g, %g"%(self.centregravite[0], self.centregravite[1])
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)
#
        i = u'longueur'
        try :
            i1 = u"%g"%(self.longueur('r'))
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)

        infos1 = [
#             u"nb updates          = %s"%self.nbupdate,
            u'%25s = '%u"ruptures"+"%s"%self.ruptures,
            u'%25s = '%u"composantes"+"%s"%[s.name for s in self.splines],
            u'%25s = '%u"nb pts ctl par comp"+"%s"%[len(s) for s in self.splines]
            ]
        return infos + infos1

        # return
#     def debug(self):
# #         return self.__str__()
#         debug(self)
        # for k, spline in enumerate(self.splines) :
        #     print "\nSpline composante numero %d"%k
        #     print 27*"="
        #     print spline

    def toDump(self):
        return {
                'classename' : className(self),#besoin de savoir quel type de spline.
                'cpoints'    : self.cpoints.tolist(),
                'role'       : self.role,
                'name'       : self.name,
                'methode'    : self.methode,
                ############################
                'nbpd'       : self.nbpd,
                'nbpe'       : self.nbpe,
                'ruptures'   : self.ruptures,
                'mode'       : self.mode
                }
################################################################
    @property
    def methode(self):
        u"""A l'initialisation, c'est dans self._methode, ensuite, c'est recalculé de manière dynamique,
        a partir des brins de la spline"""
        return [s.methode for s in self.splines]

    @property
    def mode(self):
        return [s.mode for s in self.splines]

    @mode.setter
    def mode(self, newmode):
        u"""
        Changement mode d'echantillonnage
        """
        if newmode==self.mode :
            return
        else :
            self._mode = newmode
            try : del self._epoints
            except AttributeError : pass
            try : del self._tech
            except AttributeError : pass


    def computeSpline(self, methodes):
#         #debug(self.splines)
        for k, spline in enumerate(self.splines) :
#             deb, fin = self.ruptures[k], self.ruptures[k+1]
            spline.computeSpline(methodes[k])

    def __getitem__(self,k):
        return self.cpoints[k]

    def __setitem__(self,k,value):
        u"""
        On bouge un point de contrôle existant :
        mise à jour du polygone de contrôle et _update
        """
        #On parcourt les splines, on modifie celle(s) qui est (sont) concernée(s)
        ruptures = self.ruptures
        for ks in range(len(self.splines)) :
            deb, fin = ruptures[ks], ruptures[ks+1]
            if deb <= k <= fin : #k est dans la spline ks
                debug(ks=ks, k=k-deb, value=value)
                self.splines[ks][k-deb] = asarray(value)
                # ici ne pas faire de break, car k peut être dans 2 splines
                # p.ex. k==fin de la spline 0 et k==deb de la spline 1
                # il faut modifier les deux splines 0 et 1
        self.cpoints[k] = value
        self._update()#indispensable pour mettre à jour self.cpoints et dpoints
    @property
    def ruptures(self):
        u"""A l'initialisation, c'est dans self._ruptures, ensuite, c'est recalculé de manière dynamique,
        a partir des brins de la spline"""
        if hasattr(self,'_ruptures') :
            return self._ruptures
        else :
            R = [0]
            for k, spline in  enumerate(self.splines):
                R.append(R[k]+len(spline)-1)
            return R

    def insertPoint(self, pos, i=None):
        u"""
        Si i = None : calcule dans quel segment [i, i+1] il est raisonable d'insérer le point pos
            et demande à la spline contenant ce segment de faire l'insertion
        si i != None, insertion du point pos entre i-1 et i (de sorte qu'il devienne le point i)
            c'est la spline simple contenant [i-1,i] qui fait l'insertion
        """
        cpoints = self.cpoints
        if i is None :
            i, _ = segmentPlusProche(cpoints, pos)#segment [i, i+1]
            if i is None :
                debug('Je ne sais pas ou (dans quel segment) inserer ce point %s'%str(pos))
                return
        else :
            i = i-1
        si  = self.localiser(i)#si=numero de la spline
        ks, kp = si[-1]
        #i est le kp-ieme point de la spline ks, H est inséré entre i et i+1
        self.splines[ks].insertPoint(pos, 1+kp)
        self._update()
        del self._ruptures
        return i

    def appendPoint(self,pos):
        i = self.splines[-1].appendPoint(pos)
        self._update()
        return self.ruptures[-2] + i

    def localiser(self, k):
        u"""
        :param k: int, est un numéro de point, dans l'intervalle [0, len(self)[
        :return : les numéros des splines auxquelles appartient le point k,
            ainsi que la position du point k dans ces splines, sous une des deux
            formes suivantes :
        - [(ks, i),] ks numero de spline, i = place de k dans cette spline
            si k appartient à une seule spline ou bien
        - [(ks0, i0), (ks1, i1)] : un point k peut appartenir à (au plus)
            deux splines qui sont alors consécutives. Dans cette hypothèse,
            k est la fin de l'une (ks0) et le début de la suivante (ks1=1+ks0)
            et i0=len(self.splines[ks0]), i1 = 0
        """
        numspl = []
        ruptures = self.ruptures
        for ks in range(len(self.splines)) :
            deb, fin = ruptures[ks], ruptures[ks+1]
            if deb <= k <= fin : #k est dans la spline ks
                numspl.append((ks, k-deb))
                #pas de break ici, il y a peut etre une seconde spline contenant k
        return numspl

    def removePoint(self, k):
        u"""
        - Suppression du point k
        """
        ls = self.localiser(k)#localise le point k : le num. de sa spline, et son numero local
        if len(ls) > 1 : #c'est un point de rupture, on joint les deux splines adjacentes et on supprime le point
            kr = self.ruptures.index(k)#le num du point de rupture (dans self.ruptures)
#             debug(ls=ls, kr=kr, ruptures=self.ruptures)
            kr = self.join(kr)#le num du point de rupture (dans self.cpoints)
            self.removePoint(kr)#recursif
            #On devrait pas passer ici...
            raise NotImplementedError(u'On devrait pas passer ici... On tente la suppression d\'un point anguleux (%d)'%k)
        else :
#         point = self[k]
            ks, kp = ls[0]#numero de spline, numero de point dans la spline
            pt = self.splines[ks].removePoint(kp)
            self._update()
            try : del self._ruptures
            except AttributeError : pass
            return pt

    def hardRotate(self, alfa, centre=None, unit='degres'):
        if unit == 'degres' :
            alfa = math.radians(alfa)
        if centre is None :
            centre = self.barycentre
        for spline in self.splines :
            spline.hardRotate(alfa, centre, unit='radians')
        # Les methodes des splines composantes changent (les valeurs des vecteurs dérivées aux extrémités)
        #on n'a pas le droit de changer self.methode...alors on contourne en changeant self._methode
        self._methode = [spl.methode for spl in self.splines]
        self._update()

    def hardScale(self, scale, centre=None):
        u'''
        modification in situ de self, par application d'une homothétie
        de centre 'centre' et de rapport 'scale'.
        Les points de contrôle eux même sont modifiés.
        :param centre: (float, float) les coordonées (x0, y0) du centre d'homothétie.
        :param scale: float ou (float, float)

            - si scale est un float, le rapport d'homothétie est (hx, hy)=(scale, scale)
            - si scale est un (float, float), le rapport d'homothétie est (hx, hy)=(scale[0], scale[1])

        chaque point de contrôle P=(x,y) est transformé en P'=(x',y') avec
                x' = centre[0] + hx*(x-centre[0]), et
                y' = centre[1] + hy*(y-centre[1])
        Chaque spline composante de selgf est transformée de manière adéquate,
        en transformant également les valeurs des dérivées aux extrémités
        lorsqu'elles sont précisées. Cf NSplineSimple.hardScale()
        '''
        if centre is None :
            centre = self.barycentre
#         hardScale(self.cpoints, scale, centre)
#         debug(cpoints=self.cpoints)
        for spline in self.splines :
            #Attention, il faut que les deux splines soient mises à l'echelle avec le meme centre
            spline.hardScale(scale, centre)
#             debug(spline.cpoints)
        # Les methodes des splines composantes changent (les valeurs des vecteurs dérivées aux extrémités)
        #on n'a pas le droit de changer self.methode...alors on contourne en changeant self._methode
        self._methode = [spl.methode for spl in self.splines]
        self._update()

    def symetriser(self, axe, centre=None):
        if centre is None :
            centre = self.barycentre[0][axe]
        for spline in self.splines :
            spline.symetriser(axe, centre)
        # Les methodes des splines composantes changent (les valeurs des vecteurs dérivées aux extrémités)
        #on n'a pas le droit de changer self.methode...alors on contourne en changeant self._methode
        self._methode = [spl.methode for spl in self.splines]
        self._update()

#     def copy(self):
#         u"""retourne une copie de self"""
#         dump = self.toDump()
#         return NSplineComposee(**dump)

    def translate(self,vecteur):
        for s in self.splines :
            s.translate(vecteur)
        self._update()

#     def absCurv(self, T, witherr=False):

    def join(self, k):
        u"""Supprimer le k-ieme point de rupture, i.e.
        joindre les deux splines k-1 et k en une seule,
        le nb de splines, ... diminue de 1"""
        if k in (0,-1,len(self.ruptures)-1) :
            debug(u'impossible de supprimer le premier ou dernier point de rupture')
            return None
        r = self.ruptures[k]
        ruptures  = self.ruptures
        methodes  = self.methode
        cpoints   = self.cpoints
        modes     = self.mode
        nbpds     = self.nbpd
        nbpes     = self.nbpe
        name      = self.name
        splines = self.splines
        ruptures.pop(k)
        methodes.pop(k-1)
        modes.pop(k-1)
        nbpds.pop(k-1)
        nbpes.pop(k-1)
        splines.pop(k-1)
        self.load(dict(cpoints=cpoints,
                  methode=methodes,
                  nbpd=nbpds,
                  mode=modes,
                  nbpe=nbpes,
                  ruptures=ruptures,
                  name=name))
        self._update()
        return r#le numero du point

    def split(self, n, mg=None, md=None):
        u"""Découpage de cpoints au point n, avec méthodes mg à gauche et md à droite
        méthodes mg et md sont définies/expliquées dans le fichier splineabstraite.py,
        fonction computeSpline(...)
        Par défaut, on prend les methode mg=md=self.methode """
#         debug(n=n, mg_md=(mg,md), methodes=self._methodes)
        if isinstance(n, Iterable) :
            n = list(set(n).difference(set(self.ruptures)))
            n.sort(reverse=True)
            for i in n :
                self.split(i, mg, md)
#             return
        elif isinstance(n, int) :
            ruptures  = self.ruptures
            methodes  = self.methode
            cpoints   = self.cpoints
            modes     = self.mode
            nbpds     = self.nbpd
            nbpes     = self.nbpe
            for k, r in enumerate(ruptures) :
                if r>n or r in (-1, len(self)-1):
                    break
            k -= 1
            #n est entre les ruptures k et k+1. Elle tombe dans la k-ième spline
            #On coupe en deux la spline S=self.splines[k].
            # S => Sg,Sd. deux demi splines
            # apres l'operation, Sg=self.splines[k], Sg=self.splines[k+1]
#             debug(n=n, k=k)
            mode, nbpe, prec = modes[k], nbpes[k], nbpds[k]
#             mg, md = methodes[k][:], methodes[k][:]
            if mg is None : #la methode de la demi spline de gauche Sg
                mg = methodes[k][:]#une copie
                #Si les derivées sont précisées, on transforme en (2,0,0), la condition de droite
                if mg[0] in ('cubic',) and isinstance(mg[1], (list, tuple)):
                    mg = (mg[0], (mg[1][0], (2,0,0)))

            if md is None : #la methode de la demi spline de droite Sd
                md = methodes[k][:]#une copie
                #Si les derivées sont précisées, on transforme en (2,0,0), la condition de gauche
                if md[0] in ('cubic',) and isinstance(md[1], (list, tuple)):
                    md = (md[0], ((2,0,0), md[1][1]))
    #         mg = methodes[k-1] if mg is None else mg
    #         md = methodes[k] if md is None else md
            ruptures.insert(k+1, n)
    #         nbpcs     = diff(asarray(ruptures))#le nb de points de contrôle des brins
    #         nbpcinf3  = where(nbpcs<3)[0]#les numeros des brins avec 1 ou 2 points de controle points
    #         debug(rupture=n, k=k, methodes_g_d_avant=(mg,md))
    #         if k in nbpcinf3 : #moins de 3 points, on fait une spline 'ius',1
    #             mg = 'ius',1
    #             debug('>>>>> k in nbpcinf3',k=k,nbpcinf3=nbpcinf3)
    # #         debug(k=k,n=n,ruptures=ruptures,nbpcinf3=nbpcinf3, nbpcs=nbpcs, mg=mg)
    #         debug(rupture=n, k=k, methodes_g_d_apres=(mg,md))
    #         exit()
            methodes.insert(k, mg)
            methodes[k+1] = md
            nbpds.insert(k, prec)
            modes.insert(k, mode)
            nbpes.insert(k, nbpe)
    #         self._ruptures = ruptures
    #         self._methode = methodes
#             debug(n=n, methodes=methodes)
            self.mode = modes
            self.nbpe = nbpes
#             return methodes, nbpds, modes, nbpes, ruptures
        #pour finir :
#             dd  = dict(cpoints=cpoints,
#                        methode=methodes,
#                        nbpd=nbpds,
#                        mode=modes,
#                        nbpe=nbpes,
#                        ruptures=ruptures,
#                        classname=className(self))
#             debug()
#             pprint(dd)
            self.load(dict(cpoints=cpoints,
                           methode=methodes,
                           nbpd=nbpds,
                           mode=modes,
                           nbpe=nbpes,
                           ruptures=ruptures,
                           classname=className(self),
                           name=self.name))
            for s in self.splines :
                #Pour s'assurer que le update est bien fait car on a affecté
                # s.methode avec s._methodes, sans passer par le setter methode.setter
#                 mav = s.methode
                s.methode = s.methode
#                 if s.methode != mav:
#                     debug(avant=mav,apres=s.methode)
    #         debug(self)
            self._update()
#         raise NotImplementedError

    def _update(self):
        u"""
        Suppression de tous les attributs self._xxx volatiles, ils sont
            recalculés à la demande i.e. quand on appelle self.xxx
        Est appelé à chaque modification
        - (géométrique) d'un point de contrôle de la spline
        - ou bien du PolygonF de base
        - ou bien de methode spline (cubic, IUS, US,...), ou des dérivées aux extremites
        - ou bien de mode d'échantillonage
        ultra privée, ne pas l'appeler de l'extérieur.
        """
#         super(NSplineComposee, self)._update()
        try : del self._epoints
        except AttributeError : pass
        try : del self._dpoints
        except AttributeError : pass
        try : del self._height
        except AttributeError : pass
        try : del self._width
        except AttributeError : pass
        try : del self._longueur
        except AttributeError : pass
        try : del self._dac
        except AttributeError : pass
        try : del self._sx
        except AttributeError : pass
        try : del self._sy
        except AttributeError : pass
        try : del self._tech
        except AttributeError : pass

#         try : del self._cpoints#c'est juste un intermediaire à la construction
#         except AttributeError : pass
#         try : del self._ruptures#c'est juste un intermediaire à la construction
#         except AttributeError : pass

    def plot(self, figure=None, aspect={}, titre=None, more=[], texts=[], show=True,
             buttons=True, numbers=['3p']):
        """
        :param figure: une instance de matplotlib.figure.Figure
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
        :param texts : list, [(x,y, 'text', style=, fontsize=...), ...] texts en x,y
        :param numbers: list ou set ou tuple ['3c','12e','100d'] numéroter les points
            controle, echantillon, discretisation '3c' signifie que l'on numérote
            les cpoints par pas de 3, '100d' les dpoints par pas de 100 etc ...
        :param show: bool, True pour affichage de TOUS les plots en attente, False sinon
        :return figure: la figure passée en argument ou nouvelle figure si None

        """
    #     renderer = RendererBase()
        defaultaspect = {
                    'c':'ro',#control :        red
                    'd':'b-',#discretisation : blue
                    'e':'g.',#echantillon :    green
                    'r':'k*',#rupture:         black
                  }
        defaultaspect.update(aspect)
        if figure is None :#pas de figure => on la crée
            figure = plt.figure('plot(self)')
#         debug(figure=figure)

    #     if nbpd is None : nbpd = self.nbpd
    #     if nbpe is None : nbpe = self.nbpe
    #     if mode is None : mode = self.mode
        D = self.dpoints
        C = self.cpoints
        E = self.epoints
#         debug(E.shape)
    #     exit()
        axes = figure.get_axes()
        if axes : ax = axes[0]
        else : ax = figure.subplots()
    #    debug(axes=axes, ax=ax,gca=plt.gca())

    #     ax = figure.subplots()
#         debug(ax=className(ax),figure=className(figure))
        if titre : ax.set_title(titre)
#         ax.set_title('prout')
    #         titre = self.name+str(self.methode)

        fmtc, fmtd, fmte, fmtr = (defaultaspect['c'], defaultaspect['d'],
                                  defaultaspect['e'], defaultaspect['r'])
    #     debug('fmtc=%s, fmtd=%s, fmte=%s'%(fmtc, fmtd, fmte))
        if fmtc : control,     = ax.plot(C[:,0], C[:,1], fmtc, lw=1, label=u'Contrôle')
        if fmtd : spline,      = ax.plot(D[:,0], D[:,1], fmtd, lw=1, label=u'Spline')
        if fmte : echantillon, = ax.plot(E[:,0], E[:,1], fmte, lw=1, label=u'Échantillon')
        R = asarray(self[self.ruptures])#des points
        fmtr = defaultaspect['r']
        if fmtr : ruptures,     = ax.plot(R[:,0], R[:,1], fmtr, markersize=12, label=u'Rupture')

    #     for x, y, color, name in more:
    #         _, = ax.plot(x, y, color, label=name)
        figure.legend()

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
        ax.axis('equal')

        if buttons :
    #         figure.add_subplot('112')
    #         rax = plt.axes([0.05, 0.4, 0.1, 0.15])
            rax = figure.add_axes([0.05, 0.4, 0.1, 0.15])
#             debug(rax=className(rax))#,fsubplot=className(fsubplot))
            labels = [u'control', u'spline',u'echantillon']
            values = [True, True, True]
            draws = [control, spline, echantillon]
            if isinstance(self, NSplineComposee) :
                labels.append(u'ruptures')
                values.append(True)
                draws.append(ruptures)

        #         for x, y, color, name in more:
        #             temp, = ax.plot(x, y, color)
        #             if not name : continue
        #             if buttons :
        #                 draws.append(temp)
        #                 labels.append(name)
        #                 values.append(True)
            figure.subplots_adjust(left=0.2)
            for item in draws :
                item.set_visible(True)
    #     if buttons :
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
    #         figure.subplots_adjust(left=0.2)
        axes = figure.get_axes()
#         debug(axes=axes)
        if show : plt.show()
        return figure

    def plotCourbure(self):
        from matplotlib.widgets import CheckButtons
#         D = self.dpoints
#         C = self.cpoints
#         _, ax = plt.subplots()
        titre = self.name+' courbure'
        plt.title(titre)
#         spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
#         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
#         control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        n = 1001
        ns = len(self.splines)
        T0 = linspace(0, 1, n+1)[1:-1]
        Ts = linspace(0, float(ns), ns*n+1)
        Cs = empty(ns*n+1)
        Cs[:] = nan
        for k,s in enumerate(self.splines):
#             debug(methode=s.methode)
            Cs[1+k*n:(k+1)*n] = log(abs(s.courbure(T0)))
#             Cs[1+k*n:(k+1)*n] = s.courbure(T0)
        plt.plot(Ts, Cs)
        Tnan = Ts[where(isnan(Cs))]
        Y = empty(len(Tnan))
        Y[:] = min(Cs[where(isfinite(Cs))])
        plt.plot(Tnan,Y,'ro')
        plt.subplots_adjust(left=0.2)

        plt.show()
        return plt

    def longueur(self, p='r'):
        if p=='r':#longueur vraie
            if hasattr(self, '_longueur') :
                return self._longueur
            else :
                self._longueur = sum([s.absCurv(1) for s in self.splines])
                return self._longueur
        else :
            return super(NSplineComposee, self).longueur(p)

if __name__=="__main__":
    from testsplinecomposee import testMain
    import config
    config.TEST_MODE = False#pour avoir les liens clickables
    testMain()