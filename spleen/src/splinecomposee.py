#!/usr/local/bin/python2.7
# encoding: utf-8
from utilitaires.utilitairesdivers import diff
u'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSplineComposee
Description :
@author:      puiseux
@copyright:   2016 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
import sys,os,math
import numpy as np
from numpy import (asarray, linspace, log, vstack, zeros, ndarray, empty, nan,
    abs, where, isnan, zeros_like, isfinite)
from collections import Iterable
from splinesimple import NSplineSimple
from utilitaires import (dist2, dist, XY, segmentPlusProche,debug, rdebug,
                         className)
from matplotlib import pyplot as plt

class NSplineComposee(NSplineSimple):
    # prefs = preferences.SplinePrefs()
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
    def __init__(self, **dump):#cpoints=None, parent=None, ruptures=[0,-1], methode=[('cubic',"not-a-knot"),]):
        u"""points doit contenir tous les points de toutes les splines-composantes"""
#         debug(dump=dump)
        self.splines = []
#         debug(dump=dump)
        #super(...) appelle setDefaultValues() puis load()
        super(NSplineComposee, self).__init__(**dump)
        self._update()
#         rdebug('fin init', self.splines[0], '\n', self.splines[-1])
#         les deuix splines ont 0 points...
#         self.ruptures = ruptures
#         self.methode = methode
################################################################
    def setDefaultValues(self):
        """Valeurs par defaut:"""
        self._methode    = [('cubic',"not-a-knot"),]#les methodes pour chaque brin de spline
        self._ruptures   = [0,-1]#les points de séparation des splines
        self._nbpd       = [500,]
        self._mode       = ['linear']#polyligne
        self._nbpe       = [100,]
        self._cpoints    = zeros((0,2))
        self.name        = className(self)#+'(%d)'%(id(self))
        self.role        = className(self)
        self.nbspline    = 0# nb appels à computeSpline
        self.nbdpoint    = 0# nb calculs dpoints
        self.nbech       = 0# nb echantillonnages (epoints)

    def load(self, dump):
#         try : self._ruptures = dump.pop('ruptures')
#         except KeyError : pass

#         super(NSplineComposee, self).load(dump)#load tout le rest sauf cpoints
        u"""Ce qui est commun à toutes les splines.
        Attention, la cle est 'points'  dans les anciens projets, et 'cpoints' dans les nouveaux.
        OUI=>On affecte _cpoints directement qui évite l'appel au setter de cpoints (qui appelle _update()).
        cpoints DOIT être liste, tuple, array ou ndarray, sans points doubles consécutifs
        """
        keys = dump.keys()
        if 'name'      in keys : self.name      = dump.pop('name')
        if 'gparent'   in keys : self.gparent   = dump.pop('gparent')
        if 'role'      in keys : self.role      = dump.pop('role')
        if 'ruptures'  in keys : self._ruptures = dump.pop('ruptures')
        if 'methode'   in keys : self._methode  = dump.pop('methode')#type spline
        if 'mode'      in keys : self._mode     = dump.pop('mode')#mode echantillonnage
        if 'nbpe'      in keys : self._nbpe     = dump.pop('nbpe')#nbp echantillonnage
        if 'nbpd'      in keys : self._nbpd     = dump.pop('nbpd')#nbp discretisation
        if 'precision' in keys : self._nbpd     = dump.pop('precision')#idem
        for key in ('classename', 'classname') :
            if key in keys and dump.pop(key) == 'NSplineSimple' :
                #C'est une spline simple, c'est un load après construction
                #Les valeurs sont à defaultValue(),
                #on la charge comme spline composée à une seule composante

    #         if not 'ruptures' in keys :
                #Si pas de rupture, c'est une spline simple, on la charge comme spline
                #composée à une seule composante
                self._methode = [self._methode]
                self._mode = [self._mode]
                self._nbpe = [self._nbpe]
                self._nbpd = [self._nbpd]
                break
#         debug(self.methode)

        try : del self._epoints
        except AttributeError : pass
        u"""
        En cas de load() avec plusieurs ruptures (plusieurs brins), mais avec une seule méthode, mode...
         on duplique methode, mode... pour tous les brins: p.ex.
        NSplineComposee(cpoints=...,
                        ruptures=[0,5,-1], # deux brins
                        methode=[('cubic','clamped')], # devrait être une liste a deux elements
                        nbpd=[100], # devrait être une liste a deux elements
                        mode=['linear'], # devrait être une liste a deux elements
                        nbpe=[10])       # devrait être une liste a deux elements"""
        ns = len(self.ruptures)-1#nb splines
        if len(self._methode) == 1 : self._methode = ns*[self._methode[0]]
        if len(self.nbpd) == 1 : self.nbpd = ns*[self.nbpd[0]]
        if len(self.mode) == 1 : self.mode = ns*[self.mode[0]]
        if len(self.nbpe) == 1 : self.nbpe = ns*[self.nbpe[0]]
        if len(self._methode) != len(self.ruptures)-1 :
            msg = ['Il faut preciser une methode pour chaque brin de spline.',
                   'Il doit y avoir nr-1 methodes, nr etant le nombre de ruptures.',
                   'Ici, on a %d methodes et %d ruptures'%(len(self.methode), len(self.ruptures))]
            raise ValueError('\n'.join(msg))
        methodes   = self._methode
#         debug(methodes)
        modes      = self.mode
        nbpes      = self.nbpe
        nbpds      = self.nbpd
        ok = isinstance(methodes,  (list, tuple)) and\
             isinstance(modes,     (list, tuple)) and\
             isinstance(nbpes,     (list, tuple)) and\
             isinstance(nbpds,     (list, tuple)) and\
             isinstance(self.ruptures  ,(list, tuple))
        if not ok :
            raise ValueError('Les parametres ruptures, methodes, modes, nbpes, nbpds doivent être des listes')
        ok = ok and len(methodes) == len(modes) == len(nbpes) == len(nbpds)
        if not ok :
            raise ValueError('Il manque ou il y a trop de données. Nb ruptures=%d,methodes=%d,modes=%d, nbpes=%d, nbpds=%d '%(
                len(self.ruptures),
                len(methodes),len(modes),
                len(nbpes), len(nbpds)))
        self.splines = []
        u"""
        maintenant il y a une methode, mode, ... par brin
        On charge les points de contrôle"""
        cpoints = zeros((1,2))
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
#         debug(ruptures=ruptures, shape=cpoints.shape)
        ok = self._ruptures[0] == 0 and self._ruptures[-1] in (-1,len(cpoints)-1)
        for r in self._ruptures[1:-1] :
            ok = ok and 0 <= r <= len(cpoints)-1
        if not ok :
            msg = [u'Probleme avec les points de rupture %s:'%self._ruptures,
                   u'Le premier devrait etre 0, le dernier devrait etre %d (ou -1)'%(len(cpoints)-1),
                   u'Les autres devraient etre dans l\'ordre croissant, et dans l\'intervalle ]0,%d[ '%(len(cpoints)-1),
                   ]
            raise ValueError('\n'.join(msg))
        self._ruptures[-1] = len(cpoints)-1
#         if len(cpoints) <= 1 :NON, pour un profil vide, il faut 2 splines
#             self.splines = [NSplineSimple(cpoints=zeros((1,2)))]
        if 1: #else :
            for k in range(len(methodes)) :
                deb, fin = self.ruptures[k], self.ruptures[k+1]
#                 debug('construction spline %d'%k, methode=methodes[k], deb=deb, fin=fin)
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
                                                name      = self.name+'-#%d'%k,
                                                role      = 'piece-#%d'%k))
#         for k, s in enumerate(self.splines) :
#             debug('spline %d : \n%s'%(k,s))

    def elaguer(self, eps, replace=False):
        A = []
        for spline in self.splines :
            a = spline.elaguer(eps, replace)
            A.append(a)
        self._update()
        return A

    def echantillonner(self, nbp=None, mode=None):
        if mode is None : mode = self.mode
        if nbp is None : nbp=self.nbpe
        epoints = self.splines[0].epoints
        for n, m, spline in zip(nbp, mode, self.splines[1:]) :
            spline.nbp, spline.mode = n, m
            #on ne met pas le premier epoint sinon il serait en double
            #avec celui de la spline précédente.
            epoints=vstack([epoints,spline.epoints[1:]])
        self._epoints = epoints
        return epoints

    @property
    def epoints(self):
        if hasattr(self, '_epoints') :
            return self._epoints
        else :
            self.echantillonner()
            return self._epoints
#     def echantillon(self, corde):
#         return  self.epoints*corde

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
        return self._nbpd
    @nbpd.setter
    def nbpd(self, prec):
        self._nbpd = prec
        for k, spline in enumerate(self.splines) :
            spline.nbpd = prec[k]
    @property
    def info(self):
        infos1 = [
#             u"nb updates          = %s"%self.nbupdate,
            '%20s = '%u"ruptures"+"%s"%self.ruptures,
            '%20s = '%u"composantes"+"%s"%[s.name for s in self.splines],
            '%20s = '%u"nb pts ctl par comp"+"%s"%[len(s) for s in self.splines]
            ]
        return super(NSplineComposee, self).info + infos1

        # return
#     def debug(self):
# #         return self.__str__()
#         debug(self)
        # for k, spline in enumerate(self.splines) :
        #     print "\nSpline composante numero %d"%k
        #     print 27*"="
        #     print spline

    def toDump(self):
        dump = super(NSplineComposee, self).toDump()
        dump['ruptures'] = self.ruptures
        dump['methode'] = self.methode
        return dump
################################################################
    @property
    def methode(self):
        u"""A l'initialisation, c'est dans self._methode, ensuite, c'est recalculé de manière dynamique,
        a partir des brins de la spline"""
        return [s.methode for s in self.splines]

    def computeSpline(self, methodes):
#         #debug(self.splines)
        for k, spline in enumerate(self.splines) :
#             deb, fin = self.ruptures[k], self.ruptures[k+1]
            spline.computeSpline(methodes[k])

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
                self.splines[ks][k-deb] = value
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
        return i

    def appendPoint(self,pos):
        i = self.splines[-1].appendPoint(pos)
        self._update()
        return self.ruptures[-2] + i

    def localiser(self, k):
        u"""
        k est un numéro de point, entier dans l'intervalle [0, len(self)[
        Calcule et retourne les numéros des splines auxquelles appartient le point numéro k, ainsi que
        la position du point k dans ces splines, sous une des deux formes suivantes
        - [(ks, i),]  ks numero de spline, i = place de k dans cette spline si k appartient à une seule spline ou bien
        - [(ks0, i0), (ks1, i1)] : un point k peut appartenir à (au plus deux) splines qui sont alors consécutives.
            Dans cette hypothèse, k est la fin de l'une (ks0) et le début de la suivante (ks1=1+ks0)
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
            raise NotImplementedError('On ne devrait pas passer ici... Suppression d\'un point anguleux (%d), a verifier'%k)
        else :
#         point = self[k]
            ks, kp = ls[0]#numero de spline, numero de point dans la spline
            return self.splines[ks].removePoint(kp)

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
                  ruptures=ruptures))
        self._update()
        return r#le numero du point

    def split(self, n, mg=None, md=None):
        u"""Découpage de cpoints au point n, avec méthodes mg à gauche et md à droite
        méthodes mg et md sont définies/expliquées dans le fichier splineabstraite.py,
        fonction computeSpline(...)
        Par défaut, on prend les methode mg=md=self.methode """
#         debug(n=n, mg=mg, md=md)
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
                if r>n or r==-1:
                    break
            k -= 1
            mode, nbpe, prec = modes[k-1], nbpes[k-1], nbpds[k-1]
            if mg is None : #la spline de gauche
                mg = methodes[k-1][:]#une copie
                #Si les derivées sont précisées, on transforme en (2,0,0), la condition de droite
                if mg[0] in ('cubic',) and isinstance(mg[1], (list, tuple)):
                    mg[1]=[mg[1][0], (2,0,0)]

            if md is None : #la spline de droite
                md = methodes[k][:]#une copie
                #Si les derivées sont précisées, on transforme en (2,0,0), la condition de gauche
                if mg[0] in ('cubic',) and isinstance(mg[1], (list, tuple)):
                    mg[1]=[(2,0,0), mg[1][1]]

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
            nbpds.insert(k-1, prec)
            modes.insert(k-1, mode)
            nbpes.insert(k-1, nbpe)
    #         self._ruptures = ruptures
    #         self._methode = methodes
            debug(n=n, methodes=methodes)
            self.mode = modes
            self.nbpe = nbpes
#             return methodes, nbpds, modes, nbpes, ruptures
        #pour finir :
            self.load(dict(cpoints=cpoints,
                      methode=methodes,
                      nbpd=nbpds,
                      mode=modes,
                      nbpe=nbpes,
                      ruptures=ruptures))
            for s in self.splines :
                mav = s.methode
                s.methode = s.methode
                if s.methode != mav: debug(avant=mav,apres=s.methode)
    #         debug(self)
            self._update()
#         raise NotImplementedError

    def _update(self):
        u'''Est appelé à chaque modification (géométrique) d'un point de contrôle de la spline
        ou bien du PolygonF de base
        Méthode ultra privée, ne pas l'appeler de l'extérieur.
        '''
        super(NSplineComposee, self)._update()
        try : del self._cpoints#c'est juste un intermediaire à la construction
        except AttributeError : pass
        try : del self._ruptures#c'est juste un intermediaire à la construction
        except AttributeError : pass

    def plot(self, plt, control=True, nbpd=None, nbpe=None, mode=None,
             titre=None, more=[], texts=[], show=True, buttons=True, numbers=[]):
        R = asarray(self[self.ruptures])#des points
        Rx, Ry = XY(R)
#         debug(Rx,Ry)
        more.append((Rx,Ry,'ko','ruptures'))
        super(NSplineComposee, self).plot(plt,control=control, nbpd=nbpd,
                                          nbpe=nbpe, mode=mode, titre=titre,
                                          more=more, texts=texts, show=show,
                                          buttons=buttons, numbers=numbers)
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
            debug(methode=s.methode)
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


if __name__=="__main__":
    from testsplinecomposee import testMain
    testMain()