#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016-2017-2018 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
__updated__="2018-05-21"
from utilitaires import (rstack, eliminerPointsDoublesConsecutifs, diff,pointsFrom)
from splineabstraite import absCurvReal

import sys,os,math
from array import array

import numpy as np
from numpy import log, linspace, asarray, sqrt, arange
from numpy.linalg import  norm
import scipy as sp
# from PyQt4 import QtGui
# from PyQt4 import QtCore
# from PyQt4.QtGui import (QApplication, QPolygonF)
# from PyQt4.QtCore import (Qt, QPointF, QString, QObject)
from config import VALIDATION_DIR, WORK_DIR
from scipy.optimize import newton, minimize
from pprint import pprint
# from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, UnivariateSpline
#from snapper import Snapper
# import preferences
# import cPickle
from utilitaires import (Path, segmentPlusProche, stack, debug, rdebug, dist,
                                     hardScale, absCurv,dist2,rotate,
                                     courbure,symetrieAxe
#                                      baryCentre,centreGravite,
#                                      splineInterpolation, pointsDoubles,
#                                      eliminerPointsDoublesConsecutifs, simpson
                                     )
# from gui.graphicsbase.graphicscommon import (p2s,
# #                                      p2t, pointsFrom,qpolygonFrom,
#                                      qpolygonFrom, pointsFrom, p2t)
# from inout.writerpts import writeProfil
# from preferences import SplinePrefs, NPrefs
from splineabstraite import NSplineAbstract, computeSpline, distance2PointSpline
# import pprint
# from inout.lecteurs import LecteurUniversel
from config import RUNS_DIR
from math import radians
import cPickle
# from numpy import dtype
# from matplotlib import pyplot as plt
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

class NSplineSimple(NSplineAbstract):
    class Default():
        precision = 1000
        _methode  = ('ius',2)
        mode      = 'segment'
        nbpe      = 30        

    def __init__(self, **dump):# cpoints=None, name='spline', role='spline'):
        u"""
        * dump DOIT contenir les clés :
             - aucune
        * dump PEUT contenir les clés suivantes (sinon il y a des valeurs par défaut)
            - 'points' dont la valeur est un np.ndarray, ou un QPolygon
                ou une liste de points, ou un filename ou encore None
            - 'precision' = nb de points pour le tracé de la spline
            - 'methode' : cf computeSpline() paramètre hyper important.
                détermine s'il s'agit d'une spline cubique, de degré autre,
                périodique, si les dérivées aux extrémités sont données, etc...
            - 'mode' : le mode de répartition des points en cas d'échantillonage
            - 'name' : str ou unicode pour affecter un nom
            - 'role' : str ou unicode pour affecter un role...
        """
#         debug(dump)
        #appel de setDefaultValues() puis de load(dump)
        super(NSplineSimple,self).__init__(**dump)
        if self.mode in ('rayon', 'courbure') :
            if self.methode == ('ius',1) :
                raise ValueError(u"Une spline lineaire ne peut pas etre echantillonnee avec le mode '%s'"%self.mode )
        self._update()#non appelé par SplineAbstract

    def load(self, dump):
        u"""Ce qui est commun à toutes les splines."""
        super(NSplineSimple, self).load(dump)
#         msg = [
#             "Attention, le load (le dump ?) ne se fait pas bien, les derivees ne sont pas a jour",
#             "Peut etre la mise a jour des vecteurs derivees aux extremites lors d'une rotation ou symetrie... TODO"
#             ]
#         rdebug('\n'.join(msg))
        u"""
        Attention, la cle est 'points'  dans les anciens projets, et 'cpoints' dans les nouveaux.
        OUI=>On affecte _cpoints directement qui évite l'appel au setter de cpoints (qui appelle _update()).
        cpoints DOIT être liste, tuple, array ou np.ndarray, sans points doubles consécutifs
        """
        for key in ('points', 'cpoints') :#la cle 'cpoints' est prioritaire, en cas
            if dump.has_key(key):
                cpoints = dump.pop(key)
                if isinstance(cpoints, (list, tuple, array)) :#types Python de base
                    #obligé de considerer ce cas, car pour le load d'un pkl cpoints est une liste
                    cpoints = asarray(cpoints)
                if not isinstance(cpoints, np.ndarray) :
                    msg = u'Les points de controle doivent etre de type np.ndarray((*,2)), et non pas %s.'%cpoints.__class__.__name__
                    msg = msg + u' \nUtilisez la fonction points=pointsFrom(points) pour transformer a peut pres n\'importe quoi en np.ndarray((n,2)).'
                    msg = msg + u' \nLes points de controle ne doivent pas avoir de points doubles consecutifs.'
                    msg = msg + u' \nUtilisez points = eliminerPointsDoublesConsecutifs(points) en cas de doute.'
                    raise TypeError(msg)
                if self.methode[1]=='periodic' and len(cpoints) >0 and np.any(cpoints[0] != cpoints[-1]) :
                    p1, p2 = list(cpoints[0]), list(cpoints[-1])
#                     raise ValueError
                    debug(u'Pour une spline periodique, les points de controle \ncpoints[0]=%s \net cpoints[-1]=%s doivent etre identiques. \nIci la distance entre ces deux points est d = %.2g'%(p1, p2, dist(p1,p2)))
                self.cpoints = cpoints

        try :
            del self._epoints
        except AttributeError :
            pass


    # @property
    # def qcpolygon(self):
    #     try :
    #         return self._qcpolygon
    #     except AttributeError :
    #         self._qcpolygon = qpolygonFrom(self.cpoints)
    #         return self._qcpolygon
    # @property
    # def qepolygon(self):
    #     try :
    #         return self._qepolygon
    #     except AttributeError :
    #         self._qepolygon = qpolygonFrom(self.epoints)
    #         return self._qepolygon
#     qpolygon = qepolygon

    @property
    def dpoints(self):
        u"""
        Les points discretises, sont recalculés
            - si _update() a été appelé (_dpoints a été supprimé) ou bien
            - si self.precision a changé
        Si on veut des dpoints aux abscisses T=(t1,...tn), on appelle directement X=self(T) plutôt que X=self.dpoints
        """
        try :
            if len(self._dpoints) == self.precision :
                return self._dpoints
            else :
                return self(linspace(0,1,self.precision))
        except AttributeError :#self._dpoints a été supprimé par un _update()
            if hasattr(self, 'sx') and self.sx is not None :
                T = linspace(0.0, 1.0, self.precision)
                self._dpoints = asarray(self(T))
            else :
                self._dpoints = self.cpoints#Pas de spline, on retourne les points de contrôle
        return self._dpoints

    def abscurv(self, T=None):
        u"""
        Ne fait qu'appeler le abscurv de NSplineAbstract.
        Toute la première partie est un ensemble de vérifications pour debogage
        En production, virer cette verification. (mettre verif=False)
        Les abscisses curvilignes normalisées T des points de contrôle, entre 0 et 1.
        Ce sont exactement les paramètres des splines cubiques :
        - CubicSpline (cubic): T = self.sx.x = self.sy.x
        - UnivariateSpline (us) et InterpolatedUnivariateSpline (ius):
            T = self.sx.get_knots() = self.sy.get_knots()
        :FIXME: pour ('cubic','periodic'), len(self) != len(self.sx.x) => je crois que c'est bon maintenant
            probablement self.cpoints n'est pas fermée self[0] != self[-1]
        :TODO AVEC PRÉCAUTIONS :
            - remplacer les appels à self.abscurv() par self.knots
            - et supprimer la valeur par defaut self.abscurv(T=None) => self.abscurv(T)
            - veiller à ce que self.abscurv retourne la vrai abscisse curviligne (absCurvReal)
        """
        verif = True
        if not verif :
            return super(NSplineSimple, self).abscurv(T)
        ##################################################################
        #  En production, virer cette verification.
        ##################################################################
        if not hasattr(self, 'sx') : #verification, : a-t-on loupé un update ?
            msg = u"A-t-on loupe un _update ?"
            msg += u"\nEn production, virer cette verification if hasattr(self, 'sx'):..."
            rdebug(msg)
            rstack()
            pass
        else :
            W = absCurv(self.cpoints, normalise=True)
            K = self.knots
            if len(K) != len(self) :
#                 rdebug(u'\n'+30*u'!'+u'\n    pour debogage de "%s"'%self.name,u'\n'+30*u'!')
                rdebug(u"BUG : je fais un _update de la spline %s car len(K)=%d != %d=len(self)"%(len(self.knots), self.name, len(self.cpoints)))
                self._update()
#             if self.methode[0] in ('ius', 'us', 'cubic') :
#                 k = self.methode[1]-1
#                 debug(delta=self.knots-W)
#                 ns = np.linalg.norm(W-self.knots, np.inf)
#             elif self.methode[0] == 'cubic' :
#                 if len(self.sx.x) != len(self) :
#                     self._update()
            else :
                ns = np.linalg.norm(W-self.knots, np.inf)
#             else :
#                 ns = 0.0 #tant pis, je sais pas
            if  ns > 1.0e-4:#norme du sup, tolerance 0.1 mm
                rdebug(u'\n'+30*u'!'+u'\n    pour debogage de "%s"'%self.name,u'\n'+30*u'!')
                try : pprint(self.profil.toDump(), sys.stderr)#cas d'une cloison
                except : pprint(self.toDump(), sys.stderr)
                msg = u'BUG : il faut faire un update de la spline car max|self.knots - abscurv|= %g'%ns
                raise ValueError(msg)
        ##################################################################
        #  fin verification.
        ##################################################################
        return super(NSplineSimple, self).abscurv(T)


    def setDefaultValues(self):
        """Valeurs par defaut:"""
        self.precision = self.Default.precision
        self._methode  = self.Default._methode
        self.mode      = self.Default.mode
        self.nbpe      = self.Default.nbpe
        self._cpoints  = np.zeros((1,2))
#         self._qcpolygon = QPolygonF([QPointF(0,0)])
        self.role       = self.classname
        self.name      = self.classname#+'(%d)'%(id(self))
#         self.qpolygon   = np.zeros((1,2))
        #debug("FIN")

#     @property
#     def qcpolygon(self):
#         return self._qcpolygon

#     @qcpolygon.setter
#     def qcpolygon(self, points):
#         u"""On devrait éliminer les points doubles consécutifs,
#         sinon ca plante dans computeSpline()"""
#         qp = qpolygonFrom(points)
#         if qp is not None : self._qcpolygon = qp
#         else : self._qcpolygon = QPolygonF()
#         self._update()

    @property
    def cpoints(self):
        u"""points de contrôle, sous forme np.ndarray((n,2)),
        [re-assemblé à chaque fois => NON]"""
#         try :
        return self._cpoints
#         except AttributeError :
#             self._cpoints = pointsFrom(self.qcpolygon)
#             return self._cpoints
    @cpoints.setter
    def cpoints(self, points):
        u"""points de contrôle, sous forme np.ndarray((n,2))
        - points est un np.ndarray((n,2)) et ne doit pas avoir de points doubles consécutifs."""
#         p = eliminerPointsDoublesConsecutifs(pointsFrom(points))
#         p = points
        if points is None or len(points)<1 :
            self._cpoints = np.zeros((0,2))
        else :
            self._cpoints = points
        self._update()

    @property
    def knots(self):
        if hasattr(self, 'sx') and hasattr(self, 'sy'):#Normalement ac, sx.x, sy.y sont identiques. Ce sont les knots
#             debug(self.sx.__dict__)
            try :
                return self.sx._data[0]
                return self.sx.get_coeffs()#si methode=('ius',k)

            except AttributeError :
                return self.sx.x#si methode = (cubic','xx')
        else :
            return None

#             if len(ac) != len(sxx) or norm(ac-sxx) != 0.0 :
#                 self._update()
#                 return self.knots#recursif...
#             else :
#                 return sxx
#         ac = self.abscurv()
#         if len(ac) <= 1 : #sx n'esiste pas
#             return ac
#         else :
#             return ac

    def echantillonner(self, nbp=0, mode='linear', ta=0, tb=1, onlyT=False):
#     def echantillonner(self, nbp=0, mode='segment', ta=0, tb=1, onlyT=False):
        u"""
        répartition (échantillonnage) de nbp points sur la spline self,
        entre les abscisses ta et tb, suivant le mode précisé par 'mode'
        *modifie self.epoints*
        :return

            - si onlyT=False : le tableau E=np.ndarray(shape=(n,2)), des points échantillonnés
            - si onlyT=True : les abscisses T=np.ndarray(shape=(n,1)) (je crois!) des points échantillonnés

        :param mode :

            - si mode = ``'segment'`` decoupage de chaque segment P[i], P[i+1] en autant de points
                que nécessaires pour en avoir nbp au final.
                :TODO : ce mode est *A SUPPRIMER car LES POINTS ECHANTILLONNES NE SONT PAS SUR LA SPLINE*,
                mais sur les segments reliant les points de contrôle.
                et les points de contrôle *SONT* des points d'echantillon.
                Si on a besoin de ce mode, la spline doit être un Polyligne
            - si mode = 'linear', les points d'échantillonnage sont régulièrement répartis
                tout au long de la spline
            - si mode = 'rayon' ou 'courbure' :
                la densité de points est proportionnelle à la courbure. (Approximativement)
            - si mode = 'cpoints' :
                retourne simplement les points de contrôle.
            - si mode = 'telkel' :
                ta doit être un tableau des abs curv. des points echantillonnés.
                retourne self(ta)

        :type mode: str ou unicode

        :param nbp : nombre de points d'échantillonnage.

            - Si mode='telkel' ou 'cpoints', ce paramètre est inutile.
            - Dans tous les autres cas il est indispensable, nbp>0
                (sinon raise ValueError)

        :type nbp : int

        :param ta : facultatif

            - l'abs. curv. du premier point d'échantillonnage
            - si mode = 'telkel' ou 'cpoints', ta est le tableau des abs. curv. des points échantillonnés.

        :type ta : float dans [0,1] ou np.ndarray((n,1)) comme retourné par self.absCurv()
        :param tb :

            - l'abs. curv. du dernier point d'échantillonnage
            - si mode='telkel' ou 'cpoints' ce paramètre est inutilisé.

        :type tb : float dans [0,1]
        """
        if nbp==0 : nbp = self.nbpe
        else : self.nbpe = nbp
        self.mode = mode
        if nbp==0 and mode not in ('telkel','cpoints'):
            msg = u'Echantillonnage : préciser le nb de points et le mode'
            raise ValueError, msg
#         debug(echantillonner=self.echantillonner)

#         debug(t=(ta,tb), mode=mode, methode=self.methode)
#         stack()
    #         #debug(self)
        if mode == 'telkel' :
#             debug(ta=ta)
            T = ta#en fait c'est un tableau d'abscisses
            self._epoints = asarray(self(ta))
            return self._epoints

        if mode == 'cpoints' :
            #On retourne les cpoints, simplement.
            #Utile pour conserver les shapes des anciens pj sans modification
            #Utile pour les Polylines et Polygones
            T = self.abscurv()
            self._epoints = np.copy(self.cpoints)
            if onlyT :
                return T
            else :
                return self._epoints

        if mode == 'segment' :
            msg = u"""Le mode d'échantillonnage '%s' est obsolète, car les points échantillonnés
    sont sur les segments joignant les points de contrôle, ils ne sont donc pas sur la spline.
    Pour que les points échantillonnés soient sur les segments joignant les points de contrôle,
    il suffit d'utiliser :
        - un NPolygone ou
        - un NPolyligne ou
        - une NSplineSimple(..., methode=('ius',1))
    Fixer ensuite le mode d'échantillonnage à 'linear', 'cpoints', 'telkel' ou 'rayon')"""%mode
            rdebug(msg)
            raise ValueError(u"Mode d'echantillonnage obsolete")
            #création d'une spline de degré 1 (linéaire)
            #on n'a pas le droit de changer la méthode d'une spline.
            #on en reconstruit une, lineaire.
            #TODO:A supprimer de NSplineSimple, utiliser Polygone ou PolyLine
            cpoints = self.cpoints
            if self.methode == ('cubic', 'periodic') :
                #On rajoute un point a la fin, sinon il echantillonne pas tout
                pass
    #                 cpoints = np.vstack((cpoints, cpoints[0]))#deja fait à l'init.
            slin = NSplineSimple(cpoints=cpoints, methode=('ius', 1), precision=0, name='Pour echantillonnage par segment')
            #On prend nbp-len(self) points linéairement espacés
            # auxquels on rajoute les points de contrôle
            T1 = slin.abscurv()
            if nbp == 0 : nbp = len(self.cpoints)
            if len(T1)>=nbp : #Pas assez de points (nbp) pour mettre les points de contrôle (T1).
                T = linspace(ta,tb,nbp)
    #                 debug(nbp=nbp, T1_shape=T1.shape)#, linshape=linspace(0,1,nbp-len(T1)).shape)
            else :
                T = np.hstack((linspace(ta, tb,nbp-len(T1)), T1))#/T1[-1]))
                np.sort(T)
    #             #debug(T=T)
    #             #debug(asarray(slin(T)))
            self._epoints = asarray(slin(T))
            if onlyT : return T
            else : return self._epoints
        elif mode in ('linear', 'lineaire', 'lin') :
            T = linspace(ta, tb, nbp)#y compris ta et tb
            self._epoints = asarray(self(T))
            if onlyT : return T
            else : return self._epoints
        elif mode in ('rayon', 'courbure') :
            u"""On ne touche plus à RIEN !!!
            Calcul des points d'échantillonnage, de sorte que la densité de points
            soit localement proportionnelle à la courbure
            - la "COURBURE ABSOLUE" en un point t est ca(t) = sqrt(abs(self.courbure(t)))
            On la calcule en N points d'une discrétisation fine T, elle se trouve dans CA
            - La SINUOSITÉ sur un micro-intervalle [t, t+dt] est l'intégrale de t à t+dt de ca(s)ds
                On la calcule en chaque micro-intervalle de la discrétisation fine, elle est dans S
                La SINUOSITÉ TOTALE est la sinuosité sur [a,b], c'est la somme des sinuosités des micro-intervalles
            - On partage l'intervalle [ta,tb] en nbp-1 intervalles (de largeur variable)
                ta=te[0]<te[1]<te[2]<...<te[nbp-1]=tb;
            tels que sur chaque intervalle j, la sinuosité sj soit constante=s0=SINUOSITÉ TOTALE/(nbp-1)
            ou, ce qui revient au même, tq la sinuosité de ta à te[j] soit egale à Sj=j*s0
            """
#             assert (nbp!=1)
            if nbp==1 :
                debug(nbp=nbp, ta=ta,tb=tb)
                stack()
                if onlyT : return ta
                else : return self(ta)
            N = 100001
            T = linspace(ta,tb,N)#N-1 micro intervalles
            dt = T[1]-T[0]#micro intervalles de largeur dt
            CA = np.sqrt(np.abs(self.courbure(T)))
            S = (CA[1:]+CA[:-1])*(dt/2)#Les sinuosités des N micro-intervalles = dt*moyenne des courbures aux extrémités
            s0 = sum(S)/(nbp-1)#La sinuosité (constante) de chaque intervale Te[j],Te[j+1] défini ci dessous
            Te = np.ones(nbp)*ta# les te[j] cherchés
            Te[-1] = tb
            i, Sj = 0, 0.0#i=numéro du micro-intervalle
#             dbg = False
            for j in range(1, nbp) :#j=les numeros des points d'échantillonnage
                while Sj < j*s0 :#Sj=Sinuosité de ta à t[j]
                    try :
                        Sj += S[i]
                    except IndexError as msg :
                        msg = ['IndexError :', '[ta, tb]=[%.4g, %.4g]'%(ta,tb),
                               'nbp=%d, j=%d'%(nbp,j),
                               's0 = %.4g ; j*s0-Sj=%.4g > 0'%(s0,j*s0-Sj),
                               ]
                        if j*s0-Sj < s0/1000 :
                            deb, msg = strace, [u'\nPas grave %s'%self.name]
                        else :
                            deb, msg = rdebug, [u'\nAttention']
                            msg += [
                                   'i=%d'%i,
                                   'Te=%s'%(str(Te))]
                        deb('\n    '.join(msg))
                        break
                    i += 1#i=
                Te[j] = T[i]
#             if dbg :
#                 rdebug()
            self._epoints = self(Te)
            if onlyT : return Te
            else : return self._epoints
        else :
            rdebug('mode echantillonnage inconnu : ',mode=mode,nbp=nbp, ta=ta, tb=tb, onlyT=onlyT)

    def scaled(self, echelle):
        u'''
        :return: une COPIE de self, mise à l'échelle.
        :rtype: NSplineSimple.
        :TODO: retourner le type de self plutot qu'une NSplineSimple
        Pas de modif de self.
        :attention: Ne pas faire simplement un scale de self.cpoints, il faut
        appeler spl.hardScale(echelle) qui fait les mises à jour des dérivées le cas échéant.
        :param echelle: (echx, echy) l'échelle en x et y :
        :type echelle: (float, float)
        '''
#         points=self.cpoints
#         hardScale(points,echelle)
        echx,echy=echelle[0], echelle[1]
        name='scaled(%.2g,%.2g)x'%(echx,echy)
        try : name=name+self.name
        except TypeError : name=name+'-scaled'
        dump = self.toDump()
#         dump['cpoints'] = points
        dump['name'] = name
        dump['parent'] = self.parent()
        spl = NSplineSimple(**dump)
        spl.hardScale(echelle)
        return spl

    @property
    def methode(self):
        return self._methode
    @methode.setter
    def methode(self, newmethode):
        u"""Interdiction de changer de méthode (cubic, us, ius),
        On peut changer les paramètres de la méthode (clamped, not-a-knot,...)
        Pour changer de méthode, quand même, on peut faire :
        >>> S = SplineSimple(points=..., methode=meth1)
        >>> dump = S.toDump()
        >>> dump['methode'] = newmethode
        >>> S.load(dump)"""
        if newmethode[0] == self.methode[0] :
            self._methode = newmethode
            self._update()
        else :
            dump = self.copy().toDump()
            dump['methode'] = newmethode
            self.load(dump)
#             raise RuntimeError("Interdiction de changer de méthode")

    def computeSpline(self, methode=None):
        u"""
        Calcule la spline (sx, sy) considérée comme une courbe paramétrée sx(t), sy(t).
        sx et sy sont deux splines à une seule variable au sens scipy.
        """
#         debug(name=self.name, len=len(self))
        if methode is None :
            methode=self.methode
        if methode[1] != self.methode[1] :
            debug(u'Attention, changement de parametres dans \'%s\': self.methode= %s => methode=%s'%(self.name,str(self.methode),str(methode)))
            self.methode = methode
        try :
#             debug('Avant computeSpline', self)
#             pprint(self.toDump())
            self._knots, self.sx, self.sy = computeSpline(self.cpoints, methode)
#             rdebug('computeSpline:OK?')
#             debug('Apres computeSpline', self.name)
            try : #On regarde si sx et sy existent. Sinon, on remonte l'exception
                sx_05=self.sx(0.5)
                sy_05=self.sy(0.5)
#                 rdebug(sx_05=sx_05, sy_05=sy_05)
                if np.isnan(sx_05) or np.isnan(sy_05) :# in (np.nan,np.inf) or sy_05 in (np.nan,np.inf) :
                    raise ValueError('sx(0.5) ou sy(0.5) : not a number')
            except AttributeError as msg :
                rdebug(msg)
                raise msg
            except TypeError as msg :
                #sx=None, n'existe pas encore
                pass
        except ValueError as msg :
            # le message est `x` must be strictly increasing sequence.
            # il y a des points doubles certainement, on les vire
#             rdebug(self._knots)#existe pas
#             raise
#             rdebug(msg, u'"%s", Je tente d\'eliminer les points doubles consecutifs'%self.name)
#             debug(points_sales=self.cpoints.tolist())
            cpoints, avirer = eliminerPointsDoublesConsecutifs(self._cpoints, vires=True)
#             rdebug(avirer = avirer)
            if not avirer and len(cpoints)!=0:
                debug(BUG=u"Un update doit manquer..., len(cpoints)=%d"%len(cpoints))
#             debug(abscurv=absCurv(cpoints, False))
            if len(cpoints) < len(self._cpoints) :#sinon ca boucle
                try : role=self.role
                except : role='pas de role'
                rdebug(avirer=avirer, classename=self.classname, name=self.name, role=role)
                self._cpoints = cpoints
                self.computeSpline(methode)
            # Ca fait le _update et ca rappelle self.computeSpline()
#             ac = self.abscurv()
#             dac = ac[1:]-ac[:-1]
#             avirer = []
#             for k, d in enumerate(dac) :
#                 if abs(d)<=1.0e-10 :
#                     avirer.append(k)
# #                     self.removePoint(k)
#             if avirer :
#                 msg1 = 'Les points de controle %s sont des doublons consecutifs.'%avirer
#                 msg1 += 'Suppression des doublons'
#                 rdebug(msg1,'')
# #                 rdebug(avirer=avirer, delta_abscurv=dac.tolist(), controle=sum(dac))
#                 for k in reversed(avirer) :
#                     self.removePoint(k, update=False)
# #                 rdebug('apres nettoyage : ', cpoints = self.cpoints)
#                 self.computeSpline(methode)
# #                 #le cpoint.setter appele _update() qui rappelle computeSpline()
# #                 self.cpoints = eliminerPointsDoublesConsecutifs(self.cpoints)

        except RuntimeWarning as msg :
#             rdebug(msg)
            raise(msg)
        except Exception as msg :
#             rdebug(type(msg).__name__, msg)
            raise(msg)

    def appendPoint(self,pos):
        i = self.cpoints.shape[0]#l'indice du dernier point
        if len(self) and dist2(pos, self.cpoints[-1]) == 0 :
            raise ValueError('Point double, Impossible d\'ajouter ce point:%s'%(pos))
            return #eviter points doubles, les abs curv douvent être strictement croissantes.
        else :
            self.cpoints = np.vstack((self.cpoints, pos))
            self._update()
            return i
    def index(self, t):
        """Localise t dans self.abscurv:retourne les deux entiers
        k1,k2  tels que k1 < t <= k2"""
        if t<0 or t>1 : return np.nan, np.nan
        T = self.abscurv()
        k = 0
        while T[k]<t : k+=1
        if T[k] == t : return k,k
        else : return k-1, k

    def insertPoint(self, pos, i=None):
        u"""
        segment_i = [self[i-1], self[i]]
        Si i = None :
            Calcule dans quel segment il est raisonable d'insérer le point pos
        si i != None,
            insertion du point pos entre i-1 et i (de sorte qu'il deviennet le point i)
        """
        if isinstance(pos, QPointF) :
            pos = pos.x(), pos.y()
        cpoints = self.cpoints
        if i is None :
            im1, _ = segmentPlusProche(cpoints, pos)
            i = im1 + 1
            #Les segments sont numerotes de 0 à len(cpoints)-2
            #le segment i est délimité par les points i-1 et i
#             debug('Calcul segment le plus proche de %s'%pos)
#             if im1 is None :#Ca n'arrive jamais, cf segmentPlusProche()
#                 debug('Je ne sais pas ou (dans quel segment) inserer ce point %s'%str(pos))
#                 return
        elif i == 0 :
            #si i est le premier point, on insere dans le premier segment
            im1, i = 0, 1
        else :
            im1 = i-1
        if dist2(pos, cpoints[i]) == 0\
        or dist2(pos, cpoints[im1]) == 0 :
#             debug('segment le plus proche',i=i,H=H,pos=pos)
            raise ValueError('Point double, impossible d\'inserer ce point:%s en position %d'%(str(pos),i))
            return #eviter points doubles, car les abs curv doivent être strictement croissantes.
        else:
            cpoints = np.insert(cpoints, i, (pos,), axis=0)
#             debug(cpoints=cpoints)
            self.cpoints = cpoints
#             self._update()#deja fait par le cpoint.setter
            return i#La position du point inséré

    def removePoint(self, pnum, update=True):
        u"""
        - Suppression du point pnum de qcpolygon
        """
        point = self[pnum]
        if update : #le setter de cpoints fait un update
            self.cpoints = np.delete(self.cpoints, pnum, 0)
        else :
            #le setter de cpoints fait un update, parfois c'est indésirable
            #dans ce cas, on intervient directement sur _cpoints
            self._cpoints = np.delete(self._cpoints, pnum, 0)
#         self._update()#c'est fait par le setter de cpoints
        return point

    def hardRotate(self, alfa, centre=None, unit='degres'):
        u'''
        Rotation de self, alfa est en degres par défaut
        si les dérivées aux extremités sont précisées, il faut les tourner elles aussi
        methode est de la forme : ('cubic',(D0, Dn))
        avec D0 = (d, a, b) ou (a,b) est le vecteur dérivée d-ieme au point 0
        et Dn idem pour le dernier point.
        '''
#         debug('avant:', alfa=alfa)
        if unit == 'degres' :
            alfa = math.radians(alfa)
#         debug('apres:',alfa=alfa)
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On fait tourner les dérivées premiere et seconde du point 0, puis du dernier point.
            #Si les data sont x' et y'' par exemple, il faut faire tourner x',y' puis x'', y''
            newderivees = []
            for k in range(2) : #on tourne les deux vecteurs dérivés aux extrémités
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être tournee de alfa
                d = derivees[0]#ordre de derivation
                u = derivees[1:], #u=(x', y') ou bien (x'', y'') c'est lui qui tourne
                #debug('Avant rotate:x',d=d, u=u)
                u = rotate(u, alfa, (0,0))[0]
                #idem pour la derivée k-ième
#                 #on ne garde que u[0] et v[1]
                newderivees.append((d, u[0], u[1]))

            self._methode = [self.methode[0], newderivees]

#         debug('Avant rotate', cpoints=self.cpoints, barycentre=self.barycentre, gravity=self.gravitycenter)
        if centre is None :
            centre = self.barycentre
        points = rotate(self.cpoints, alfa, centre)
#         debug('apres rotate', cpoints=self.cpoints, barycentre=self.barycentre, gravity=self.gravitycenter)
        self.cpoints = points
#         self._update()c'est fait dans le setter cpoints()

    def symetriser(self, axe, centre=None):
        u'''
        Modif de la spline elle meme. self est modifié.
        C'est à dire
        >>> self[i,axe] = 2*centre -self[i,axe]
        si centre == None on prend l'isobarycentre[axe]
        '''
        if len(self)==0 : return
        if centre is None :
            centre = self.barycentre[0][axe]
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On fait tourner les dérivées du point p0, puis du dernier point.
            #(dérivées premiere ou seconde suivant methode)
            newderivees = []
            for k in range(2) :
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être symétrisée
                d = derivees[0]#ordre de derivation
                u = asarray(derivees[1:], dtype=float).reshape((1,2))
#                 debug('Avant sym:',u=u, axe=axe, centre=centre)
                u = symetrieAxe(u,axe,0.0)[0]
#                 debug('Apres sym:',u=u, axe=axe, centre=centre)
                newderivees.append((d, u[0], u[1]))

            self._methode = (self.methode[0], tuple(newderivees))
        self.cpoints = symetrieAxe(self.cpoints,axe,centre)
#         self._update()c'est fait par le setter de cpoints
#         self.qcpolygon=points#qpolygonFrom(points)
#         self.update()C'est fait dans le setter qcpolygon()

    def hardScale(self, scale=(1.0,1.0), centre=None):
        u'''
        Modif de la spline elle même. self est modifié.
        Mise à l'échelle d'un scale=(fx,fy), centrée sur centre. (en fait, une homothétie)
        C'est à dire

            - self[i][0] = centre[0] + scale[0]*(self[i]-centre)[0]
            - self[i][1] = centre[1] + scale[1]*(self[i]-centre)[1]

        Le cas échéant, les valeurs des dérivées aux extrémités de la spline sont
        mises à l'échelle, de sorte que self transformée soit réellement homothétique de self.
        En modifiant self.cpoints, le _update est appelé, tous les attributs de self
        (_epoints, _absCurv, _dpoints, etc...) sont mis à jour.
        :param centre: (float, float) les coordonées (x0, y0) du centre d'homothétie.
            si centre==None on prend l'isobarycentre des points de contrôle.
        :param scale: (float, float)
        :return: None

        '''
        if len(self) == 0 : return
        if centre is None :
            centre = self.barycentre
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On scale les dérivées du point p0, puis du dernier point.
            #(dérivées premiere ou seconde suivant methode)
            newderivees = []
            for k in range(2) :
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être mise a l'echelle
                d = derivees[0]#ordre de derivation
                #u = (x'(t), y'(t)) ou (x"(t),y"(t)) en t=0 ou 1
                u = asarray(derivees[1:], dtype=float)
                u[0] *= scale[0]
                u[1] *= scale[1]
                newderivees.append((d, u[0], u[1]))

            self._methode = (self.methode[0], tuple(newderivees))
#         debug(Avant=self.cpoints)
        self.cpoints = hardScale(self.cpoints, scale, centre)#ca fait le _update
#         debug(APRES=self.cpoints)
        return

#     def copy(self):
#         """retourne une copie de self"""
#         dump = self.toDump()
#         return NSplineSimple(**dump)

    def translate(self,vecteur):
        self.cpoints = self.cpoints+vecteur
#         self._update()#C'est fait dans la property cpoints
    hardMove = translate

    def bidouille(self, alfa):
        """Je recupère les coefficients des splines, et je leur applique la rotation"""
        cx = self.sx.c
        cy = self.sy.c
        for k, c in enumerate(cx) :
            debug('==> cx[%d] = '%k, c)
#         debug (cx=cx)
        ca, sa = math.cos(alfa), math.sin(alfa)
        R = np.matrix([[ca, -sa],[sa, ca]])
        Rcx, Rcy = np.zeros(cx.shape), np.zeros(cy.shape)
        for k in range(4) :
            Rcx[k], Rcy[k] = R*[cx[k], cy[k]]
        self.sx.c = Rcx
        self.sy.c = Rcy
        del self._dpoints, self._epoints#, self._cpoints
#     @property
#     def longueur(self):
#         u"""longueur = diagonale du rectangle d'encombrement????"""
#         return absCurv(self.dpoints, normalise=False)[-1]
# #         return self.abscurv()
#         return np.linalg.norm([self.height, self.width])

    def pourMille(self, longueur):
        u"""longueur convertie en ‰ de la longueur de self"""
        return 1000*longueur/self.longueur
#         a creuser
#         self.plot()
    def plotCourbure(self):
        from matplotlib import pyplot as plt
        from matplotlib.widgets import CheckButtons
#         nbpd = self.precision
        nbpe = self.nbpe
        self.echantillonner(nbpe)
        D = self.dpoints
        C = self.cpoints
        E = self.epoints
        _, ax = plt.subplots()
        titre = self.name+' courbure'
        plt.title(titre)
        T = linspace(0,1,100)
        spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
#         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
        control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        courbure = self.courbure(T)
        courbure += (1.0 + abs(min(courbure)))
        courbure = log(courbure)
        courbure /= max(abs(courbure))
        courbure, = ax.plot(T, courbure)
        buttons = ['points controle','courbure','spline',]
        values = [True, True, True, True]
        draws = [control, courbure, spline]
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        rax = plt.axes([0.05, 0.4, 0.1, 0.15])
        check = CheckButtons(rax, buttons, values)

        def func(label):
            if label == 'spline':
                spline.set_visible(not spline.get_visible())
            elif label == 'points controle':
                control.set_visible(not control.get_visible())
            elif label == 'courbure':
                courbure.set_visible(not courbure.get_visible())
            else :
                draw = draws[buttons.index(label)]
                draw.set_visible(not draw.get_visible())
            plt.draw()
        check.on_clicked(func)
        plt.show()
        return plt

    def ajuster(self, eps=0.5, replace=False):
        u"""
        On cherche une spline s1 avec un minimum de points de contrôle
            et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)
        On utilise une spline d'AJUSTEMENT scipy
        retourne la spline, s1
        En général, self.elaguer fonctionne mieux que self.ajuster
        """
        def refine(n=3):
            res = minimize(f1, np.ravel(s1.cpoints), (s0, Td,), 'BFGS', options={'maxiter':n})#, 'disp':True})
            c1 = res.x
            c1.shape = (-1,2)
            #recalage des points 0 et -1
            c1[0] = c0[0]
            c1[-1] = c0[-1]
            s1.cpoints = c1
        s0 = self
        Td = linspace(0,1,1000)
        T0 = self.abscurv()
        c0 = s0.cpoints.copy()
        d0 = asarray(s0(Td))

#         A, B = c0[0], c0[-1]
#         c1 = pointsFrom((A,B))
        w = np.ones(len(c0))
        w = abs(s0.courbure(T0))
        w[0]*=10000
        w[-1]*=10000
        parametres = {'w':w, 'k':3, 's':1.0e-6, 'ext':2, 'check_finite':True}
        T, sx, sy = computeSpline(c0, methode=('us', parametres))
        Tx = sx.get_knots()
        Ty = sy.get_knots()
        T1 = sorted(set(np.hstack((Tx,Ty))))
        c1 = np.asarray(zip(sx(T1), sy(T1)))
        c1[0] = c0[0]
        c1[-1] = c1[-1]
#         debug(c1.shape)
        s1 = NSplineSimple(cpoints=c1,
                           methode=s0.methode,
                           precision=s0.precision,
                           name='adjust(%s)'%s0.name,
                           mode=s0.mode)
        d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
        debug('AJUSTEMENT avant optimisation : dist = %.2g mm/m '%(self.pourMille(d)))
        for k in range(1) : refine()
        d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
        #courbure totale
        courburetotale0 = self.integraleCourbure(a=0.01, b=0.99)
        courburetotale1 = s1.integraleCourbure(a=0.01, b=0.99)
        debug('Apres AJUSTEMENT : dist = %.2g mm/m ; courbure totale :c1/c0 %f '%(self.pourMille(d),courburetotale1/courburetotale0))

#         debug('Apres ajustement : courbure totale :c1/c0 %f '%(courburetotale1/courburetotale0))
        n0 = len(self)
        n1 = len(s1)
        if n0 == n1 :
            s1.cpoints = self.cpoints.copy()
        elif replace :
            self.cpoints = s1.cpoints
        return s1, self.pourMille(d),(n0,n1)

    def elaguer(self, eps=0.5, replace=False):
        u"""
        On cherche une spline s1 avec un minimum de points de contrôle
        et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)

        La distance(self.cpoints, s1) est le plus grand écart entre
        la spline calculée et les points de contrôle de la spline self.
        autrement dit le max des distances d(self.cpoints[k], s1), k=0,1,...
        où d(P, s1) est le min de la fonction t -> norme(P-s1(t)).
        Voir la fonction distance2PointSpline().
        On discrétise finement self (une seule fois) -> tableau de points D0
        On initialise s1 partant des quelques points de contrôle dont
        self.cpoints[0] et self.cpoints[-1].
        On discrétise finement s1 (à chaque ajout de point) -> tableau de points D1
        puis on rajoute des points de contrôle à s1 jusqu'a obtention de la précision désirée (eps)
        Pour rajouter un point, on repère la distance point à point de D0 et D1,
        disons qu'elle est au i-eme point et on insère dans s1 un point de contrôle
        que l'on positionne exactement en D0[i].

        [ne marche pas =>] Quand la précision désirée est atteinte, on met
        un petit coup d'optimisation pour améliorer la position des points
        de contrôle de s1.

        :param eps: float, la spline resultante est à une distance eps au maximum de self.
        :param replace: bool, si True, la spline élaguée remplace self.
        :return: s1, pm, (n0,n1)

            - s1 est la spline NSplineSimple élaguée
            - pm = self.pourMille(d) est la précisison en valeur relative :
                la spline s1 ne s'éloigne pas à plus de pm ‰ de self.
            - n0, n1 = nb de points de contrôle avant et après élagage

        """
#         eps = self.pourMille(eps)
#         debug(eps=eps)
        if len(self) < 10 :
            rdebug(u"Élagage inutile, %d points de contrôle seulement"%len(self))
            return self, self.pourMille(0),(len(self),len(self))
#         from matplotlib import pyplot as plt
        s0 = self
        nd = 5000
        Td = linspace(0,1,nd)
        T0 = self.abscurv()
        c0 = s0.cpoints.copy()
        d0 = asarray(s0(Td))
        t, m = self.methode
        if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
            init = (c0[0], self(0.33), self(0.66), c0[0])
            init = (c0[0], self(0.25), self(0.5), self(0.75), c0[0])
#             A, B, C, D = c0[0], c0[1], c0[-2], c0[0]
            c1 = pointsFrom(init)
        elif dist2(c0[0], c0[-1]) < 1.0e-6 :
            c1 = pointsFrom((c0[0], self(0.5), c0[0]))
        else:
            c1 = pointsFrom((c0[0], c0[-1]))
        s1 = NSplineSimple(cpoints=c1,
                           methode=s0.methode,
                           precision=s0.precision,
                           name='%s-elaguee'%s0.name,
                           mode=s0.mode)
        d1 = asarray(s1(Td))
#         rdebug(s1.longueur)
        for k in range(len(c0)) :
            ld2 = asarray([distance2PointSpline(c0[i], s1, t0=t).fun for i, t in enumerate(T0)])
            ld2[0] = ld2[-1] = 0.0
            imaxd = (ld2 == max(ld2)).nonzero()#indices du-des point-s ou la distance entre les deux splines est max
            idx = imaxd[0][0]
#             debug(imaxd=imaxd)
            d = np.sqrt(max(ld2))
            d = self.pourMille(d)
            pos0 = c0[idx]#position(x,y)
#             debug(u'idx=%d, t=%.4f, d=%.2g ‰'%(idx, t, d))#, self.pourMille(d)))
#             print(u'd%d = np.asarray(%s)'%(k, str(np.sqrt(ld2).tolist())))
#             print(u'd0%d = np.asarray(%s)'%(k, str(d0.tolist())))
#             print(u'd1%d = np.asarray(%s)'%(k, str(d1.tolist())))
#             print(u'insert_pos=%s'%pos0)
            try :
                res = s1.insertPoint(pos0)
#                 debug(res_insert=res)
#                 s1.plot(plt)
            except ValueError as msg :#Impossible d'inserer le point (point double ?)
                debug(u'je tente autre chose', msg, pos0=pos0)
#                 debug(s0=s0.cpoints)
#                 debug(s1=s1.cpoints)
#                 rdebug()
#                 nd += 10
                Td = linspace(0,1,nd)
                d0 = asarray(s0(Td))
                d1 = asarray(s1(Td))
                ad = norm(d0-d1, 2, axis=1)#[1:-1]#ecart point à point des deux splines discrétisées
                mad = (ad == max(ad)).nonzero()#indice des points ou la distance entre les deux splines est max
    #             debug(mad=len(mad), dist=ad[mad])
                idx = mad[0]
                t = Td[idx][0]
                pos0 = d0[idx][0]
                try :
                    s1.insertPoint(pos0)
                except ValueError as msg :
                    rdebug(msg, pos0=pos0)
                    rdebug(u'Precision non atteinte, iteration %d : dist = %.2e '%(k,d))
                    break
#                     c1 = s1.cpoints.copy()
            d = [distance2PointSpline(c0[i], s1, t0=t).fun for i, t in enumerate(T0)]
            d = np.sqrt(max(d))
            d = self.pourMille(d)
            debug(u'dist-%d = %.2g‰ '%(k,d))
            d1 = asarray(s1(Td))
            c1 = s1.cpoints.copy()
            if d<eps : break
        if len(s1) == len(self) :#même nb de points : on garde la spline initiale.
            s1.cpoints = self.cpoints.copy()
            n0 = len(self)
            n1 = len(s1)
            return s1, d,(n0,n1)
        #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de s1
        #ca ne marche pas du tout !!!
        def refine(n=3):
            res = minimize(f1, np.ravel(s1.cpoints), (s0, Td,), 'BFGS', options={'maxiter':n})#, 'disp':True})
    #         print res.fun
            c1 = res.x
            c1.shape = (-1,2)
            #recalage des points 0 et -1
            c1[0] = c0[0]
            c1[-1] = c0[-1]
            s1.cpoints = c1
#         debug('Avant optimisation : dist = %.2g mm/m ... wait...'%(d))
#         for k in range(3) : refine()
#         d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
#         d = self.pourMille(d)
#         debug('Apres optimisation : dist = %.2g mm/m '%(d))
        #courbure totale
#         courburetotale0 = self.integraleCourbure(a=0.01, b=0.99)
#         courburetotale1 = s1.integraleCourbure(a=0.01, b=0.99)
#         ratio = courburetotale1/courburetotale0
        n0 = len(self)
        n1 = len(s1)
        debug('Apres ELAGAGE : dist = %.2g mm/m ; %d => %d '%(d,n0,n1))
#         debug('Apres ELAGAGE, courbure totale c1/c0 = %f '%(courburetotale1/courburetotale0))

        if replace :
            self.cpoints = s1.cpoints
        return s1, d,(n0,n1)

################################################################################

    def elaguerOld(self, eps=0.5, replace=False):
        u"""
        On cherche une spline s1 avec un minimum de points de contrôle
            et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)
        La distance(self.cpoints, s1) est le max des distances d(self.cpoints[k], s1), k=0,1,...
        où d(P, s1) est le min de la fonction t -> norme(P-s1(t)). Voir la fonction distance2PointSpline.
        On discrétise finement self (une seule fois) -> tableau de points D0
        On initialise s1 partant des deux points de contrôle self.cpoints[0] et self.cpoints[-1]
        On discrétise finement s1 (à chaque ajout de point) -> tableau de points D1
        puis on rajoute des points de contrôle à s1 jusqu'a obtention de la précision désirée (eps)
        Pour rajouter un point, on repère la distance point à point de D0 et D1, disons qu'elle est au i-eme point
        et on insère dans s1 un point de contrôle que l'on positionne exactement en D0[i].
        Quand la précision désirée est atteinte, on met un petit coup d'optimisation pour améliorer
        la position des points de contrôle de s1.
        retourne la spline, s1
        """
#         eps = self.pourMille(eps)
#         debug(eps=eps)
        s0 = self
        Td = linspace(0,1,1000)
        T0 = self.abscurv()
        c0 = s0.cpoints.copy()
        d0 = asarray(s0(Td))
        t, m = self.methode
        if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
            A, B, C, D = c0[0], self(0.33), self(0.66), c0[0]
            A, B, C, D, E = c0[0], self(0.25), self(0.5), self(0.75), c0[0]
#             A, B, C, D = c0[0], c0[1], c0[-2], c0[0]
            c1 = pointsFrom((A,B,C,D,E))
        else :
            A, B = c0[0], c0[-1]
            c1 = pointsFrom((A,B))
        s1 = NSplineSimple(cpoints=c1,
                           methode=s0.methode,
                           precision=s0.precision,
                           name='%s-elaguee'%s0.name,
                           mode=s0.mode)
        d1 = asarray(s1(Td))
    #     debug(s1)
        for k in range(len(c0)) :
    #         c1 = s1.cpoints#.copy()
            u"""
            discrétisation fine de s0 et s1 (=> d0 et d1), et calcul de l'écart point à point
            entre ces deux splines discrétisées
            Là ou l'écart est maximum, on rajoute un point de contrôle P à s1,
            et on positionne P sur s0"""
            ad = norm(d0-d1, 2, axis=1)#[1:-1]#ecart point à point des deux splines discrétisées
            mad = (ad == max(ad)).nonzero()#indice des points ou la distance entre les deux splines est max
#             debug(mad=len(mad), dist=ad[mad])
            idx = mad[0]
            t = Td[idx][0]
            pos0 = d0[idx][0]
            try :
                s1.insertPoint(pos0)
            except ValueError as msg :#Impossible d'inserer le point (point double ?)
                rdebug(msg, pos0=pos0)
                debug(s0=s0.cpoints)
                debug(s1=s1.cpoints)
#                 rdebug()
                d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
                rdebug('Precision non atteinte, iteration %d : dist = %.2e '%(k,self.pourMille(d)))
                d1 = asarray(s1(Td))
                c1 = s1.cpoints.copy()
                break#TODO : pour le moment on break car sinon, ca marche pas
#                 debug(s0=s0.cpoints)
#                 debug(s1=s1.cpoints)
#                 exit()
#                 k -= 1# on revient sur le point problematique
#                 Td = linspace(0,1,len(Td)+len(Td)/10)
#                 s0.precision=len(Td)
#                 d0 = asarray(s0(Td))
            d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
            debug('iteration %d : dist = %.2e '%(k,self.pourMille(d)))
            d1 = asarray(s1(Td))
            c1 = s1.cpoints.copy()
            if self.pourMille(d)<eps : break
        if len(s1) == len(self) :#même nb de points : on garde la spline initiale.
            s1.cpoints = self.cpoints.copy()
            n0 = len(self)
            n1 = len(s1)
            return s1, self.pourMille(d),(n0,n1)
        #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de s1
        def refine(n=3):
            res = minimize(f1, np.ravel(s1.cpoints), (s0, Td,), 'BFGS', options={'maxiter':n})#, 'disp':True})
    #         print res.fun
            c1 = res.x
            c1.shape = (-1,2)
            #recalage des points 0 et -1
            c1[0] = c0[0]
            c1[-1] = c0[-1]
            s1.cpoints = c1
#         debug('Avant optimisation : dist = %.2g mm/m ... wait...'%(self.pourMille(d)))
#         for k in range(3) : refine()
#         d = np.sqrt(max(distance2PointSpline(c0[k], s1, t0=t).fun for k, t in enumerate(T0)))
#         debug('Apres optimisation : dist = %.2g mm/m '%(self.pourMille(d)))
        #courbure totale
        courburetotale0 = self.integraleCourbure(a=0.01, b=0.99)
        courburetotale1 = s1.integraleCourbure(a=0.01, b=0.99)
        ratio = courburetotale1/courburetotale0
        n0 = len(self)
        n1 = len(s1)
        debug('Apres ELAGAGE : dist = %.2g mm/m ; courbure totale :c1/c0 %f ; %d => %d '%(self.pourMille(d),courburetotale1/courburetotale0,n0,n1))
#         debug('Apres ELAGAGE, courbure totale c1/c0 = %f '%(courburetotale1/courburetotale0))

        if replace :
            self.cpoints = s1.cpoints
        return s1, self.pourMille(d),(n0,n1)

################################################################################

def adet(u,v):
    u"""Retourne 2*surface du triangle défini par les vecteurs u et v"""
    return abs(u[0]*v[1]-u[1]*v[0])

def dist2s(P0, P):
    u"""
        - distance entre les polygones et P0 et P au sens suivant :
        la surface en valeur absolue entre les deux polygones P0 et P,
        - P et P0 ont le meme nombre de points.
        - La surface est calculée comme la somme des surfaces (en val. abs.)
        des n-1 quadrilatères P0[i-1], P0[i], P[i], P[i-1]"""
    if len(P0) != len(P) :
        raise ValueError("len(P0)=%d doit etre egal a len(P)=%d"%(len(P0), len(P)))
    AB = zip(P0, P)[:-1]
    surf = 0.0
    for (a0,b0), (a1,b1) in zip(AB[1:], AB[:-1]) :
#         a0a1, a1b0, b0a1, a1b1 = a1-a0, b0-a1, a1-b0, b1-a1
        surf += adet(a1-a0, b0-a1) + adet(a1-b0, b1-a1)
    return surf

def distance(s1, s2, precision):
    u"""
    distance separant deux splines, au sens :
    valeur absolue de la surface de s1-s2
    """
#     s2.precision = s1.precision = precision
    T = linspace(0, 1, precision)
    d = s1(T) - s2(T)
    return np.linalg.norm(d)
# #############################

def elaguer(p0, eps=1.0e-3):
    Np = p0.precision#nb de points de discretisation
    D0 = p0.dpoints#les points discretisés de la spline de reference
    C0 = p0.cpoints#les points de contrôle de la spline de reference
#     L0 = p0.longueur()
    T = linspace(0, 1, Np)
    methode = p0.methode
#     D0courb = courbure(D0)
#     C0courb = courbure(C0)
    global nbiter
    nbiter = 0
    def f(C, T, methode):
        u"""Fonction à minimiser : distance de la spline S0 à la spline d'interpolation de C.
        C = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(C))"""
#         print T
        global nbiter
        nbiter += 1
        oldshape = C.shape#tableau a plat (np.ravel(C))
        C.shape = (len(C)/2,2)
#         Ts, Sx, Sy = splineInterpolation(C, 'c cubic')#la spline d'interpolation de C
        S = NSplineSimple(cpoints=C, methode=methode)
#         print Ts
        C.shape = oldshape
        D = S(T)#discretisation de la spline S
#         D = asarray(zip(Sx(T), Sy(T)))
#         return np.sqrt(np.linalg.norm(D0-D))
        E = D-D0
        E[ 0] *= 1.0e10#pénalisation des extrémités qui bouge (on ne veut pas qu'elles bougent)
        E[-1] *= 1.0e10

        norm = np.linalg.norm(E)
        if nbiter%100 == 0 : print 'iteration %d : %.2g'%(nbiter, norm)
        return norm

#     def fold(C, T):
#         u"""Fonction à minimiser : distance de la spline S0 à la spline d'interpolation de C.
#         C = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(C))"""
# #         print T
#         global nbiter
#         nbiter += 1
#         oldshape = C.shape#tableau a plat (np.ravel(C))
#         C.shape = (len(C)/2,2)
# #         Ts, Sx, Sy = splineInterpolation(C, 'c cubic')#la spline d'interpolation de C
#
#         S = NSplineSimple(cpoints=C, methode=methode)
# #         print Ts
#         C.shape = oldshape
#         D = S(T)#discretisation de la spline S
# #         D = asarray(zip(Sx(T), Sy(T)))
# #         return np.sqrt(np.linalg.norm(D0-D))
#         norm = np.linalg.norm(D0-D)
#         if nbiter%10 == 0 : print 'iteration %d : %.2g'%(nbiter, norm)
#         return norm

    def numPointsToKeep(C, lref, genre='shark'):
#         T = absCurv(C)
#         T, sx, sy = splineInterpolation(C, 'c cubic')
#         S = NSplineSimple(cpoints=C)
#         T = S.sx
#         T/=T[-1]
#         curv = scourbure((sx,sy), T)
        curv = courbure(C)*lref
        dcurv = np.abs(curv[:-1]-curv[1:])
#         debug(dcurv=dcurv)
        tokeep = dcurv>0.1#pour sharknose : 0.12
#         tokeep = dcurv>0.15#sans sharknose : 0.15
        #Dans tous les cas on en supprime un peu moins que la moitié
        nptk = set([idx for idx, tk in enumerate(tokeep) if tk])
        nptkp = set([idx+1 for idx in tokeep if idx+1 < len(curv)]) - nptk
        nptkm = set([idx-1 for idx in tokeep if idx-1 >= 0 ]) - nptk
#         debug(tokeep=nptk)
#         debug(nptkm=nptkm, nptkp=nptkp)
        nptk = nptk.union(nptkp.union(nptkm))
        return sorted(list(nptk))
#
#         encore = True
#         while (encore) :
#             encore = False
#             for k1, k2 in zip(nptk[:-1], nptk[1:]) :
#                 if np.linalg.norm(C[k1]-C[k2]) > 0.2 or k2-k1 > 10 :
#                     encore = True
#                     nptk.append((k1+k2)/2)
#             nptk.sort()
#             debug(nptk=nptk)
#
#         return nptk
#     n = 2 # On garde 1 point sur n
#     nba = len(C0)-1
    A = C0[ 0]
    B = C0[-1]
    corde = np.linalg.norm(B-A)
#     debug(corde)
#     numeros = list(set(range(0, len(p0), n) + [nba, len(p0)-1,]))
    numeros = list(set(numPointsToKeep(C0, corde) + [0,len(p0)-1,]))
    numeros.sort()
#     nnba = numeros.index(nba)#On conserve l'emplacement du BA pour le modifier apres optimisation.
#     numeros.remove(nba)
    fuera = list(set(range(len(p0))) - set(numeros))
    fuera.sort()
    print u'numeros des points conservé=%s'%str(numeros)
    print u'numéros des points supprimés=%s'%str(fuera)
    C0 = p0.cpoints[numeros].copy()
#     print 'C0=',C0.tolist()
    print 'f(C0,T) = %.2e'%f(np.ravel(C0),T,methode)
    res = minimize(f, np.ravel(C0), (T,methode), 'BFGS', options={'maxiter':0,'gtol':eps})#, 'disp':True})
    print 'res.fun=%.2e'%res.fun
    C1 = res.x.reshape((-1,2))
#     C1[ 0]   = A
#     C1[-1]   = B
    print 'x-C0', C1-C0
    p1 = NSplineSimple(cpoints=C1, precision=10*Np, methode=p0.methode, mode='courbure') #mode=p0.mode)#
    p0.precision = 10*Np
    p0._update()
    D0 = p0.dpoints
    D1 = p1.dpoints
    L0 = p0.dlongueur
    L1 = p1.dlongueur
    C1 = p1.cpoints
#     debug(len_C1=len(C1), len_T=len(T))

    print u'  longueur  L1          = %.2e '%abs(L1)
    print u'  δ-longueur            = %.2e '%abs(L1-L0)
    print u'  norme frobenius D0-D1 = %.2e '%norm(D0-D1, 'fro')
    print u'  max norme 2     D0-D1 = %.2e <=='%max(norm(D0-D1, 2, axis=1))
    print u'  distance(p0,p1)/L0    = %.2e '%(distance(p0, p1, Np)/L0)
    print u'  f(p1)                 = %.2e '%f(np.ravel(C1), linspace(0, 1, 10*Np),methode)
    print u'  nb de points éliminés = %d / %d'%(len(fuera), len(p0))

#     from matplotlib import pyplot
# #     C0 = p0.cpoints[numeros]
#     minx, maxx = np.min(C0[:,0])-0.1, np.max(C0[:,0])+0.1
#     miny, maxy = np.min(C0[:,1])-0.1, np.max(C0[:,1])+0.1
#     pyplot.plot(
#                 D0[:,0], D0[:,1],'r-',
#                 C0[:,0], C0[:,1],'ro',
#                 D1[:,0], D1[:,1] ,'b-',
#                 C1[:,0], C1[:,1] ,'b>',
#                 (minx,maxx,0.0,0.0), (0.0, 0.0,miny, maxy), 'w.',#pour cadre large
#                 )
#     pyplot.show()
    return p1
######################
global nbiter
nbfunc=0

def f1(c, s0, Td):
    u"""Fonction à minimiser :
    distance des points de contrôle s0.cpoints à la spline d'interpolation de c.
    c = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(c))
    s0 = spline de reference
    Td = les abscurv dans[0,1] des points de discrétisation pour la comparaison"""
#         print T
    global nbfunc
    nbfunc += 1
    oldshape = c.shape#tableau a plat (np.ravel(C))
#     debug (C.shape, len(C)/2,2)
    c.shape = (len(c)/2,2)
    s = NSplineSimple(cpoints=c, methode=s0.methode)
#     E = S0(T)-S(T)
    E = [distance2PointSpline(s0.cpoints[k], s, t0=t).fun for k, t in enumerate(s0.sx.x)]
    E[0]*=1.0e3
    E[-1]*=1.0e3
#     norme = np.sqrt(sum(E[1:-1]))
    norme = np.sqrt(sum(E))
    c.shape = oldshape
#     norm = np.linalg.norm(E)
#     if nbfunc%100 == 0 :
#         print 'nbfunc %d : %.2g'%(nbfunc, norme)
    return norme
def f(C, S0, T):
    u"""Fonction à minimiser : distance de la spline S0 à la spline d'interpolation de C.
    C = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(C))
    T = les abscurv dans[0,1] des points de discrétisation pour la comparaison"""
#         print T
    global nbiter
    nbiter += 1
    oldshape = C.shape#tableau a plat (np.ravel(C))
#     debug (C.shape, len(C)/2,2)
    C.shape = (len(C)/2,2)
    S = NSplineSimple(cpoints=C, methode=S0.methode)
    E = S0(T)-S(T)
    E[0]*=1.0e10
    E[-1]*=1.0e10
    C.shape = oldshape
    norm = np.linalg.norm(E)
    if nbiter%100 == 0 :
        print 'iteration %d : %.2g'%(nbiter, norm)
    return norm

def testElaguer(filename, trados, fonction='elaguer'):
    from matplotlib.widgets import CheckButtons
    from matplotlib import pyplot as plt
    if fonction=='ajuster' :
        msg=u"""L'élagage par "ajustement" ne marche pas, et il est très long.
        Je le supprime des tests.
        :TODO: voir au besoin si on peut l'améliorer"""
        rdebug(msg)
        return

    def compareCourbures(s0,s1, corde, trados):
        u"""s0 et s1 sont des intrados (ou des extrados), s0:original, s1:elagué
        On trace les deux splines et les log des abs(courbures), mises à l'échelle pour y voir qque chose"""
#         from matplotlib import pyplot as plt
#         from matplotlib.widgets import CheckButtons
    #         nbpd = self.precision
        s1.nbpe = nbpe = s0.nbpe
        titre = ' '.join(('courbures', s0.name, s1.name))
        plt.title(titre)
    #     E = s0.epoints
        _, ax = plt.subplots()
        T = linspace(0,1,100)

#         s0.echantillonner(nbpe)
        D0 = s0.dpoints
        C0 = s0.cpoints
        spline0, = ax.plot(D0[:,0], D0[:,1], 'r-', lw=1)
    #         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
        control0, = ax.plot(C0[:,0], C0[:,1], 'ro', lw=1, label=u'%d points contrôle'%len(C0))
        courbure0 = s0.courbure(T)
        courbure1 = s1.courbure(T)
        minc = min(min(courbure0), min(courbure1))
        maxc = max(max(abs(courbure0)), max(abs(courbure1)))
        courbure0 += (1.0 + abs(minc))
        courbure0 = log(courbure0)
        courbure0 /= maxc
        if trados == 'e' :
            courbure0 = courbure0[::-1]
        courbure0, = ax.plot(corde*T[1:], corde*courbure0[1:])
        D1 = s1.dpoints
        C1 = s1.cpoints
        spline1, = ax.plot(D1[:,0], D1[:,1], 'b-', lw=1)
        control1, = ax.plot(C1[:,0], C1[:,1], 'bo', lw=1, label=u'%d points contrôle'%len(C1))

        courbure1 += (1.0 + abs(minc))
        courbure1 = log(courbure1)
        courbure1 /= maxc
        if trados == 'e' :
            courbure1 = courbure1[::-1]
        courbure1, = ax.plot(corde*T[1:], corde*courbure1[1:])
        ax.legend(loc='upper right')
        buttons = ['control0', 'control1', 'courbure0', 'courbure1', 'spline0', 'spline1']
        values = [True, True, True, True, True, True]
        draws = [control0, control1, courbure0, courbure1, spline0, spline1]
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        rax = plt.axes([0.05, 0.4, 0.1, 0.15])
        check = CheckButtons(rax, buttons, values)

        def func(label):
            if label == 'spline0':
                spline0.set_visible(not spline0.get_visible())
            elif label == 'spline1':
                spline1.set_visible(not spline1.get_visible())
            elif label == 'control0':
                control0.set_visible(not control0.get_visible())
            elif label == 'control1':
                control1.set_visible(not control1.get_visible())
            elif label == 'courbure0':
                courbure0.set_visible(not courbure0.get_visible())
            elif label == 'courbure1':
                courbure1.set_visible(not courbure1.get_visible())
            else :
                draw = draws[buttons.index(label)]
                draw.set_visible(not draw.get_visible())
            plt.draw()
        check.on_clicked(func)
        plt.show()
        return plt

    def comparePointsSpline(P,s):
        u"""comparaison graphique, P famille de points s spline
        On trace les projections de P sur la spline s"""
        projections = []
        for p in P :
            optres = distance2PointSpline(p, s)
            #t=abscisse scurviligne de la pj p' de p sur s, et d=distance d(p,p')
            t, d = optres.x, optres.fun
            projections.append(s(t))
        projections = asarray(projections)
        s.plot(plt, more=([P[:,0], P[:,1],'g.', 'C0'], [projections[:,0], projections[:,1],'y.','projections']))
        plt.show()

    u"""Extraction intrados et extrados, normalisation"""
#     if '.spl' in filename :
#         with open(filename,'r') as f :
#             lines = f.readlines()
#         for line in lines :
#             dump = eval(line)
#             points = dump['cpoints']
#         print S

    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points*=(1000/corde)#en mm
    corde = 1000
    if trados == 'e' :#Extrados
        methode = ('cubic',((2, 0, 0), (1, 0, -corde/2))) #extrados
        c0 = points[:1+nba]
    elif trados == 'i':#intrados
        methode = ('cubic',((1, 0, -corde/4), (2, 0, 0)))#intrados
        c0 = points[nba:]
    u"""Construction de la spline d'interpolation du trados"""
    precision = 1000
    T0 = linspace(0,1,precision)
    s0 = NSplineSimple(cpoints=c0, precision=precision, methode=methode, mode='courbure')
    u"""La spline elaguée"""
    if fonction == 'elaguer' :
        s1, d, (n0,n1) = s0.elaguer(1)
    elif fonction == 'ajuster' :
        s1, d, (n0,n1) = s0.ajuster(1)
    d0 = s0(T0)
    d1 = s1(T0)
    print u'  max dist(C0, s1)      = %.2e max des distances aux points de contrôle (‰) <  ################'%(d)
    print u'  max norme 2     D0-D1 = %.2e max des distances point à point des splines discrétisées'%max(norm(d0-d1, 2, axis=1))
    print u'  norme frobenius D0-D1 = %.2e somme des distances point à point des splines discrétisées'%norm(d0-d1, 'fro')
    print u'  norme frobenius D0-D1 = %.2e moyenne des distances point à point des splines discrétisées'%(norm(d0-d1, 'fro')/len(d0))
    print u'  nb pts contrôle       : %d => %d'%(len(s0),len(s1))
    compareCourbures(s0, s1, corde, trados)
#     comparePointsSpline(s0.cpoints, s1)
    return


def testElagage1(filename):
    u"""Elagage d'une spline simple=suppression de points de contrôle.
    filename doit contenir un profil.
    1-On détermine son extrados, son intrados, le ba, le bf, la corde...
    2-On construit la spline de l'intrados
    3-On l'élague
    3-On affiche"""
    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points/=corde
    corde = 1.0
#     return
# #     #Extrados
#     methode = ('cubic',((2, 0, 0), (1, 0, -corde/2))) #extrados
#     s0 = NSplineSimple(points=points[:nba+1], precision=100, methode=methode)
    #intrados
    methode = ('cubic',((1, 0, -corde/2), (2, 0, 0)))#intrados
    s0 = NSplineSimple(points=points[nba:], precision=100, methode=methode, mode='courbure')
    print s0
    s1 = elaguer(s0, eps=1.0e-4)
#     s1 = elaguer(s1, eps=1.0e-4)
#     s1 = elaguer(s1, eps=1.0e-4)

    Td = linspace(0,1,1000)
    c0 = s0.cpoints.copy()*1000#en mm
    d0 = asarray(s0(Td))*1000
#     debug(d0_shape=d0.shape)
    c1 = s1.cpoints.copy()*1000
    d1 = asarray(s1(Td))*1000
#     debug(d1_shape=d1.shape)
#     d1 = s1.dpoints.copy()
#     plt = s0.plot()
    from matplotlib.widgets import CheckButtons
    from matplotlib import pyplot as plt
    _, ax = plt.subplots()
    spline0, = ax.plot(d0[:,0], d0[:,1], 'b-', lw=1)
#     echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
    control0, = ax.plot(c0[:,0], c0[:,1], 'b.', lw=1)
    spline1, = ax.plot(d1[:,0], d1[:,1], 'r-', lw=1)
#     echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
    control1, = ax.plot(c1[:,0], c1[:,1], 'r.', lw=1)

    plt.subplots_adjust(left=0.2)
    plt.axis('equal')
    rax = plt.axes([0.05, 0.4, 0.1, 0.15])
    check = CheckButtons(rax, ('control0', 'control1', 'spline0', 'spline1'), (True, True, True, True))

    def func(label):
        if label == 'spline0':
            spline0.set_visible(not spline0.get_visible())
        elif label == 'control0':
            control0.set_visible(not control0.get_visible())
        elif label == 'spline1':
            spline1.set_visible(not spline1.get_visible())
        elif label == 'control1':
            control1.set_visible(not control1.get_visible())
        plt.draw()
    check.on_clicked(func)

    plt.show()
#     for precision in range(1000, 1001, 1000):
#         print 'precision=%d, distance(s0, s1)=%.2e'%(precision, distance(s0, s1, precision))
#     p1 = elaguer2(p0)

def testConstructeurs(filename) :
    from matplotlib import pyplot as plt
    p = NSplineSimple()
    print p
    p = NSplineSimple(cpoints=np.zeros((0,2)),
                      parent=None,
                      methode=('cubic','periodic'),
                      mode='rayon',
                      name='SPLINE1')
    print p
    for point in p.qcpolygon :
        print point
    for point in p.qdpolygon :
        print point
    debug(filename)

#     return
#     p = NSplineSimple()
    p.cpoints = pointsFrom(filename)
    p.plot(plt, titre='cpoints=pointsFrom()')
    S0 = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic','periodic'),
                       mode='courbure',name=filename.name)
    debug(S0)

    S0.plot(plt, titre=filename.name)

    filename = Path(RUNS_DIR, 'spline.pkl')
    S0.save(filename=filename)
    S0.open(filename)
    debug(S0)
    S0.plot(plt, titre='open pkl')
    plt.show()
    S1 = S0.copy()
    debug(S1)
    S1.plot(plt, titre='copy')

#     return

def testMethodes(filename):
    from matplotlib import pyplot as plt
    k = 0
    for methode in (
                    ('cubic','not-a-knot'),
                    ('cubic','periodic'),
                    ('cubic','natural'),
                    ('cubic',((1, 0, -5), (1, 0, -5))),
                    ('cubic',((2, 0, 0), (1, 0, -5))), #extrados
                    ('cubic',((1, 0, -5), (2, 0, 0))), #intrados
                    ):
        k +=1
        S = NSplineSimple(points=pointsFrom(filename),
                          methode=methode,
                          mode='courbure',
                          name='SPLINE%d'%(k+1))
        S.plot(plt, titre=filename.name+str(methode))
        S.plotCourbure()

#     filename=Path(VALIDATION_DIR,'simple','simple2d2.gnu')
    for methode, mode in (
#                     (('ius',1),'linear'),
#                     (('ius',2),'linear'),
#                     (('ius',3),'courbure'),
                    (('ius',4),'courbure'),
                    (('ius',7),'courbure'),
                    (('us',1),'linear'),
                    (('us',2),'linear'),
                    (('us',3),'courbure'),
                    (('us',4),'courbure'),
                    (('us',7),'courbure'),
                    ):
        try :
            S = NSplineSimple(points=pointsFrom(filename), methode=methode, mode=mode)
            print S
            rdebug(u'\n'+70*u'*'+u'\n'+ u' TODO : methode=("us",k) ne marche pas. Va y avoir plein d\'exceptions'+u'\n'+70*u'*'+u'\n')

            S.plot(plt, titre=filename.name+str(methode))
        except Exception as msg:
            rdebug('(methode=%s; mode=%s) => pas de spline'%(str(methode),str(mode)),
                   '%s : %s'%(type(msg).__name__,unicode (msg)))
    plt.show()

def testModifGlobales(filename) :
    from matplotlib import pyplot as plt
#     numfig = 0
#     plt.figure(numfig)
    S0 = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                       precision=1000,
                       mode='courbure',
                       name='SPLINE01')
    S0.plot(plt, titre='NSplineSimple(points=pointsFrom(filename),...')
#     numfig += 1
#     plt.figure(numfig)
    S0.cpoints = pointsFrom(filename)
    S0.plot(plt, titre=filename.name+'cpoints = pointsFrom()')

#     numfig += 1
#     plt.figure(numfig)
    S0.hardScale((2,0.5))
    S0.plot(plt, titre='hardScale((2,0.5))')

#     numfig += 1
#     plt.figure(numfig)
    S0.symetriser(0)
    S0.plot(plt, titre='symetriser(0)')

#     numfig += 1
#     plt.figure(numfig)
    S0.symetriser(1)
    S0.plot(plt, titre='symetriser(1)')

#     numfig += 1
#     plt.figure(numfig)
    S0.hardRotate(30)#en degres par defaut
    S0.plot(plt, titre='hardRotate(30)')

#     numfig += 1
#     plt.figure(numfig)
    S0.translate((2,3))
    S0.plot(plt, titre='translate((2,3))')

    plt.show()

def testModifLocales(filename)  :
    from matplotlib import pyplot as plt
    S = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                       mode='courbure',
                       name='SPLINE01')

    S.plot(plt, titre=filename.name)

    try : S.insertPoint((0,0))
    except ValueError as msg : print msg

    S.appendPoint((2,1.5))
    S.plot(plt, titre='appendPoint((2,1.5))')

    S.insertPoint((2,-2))
    S.plot(plt, titre='insertPoint((2,-2))')

    try :
        S.insertPoint((6,0))
        S.plot(plt, titre='insertPoint((6,0))')
    except ValueError as msg :
        rdebug(msg)
    S.removePoint(3)
    S.plot(plt, titre='removePoint(3)')

    S[1] = (4,-0.1)
    S.plot(plt, titre='S[1] = (4,-0.1)')
    print S[1]
#     S.plot(plt, )
    try : S._update()
    except AttributeError as msg : rdebug(msg)
    pprint(S.__dict__)
    dump = S.toDump()
    debug()
    pprint(dump)
    dump['methode'] = ('cubic','natural')
    S.load(dump)
    S.name='SPLINE3'
    S.plot(plt, titre='load(toDump())')
    print S
    plt.show()

def testSaveRead(filename):
    u"""Diverses manières de lire, construire, sauvegarder une NSplineSimple"""
    # from matplotlib import pyplot as plt
    print "    ==> NSplineSimple::__init__()"
    S = NSplineSimple(points=pointsFrom(filename),
                      methode=('cubic',((2, 0, 0), (1, 0, -5))),
                      mode='courbure',
                      name='SPLINE01')
    S4 = NSplineSimple()
    S4.load(S.toDump())
    dump4 = S4.toDump()
    fname = Path(filename.dirname, filename.namebase+'#.pkl')
    print fname
    print "    ==> cPickle::dump()/load()"
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = NSplineSimple(**dump)
#     S1.plot(plt)
#     plt.show()
    dump  = S.toDump()
    dump1 = S1.toDump()
    print 'dump==dump1 ?', dump==dump1

    print "    ==> file::read()/write() "
    fname = Path(filename.dirname, filename.namebase+'#.spl')
    S2 = NSplineSimple(**dump)
    dump2 = S2.toDump()
    with open(fname,'w') as f :
        f.write(str(S2.toDump()))
    with open(fname,'r') as f :
        dump = eval(f.read())
    print "    ==> NSplineSimple::load()"
    S3 = NSplineSimple(**dump)
    dump3 = S3.toDump()
#     dump3_1 = S3.toDump()
    print 'dump==dump1=dump2=dump3=dump4 ?', dump==dump1==dump2==dump3==dump4
    debug(dump=dump)
    debug(dump1=dump1)

    debug(dump2=dump2)
    debug(dump3=dump3)
    debug(dump4=dump4)
#     exit()

def testEchantillonnage(filename, trados):
    u"""echantillonnage entre ta et tb. [ta,tb]=[0,1] => echantillonnage complet"""
#     s = NSplineSimple(points=pointsFrom(filename))
    u"""Extraction intrados et extrados, normalisation"""
    from matplotlib import pyplot as plt

    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points*=(1000/corde) #en mm
    corde = 1000
    if trados == 'e' : #Extrados
        methode = ('cubic',((2, 0, 0), (1, 0, -corde/2)))#extrados
        c0 = points[:1+nba]
    elif trados == 'i': #Intrados
        methode = ('cubic',((1, 0, -corde/4), (2, 0, 0)))#intrados
        c0 = points[nba:]
    u"""Construction et elagage de la spline d'interpolation du trados"""
    precision = 1000
    T0 = linspace(0,1,precision)
    s0 = NSplineSimple(cpoints=c0, precision=precision, methode=methode, mode='courbure')
    s0.elaguer(1, replace=True)
    ne, ta, tb = 30, 0.2, 0.8
    E = s0.echantillonner(nbp=ne, ta=ta, tb=tb, mode='courbure')
    s0.plot(plt, more=[(E[:,0],E[:,1],'y^','echantillon \nn=%d, %.2f=>%.2f'%(ne,ta,tb)),])
#     plt.plot(E[:,0],E[:,1],'y^','echantillon')
    plt.show()

def testDivers(filename):
    msg=u"""Dans cet exemple, on a des points doubles p4=p5=p6 et p8=p9,
    donc la spline x(t) ne peut pas être construite.
    Que faire ?
    - supprimer p[6] <=== Oui c'est la solution adoptée ici
    - déplacer légèrement le point
    - faire une spline composée, splittée en 6
    """
    from matplotlib import pyplot as plt
    cpoints = [
               [ 1.00827138,  0.08545793],
               [ 1.0063704 ,  0.06720821],
               [ 1.00541204,  0.04849635],
               [ 1.00541204,  0.02962959],
               [ 1.0063704 ,  0.01091773],
               [ 1.0063704 ,  0.01091773],
               [ 1.0063704 ,  0.01091773],
               [ 1.00827138, -0.00733198],
               [ 1.12226291,  0.03906297],
               [ 1.12226291,  0.03906297],
               [ 1.12178274,  0.02023498],
               [ 1.00827138,  0.08545793],
               [ 1.00827138,  0.08545793],
               ]
    cpoints = pointsFrom(filename)
    S0 = NSplineSimple(cpoints=cpoints, name='points doubles')
    debug(S0)
    S0.plot(plt, show=True)
    ########################################
    msg=u"""S est une NSplineSimple. Teste les méthodes :
    - S.knots() qui retourne les knots de la spline.
        Normalement on doit avoir knot==abscurv()==sx.x==sy.x sauf à l'initialisation
    - S.abscurv(T) qui retourne les abscisses curvilignes VRAIES des points S(T)
    - S.abscurv(None) qui retourne les abscisses curvilignes NORMALISÉES des points S(T)
    - S.longueur par intégration numérique (voir la méthode S.abscurv(T))
    - S.distance2To(p) et affiche le résultat pour 3 points
        Calcule la distance du point 'p' à la spline 'S' et le point 'h' de la spline
        le plus proche de 'p'. Retourne également la precision et le nb d'appels à fonction
        Utilise :
        "https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html
    """
    dump = {'classename': 'NSplineSimple',
            'cpoints': cpoints,
            'methode': ('cubic', 'not-a-knot'),
#             'mode': 'segment',
#             'methode': ('ius', 1),
            'mode': 'rayon',
#             'mode': 'cpoints',
            'name': u'profil',
            'nbpe': 30,
            'precision': 1000,
            'role': 'NSplineSimple'}
    dump = {'cpoints': [[0.5279636979103088, 0.07332829385995865], [0.7137287259101868, 0.3275330364704132], [1.268811468596966, 0.3365727279966774], [1.1390328407287598, 0.07332829385995865], [1.2571670716747245, 0.2148051133421408], [1.2206038430660453, -0.0507648238639925], [1.5545598268508911, -0.0048885527066886425]],
#             'methode': ('cubic', 'not-a-knot'),
            'methode': ('ius', 3),
            'name': u'test',
            'nbpe': 30,
            'precision': 1000,
            'mode': 'lineaire'}
    debug(msg)
    S = NSplineSimple(**dump)
    S.plot(plt)
    knots = S.knots
    acn = S.abscurv(None)#T Normalisé == knots
    acr = S.abscurv(acn)
    debug('S.abscurv(None)',acn.tolist())
    debug('S.knots',knots.tolist())
    l = S.longueur
    debug('(acr-acn)*longueur',(acr-acn*l).tolist())
    debug('longueur=%.10g'%l)
    lr = acr[-1]
    debug('vraie abscurv[-1]=%.10g'%lr)
    ln = acn[-1]
    debug('abscurv[-1]=%.10g'%ln)
    S.precision=1000
    ld3 = S.dlongueur
    debug('S.dlongueur=%g'%ld3+u"(calculé sur %d points)"%S.precision)
    S.precision=10000
    ld4 = S.dlongueur
    debug('S.dlongueur=%g'%ld4+u"(calculé sur %d points)"%S.precision)
    S.precision=100000
    ld5 = S.dlongueur
    debug('S.dlongueur=%g'%ld5+u"(calculé sur %d points)"%S.precision)
    le = S.elongueur
    debug('S.elongueur=%g'%le)
    lc = S.clongueur
    debug('S.clongueur=%g'%lc)
    diff = [l-lr,l-ln,l-ld3,l-ld4,l-ld5,l-le,l-lc]
    plt.bar(range(len(diff)),diff)
    plt.xticks(range(len(diff)),'lr,ln,ld3,ld4,ld5,le,lc'.split(','))
    plt.title('differentes longueurs calculees')

    plt.show()
    p1, p2, p3 = [0.2,0.01], [1.1,0.2], [-0.1,0.3]
    res1, res2, res3 = S.distance2To(p1), S.distance2To(p2), S.distance2To(p3)
    pj1, pj2, pj3 = S(res1.x), S(res2.x), S(res3.x)
#     rdebug(str(res))
    debug("****  DISTANCE POINT-SPLINE  ****")
    debug('message="%s"'%res1.message,'d(S,p1)=%g'%res1.fun, 't=%g'%res1.x, 'nb appels a fct=%d'%res1.nfev)
    debug('message="%s"'%res2.message,'d(S,p2)=%g'%res2.fun, 't=%g'%res2.x, 'nb appels a fct=%d'%res2.nfev)
    debug('message="%s"'%res3.message,'d(S,p3)=%g'%res3.fun, 't=%g'%res3.x, 'nb appels a fct=%d'%res3.nfev)
    more = [([p1[0],pj1[0]],[p1[1],pj1[1]],'g-o', 'd(S,p1) : %.2g'%sqrt(res1.fun)),
            ([p2[0],pj2[0]],[p2[1],pj2[1]],'g-o', 'd(S,p2) : %.2g'%sqrt(res2.fun)),
            ([p3[0],pj3[0]],[p3[1],pj3[1]],'g-o', 'd(S,p3) : %.2g'%sqrt(res3.fun))]
    projections = [(pj1[0],pj1[1],'y.',''), (pj2[0],pj2[1],'y.', ''), (pj3[0],pj3[1],'y.', '')]
    more += projections
    S.plot(plt, more=more, titre='distance point-spline', show=True)
#     plt.plot(p1[0],p1[1],'go', label='p1 : %g'%res1.fun)
#     plt.legend()
    plt.show()
#     return
    debug('ORIGINAL', S=S)
    debug(rect=S.boundingRect())
#     debug('CLONE', clone=S.clone())
    debug('COPY', clone=S.copy())
    S.hardScale((0.5,0.5), S.centre)
    debug('APRES HARDSCALE',S=S)
    debug(rect=S.boundingRect())
    S.hardMove((10,10))
    debug('APRES HARDMOVE',S=S)
    debug(rect=S.boundingRect())
#     debug(top='%g'%rect.top(), left='%g'%rect.left(), width='%g'%rect.width(), height='%g'%rect.height(), )
#     debug(bottom='%g'%rect.bottom(), top='%g'%rect.top())
    S.plot(plt, titre="**dump")
#     S.elaguer(1, replace=True)
#     S.plot(plt, titre="elaguer")
def testPlacementRM():
    ##############################################
    ### debut construction exemple
    ##############################################
    rdebug()
    B1 = asarray([[0,0],[0.5,0.25],[0.5,0.8],[0,1.1],[0,1.7],[0.5,2]], float) #un S
    B2 = asarray([[2,0],[2.5,1],[2,2]], float)#un C à l'envers
#     delta = 0.0186#distance entre reperes de montage (en m)
    from matplotlib import pyplot as plt
    s1 = NSplineSimple(cpoints=B1, parent=None,
#                        methode=('ius',1),
                        methode=('cubic','natural'),
                       name='B1',
                       )
    s2 = NSplineSimple(cpoints=B2, parent=None,
#                        methode=('ius',1),
                        methode=('cubic','natural'),
                       name='B2',
                       )

#     ac1 = absCurv(s1.cpoints, normalise=False)
#     ac2 = absCurv(s2.cpoints, normalise=False)
    l1, l2 = s1.longueur, s2.longueur
    debug(u'avant hardScale',l1=l1, l2=l2)
    #Je les ramène à la même longueur
    s2.hardScale((l1/l2, l1/l2))
    l1, l2 = s1.longueur, s2.longueur
    lref = min(l1,l2)
    debug(u'après hardScale',l1=l1, l2=l2)
    B1, B2 = s1.epoints, s2.epoints#le tableaux echantillonnés que va me passer ML
    ##############################################
    ### fin construction exemple
    ##############################################
    delta = 0.2#distance entre reperes de montage (en m)
    s1, s2, T = placementReperesMontage(B1,B2,delta,)
    n=len(T)
    rm1, rm2 = s1(T), s2(T)
    ac1, ac2 = s1.abscurv(T), s2.abscurv(T)#les abs curv. vraies des reperes
    delta1 = ac1[1:]-ac1[:-1]
    delta2 = ac2[1:]-ac2[:-1]
    dd = np.max(np.abs(delta1-delta2))
    debug("Erreur : %.1e mm"%(1000*dd/lref))

    plt.plot(s1.dpoints[:,0], s1.dpoints[:,1], 'b-')
    plt.plot(s1.cpoints[:,0], s1.cpoints[:,1], 'bx', label=u'Points de contrôle des splines')
    plt.plot(s2.dpoints[:,0]-1, s2.dpoints[:,1], 'b-')
    plt.plot(s2.cpoints[:,0]-1, s2.cpoints[:,1], 'bx')
    plt.plot(rm1[:,0], rm1[:,1], 'ro', label=u'Repères de montage')
    plt.plot(rm2[:,0]-1, rm2[:,1], 'ro')
    plt.axis('equal')
    plt.legend()
    msg = u"""Deux bords de pièces à assembler, de longueur $l_1$=%.4g m $\simeq$ $l_2$=%.4g m
    Placements de %d repères de montage espacés de %.2g m.
    L'erreur est de %.2f mm soit %.2f‰"""%(l1,l2,n,delta,1000*dd,1000*dd/lref)
    plt.title(msg)
    plt.show()

def placementReperesMontage(B1,B2,delta,):
    u"""Besoin:
    ----------
    Calculer une correction sur les repères de montage des cloisons, intrados, extrados...
    pour éviter qu'ils se trouvent légèrement décalés à cause du décalage
    de montage et de la courbure des pièces
    Actuellement les repères sont décalés perpendiculairement au bords de pièce,
    du coup ils s'écartent si le bord est convexe et se rapprochent si le bord est concave.
    L'idée est de recalculer les repères de montage avec respect des distances
    des points initiaux des bords de pièces,
    ainsi les repères de toutes les pièces vont s'aligner parfaitement.

    Algo correction:
    -------------------
    Soit B un tableau de points initiaux 2D correspondant à un bord de pièce à plat.
    Soit RM1 le tableau des repères de montage calculés avec une fréquence d'échantillonage de 1
    (donc pas sous échatillonés, à un point de B correspond un point de RM1,
    je précise car après je les réchantillone).
    1. Calculer le tableau D des distances successives des points de B
    2. Calculer une spline SRM1 qui passe par les points de RM1
    3. Calculer un tableau de points RM2 sur la spline SRM1
       tel que les distances successives D soient respectées

    => le plus simple pour moi serait une fonction qui prenne en entrée B et RM1
    (donc deux tableaux de points avec le même nombre de points 2D qui
    représentent des courbes ayant à peu près la même longueur géométrique...
    mais la 2e peut être plus ou moins longue que la 1ère,
    je ne sais pas comment tu gères si les points dépassent au bout ou pas...)
    et qui sorte RM2

    Réponse : la fonction placementReperesMontage(B1, B2, delta)
    -------
    :param B1, B2: les tableaux de points des deux bords à assembler
        (quelconques mais de même longueur, à peu près)
    :type B1, B2 : np.ndarray((n,2))
    :param delta: (float) le  pas des repères
        (distance entre deux repères consécutifs, constant pour le moment mais
        ça peut se changer)
    :return S1, S2, T: deux splines NSplinesSimple, et un tableau T de paramètres t_i,
        de type ndarray((n,1)).
    :utilisation:
        - n=len(T) est le nombre de repères calculé, les reperes eux mêmes sont obtenus par
        - S1(T): un ndarray((n,2)) avec les coordonnées des repères sur B1
        - S2(T): un ndarray((n,2)) avec les coordonnées des repères sur B2
    :remarque: sur le cas test, l'erreur est de 0,4 mm pour des bords de longueur 2 m environ.
        Soit 2 pour 10 000.
    """
    ##############################################
    ### construction des repères de montage
    ##############################################
    s1 = NSplineSimple(cpoints=B1, parent=None,
#                         methode=('cubic','natural'),
                        methode=('ius',1),
                       name='s1')
    s2 = NSplineSimple(cpoints=B2, parent=None,
#                         methode=('cubic','natural'),
                        methode=('ius',1),
                       name='s2')
    l1, l2 = s1.longueur, s2.longueur
    lref = min(l1,l2)
    debug(longueur1=l1, longueur2=l2)
    n, r = divmod(lref, delta)#nb de reperes de montage
    debug(u'nb reperes=%d, reste=%.2g'%(n,r))
    n = int(n)
    T = arange(1,n+1)*(delta/lref)
    return s1, s2, T

def testCorrectionRM():
    ##############################################
    ### debut construction exemple
    ##############################################
    from matplotlib import pyplot as plt
    rdebug()
    vc = 0.25
    B = asarray([[0,0],[0.5,0.25],[0.5,0.8],[0,1.1],[0,1.7],[0.5,2],
                 [0.6,1.8],[1.5,1.0]], float) #un S
    # construction valeurs couture
    SB1 = NSplineSimple(cpoints=B, parent=None,
#                         methode=('ius',1),
                        methode=('cubic','natural'),
                        name='B',
                        )
#     SB1.plot(plt)
    knots = SB1.knots#==B
    tangentes = SB1(knots,1)#dérivées premières
#     debug(tangentes=tangentes)
    normales = asarray(zip(-tangentes[:,1],tangentes[:,0]))
    for k, normale in enumerate(normales):
        normales[k] = normale/np.linalg.norm(normale)
#     u = normales[0]
#     debug (u[0]**2+u[1]**2)
#     debug(normales=normales)
    #points couture
    Bc = B + vc*normales
#     debug(Bc=Bc)
    #spline coutures/reperes
    C = NSplineSimple(cpoints=Bc,
                      methode=('cubic','natural'),
                      name='v.coutures')
#     debug(C)
#     C.plot(plt)
    DC = C.dpoints
    DB = SB1.dpoints
    ####################
    ## reperes mal placés
    ####################
    nrep = len(B)#nb reperes
    TR = absCurv(B, normalise=True)
    R = C(TR)#repere sur spline couture, mal placés
    ##############################################
    ### fin construction exemple
    ### les distances dans B ne sont pas les mêmes
    ###    que les distances dans R
    ##############################################
#     B,R = R,B#permutation B et R
#     DB,DC = DC, DB#permutation tracés
    ###################
    lR = absCurv(R, normalise=False)[-1]
    lB = absCurv(B, normalise=False)[-1]
    ##############################################
    #  correction
    ##############################################

    SCR, T = correctionReperesMontage(B, R, 'test')
    CR = SCR(T)

    ##############################################
    #tracé
    ##############################################
    plt.subplot(1,2,1)
    plt.plot(DC[:,0], DC[:,1],'g,')
#     plt.plot(CC[:,0], CC[:,1],'r.', label='spline coutures')
    plt.plot(DB[:,0], DB[:,1],'b,')
#     plt.plot(CB[:,0], CB[:,1],'b.', label='spline bord')
    plt.suptitle(u'Cas ou longueur(R) = %.5g >  %.5g = longueur(B)'%(lR, lB), fontsize=16, fontname='serif')
    plt.plot(B[:,0], B[:,1],'b-o',label=u'Points du bord (%d points)'%nrep)
    plt.plot(R[:,0], R[:,1],'g-o',label=u"Initialement %d repères mal placés"%nrep)
    plt.plot(CR[:,0], CR[:,1],'ro',label=u'Correction : %d repères '%nrep)
    plt.legend()
    plt.axis('equal')
#     plt.show()
    def deltaLongueur(B, R, sCR, T):
        Bac, Rac, CRac = absCurv(B), absCurv(R), sCR.abscurv(T)
        lB, lR, lCR = diff(Bac), diff(Rac), diff(CRac)
        return lB-lR, lB-lCR
    lBR, lBCR = deltaLongueur(B, R, SCR, T)
    lB, lR, lCR = absCurv(B)[-1], absCurv(R)[-1], SCR.abscurv(T[-1])
    plt.subplot(1,2,2)
    plt.title(u"$\Delta$longueurs. (Longueur totale R Corrigé=%.5g)"%(lCR))
    plt.plot(lBR,'b.-', label=u'initial : max($\delta l)=$%.2g'%max(abs(lBR)))
    plt.plot(lBCR,'r.-', label=u'corrigé : max($\delta l)=$%.2g'%max(abs(lBCR)))
    plt.plot(lBR,'b.')
    plt.plot(lBCR,'r.')

    plt.legend()
    plt.show()
#     P7 = asarray(C(T[7]))
#     C.precision*=10
#     pts = C.dpoints-P7
#     debug([(i, norm(p)) for i, p in enumerate(pts[5261:5565])])
#     debug(C)


#     msg = u"""Deux bords de pièces à assembler, de longueur $l_1$=%.4g m $\simeq$ $l_2$=%.4g m
#     Placements de %d repères de montage espacés de %.2g m.
#     L'erreur est de %.2f mm soit %.2f‰"""%(l1,l2,n,delta,1000*dd,1000*dd/lref)
#     plt.title(msg)
#     plt.show()


def correctionReperesMontage(B,R, mode='production'):
    u"""Besoin:
    ----------
    Calculer une correction sur les repères de montage des cloisons, intrados, extrados...
    pour éviter qu'ils se trouvent légèrement décalés à cause du décalage
    de montage et de la courbure des pièces
    Actuellement les repères sont décalés perpendiculairement au bords de pièce,
    du coup ils s'écartent si le bord est convexe et se rapprochent si le bord est concave.
    L'idée est de recalculer les repères de montage avec respect des distances
    des points initiaux des bords de pièces,
    ainsi les repères de toutes les pièces vont s'aligner parfaitement.

    Algo correction:
    -------------------
    Soit B un tableau de points initiaux 2D correspondant à un bord de pièce à plat.
    Soit R le tableau des repères de montage calculés avec une fréquence d'échantillonage de 1
    (donc pas sous échatillonés, à un point de B correspond un point de R,
    je précise car après je les réchantillone).
    1. Calculer le tableau dB des distances successives des points de B
    2. Calculer une spline sR qui passe par les points de R
    3. Calculer un tableau de points R2 sur la spline sR
       tel que les distances successives dB soient respectées

    => le plus simple pour moi serait une fonction qui prenne en entrée B et R
    (donc deux tableaux de points avec le même nombre de points 2D qui
    représentent des courbes ayant à peu près la même longueur géométrique...
    mais la 2e peut être plus ou moins longue que la 1ère,
    je ne sais pas comment tu gères si les points dépassent au bout ou pas...)
    et qui sorte sR

    Réponse : la fonction correctionReperesMontage(B, R)
    -------
    qui fournit un tableau RC de n=len(B) repères de montage, répartis sur
    le **chemin polygonal** défini par R, de sorte que les distances entre points de RC
    soient les mêmes que les distances entre les points de B
    **Les distances entre points de RC sont mesurées le long du chemin R**
    (et non pas de point à point)

    :param B, R: points du bord et repères de montage d'origine (chemin R)
    :type B : np.ndarray((n,2))
    :type R : np.ndarray((m,2)) (n et m n'ont pas besoin d'être égaux)
    :return:
        - si mode = 'test' : (sR, T) = la spline polygonale passant par les
            repères d'origine, que j'appelle chemin polygonal ci-dessus,
            et les valeurs du paramètre t de sorte que les n repères de montage
            corrigés sont RC=sR(T).
        - si mode = 'production' :les repères de montage corrigés,
            sous forme d'un np.ndarray((n,2))
    :remarque: sur le cas test, construit par la fonction testCorrectionRM()
    la longueur du bord B est 4.303, la longueur du chemin R est 4.3365
    l'erreur max initiale sur les longueurs est 0.033 soit 8‰
    l'erreur max sur les longueurs corrigées est 4.5e-16 soit
        moins de 1.0e-15 (en valeur relative)
    """
    #spline repères
    sR = NSplineSimple(cpoints=R, methode=('ius',1))
    # distances à respecter (distances à B[0])
    distB = absCurv(B)
    # Comme sR est de degré 1, sur sR le paramètre t EST EGAL à
    # dist(sR(0.),sR(t))/sR.longueur == t
    #rdebug(u'verif dist(sR(0.),sR(t)) == t', longueur=sR.longueur)
    #l=sR.longueur
    #for t in linspace(0, 1., 20) :
    #    print (sR.abscurv(t)/(t*l), dist(sR(0),sR(t)))
    # OK : Sur sR : dist(sR(0.),sR(t)) == t*sR.longueur
    # paramètres T correspondant aux distances à respecter
    T = distB/sR.longueur
    if mode == 'test' :
        return sR, T#pour tests
    elif mode == 'production' :
        return sR(T)#en production

if __name__=="__main__":
#     placementReperesMontage(T=np.asarray([[2,-2],[1,-1],[0,0],[1,1],[2,2]]),
#                              TR=np.asarray([[1.5,-2.5],[0.5,-1.5],[-math.sqrt(2)*0.5,0],[0.5,1.5],[1.5,2.5]]))
#     exit()
    # app = QApplication(sys.argv)
    s = NSplineSimple()
    debug('constructeur vide',s)
#     testElagage2(None)
#     exit()
    filename = Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename = Path(VALIDATION_DIR,'1.gnu')
    filename = Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    filename = Path(VALIDATION_DIR,'unenervure2d.gnu')
#                         filename = Path(VALIDATION_DIR,'faial2.dxf')#marche pas
    filename = Path(VALIDATION_DIR,'simple','anguleux.gnu')
    filename = Path(VALIDATION_DIR,'simple','intrados.gnu')
    filename = Path(VALIDATION_DIR,'simple','extrados.gnu')
    filename = Path(VALIDATION_DIR,'simple','profil.gnu')
    filename = Path(VALIDATION_DIR,'simple','demi-cercle.gnu')
    filename = Path(VALIDATION_DIR,'reference.pts')
    filename = Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename = Path(VALIDATION_DIR,'shark.pts')
    filename = Path(WORK_DIR,'profils','D2v5p2.pts')
    filename = Path(VALIDATION_DIR,'diamirM.pts')
    filename = Path(VALIDATION_DIR,'P0.spl')
    with open(Path(VALIDATION_DIR,'blocjonc.spl'),'r') as f :
        dmodel = eval(f.read())['dmodel']
#         debug(s)
    model = NSplineSimple(**dmodel)
    debug(model)
#     exit()
    #
    rdebug(u'\n*************************************\n TODO : methode=("us",k) ne marche pas\n*************************************\n ')
    if 1 :
        testCorrectionRM()
        debug('fin testCorrectionRM')
    if 1 :
        testPlacementRM()
        debug('fin testPlacementRM')
    if 1:
        testSaveRead(filename)
        debug('fin testSaveRead')
    if 1 :
        try :
            testDivers(filename)
        except AttributeError as msg:
#             stack()
            rdebug('AttributeError :'+str(msg))
    #         raise
        debug('fin testDivers')
#     sys.exit(app.exec_())
    if 1:
        testEchantillonnage(filename, trados='e')
        testEchantillonnage(filename, trados='i')
        debug('fin testEchantillonnage')
    if 1:
        testConstructeurs(filename)
        debug('fin testConstructeurs')
    if 1:
        testMethodes(filename)
        debug('fin testMethodes')
    if 1:
        testModifGlobales(filename)
        debug('fin testModifGlobales')
    if 0:
        testModifLocales(filename)
        debug('fin testModifLocales')
    if 0:
        testElaguer(filename, trados='e', fonction='ajuster')
        debug('fin testElaguer extrados')
    if 0:
        testElaguer(filename, trados='i', fonction='ajuster')
        debug('fin testElaguer intrados')
    if 0:
        testElaguer(filename, trados='e', fonction='elaguer')
        debug('fin testElaguer extrados')
    if 0:
        testElaguer(filename, trados='i', fonction='elaguer')
        debug('fin testElaguer intrados')
    if 0:
        testElagage1(filename)
        debug('fin testElagage1')
    rdebug('fin tests')
    sys.exit(app.exec_())

