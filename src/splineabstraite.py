#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
import logging
from numbers import Number
# logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
import sys,os,math
import numpy as np
# from array import array

# from PyQt4.QtGui import (QApplication, QPolygonF)
# from PyQt4.QtCore import (Qt, QPointF, QString, QObject, QRectF)
# from config import VALIDATION_DIR, WORK_DIR
from scipy.optimize import newton, minimize, minimize_scalar
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, UnivariateSpline
import cPickle
from utilitaires import (Path, whoami, debug, rdebug,stack,rstack,
                        absCurv,dist2, baryCentre, centreGravite, simpson)
# from gui.graphicsbase.graphicscommon import (qpolygonFrom, pointsFrom)
# from inout.writerpts import writeProfil
# from preferences import SplinePrefs, NPrefs
# from pprint import pprint
# from inout.lecteurs import LecteurUniversel
# from config import RUNS_DIR
from numpy.core.defchararray import center
from numpy import asarray, sqrt
from scipy.integrate import quad


def absCurvReal(S, T):
    u"""
    Calcule et retourne l'abscisse curviligne réelle des points S(T) sur la spline S.
    L'abscisse curviligne d'un point de paramètre t dans [0,1] est
    l'intégrale de 0 à t de phi(t) = sqrt(S.sx(t)**2+S.sy(t)**2)
    L'intégration est assurée par scipy.integrate.quad()
    Si la spline S a trop de points de contrôle, ca rame et l'erreur est importante
    err
    :param S: une NSplineAbstract
    :param T: les n paramètres t des points dont on veut l'abscisse curviligne.
        Ces n paramètres doivent être dans l'intervalle [0,1]
    :type T: au choix
        - un ndarray de shape (n,1) à valeurs réelles dans [0,1],
        - une liste de n valeurs réelles dans [0,1],
        - un tuple de n valeurs réelles dans [0,1],
        - un réel unique t dans [0,1]
    :return ac, err:
        - deux ndarray de n réels avec ac = abscisses curvilignes, err erreur estimée
        - ou bien deux réels, si T est réel
        (voir scipy.integrate)
    """
    if isinstance(T, Number) :
        #On integre sur les sous intervalles de S délimités par les knots
        phi = lambda s : sqrt(S.sx(s,1)**2+S.sy(s,1)**2) #la fonction à integrer
        bornes = [tk for tk in S.knots if tk<T]+[T]#bornes de sous intervalles
        intervalles = zip(bornes[:-1], bornes[1:])#les intervalles
        ac, err = 0, 0
        for (t1, t2) in intervalles :
            int_t1_t2, err12 = quad(phi, t1, t2)#integration intervalle [t1, t2]
            ac += int_t1_t2
            err = max(err,err12)
        return ac, err

#             return ac1+ac2, max(err1, err2)
    else :
        res = asarray([absCurvReal(S, t) for t in T])
#         res = asarray([quad(lambda s : sqrt(S.sx(s,1)**2+S.sy(s,1)**2), 0.0, t) for t in T])
        return res[:,0],  res[:,1]



def distance2PointSpline(p, S, t0=0.5, tol=1.0e-9):
    u"""
    Comme son nom l'indique, calcule le carré de la plus courte distance d'un point p à
    une spline paramétrée S(t) = x(t), y(t), 0<=t<=1
    :param S: la spline paramétrique
    :type S: NSplineAbstract.
    :param p: (x,y)=(float, float) le point.
    :param t0: float, initialisation des itérations
    :param tol: float, tolérance passée à scipy.optimize.minimize_scalar
    :return: le resultat res,
    voir https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
    en particulier :

        - res.x : float, valeur de t réalisant cette distance
        - res.nfev : int, nb evaluation de phi(t)
        - res.fun : float, valeur finale de phi(t)

    local : la fonction phi(t) à optimiser est la distance de p à S(t)
    """
    a, b = p[0], p[1]
    x, y = S.sx, S.sy

#     def phi(t):
#         return (a-x(t))**2 + (b-y(t))**2
#     def dphi(t):
#         d1 = 2*(x(t,1)*(x(t)-a) + y(t,1)*(y(t)-b))
#         return d1
#
#     def d2phi(t) :
#         d2 = 2*(x(t,2)*(x(t)-a) + x(t,1)**2 +
#                 y(t,2)*(y(t)-a) + y(t,1)**2)
#         return d2
#     res = minimize_scalar(phi, bounds=(0.0, 1.0), method='bounded', tol=1.0e-9)
    res = minimize_scalar(lambda t: (a-x(t))**2 + (b-y(t))**2,
                                bounds=(0.0, 1.0),
                                method='bounded',
                                options={'xatol':1.0e-9})
#     res = newton(dphi, t0, fprime=d2phi,
#                 tol=1.0e-3, maxiter=50)
    return res

def computeSpline(cpoints, methode):
    u"""
    Calcule la spline (sx, sy) considérée comme une courbe paramétrée sx(t), sy(t).
    sx et sy sont deux splines à une seule variable au sens scipy.
    Si cpoints contient des points double consécutifs, la construction échoue.

    Retourne: T, sx, sy
    --------
    - T np.ndarray(n) = abscisses curvilignes des points de cpoints
    - sx et sy deux splines d'interpolation de cpoints, i.e. vérifiant [sx(T[i]), sy(T[i])] = cpoints[i]
        sx, sy = None, None si cpoints.shape = (0,2) ou cpoints = None
        si cpoints contient deux points consécutifs identique, alors
        l'exception "ValueError: `x` must be strictly increasing sequence." est levée
    Parametres :
    ----------
    - cpoints = np.ndarray((n,2)) les points de contrôle.
    - methode : peut prendre differentes valeurs. Il est de la forme :
        methode = (type,parametres)
        avec :
    # type est une chaîne de caractères parmi ('cubic', 'ius', 'us')
        c'est le type de spline (suivant la terminologie scipy.interpolate).
    # paramètres sont les paramètres associés au type de spline
    # si type = 'cubic' c'est une spline d'INTERPOLATION
        - parametres = 'periodic' pour splines fermées, sans point anguleux (trous par exemple)
        - parametres = 'clamped' fixe à zéro les derivées premières de x et y aux extrémités
                équivaut à parametres=((1, 0, 1, 0), (1, 0, 1, 0)) (un extrados p.ex.)
        - parametres = 'natural' fixe à zéro les derivées secondes de x et y aux extrémités
                équivaut à parametres=((2, 0, 2, 0), (2, 0, 2, 0))
        - parametres = 'not-a-knot' : The first and second segment at a curve end are the same polynomial.
            It is a good default when there is no information on boundary conditions.
        - [obsolete, simplification, cf ci-dessous]
            parametres = ((dx0, vx0, dy0, vy0), (dxn, vxn, dyn, vyn)) permet de fixer les dérivées aux extrémités
            * le premier tuple concerne le premier point de la spline
                - dx0, vx0 : ordre de dérivation en x et valeur de la dérivée dx0-ième.
                   P.ex. si dx0, vx0 = 1, 3.14 on demande que x'(0)=3.14.
                   Usuellement, dx0 et dy0 valent 0,1 ou 2
                - dy0, vy0 : analogue pour y(0)
            * le deuxième tuple définit de la même manière les dérivées de x(t) et y(t) en t=1,
                 donc pour le dernier point de la spline
            [fin obsolete]
        - parametres = ((d0, dx0, dy0), (dn, dxn, dyn)) permet de fixer les dérivées aux extrémités
            * le premier tuple concerne le premier point de la spline
                - d0 = ordre de dérivation et
                - dx0,dy0 = le vecteur dérivée d0-ième.
                   P.ex. si d0, dx0, dy0 = 1, 3.14, 7 on demande que x'(0),y'(0) = 3.14, 7
                   Usuellement, d0 vaut 0,1 ou 2
            * le deuxième tuple définit de la même manière le vecteur dérivé-dn-ieme en t=1,
                 donc pour le dernier point de la spline
            * NB, il est possible de fixer plus finement les valeurs des dérivées aux extrémités
             (cf § obsolete ci-dessus, on peut fixer (par exemple) x'(0) et y"(0) indépendemment,
             alors qu'actuellement on fixe v'(0)=(x'(0),y'(0)) ou bien v"(0)=(x"(0),y"(0))
    # si type = 'ius' ou 'interpolatedunivariatespline', il s'agit de spline d'INTERPOLATION,
        (utilise la classe InterpolatedUnivariateSpline de scipy, fille de UnivariateSpline, avec s=0)
         parametres est alors un entier, le degré de la spline.
         Un polyligne est de type 'ius', avec degré=1.
         Les autre degrés ne sont pas forcément utilisés
         je crois que ('ius', 3) fait double emploi avec ('cubic', 'not-a-knot')
    # si type = 'us' ou 'univariatespline', c'est une spline d'AJUSTEMENT. Dans ce cas, les paramètres sont
        un dictionnaire :
        parametres = {'w':None, 'bbox':[None, None], 'k':3, 's':None, 'ext':0, 'check_finite':False}
        voir  UnivariateSpline(x, y, w=None, bbox=[None, None], k=3, s=None, ext=0, check_finite=False)
        dans la doc scipy.

    """
    type_, parametres = methode
#     debug(methode=(type_,parametres),shape=cpoints.shape)

    if len(cpoints)<2 :#moins de deux points de controle
        _knots, sx, sy = np.zeros((2,)), None, None
        return _knots, sx, sy

    if type_ == 'cubic' :
        if parametres in ('periodic', 'per', 'periodique') :
            #il faut au moins 3 points, avec P[0]==P[-1] et un des points intermediaires P[k]!=P[0]
            if len(cpoints)<3 :
                _knots, sx, sy = np.zeros((2,)), None, None
                return _knots, sx, sy
#             debug(cpoints_shape=cpoints.shape, p0=cpoints[0],pn=cpoints[-1])
#             if all(cpoints[0] == cpoints[-1]) : pass
            if np.linalg.norm(cpoints[0] - cpoints[-1])<1.0e-10 :
                cpoints[-1]=cpoints[0]
            else : #On rajoute le premier point a la fin
#                 debug('cpoints[0]-cpoints[1]=%s'%(cpoints[0]-cpoints[1]))
#                 raise ValueError('Spline periodique, %s != %s'%(cpoints[0],cpoints[-1]))
                cpoints = np.vstack((cpoints, cpoints[0]))

    N = len(cpoints)
#     debug(shape=cpoints.shape)
    T = absCurv(cpoints, normalise=True)
#     debug(abscurv_cpoints=T)
#     _knots = T
    X = cpoints[:,0]
    Y = cpoints[:,1]

    if type_ == 'cubic' :
#         debug(parametres=parametres)
        if isinstance(parametres, (str, unicode)):
            sx = CubicSpline(T, X, bc_type=parametres)
            sy = CubicSpline(T, Y, bc_type=parametres)
#             debug(sx,sy)
        else :
            (d0, vx0, vy0), (dn, vxn, vyn) = parametres
            sx = CubicSpline(T, X, bc_type=((d0, vx0),(dn,vxn)))
            sy = CubicSpline(T, Y, bc_type=((d0, vy0),(dn,vyn)))
    elif type_ == 'ius':
#         try :
            sx = InterpolatedUnivariateSpline(T, X, k=parametres)#s=la précision de l'ajustement s=0 <=> interpolation
            sy = InterpolatedUnivariateSpline(T, Y, k=parametres)
#         except Exception as msg:
#             rdebug('pas de spline (degre trop eleve ?)')
#             print unicode (msg)
#             sx = sy = None

    elif type_ == 'us' :
        weights = np.ones(N)
        W = 1000.0
        eps = 1.0e-5#NPrefs.EPS
        # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
        # le choix de s suivant implique
        # abs(xi-f(ti))<eps et
        # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
        weights[0] = weights[-1] = W
        weights /= np.sum(weights)
        s = eps/(N*W)
#         debug(len(T), parametres)
        sx = UnivariateSpline(T, X, **parametres)#s=la précision de l'ajustement s=0 <=> interpolation
        sy = UnivariateSpline(T, Y, **parametres)
#         sx = UnivariateSpline(T, X, w=weights, k=parametres, s=s)#s=la précision de l'ajustement s=0 <=> interpolation
#         sy = UnivariateSpline(T, Y, w=weights, k=parametres, s=s)
#         weights = np.ones(N)
#         W = 1000.0
#         eps = NPrefs.EPS
#         # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
#         # le choix de s suivant implique
#         # abs(xi-f(ti))<eps et
#         # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
#         weights[0] = weights[-1] = W
#         weights /= np.sum(weights)
#         s = eps/(N*W)
#
#         sx = UnivariateSpline(T, X, w=weights, k=parametres, s=s)#s=la précision de l'ajustement s=0 <=> interpolation
#         sy = UnivariateSpline(T, Y, w=weights, k=parametres, s=s)
    return T, sx, sy

class NSplineAbstract(object):
#NSplineAbstract a besoin d'heriter de QObject pour pouvoir emettre un SIGNAL
    u"""
    TODO :
    - l'echantillonnage ne doit plus être fixé à l'__init__ (suppression de nbpe, mode)
    - self.epoints ne doit plus exister : il devrait etre calculé à la demande
        epoints = self.echantillonner(nbp, mode, ta, tb) avec tous les parametres obligatoires
        les epoints n'ont pas à être mis à jour, car dans le l'interface graphique,
        faire suivre les points échantillonnés devient trop lourd
    FIN_TODO
    NSplineAbstract est une interface commune à toutes les splines
    (simple, composées, profil...)
    Ne peut pas être instancié, car elle a des méthodes virtuelles pures.
    ========================================================================
    = La spline est auto-suffisante, le _update() ne doit pas être appelé  =
    = par d'autres objets. Sauf éventuellement dans les classes filles.    =
    ========================================================================
    Une NSpline est constitué de
    - self.cpoint (property) les points de contrôle sous forme np.ndarray((n,2))
        C'est lui qui fait référence => et c'est une connerie. Les calculs se font avec points
    - self.cpoints(points) = le setter de self.cpoints, appelle self._update()
    - self.sx, self.sy : une spline parametrique d'interpolation calculées par scipy
    - self.dpoints (property) les points de la spline discrétisée pour visualisation,
        ce tableau contient self.precision points (e.g. self.precision est grand:1000)
    - self.epoints (property) : les points echantillonnés de la spline, suivant
        le mode de répartition défini par self.mode et self.nbpe.
        S'ils n'existent pas, la property epoints fait appel à self.echantillonner()
    - self.mode : mode d'echantillonnage
    - self.nbpe : nb points d'echantillonnage.
    - self.name : le nom de la spline
    - self.role : chaine de caractères, le rôle.
    - self.qcpolygon de type QPolygonF () : les N points de contrôle.
        => ATTENTION, à chaque appel l'intégralité du self.cpoints est transformé en QPolygonF.
    - self.qdpolygon de type QPolygonF () : self.dpoints transformé en QPolygonF
        => ATTENTION, à chaque appel l'intégralité du self.dpoints est transformé en QPolygonF.
    Méthodes :
    --------
    - self.abscurv() : recalcule les abscisses curvilignes des cpoints, stocké dans _knots.
        NON => Déclenche un recalcul des splines sx et sy <= NON
        Ces abscisses sont les paramètres sx.x et sy.x des deux splines
        en principe ce sont les trois mêmes tableaux : _knots==sx.x==sy.x
    - self.computeSpline() calcul de la spline self.sx, self.sy.
        Normalement c'est fait automatiquement dès que cpoints est modifié
    - self.save(filename) sauvegarde texte simple, un point (=deux coordonnées) par ligne.
    - self.echantillonner() : retourne la spline echantillonnée et renseigne self.epoints.
    - self.isClosed(eps=1.0e-8) return True si self[0] == self[-1] à eps près
    - self.load(dump) : permet de loader une spline sauvegardée sous la forme retournée par self.toDump()
    - self.toDump() retourne un dictionnaire en clair permettant la reconstruction de la spline
    - self.__call__(T) équivaut à P = self(T)
        où T est un tableau d'abscisses curvilignes quelconques entre 0 et 1
        retourne les points (sx(ti), sy(ti) pour ti parcourant T.
    - self.plot() pour tracer avec matplotlib
    - self.__getitem__(i), équivaut à p = self[i]
        i entier, retourne le i-eme point de contrôle.
    - self.__setitem__(i, p), équivaut à self[i] = p,
        i entier, p de type (x,y) ou QPointF : affecte p au i-eme point de contrôle.
        puis déclenche un self._update()
    - len(self) = le nb de points de contrôle.
    Les méthodes suivantes déclenchent le _update() i.e. recalcul complet de la spline (sx, sy)
    et mise à jour des dpoints,...
    - self.[insert,append]Point(p) comme leur nom l'indique, lèvent une exception en cas de point double
    - self.removePoint(p) comme son nom l'indique
    - self.hardScale, self.hardRotate, self.translate : comme leur nom l'indique.
    """
    # defaultprefs = SplinePrefs()
    SUPPORTED_FILES = ('*.gnu', '*.dxf', '*.pts', '*.pkl','*.spl','*.npkl')
    def __init__(self, **dump):
        try : parent = dump.pop('parent')
        except KeyError : parent = None
        super(NSplineAbstract, self).__init__()
        self.nbupdate = 0
        self.gparent = parent
        self.setDefaultValues()
        self.load(dump)
#         self._update()Doit etre appelé explicitement par les héritiers au bon moment
    def __call__(self, T, diff=0):
        u"""
        # Paramètre :
            * T peut être
                - un réel entre 0 et 1
                - ou un np.ndarray de n réels dans [0,1] (de shape (n,1)? ou (1,n)?)
                - ou une liste de n réels entre 0 et 1
            * diff est l'ordre de dérivation
        # Retourne :
        les 2 ndarrays de shape (n,1), zippés, contenant les points (sx(ti, diff), sy(ti, diff)) ou ti parcourt T
            * si diff=0, retourne les coordonnées du point x(t), y(t),
            * si diff=1, retourne les dérivées x'(t), y'(t)
            * si diff=2, retourne les dérivées secondes x"(t),y"(t)
        # Utilisation :
            >>> S = NSpline(....)
            >>> T = np.linspace(0.0, 1.0, 11)
            >>> S(T)
            ... : les 11 points [S(0), S(0.1), ...,S(0.9), S(1)] de la spline
            N.B. : la fonction zip fait ceci:
            >>> z = zip(range(3),range(3))
            >>> z
            ... [(0, 0), (1, 1), (2, 2)]
        """
#         #debug(self.sx(T))
#         #debug(self.sy(T))
        try :
            X, Y = self.sx(T, diff), self.sy(T, diff)
#             debug(X)
            try :
                return np.asarray(zip(X, Y))
            except TypeError :
                return [X,Y]
        except TypeError as msg :
            return
            debug(msg)
#         return np.hstack((self.sx(T), self.sy(T)))

    u"""methodes virtuelles pures"""
################################################################
    def setDefaultValues(self):
        u"""Valeurs par defaut:"""
        raise NotImplementedError(u"la methode() %s doit etre implemente"%(whoami(self)[:-2]))

    def computeSpline(self, *args, **kargs):
        u"""renseigne self.sx, self.sy"""
        raise NotImplementedError(u"la methode() %s doit etre implemente"%(whoami(self)[:-2]))

    def echantillonner(self, *args, **kargs):
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    @property
    def dpoints(self):
        raise NotImplementedError(u"la property %s doit etre implemente"%(whoami(self)[:-2]))

    @property
    def cpoints(self):
        u"""points de contrôle, sous forme np.ndarray((n,2))"""
        raise NotImplementedError(u"la property %s doit etre implemente"%(whoami(self)[:-2]))

    @cpoints.setter
    def cpoints(self, points):
        u"""points de contrôle, sous forme np.ndarray((n,2))"""
        raise NotImplementedError(u"la property %s doit etre implemente"%(whoami(self)[:-2]))

    @property
    def qcpolygon(self):
        raise NotImplementedError(u"la property %s doit etre implemente"%(whoami(self)[:-2]))

    @property
    def qepolygon(self):
        raise NotImplementedError(u"la property %s doit etre implemente"%(whoami(self)[:-2]))
#
    def hardRotate(self, alfa, centre=None):
        u'''Rotation in-situ, self est modifié'''
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def hardScale(self, scale, centre=None):
        u'''
        Modif de la spline in-situ, i.e.self est modifié.
        Mise à l'échelle scale, centrée sur centre. (en fait, une homothétie)
        C'est à dire
        >>> self[i] = centre + scale*(self[i]-centre)
        si centre == None on prend l'isobarycentre
        '''
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def symetriser(self, axe, centre=None):
        u'''
        Modif de la spline elle meme. self est modifié.
        C'est à dire
        >>> self[i,axe] = 2*centre -self[i,axe]
        si centre == None on prend l'isobarycentre[axe]
        '''
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def translate(self,vecteur):
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def appendPoint(self,pos):
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))
#
    def insertPoint(self, pos,k=None):
        u"""Calcule dans quel segment il est raisonable d'insérer le point si k=None
        ou bien insère le point en position k"""
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def removePoint(self, pnum):
        u"""Suppression du point pnum de cpoints"""
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def copy(self):
        u"""retourne une copie de self"""
        dump = self.toDump()
        return type(self)(**dump)
#         raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def distance2To(self, p):
        u"""
        Calcule et retourne la distance de la spline self à un point p
        :param p: tuple, liste ou QPointF le point
        :return res le résultat, de type scipy.optimize.OptimizeResult
        voir https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html
            en particulier :

        - res.nfev : int, nb evaluation de
            phi(t)  = (a-x(t))**2 + (b-y(t))**2
                    = (norme(self(t)-p))**2
        - res.x : float, valeur de t réalisant cette distance
        - res.fun : float, est le carré de la distance euclidienne trouvée
            i.e. phi(t)
        """
        res = distance2PointSpline(p, self, t0=0.5, tol=1.0e-9)
        return res
################################################################
    @property
    def msg(self):
        return u"%s doit etre implementer la methode %s"%(self.classname, whoami(self))
    @property
    def classname(self):
        return self.__class__.__name__
#     @property
#     def methode(self):
#         return self._methode
#     @methode.setter
#     def methode(self, newmethode):
#         u"""Interdiction de changer de méthode (cubic, us, ius),
#         On peut changer les paramètres de la méthode (clamped, not-a-knot,...)
#         Pour changer de méthode, quand même, on peut faire :
#         >>> S = SplineSimple(points=..., methode=meth1)
#         >>> dump = S.toDump()
#         >>> dump['methode'] = newmethode
#         >>> S.load(dump)"""
#         if newmethode[0] == self.methode[0] :
#             self._methode = newmethode
#             self._update()
#         else :
#             raise RuntimeError("Interdiction de changer de méthode")
    @property
    def precision(self):
        return self._precision
    @precision.setter
    def precision(self, prec):
        self._precision = prec
        try : self._dpoints = self(np.linspace(0,1,prec))
        except (AttributeError, TypeError) :
#             rdebug(u'ouuuups!')#, len=len(self._cpoints))
            pass
#     @property
#     def cpoints(self):
#         u"""points de contrôle, sous forme np.ndarray((n,2)),
#         [re-assemblé à chaque fois => NON]"""
#         try :
#             return self._cpoints
#         except AttributeError :
#             self._cpoints = pointsFrom(self.qcpolygon)
#             return self._cpoints
    @property
    def epoints(self):
        u"""points échantillonnés, sous forme np.ndarray((n,2))"""
        try :
            return self._epoints
        except AttributeError :
#             debug('*****')
#             rdebug("*** (%s) Pas de _epoints, il faut echantillonner"%(self.name))#,self.nbpe,self.mode))
#             raise ValueError
#             rstack()
            if len(self)<=1 :#pas encore de cpoints, ou une longueur nulle
                self._epoints = np.zeros((1,2))
                return self._epoints
            try :
                self.echantillonner()#nbp=self.nbpe, mode=self.mode)
            except TypeError as msg :
                rdebug(u'TypeError, Échantillonnage impossible', msg)
                self._epoints = np.zeros((1,2))
            except ValueError as msg :
                rdebug(u'ValueError', msg)
                self._epoints = np.zeros((1,2))
#             rdebug("***fin echantillonner")
            return self.epoints#recursif

    # @property
    # def qcontrol(self):#alias
    #     u"""le QPolygonF des points de contrôle"""
    #     return self.qcpolygon
    @property
    def longueur(self):
        u"""Longueur vraie de la spline défini par self.sx, self.sy"""
        try :
            return self._longueur
        except AttributeError :
            self._longueur = self.abscurv(1)#absCurvReal(self,1.0)[0]
        return self._longueur

    @property
    def clongueur(self):
        u"""Longueur du polyligne défini par les points de contrôle"""
        return absCurv(self.cpoints)[-1]
    @property
    def dlongueur(self):
        u"""Longueur du polyligne défini par les points de discrétistion"""
        return absCurv(self.dpoints)[-1]
    @property
    def elongueur(self):
        u"""Longueur du polyligne défini par les points échantillonés"""
        return absCurv(self.epoints)[-1]
    @property
    def height(self):
        try :
            return self._height
        except AttributeError :
            try :
                self._height = max(self.dpoints[:,1]) - min(self.dpoints[:,1])
                return self._height
            except (TypeError, AttributeError, ValueError) :
                return 0.0
    @property
    def width(self):
        try :
            return self._width
        except AttributeError :
            try :
                self._width = max(self.dpoints[:,0]) - min(self.dpoints[:,0])
                return self._width
            except (TypeError, AttributeError, ValueError) :
                return 0.0

    def __str__(self):
        return '\n ==> '.join(self.info)

    __repr__ = __str__
    
    @property
    def methode(self):
        return self._methode

    @property
    def info(self):
        infos=[
                u"<%s>"%self.classname,
                u'%20s = '%u'name'+u'%s'%self.name,
                u'%20s = '%u'role'+u'%s'%self.role,
                u'%20s = '%u'closed'+u'%s'%self.isClosed(),
                u'%20s = '%u'nb pts controle'+u"%d"%len(self.cpoints),
                u'%20s = '%u'methode'+u"%s"%str(self.methode),
#                 u'%20s = '%'precision affichage'+"%s"%str(self.precision),
                u'%20s = '%u'mode echantillonage'+u"%s"%self.mode,
                u'%20s = '%u'nb pts echantillon'+u"%s"%self.nbpe,
                u'%20s = '%u'nb updates'+u"%s"%self.nbupdate,
                ]
        i = u'largeur, hauteur'
        try :
            i1 = u"%g, %g"%(self.width, self.height)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(self.classname, msg)
        infos.append(u'%20s = '%i+i1)

        i = u'position cg'
        try :
            i1 = u"%g, %g"%(self.centre[0], self.centre[1])
        except Exception as msg :
            i1 = u"? (%s, %s)"%(self.classname, msg)
        infos.append(u'%20s = '%i+i1)

        i = u'longueur'
        try :
            i1 = u"%g"%(self.dlongueur)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(self.classname, msg)
        infos.append(u'%20s = '%i+i1)

        return infos

#     def fromFile(self, filename):
#         u"""Destructif, charge une spline depuis un fichier
#         >>> spl = NSplineSimple()
#         >>> spl.fromFile('toto.spl')"""
#         with open(filename,'r') as f :
#             dump = eval(f.read())
#         self.load(dump)

    def load(self, dump):
        u"""
        :param dump: un dictionnaire
        Ce qui est commun à toutes les splines"""
        try :self.name = dump.pop('name')
        except KeyError : pass

        try :self.gparent = dump.pop('gparent')
        except KeyError : pass

        try : self.role = dump.pop('role')
        except KeyError : pass
        try :
            self._methode = dump.pop('methode')
        except KeyError :
            pass
        try :
            self.precision = dump.pop('precision')
        except KeyError :
            pass
        try :
            self.nbpe = dump.pop('nbpe')
        except KeyError :
            pass
        try :
            self.mode = dump.pop('mode')
        except KeyError :
            pass

    def toDump(self, format_='new'):
        u"""Ce qui est commun à toutes les splines"""
#         raise NotImplementedError(u"la methode %s() doit etre implemente"%(whoami(self)[:-2]))
        if format_=='new':
            return {
                    'classename' : self.classname,#besoin de savoir quel type de spline.
                    'cpoints'    : self.cpoints.tolist(),
                    'role'       : self.role,
                    'name'       : self.name,
                    'methode'    : self.methode,
                    'precision'  : self.precision,
                    'mode'       : self.mode,
                    'nbpe'       : self.nbpe
                    }
        elif format_=='old' :
            #pour compatibilité ascendante
            return {
                    'points'     : self.epoints.tolist(),
                    'name'       : self.name,
                    }

    def save(self, filename):
        filename = Path(filename)
        ext = filename.ext
#         debug(filename=filename, ext=ext)
        try :
            if ext in (".gnu", '.txt'):
                #seulement les points échantillonnés
                np.savetxt(filename, self.epoints)
            # elif ext==".pts":
            #     #seulement les points échantillonnés
            #     writeProfil(filename, self.epoints)
            # elif ext=='.npkl':
                #Spline complète, pickle
                cPickle.dump(self.toDump('new'),open(filename,'w'))
            elif ext=='.pkl':
                #seulement les points échantillonnés
                cPickle.dump(self.toDump('old'),open(filename,'w'))
            elif ext=='.spl':
                #Spline complète, dict
                with open(filename,'w') as f :
                    f.write(str(self.toDump()))
            elif ext=='.dxf':
                raise NotImplementedError('Format dxf')
            debug('Sauvegarde %s : OK'%filename)
        except Exception as msg:
            rdebug(msg, 'Sauvegarde %s : impossible'%filename)

    def open(self, filename):
        filename = Path(filename)
        ext = filename.ext
        #debug(filename=filename)
        if ext in (".gnu", '.txt'):
            dump = self.toDump()
            dump['cpoints'] = np.loadtxt(filename)
        elif ext==".pts":
            dump = self.toDump()
            dump['cpoints'] = LecteurUniversel(filename).points
        elif ext=='.pkl':
            dump = cPickle.load(open(filename,'r'))
            for key in ('points', 'cpoints') :
                if dump.has_key(key) :
                    dump[key] = pointsFrom(dump.pop(key))
#         #debug(dump=dump)
        self.load(dump)

    def close_(self) :
        u'''Si le qcpolygon est EXACTEMENT fermé (eps=0.0), on ne rajoute pas de point.'''
#        self.isclosed = True
        if self.isClosed(0.0) :
            return False
        else :
            self.appendPoint(self[0])
            return True

    def isClosed(self,eps=1.0e-8):
        u'''self[0] == self[-1] a eps près'''
        cp = self.cpoints
        if len(cp)>0 : return len(self)>1 and dist2(cp[0],cp[-1])<=eps*eps# np.all(self.cpoints[0] == points[-1]):

    def __len__(self):
        return len(self.cpoints)
#         return len(self.epoints)

#     def __call__(self,k):
#         u"""
#         Pour accéder aux points du polygone echantillonné (en lecture)
#         ne traite pas les slices
#         Retourne un tuple (x,y)
#         >>> S = NSplineSimple(points=....)
#         >>> S[3]
#         [10.0,3.14] le point points[3]
#         ne gère pas les slices
#         >>> S[1:3]
#         ... AttributeError: 'QPolygonF' object has no attribute 'x'
#         """
#         return self.epoints[k]

    def __getitem__(self,k):
        u"""
        Pour accéder aux points du polygone de controle (en lecture)
        ne traite pas les slices
        Retourne un tuple (x,y)
        >>> S = NSplineSimple(points=....)
        >>> S[3]
        [10.0,3.14] le point points[3]
        ne gère pas les slices
        >>> S[1:3]
        ... AttributeError: 'QPolygonF' object has no attribute 'x'
        """
        return self.cpoints[k]
#         return self.qcpolygon[k].x(),self.qcpolygon[k].y()

    def __setitem__(self,k,value):
        u"""
        mise à jour du polygone de contrôle et du polygone seulement
        Typiquement, on pourra ecrire
        >>> S = NSplineSimple(cpoints=[(1,1),(2,2),(3,3)])
        >>> S[1] = 20,30
        >>> print S[1]
        [20.0, 30.0]
        """
        if np.all(self._cpoints[k]==value) : #inutile de déclencher un _update()
            return
        self.cpoints[k] = value
        self._update()#indispensable pour mettre à jour self.cpoints et dpoints

#
#
#     @property
#     def width(self):
#         return self.qcpolygon.boundingRect().width()
#
#     @property
#     def height(self):
#         return self.qcpolygon.boundingRect().height()

#     def eabscurv(self):
#         u"""Les abscisses curvilignes normalisées des points échantillonnés,
#         entre 0 et 1.
#             recalculées si besoin. Ce sont aussi les paramètres t de la spline sx et sy
#         """
#         if hasattr(self, '_eabscurv') and hasattr(self, 'sx'): #normalement s'il y a l'un il y a l'autre
#             return self._eabscurv
#         else :
#             self._eabscurv = absCurv(self.epoints, normalise=True)
#             return self._eabscurv


    def abscurv(self, T=None):
#         raise NotImplementedError(self.msg)
        u"""
        :return :

        - si T != None les abscisses curvilignes NON NORMALISÉES  des points S(T)
        - si T == None, les abscisses curvilignes NORMALISÉES (entre 0 et 1) des points de contrôle,
            recalculées si besoin. Ce sont aussi les knots sx.x et sy.y des deux splines sx et sy
        :TODO: remplacer les appels à self.abscurv() par self.knots et
            supprimer la valeur par defaut T=None
        """
        if T is not None :
            return absCurvReal(self, T)[0]#plus cher, plus précis ?
            ##############################################################
            #return absCurv(self(T), normalise=False)
            #Faux la vraie abscisse curviligne est absCurvReal(self, T)[0]
            ##############################################################
        else :
            ##############################################################
            # Ci dessous, _knots est exactement self.knot.
            # Inutile de le recalculer
            ##############################################################
            if hasattr(self, '_knots') and hasattr(self, 'sx'):
                #normalement s'il y a l'un il y a l'autre
                return self._knots
            else :
                self._knots = absCurv(self.cpoints, normalise=True)
                return self._knots
    @property
    def knots(self):
        u"""retourne les knots de la spline, c'est a dire les T tesls que S(T) = self.cpoints"""
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

    def courbure(self, T) :
        dx,  dy  = self.sx(T, 1), self.sy(T, 1)
        d2x, d2y = self.sx(T, 2), self.sy(T, 2)
        norm3_d2 = np.sqrt(dx**2+dy**2)**3
        sc = (dx*d2y-dy*d2x)/(norm3_d2)
        # si norm_d2=0, x"(t)=y"(t)=0, c'est une droite, courbure nulle
#         sc[np.where(norm3_d2 < 1.0e-12)] = 0.0
        return sc

    def integraleCourbure(self, a=0, b=1, n=100):
        """Integrale de la valeur absolue de la courbure, caractérise la régularité de la courbe"""
# def simpson(f, a, b, n=10):#n doit être pair, integration precise ordre 3
#     u"""Integrale de f sur [a,b], méthode de Simpson composite. (ordre 3)
#     n DOIT être pair"""
        h = float(b-a)/n
        T = np.linspace(a, b, n+1)
        C = abs(self.courbure(T))
        A1 = C[0] + C[-1]
        A2 = 2*sum(C[i] for i in range(2,n) if i%2==0)
        A4 = 4*sum(C[i] for i in range(1,n) if i%2==1)
    #         debug (h, A1, A2, A4, (h/3)*(A1 + A2 + A4))
        return (h/3)*(A1 + A2 + A4)



    @property
    def center(self):
        '''isoBayrcentre'''
        if self.isClosed() :
            return baryCentre(self.cpoints[:-1]).flatten()
        else :
            return baryCentre(self.cpoints).flatten()
    centre=center

    @property
    def role(self):
        return self._role
    @role.setter
    def role(self,role):
        self._role=role

    @property
    def shape(self):
        u"""La shape au sens NumPy"""
        return self.cpoints.shape

#     @property
#     def qdpolygon(self):
#         u"""Les points discrétisés de la spline"""
# #         #debug ()
#         return qpolygonFrom(self.dpoints)

    def boundingRect(self):
        dpts = self.dpoints
        xM, xm, yM, ym = dpts[:,0].max(), dpts[:,0].min(), dpts[:,1].max(), dpts[:,1].min()
        return xm,ym,xM,yM

    # def controlBoundingRect(self):
    #     return self.qcpolygon.boundingRect()

    @property
    def gravitycenter(self):
        u"""Le centre de gravité de la surface délimitée par le polygone fermé."""
        return centreGravite(self.dpoints)

    @property
    def barycentre(self):
        u"""Le centre de gravité de la surface délimitée par le polygone fermé."""
        return baryCentre(self.cpoints)

    def _update(self):
        u'''
        Est appelé à chaque modification
        - (géométrique) d'un point de contrôle de la spline
        - ou bien du PolygonF de base
        - ou bien de methode spline (cubic, IUS, US,...), ou des dérivées aux extremites
        - ou bien de mode d'échantillonage
        ultra privée, ne pas l'appeler de l'extérieur.
        Supprime et recalcule tout, en particulier les splines sx et sy
        '''
        self.nbupdate +=1
#         debug(nbupdate=self.nbupdate)
        try : del self._qcpolygon
        except AttributeError : pass
        try : del self._qepolygon
        except AttributeError : pass
        try : del self._epoints
        except AttributeError : pass
        try : del self._dpoints
        except AttributeError : pass
        try : del self._knots
        except AttributeError : pass
        try : del self._height
        except AttributeError : pass
        try : del self._width
        except AttributeError : pass
        try : del self._longueur
        except AttributeError : pass

        try : #ne pas supprimer, on teste si sx existe dans abscurv
            del self.sx,
            del self.sy
        except AttributeError :
            pass
        try :
#             #debug(methode=self.methode)
            self.computeSpline(self.methode)
        except TypeError as msg :
#             rdebug(self.name, msg)
            raise(msg)
        except AttributeError as msg :
#             rdebug(self.name, msg)
            raise(msg)

    def plot(self, plt, control=True, nbpd=None, nbpe=None, mode=None,
             titre=None, more=[], show=True, buttons=True):
        """
        :param plt : une instance de pyplot, obtenue en amont par :
              >>> from matplotlib import pyplot as plt
        :type plt : cf ci-dessus
        :param control : affichage des points de controle True/False
        :type control : bool
        :param nbpd : nb points de discretisation
        :type nbpd : int
        :param nbpe : nb de points d'echantillonage
        :type nbpe : int
        :param titre : titre
        :type plt : str ou unicode
        :param more : [(X,Y, couleur,'nom'), ...] tracé supplementaire
        :type more : list
          """
        from matplotlib.widgets import CheckButtons
#         plt.figure(numfig)
#         rdebug('***********')
        if nbpd is None : nbpd = self.precision
        if nbpe is None : nbpe = self.nbpe
        if mode is None : mode = self.mode
#         debug('appel echantillonnage', type(self))

        D = self.dpoints
        C = self.cpoints
#         E = self.epoints
#         if len(E)==0 :
#         rdebug()
        E = self.echantillonner(nbp=nbpe, mode=mode)
#         debug(E=E)
        _, ax = plt.subplots()
        if titre is None : titre = self.name+str(self.methode)
        plt.title(titre)
        spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
        echantillon, = ax.plot(E[:,0], E[:,1], 'g.', lw=1)
        control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        if buttons :
            butt = ['control','echantillon', 'spline',]
            values = [True, True, True]
            draws = [spline, control, echantillon]
        for x, y, color, name in more:
            temp, = ax.plot(x, y, color)
            if not name : continue
            if buttons :
                draws.append(temp)
                butt.append(name)
                values.append(True)
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        if buttons :
            rax = plt.axes([0.05, 0.4, 0.1, 0.15])
            check = CheckButtons(rax, butt, values)

            def func(label):
                if label == 'spline':
                    spline.set_visible(not spline.get_visible())
                elif label == 'control':
                    control.set_visible(not control.get_visible())
                elif label == 'echantillon':
                    echantillon.set_visible(not echantillon.get_visible())
                else :
                    draw = draws[butt.index(label)]
                    draw.set_visible(not draw.get_visible())
                plt.draw()
            check.on_clicked(func)
        if show : plt.show()
        return plt

if __name__=="__main__":
    debug('Rien a faire')

