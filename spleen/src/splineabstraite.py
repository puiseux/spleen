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
from numbers import Number
import numpy as np
from numpy import asarray as array, linspace

from scipy.optimize import newton, minimize, minimize_scalar
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline, UnivariateSpline
import cPickle
from utilitaires.utilitaires import (Path, whoami, debug, rdebug,stack,rstack,
                        absCurv,dist2, baryCentre, centreGravite, simpson)
from utilitaires.lecteurs import pointsFrom
from utilitaires.lecteurs import LecteurUniversel
# from numpy.core.defchararray import center
from numpy import sqrt
from scipy.integrate import quad
from scipy.interpolate.fitpack2 import LSQUnivariateSpline

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
        res = array([absCurvReal(S, t) for t in T])
#         res = asarray([quad(lambda s : sqrt(S.sx(s,1)**2+S.sy(s,1)**2), 0.0, t) for t in T])
        return res[:,0],  res[:,1]

def distancePointSpline(p, S, t0=0.5, tol=1.0e-9):
    u"""
    Comme son nom l'indique, calcule le carré de la plus courte distance euclidienne 
    d'un point p à une spline paramétrée S(t) = x(t), y(t), 0<=t<=1
    :param S: NSplineAbstract, la spline paramétrique
    :param p: (x,y)=(float, float) le point.
    :param t0: float, initialisation des itérations
    :param tol: float, tolérance passée à scipy.optimize.minimize_scalar
    :return: le resultat res,
    voir https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
    en particulier :
        - res.x : float, valeur de t réalisant cette distance
        - res.nfev : int, nb evaluation de phi(t)
        - res.fun : float, valeur finale de phi(t)

    local : la fonction phi(t) à optimiser est le carré de la distance de p à S(t)
    """
    a, b = p[0], p[1]
    x, y = S.sx, S.sy
    phi = lambda t: (a-x(t))**2 + (b-y(t))**2
    res = minimize_scalar(phi,
                          bounds=(0.0, 1.0),
                          method='bounded',
                          options={'xatol':1.0e-9})
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
        #TODO les knots devraient être de shape (0,)??
        _knots, sx, sy = np.zeros((2,)), None, None
        return _knots, sx, sy

    if type_ == 'cubic' :
        if parametres in ('periodic', 'per', 'periodique') :
            #il faut au moins 3 points, avec P[0]==P[-1]
            #et un des points intermediaires P[k]!=P[0]
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
        #UnivariateSpline(x, y, w=None, bbox=[None, None], k=3,
        #                 s=None, ext=0,check_finite=False)
        weights = np.ones(N)
        W = 1000.0
        # eps = 1.0e-5#NPrefs.EPS
        # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
        # le choix de s suivant implique
        # abs(xi-f(ti))<eps et
        # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
        weights[0] = weights[-1] = W
        weights /= np.sum(weights)

        # s = eps/(N*W)
#         debug(len(T), parametres)
        sx = UnivariateSpline(T, X, k=parametres, w=None, s=1.0e-10)#s=la précision de l'ajustement s=0 <=> interpolation
        sy = UnivariateSpline(T, Y, k=parametres, w=None, s=1.0e-10)
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
    elif type_ == 'lsqus' :
        raise NotImplementedError ('LSQUnivariateSpline non implemente')
        #LSQUnivariateSpline(x, y, t, w=None, bbox=[None, None], k=3, ext=0,
        #                       check_finite=False)
        weights = np.ones(N)
        W = 1000.0
        # eps = 1.0e-5#NPrefs.EPS
        # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
        # le choix de s suivant implique
        # abs(xi-f(ti))<eps et
        # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
        weights[0] = weights[-1] = W
        weights /= np.sum(weights)
        knots = linspace(0,1,20)
        sx = LSQUnivariateSpline(T, X, knots, k=parametres, w=None)#s=la précision de l'ajustement s=0 <=> interpolation
        sy = LSQUnivariateSpline(T, Y, knots, k=parametres, w=None)
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
        C'est lui qui fait référence
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
    - self.isClosed(eps=1.0e-8) return True si dist(self[0],self[-1])< eps (distance norme 1)
    - self.load(dump) : permet de loader une spline sauvegardée sous la forme retournée par self.toDump()
    - self.toDump() retourne un dictionnaire en clair permettant la reconstruction de la spline
    - self.__call__(T) équivaut à P = self(T)
        où T est un tableau d'abscisses curvilignes quelconques entre 0 et 1
        retourne les points (sx(ti), sy(ti) pour ti parcourant T.
    - self.plot() pour tracer avec matplotlib
    - self.__getitem__(i), équivaut à p = self[i]
        i entier, retourne le i-eme point de contrôle.
    - self.__setitem__(i, p), équivaut à self[i] = p,
        i entier, p de type tuple ou liste (x,y) : affecte p au i-eme point de contrôle.
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
        super(NSplineAbstract, self).__init__()
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
                return array(zip(X, Y))
            except TypeError :
                return [X,Y]
        except TypeError :
            return

    u"""methodes virtuelles pures"""
################################################################

    def setDefaultValues(self):
        u"""Valeurs par defaut:"""
        raise NotImplementedError(u"la methode() %s doit etre implemente"%(whoami(self)[:-2]))

    def load(self):
        u"""Valeurs par defaut:"""
        raise NotImplementedError(u"la methode() %s doit etre implemente"%(whoami(self)[:-2]))

    def computeSpline(self, *args, **kargs):
        u"""renseigne self.sx, self.sy"""
        raise NotImplementedError(u"la methode() %s doit etre implemente"%(whoami(self)[:-2]))

    def echantillonner(self, *args, **kargs):
        raise NotImplementedError(u"la methode %s doit etre implemente"%(whoami(self)[:-2]))

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
#         self.nbupdate +=1
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

    def plot(self, plt, control=True, nbpd=None, nbpe=None, mode=None,
             titre=None, more=[], show=True, buttons=True):
        """
        :param plt : une instance de pyplot, obtenue en amont par :
              >>> from matplotlib import pyplot as plt
        :param control : bool, affichage des points de controle True/False
        :param nbpd : int, nb points de discretisation
        :param nbpe : int, nb de points d'echantillonage
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
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

