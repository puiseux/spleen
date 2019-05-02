#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Created on 11 mai 2012

@author: puiseux
__updated__ = '2019-02-11'
'''
import datetime
import gc  # garbage colector
import math
import pickle
import string
import sys
# import numpy
# from lecteurdxf import LecteurDXFNervures, LecteurDXF
# from lecteur import LecteurGnuplot
# from lecteurs import pointsFromFile
# from math import acos, atan, cos, sin, tan
from pprint import pprint
# from random import choice, random
from numpy.random import random#, choice
# from lecteurs import pointFromFile
# import numpy as np
# import scipy as sp
from numpy import (asarray, zeros, vstack, ndarray, matrix, absolute, copy,
    nan, sqrt, where, sin, cos, logical_and, dot, argmin, arange, delete, inf, linspace,pi,
    radians)
# from numpy.matlib import rand
from path import Path
from scipy import interpolate, ones, prod, sum
from numpy.linalg.linalg import norm
# print sp.version.version
from scipy.interpolate import (CubicSpline, InterpolatedUnivariateSpline,
                               LSQUnivariateSpline, Rbf, UnivariateSpline,
                               interp1d, splev, splrep)
from scipy.optimize import minimize_scalar
from shapely.geometry import LinearRing, Point, Polygon
from matplotlib import pylab
from collections import OrderedDict
import spleenconfig
# import logging
# logger = logging.getLogger('')
# logger.setLevel(logging.DEBUG)

def souligne(msg, c=u'-', shift=0):
    sl = shift*u' '+len(msg)*c
    return shift*u' '+msg+u'\n'+sl

def sursouligne(msg, c=u'-', shift=0):
    sl = shift*u' '+len(msg)*c

    return sl+u'\n'+shift*u' ' + msg+u'\n'+sl

def splineInterpolation(points, methode='c cubic', tension=5, degre=3):
    u"""
    Une paire de spline cubiques paramétriques qui interpole ou ajuste le polygone points
    sx(t), sy(t) sont deux splines.
    Voir la doc Scipy.interpolate...
    - si methode dans {'x cubic', 'ius','periodic','interp1d',} c'est de l'interpolation
    - si methode est {'us','univariatespline'} c'est de l'ajustement, le poids est 1 pour tous les points
    Retourne:
    --------
    T = les valeurs du paramètre t, abscisse curviligne NORMALISEE, entre 0 et 1.
    sx, sy = les deux splines
    """
    if methode in ('periodic', 'p cubic', ) :
        if all(points[0] == points[-1]) : pass
        else : #On rajoute le premier point a la fin
            points = vstack((points, points[0]))
#             points.resize((1+len(points),2)) #moins cher que vstack mais marche pas !!
#             points[-1] = points[0]
    if tension == 0.0 :
        eps = 0.0
    else :
        eps = 10.0**(-tension)

    N = len(points)
    T = absCurv(points)
    if len(points)<2 : return T, None, None
    T /= T[-1]
    X = points[:,0]
    Y = points[:,1]
    try : methode = methode.lower()
    except AttributeError : pass
#     debug(None, methode=methode, tension=tension, degre=degre)
    if methode in ('ius','interpolatedunivariatespline') :
        try :
            sx = InterpolatedUnivariateSpline(T, X, k=degre)#s=la précision de l'ajustement s=0 <=> interpolation
            sy = InterpolatedUnivariateSpline(T, Y, k=degre)
        except Exception as msg:
            debug(None)
            print unicode (msg)
#             debug(None,u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
            sx = sy = None
    elif methode in ('us','univariatespline') :
        try :
            weights = ones(N)
            W = 1000.0
            # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
            # le choix de s suivant implique
            # abs(xi-f(ti))<eps et
            # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
#             eps = 10.0**(-tension)
            weights[0] = weights[-1] = W
            weights /= sum(weights)
            s = eps/(N*W)
            sx = UnivariateSpline(T, X, w=weights, k=degre, s=s)#s=la précision de l'ajustement s=0 <=> interpolation
            sy = UnivariateSpline(T, Y, w=weights, k=degre, s=s)
        except Exception as msg:
            debug(None)
            print unicode(msg)
#             debug(None,u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
            sx = sy = None
    elif methode in ('interp1d',) :
        try :
            sx = interp1d(T, X, kind=degre)
            sy = interp1d(T, Y, kind=degre)
        except ValueError as msg:
            debug(None)
            print unicode(msg)
            sx = sy = None
#     elif methode in ('periodic',) :
#         try :
#             sx = PeriodicSpline(T, X, k=degre, s=eps)
#             sy = PeriodicSpline(T, Y, k=degre, s=eps)
#         except ValueError as msg:
#             debug(None)
#             print unicode(msg)
#             sx = sy = None
    elif 'cubic' in methode :#or isinstance(methode, (tuple, list, np.ndarray)):
        if methode == 'p cubic' : bc_type='periodic'
        elif methode == 'c cubic' : bc_type='clamped'
        elif methode == 'n cubic' : bc_type='natural'
        else : bc_type = 'not-a-knot'

        try :
#             debug(None, T)
            sx = CubicSpline(T, X, bc_type=bc_type)
            sy = CubicSpline(T, Y, bc_type=bc_type)
        except ValueError as msg:
            print unicode(msg)
            sx = sy = None

    elif isinstance(methode, (tuple, list, ndarray)):
        bc_type = methode
        try :
#             debug(None, T)
            sx = CubicSpline(T, X, bc_type=bc_type)
            sy = CubicSpline(T, Y, bc_type=bc_type)
        except ValueError as msg:
            print unicode(msg)
            sx = sy = None
    return T, sx, sy

def rotate(points, alfa, centre):
    u'''
    :param alfa: float, l'angle de rotation en radians
    :param centre: ndarray((1,2)) ou list [float,float]
    :param keep: les numeros des points à laisser en l'etat.
        On peut avoir besoin de ne pas tourner le premier ou le dernier point
        (pour les NSplineComposee)
    :return Xt: une COPIE de points rotationnée(!).
    :param points: ndarray((n,2),dtype=float) est stocké par ligne,
        chaque point est de shape (1,2), il les faudrait en colonne (shape=(2,1))
        pour faire le produit matriciel.
        Donc on transpose tout et on ecrit Xi' = C' + (Xi'-C')*A' au lieu de
        Xi = C + A*(Xi-C), pour i= 0, 1,...
    '''
#     debug(points.shape, keep=keep)
    Ct = asarray(centre).reshape((1,2))
    cosa, sina = cos(alfa), sin(alfa)
    At = matrix([[cosa,-sina], [sina,cosa]]).transpose()
    Xt = points - Ct
    Xt = Xt*At + Ct
    return asarray(Xt)

def symetriser(points, sym):
    u"""
    Symetrie in situ de 'points' par rapport à un des plan Oxy, Oyz ou Oxz
    Parametres:
    ----------
        - points = ndarray de shape (m,n,...,3) ou (m,n,...,2)
        - sym = un chaine de caracteres pour indiquer l'orientation.
            sym contient au plus 3 caracteres parmi {'x','y','z'}
            si sym contient 'x' les points sont transformés [x,y,z] -> [-x,y,z] (points est symétrisé par rapport au plan Oyz)
            si sym contient 'y' les points sont transformés [x,y,z] -> [x,-y,z] (points est symétrisé par rapport au plan Oxz)
            si sym contient 'z' les points sont transformés [x,y,z] -> [x,y,-z] (points est symétrisé par rapport au plan Oxy)
    """

    shape = points.shape
    dim = shape[-1]#2d ou 3d
    if not dim in (2,3) : return points
    vpoints = points.view()
    vpoints.shape = -1,dim # on remet a plat sauf sur la dimension. vpoints est une liste (simple) de points 2D ou 3D
    sym = sym.lower()
    if 'x' in sym : vpoints[:,0] *= -1
    if 'y' in sym : vpoints[:,1] *= -1
    if dim==3 and 'z' in sym : vpoints[:,2] *= -1
    return points


def symetrieAxe(points,axe,centre=0.0):
    '''
    symetrie
    - d'axe 'axe' = 0 (vertical) ou 1 (horizontal)
    - de centre 'centre', réel
    X = 2*centre - X pour axe vertical
    Y = 2*centre - Y pour axe horizontal
    points est modifié in situ.
    '''
#     debug(points.shape, axe, centre)
    points.shape=(-1,2)
    points[:,axe]*=-1
    if centre !=0.0 : points[:,axe]+=2*centre
    return points

def hardScale(points, echelle, centre=None, translate=False):
    '''
    homothétie de centre 'centre': X = C + ech*(X-C) = (1-ech)*C + ech*X
    points est modifié in situ.
    '''
#     if centre is None :
#         centre=[0,0]
#     centre=asarray(centre).reshape((1,2))
    points *= asarray([echelle[0],echelle[1]])
    if centre is not None :
        centre = asarray(centre)
#         debug(centre)
        centre.shape = (1,2)
        points += [1.0-echelle[0], 1.0-echelle[1]]*centre
    return points

def aire2d(points,absolute=True):
    '''Aire algébrique ou absolue délimitée par les points de points.
    = 0.5*det(AP(i), AP(i+1)) 1<=i<n avec points = {points[i], 0<=i<=n}, et A=P[0]'''
#    if isinstance(points, QPolygonF) :
#        points = pointsFromPolygon(points)
    if len(points)<=2 :
        return 0.0
    points=copy(points)#sinon effet de bord difficile à identifier !!!
    S=0.0
    A=points[0]
    points-=A
    if absolute :
        for b,c in zip(points[1:-2],points[2:-1]) :
            d=abs(b[0]*c[1]-b[1]*c[0])
            S+=d
    else :
        for b,c in zip(points[1:-2],points[2:-1]) :
            d=b[0]*c[1]-b[1]*c[0]
            S+=d
    return 0.5*S

def moyenneMobileClosed(points,molecule=[1.,2.,1.],nnodes=None):
    '''
    Lissage moyenne mobile pour les polygones fermés self[-1]==self[0]
    nnodes est la liste des numeros de points à lisser. On suppose que 0 n'y est pas
    '''
    molecule=asarray(molecule)#asarray([1., 2., 1.])/4.
    molecule/=sum(absolute(molecule))
    new=points
    if nnodes is None :
        nnodes=range(len(points))
    else :
        nnodes.sort()
        if nnodes[0]==0 and nnodes[-1]==len(points) :
            trace('',u'non implémenté')
            return new
        else : pass
    old=new[:-1]
    n=len(old)
#         deb, fin = 0, n
    deb,fin=nnodes[0],nnodes[-1]
    for k in range(deb,fin) :
        pm=old[k-1]
        p=old[k]
        if k==n-1 :
            pp=old[0]
        else:
            pp=old[k+1]
        new[k]=molecule[0]*pm+molecule[1]*p+molecule[2]*pp
    new[-1]=new[0]
    return new

def moyenneMobile(X,n=1):#,molecule=[1.,2.,1.]):
    '''
    Lissage moyenne mobile pour fonction simple i -> X[i]
    Y[i] = (X[i-1]+2X[i]+X[i+1])/4, 0<i<n
    Y[0] = (       2X[0]+X[1])/3, i=0
    Y[n] = (X[n-1]+2X[n]       )/3, i=n
    '''
    Y = X.copy()
    for k in range(1, len(X)-1):
        Y[k] = 0.25*(X[k-1] + 2*X[k] + X[k+1])
    Y[0]  = (        2*X[0] + X[1])/3.0
    Y[-1] = (X[-2] + 2*X[-1]      )/3.0
    if n==1 :
        return Y
    else :
        return moyenneMobile(Y, n-1)

def moyenneMobile1(X, n=1):#,molecule=[1.,2.,1.]):
    '''
    n lissages par moyenne mobile pour fonction simple i -> X[i]
    Y[i] = (X[i-1]+2X[i]+X[i+1])/4, 0<i<n
    Y[0] = X[0]                   , i=0
    Y[n] = X[n]                   , i=n
    '''
    Y = X.copy()
    for k in range(1, len(X)-1):
        Y[k] = 0.25*(X[k-1] + 2*X[k] + X[k+1])
    Y[0]  = X[0]
    Y[-1] = X[-1]
    if n==1 :
        return Y
    else :
        return moyenneMobile1(Y, n-1)

def centreGravite(points, surface=False):
    u'''
    Centre de gravité du polygone délimité par les points de 'points'
    Le polygone EST une plaque plane (de densité surfacique constante),
    Le polygone N'EST PAS un ensemble de cotés de masse linéique constante,
    Le polygone N'EST PAS un ensemble de masses constantes disposées aux sommets.
    La méthode utilisée est
    -de fixer un point O quelconque,
    -de remplacer chaque triangle t(i)=OP(i)P(i+1) par une masse ponctuelle proportionnelle à sa surface S(t(i))
        située aux centre de gravité du triangle.
    -puis de faire le barycentre de ces masses affectées de la surface du triangle.

    Que le polygone soit fermé ou non ne change pas le résultat.
    On obtient :
    S*OG = som(S(ti)*OG(ti), pout ti parcourant les triangles ti=OP(i)P(i+1) 0<=i<n
    S est l'aire algebrique
    S(ti) est l'aire algébrique du triangle ti , i.e. 0.5*det(OP(i), OP(i+1))
    '''
    if len(points) == 0 :
        if surface : return asarray([nan, nan]), 0.0
        else : return asarray([nan, nan])
    elif len(points) == 1 :
        if surface : return asarray(points[0]), 0.0
        else : return asarray(points[0])

    points=asarray(points)
    Sa = 0.0
    xG, yG = 0.0, 0.0
    G = asarray([xG, yG])
    A = points[0]
    T = zip(points[:-1], points[1:])+[(points[-1],points[0])]
    for b,c in T :
#    for b, c in zip(P[1:-2], P[2:-1]) :?????
        Gt = (A + b + c)/3.0
        Sat = (b[0] - A[0])*(c[1] - A[1]) - (b[1] - A[1])*(c[0] - A[0])
        G += Sat*Gt
#         debug( "pouet, a faire")
        Sa += Sat
    if Sa == 0.0 :
        if surface :  return asarray((nan, nan)), 0
        else : return asarray((nan, nan))
    else :
        if surface : return G / Sa, Sa*0.5
        else : return G / Sa

def aire(points):
    '''
    Calcul de l'aire algébrique d'un polygone.
    :ATTENTION:
    - si points est un polyligne (non fermé) il est d'abord fermé
        l'aire d'un polyligne non fermé devrait être 0.0, ça n'est pas le cas ici
    - l'aire est algébrique. Donc si le polyligne se recoupe (papillote) l'aire
        peut être nulle
    Si le polygone ne se recoupe pas, (intérieur connexe), alors l'aire donne le
    sens de rotation :
    si elle est positive, ses points tournent dans le sens trigonométrique, sinon,
    le sens des aiguilles d'une montre.
    La méthode utilisée est
    - fermer le polygone : P(n)=P(0)
    - de fixer un point 'O' quelconque, ici on prend O=P[0]
    - sommer les aires des triangles A = som(A(ti)) 0<=i<n) où
    A(ti) = (1/2).OP(i)^OP(i+1) est l'aire algébrique du triangle ti=OP(i)P(i+1)
    '''
    if len(points) <= 2 : return 0.0
    points=asarray(points)
    a = points[0]
    aire = 0.0
    T = zip(points[:-1], points[1:])+[(points[-1],points[0])]
    for b, c in T :
        At = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0])
        aire += At
    return aire*0.5

def baryCentre(points,masses=None):
    '''Retourne le barycentre du nuage de points 'points' affectés des masses 'masses' 2d'''
    if len(points)==0 :
        return None
    points=asarray(points)
#     trace('', points.shape, len(points))
    N=len(points)
    if masses is None :
        X,Y = points[:,0], points[:,1]
#         debug (X=X, Y=Y)
    else :
        masses=masses/sum(masses)
        X,Y=prod([masses,points[:,0]],axis=0),prod([masses,points[:,1]],axis=0)
    bc = [[sum(X)/N, sum(Y)/N]]
    return asarray(bc)


def p2a(qpoint):
    """
    retourne les coord d'un QPointF sous forme de np.ndarray((1,2))
    => [[x,y]] cf p2a2()
    """
#     print qpoint
    point = asarray((qpoint.x(),qpoint.y())).reshape((1,2))
#     trace('', point)
    return point

p2a12=p2a

def dist2(p1,p2,n=2):
    '''retourne le carré de la distance de p1 à p2 en norme n=2'''
    try :
        x1,y1=p1[0],p1[1]
        x2,y2=p2[0],p2[1]
        return (x2-x1)**2+(y2-y1)**2
    except TypeError : # Si p1 ou p2=None
        return nan

def dist(p1,p2,n=2):
    '''retourne la distance de p1 à p2 en norme n=2'''
    return math.sqrt(dist2(p1,p2))

def absCurv(points, normalise=False):
    """Abscisse curviligne des points de points, considéré comme polyligne"""
#     debug(points)
#     stack()
    l = len(points)
    if l==0 : return []
    elif l==1 : return [0]
    T = zeros(l)
    for k in range(l-1) :
        T[k+1] = T[k] + dist(points[k],points[k+1])
    if normalise and T[-1] != 0.0 :
        T /= T[-1]
        #sinon T=[0,0,...] On peut très bien entrer ici avec 3 points = 0
    return T


def splineCubique(points,weights=None,bbox=((0,-1),(0,-1)),eps=1.0e-2,degre=(2,3)):
    """
    bientot Obsolete. On utilise splineInterpolation
    Ajustement d'une spline cubiques paramétriques au polygone points
    sx(t), sy(t) sont deux splines.
    Retourne:
    --------
    T = les valeurs du paramètre t, abscisse curviligne NORMALISEE, entre 0 et 1.
    sx, sy = les deux splines

    TODO: parametres de précision et d'ajustement
    - ajust : proportionnel à la surface (au lieu de T[-1]*T[-1])
    - s : proportionnel à deltax et deltay suivant le cas (au lieu de T[-1])
    """
    N=len(points)
    if len(points)<3 :
        raise ValueError('Pas assez de points (%d<3), lissage impossible'%N)
    T=absCurv(points)
    T/=T[-1]
    #print whoami(), 'len(points), T[-1]', len(points), T[-1]
    ajust=eps*T[-1]*T[-1]
    for k,weight in enumerate(weights) :
        weights[k]=weight if weight>=0 else len(weights)
    weights*=N/sum(weights)#Charge totale = N
    bbx=[T[bbox[0][0]],T[bbox[0][1]]]
    bby=[T[bbox[1][0]],T[bbox[1][1]]]
    try :
        sx=interpolate.UnivariateSpline(T,points[:,0],w=weights,bbox=bbx,k=degre[0],s=ajust*T[-1])#s=la précision de l'ajustement s=0 <=> interpolation
        sy=interpolate.UnivariateSpline(T,points[:,1],w=weights,bbox=bby,k=degre[1],s=ajust*T[-1])
    except :
        raise RuntimeError("%s : impossible de calculer la spline"%whoami())
    return T,sx,sy

def aireTriangle(a,b,c):
    """
    Aire du triangle abc dans l'espace.
    C'est la moitié de la norme du produit vectoriel ab vect ac
    """
    u,v=b-a,c-a
    r=u[2]*v[0]-u[0]*v[2]
    s=u[0]*v[1]-u[1]*v[0]
    t=u[1]*v[2]-u[2]*v[1]
    return 0.5*sqrt(r*r+s*s+t*t)

def intersectionSegments(seg1, seg2):
    u'''http://www.exaflop.org/docs/cgafaq/cga1.html :
    Subject 1.03: How do I find intersections of 2 2D line segments?

This problem can be extremely easy or extremely difficult depends on your applications.
If all you want is the intersection point, the following should work:

Let A,B,C,D be 2-space position vectors.  Then the directed line segments AB & CD are given by:

        AB=A+r(B-A), r in [0,1]
        CD=C+s(D-C), s in [0,1]

If AB & CD intersect, then

        A+r(B-A)=C+s(D-C), or

        Ax+r(Bx-Ax)=Cx+s(Dx-Cx)
        Ay+r(By-Ay)=Cy+s(Dy-Cy)  for some r,s in [0,1]
Solving the above for r and s yields

            (Ay-Cy)(Dx-Cx)-(Ax-Cx)(Dy-Cy)
        r = -----------------------------  (eqn 1)
            (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------  (eqn 2)
            (Bx-Ax)(Dy-Cy)-(By-Ay)(Dx-Cx)
Let P be the position vector of the intersection point, then

        P=A+r(B-A) or

        Px=Ax+r(Bx-Ax)
        Py=Ay+r(By-Ay)
By examining the values of r & s, you can also determine some other limiting conditions:

        If 0<=r<=1 & 0<=s<=1, intersection exists
            r<0 or r>1 or s<0 or s>1 line segments do not intersect
If the denominator in eqn 1 is zero, AB & CD are parallel

If the numerator in eqn 1 is also zero, AB & CD are coincident

If the intersection point of the 2 lines are needed (lines in this context mean infinite lines) regardless whether the two line segments intersect, then

If r>1, P is located on extension of AB

If r<0, P is located on extension of BA

If s>1, P is located on extension of CD

If s<0, P is located on extension of DC

Also note that the denominators of eqn 1 & 2 are identical.
    '''
#     print seg1, seg2
    A,B=seg1# np.ndarray((2,))
#     trace('',A.shape)
    u=B-A
    C,D=seg2
    v=D-C

    w=C-A
    denominateur=u[0]*v[1]-u[1]*v[0]
    if denominateur==0.0 :
        return nan,nan,asarray((nan,nan))
    else  :
        r=(w[0]*v[1]-w[1]*v[0])/denominateur
        s=(w[0]*u[1]-w[1]*u[0])/denominateur
        return r,s,A+r*u

def segmentPlusProche(points,P):
        u'''
        Paramètres:
        ----------
            - points : ndarray((n,2)) est un polyligne dans lequel on cherche le segment
            - P : ndarray((1,2)) est un point quelconque
        Fonctionnement:
        --------------
            On recherche le segment S[i]=(points[i], points[i+1]) le plus proche de P au sens suivant :
            pour tous les segments (A,B)=(points[i], points[i+1]), P' désigne la projection orthogonale
            de P sur la droite AB. Si P' est INTERIEURE au segment [A,B], le segment est candidat
            Parmi tous les segments candidats, on retient celui qui réalise la plus courte distance PP'.
        Retourne :
        --------
            - (i, P') i=numéro du segment et P'=le projeté de P.
            - (None, None) s'il n'y a pas de candidat.

 Voir sur http://www.exaflop.org/docs/cgafaq/cga1.html#Subject%201.02:%20How%20do%20I%20find%20the%20distance%20from%20a%20point%20to%20a%20line?

    Subject 1.02: How do I find the distance from a point to a line?

Let the point be C (Cx,Cy) and the line be AB (Ax,Ay) to (Bx,By).    The length of the line segment AB is L:

        L= sqrt( (Bx-Ax)^2 + (By-Ay)^2 ) .

Let P be the point of perpendicular projection of C onto AB. Let r be a parameter to indicate P's location along the line containing AB, with the following meaning:

          r=0      P = A
          r=1      P = B
          r<0      P is on the backward extension of AB
          r>1      P is on the forward extension of AB
          0<r<1    P is interior to AB

Compute r with this:

            (Ay-Cy)(Ay-By)-(Ax-Cx)(Bx-Ax)
        r = -----------------------------
                        L^2

The point P can then be found:

        Px = Ax + r(Bx-Ax)
        Py = Ay + r(By-Ay)

And the distance from A to P = r*L.

Use another parameter s to indicate the location along PC, with the following meaning:

           s<0      C is left of AB
           s>0      C is right of AB
           s=0      C is on AB

Compute s as follows:

            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------
                        L^2

Then the distance from C to P = s*L.
        '''

        # if isinstance(P,(QPointF, QPoint,)):
        #     P = P.x(), P.y()# Point argument
        points = asarray(points)
        u = -(points-P)[:-1]
        v = points[1:]-points[:-1]
#         debug(v=v)
#         distances = [dist(point,(X,Y)) for point in points]
#         debug(points-P)
#         debug(distances=distances)
#         debug(P, len(points), distances)
#         distances = linalg.norm(u, axis=1)
        longueurs = norm(v, axis=1)
#         debug(distances_P_points=distances)
#         debug(longueurs_segments=longueurs)
        ps = [dot(ui, vi) for (ui, vi) in zip(u,v)]
        r = ps/(longueurs*longueurs)
#         debug(r_entre_0_et_1=r)
#         psn = [dot(ui, vi)/dot(vi,vi) for (ui, vi) in zip(u,v)]
#         debug(psn=psn)
#         distances = r*longueurs
#         debug(distances_P_Segment=distances)
        candidats, = where(logical_and(0.0<=r,r<=1.0))
#         debug(len(points), len(r), len(v))
#         if len(candidats)>0 :
#             projetes = points[candidats]+r[candidats]*v[candidats]
        distances = []
        projetes = []
        for candidat in candidats :
#                 debug(candidat=candidat)
#                 debug(r_candidat=r[candidat])
#                 debug(p_candidat=points[candidat])
#                 debug(v_candidat=v[candidat])
            H = points[candidat]+r[candidat]*v[candidat]
            projetes.append(H)
            distances.append(dist(P,H))

        if len(candidats) == 0 :
            pp, distances = pointLePlusProche(points, P, return_distances_array=True)
            #parcequ'il faut bien retourner quelque chose
#             if pp == 0 : pp=1#segment p[0], p[1] <= bug
            if pp == len(points)-1 : pp = len(points)-2#segment p[-2], p[-1]
            return pp, None
        else :
            winner = argmin(distances)
#             debug(winner=winner)
#             debug(candidats_winner=candidats[winner])
            i = candidats[winner]
            projete = projetes[winner]
#         debug(i=i, projete=projete)
            return i, projete

def pointLePlusProche(points,P,return_distances_array=False):
    '''
    :param points: ndarray, tableau de n points 2d de shape (n,2)
    :param P: un point 2d de type QPointF, QPoint ou tuple ou liste (x,y)
    :param return_distances_array: bool, retourne ou non le tableau des distances.

    :return: (i, dist) avec
        - i : l'indice dans 'points' du point qui réalise le min
        - dist : la distance de P à points si return_distances_array est False
        - dist : le tableau des distances de P à points si return_distances_array est True
    '''
    X,Y=P[0],P[1]
    distances=[dist(point,(X,Y)) for point in points]
    index=argmin(distances)
    if return_distances_array :
        return index,distances
    else :
        return index,distances[index]


def maintenant():
    u'''La date et l'heure formatées "human readable"'''
    return unicode(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))

def doublons(points,voisinages=[],eps=1.0e-8,sym=True):
    u"""
    Verifier qu'il n'y a pas de point double dans points.
    :param points : ndarray(nbp,3) nbp points 3d. nbp peut être une
        shape (nl, nc, 3) et dans ce cas nbp = nl*nc
    :param voisinages : ndarray(shape=(np,1+nv), dtype=int) np numéros de points,
        faisant référence au tableau 'points', chacun des np points à
        vérifier possède nv voisins. voisinage est donc de la forme
        [
         [k_0,  v_0_1,  v_0_2, ...,  v_0_nv],
         [k_1,  v_1_1,  v_1_2, ...,  v_1_nv],
         ...
         [k_np, v_np_1, v_np_2, ..., v_np_nv]
        ]
        les valeurs k_i, v_i_j sont des entiers entre 0 et nbp-1
    :param eps : réel. deux points p,q sont considérée comme confondus si dist(p,q)<eps
    :param sym : bool. True si i voisin de j <==> j voisin de i. Divise par deux
        le travail.

    :return : une liste de doublons (d_0,dd_0), ... (d_n, dd_n)
    """
    def cond(v,k):
        if sym : return v>k
        else : return v!=k

    qoints=points.view()
    qoints.shape=(qoints.size/3,3)
    nbp=qoints.shape[0]
    if voisinages!=[] :
        voisinages=asarray(voisinages)
#        debog(whoami(), voisinages.shape)
        nbp=voisinages.shape[-1]
        voisinages.shape=(-1,nbp)
    else :
        if sym :
            voisinages=asarray([[i]+range(i,nbp) for i in range(nbp)])
        else :
            voisinages=asarray([[i]+range(nbp) for i in range(nbp)])

    doublons=[]
    for voisinage in voisinages:
        k,voisins=voisinage[0],voisinage[1:]
        point=qoints[k]
        for voisin in voisins :
            if cond(voisin,k) :
                if norm(point-qoints[voisin])<eps :
                    doublons.append((k,voisin))
    return doublons

u"""Quelques fonctions de débogage qui ont évolué au cours du temps.

Ecriture dans le fichier log 'AXile.log'
----------------------------------------

1- trace et alert

    >>> debug(x, y=yzt) écrit dans le fichier log niveau INFO
    >>> rdebug(x, y=yzt) écrit dans le fichier log niveau WARNING

    par exemple supposant que nums est une variable valant (1,2,3)
    pour l'appel suivant,
    >>> trace('voila:', numeros=nums)
    la sortie ressemble  à

    GuiCloison::backup [glideobject.py, line 327]  :
        voila ; numeros = (1,2,3)
    où
    - GuiCloison::backup sont la classe et la méthode, d'où est appelé trace (alert)
    - [glideobject.py, line 327], sont le nom du fichier, le numéro de ligne,
    -     voila ; numeros = (1,2,3) sont les variables

2- strace et salert

    >>> debug(x, y='toto') écrit dans le fichier log niveau INFO
    >>> rdebug(x, y='toto') écrit dans le fichier log niveau WARNING

    la sortie est analogue à trace et alert, un peu remaniée,
    sans le nom de classe et ressemble à ceci :

    [glideobject.py:327 (backup)] valeur_de_x ; y='toto'

    où
        - glideobject.py:327 sont le nom du fichier, le numéro de ligne,
        - (backup) est la fonction d'où est appelé strace (salert)

Ecriture sur la console
-----------------------

1-debug, rdebug

    >>> debug(x, y='toto', z='tutu')
        fait la même chose que strace, sur la console stdout
    >>> rdebug(x, y='toto', z='tutu')
        fait la même chose que salert, sur la console stderr (rouge dans eclipse)

2- stack et rstack

    >>> stack(commentaire)
        écrit le commentaire et la pile des appels sur stdout
    >>> rstack(commentaire)
        écrit le commentaire et la pile des appels sur stderr (rouge dans eclipse)

    la sortie ressemble à ceci :

    ========================
    === Pile des appels ====

        le commentaire

        la pile des appels
        ...

    ======  Fin pile  ======
    ========================

"""
def sexplore(objet):
    '''retourne une chaine de caracteres avec les attributs de objet et leur valeur'''
    s=[]
    dexp=explore(objet)
    for key in sorted(dexp.keys()) :
        s.append('%20s : %s'%(key,dexp[key]))
    return '\n'.join(s)

def explore(objet):
    '''retourne un dictionnaire avec les attributs de objet et leur valeur'''
    d={}
    for key,value in objet.__dict__.iteritems() :
        if type(value) in (list,tuple) :
            try : d[key]="[%s]"%(type(value[0]).__name__)
            except IndexError: d[key]=("[]")
        else :
            d[key]='%s'%(type(value).__name__)
    return d

def className(obj):
    try : return obj.__class__.__name__
    except : return 'Inconnu'

# def phrase():
#     from config import DATA_DIR
# #     try :
#     f=open(Path(DATA_DIR,'dico'))
# #     except IOError:
# #         return
#     words=qui,quoi,comment=[],[],[]
#     n=0
#     for line in f.readlines():
#         line=line.strip()
#         if line : words[n].append(line)
#         else : n+=1
#     return ' '.join((3*'*',choice(qui),choice(quoi),choice(comment)))

'''Ecriture dans fichiers log'''
def trace(*args,**kargs):
    return _trace(sys.stdout, *args, **kargs)

def alert(*args,**kargs):
    return _trace(sys.stderr, *args, **kargs)

def _trace(output, *args,**kargs) :
    try : val = args[0] in (None,'')
    except Exception : val = False
    if val : msg = u''
    else:
        try : msg = args[0].__class__.__name__+u'::'
        except AttributeError : msg = u''
    frame = sys._getframe(2)
    fcode = frame.f_code
    filename = Path(fcode.co_filename).name
    msg0 = u'\n'+unicode(msg + fcode.co_name+u' [%s, line %d] '%(filename, frame.f_lineno))+' : \n    '
    lmsg = [unicode(arg) for arg in args[1:]]+\
           [unicode(key)+u' = '+unicode(value) for key,value in kargs.iteritems()]
    if output is sys.stdout :
        print msg0 + u' ; '.join(lmsg)
    #     logger.info(msg0 + u' ; '.join(lmsg))
    elif output is sys.stderr :
    #     logger.warning(msg0+u' ; '.join(lmsg))
        print>>sys.stderr, msg0 + u' ; '.join(lmsg)

def rstack(commentaire=''):
    _stack(sys.stderr, commentaire)
def stack(commentaire='') :
    _stack(sys.stdout, commentaire)

def _stack(output, commentaire=''):
    u"""Impression de la pile des appels sur output"""
#     print>>output,u"========================"
#     print>>output,u"=== Pile des appels ===="
    print>>output,40*u"="
    print>>output,10*u"="+u"   Pile des appels  "+10*u"="
    if commentaire :
        print>>output, commentaire
    k = 2
#     frame = sys._getframe(2)

    while 1 :
        try :
            frame = sys._getframe(k)
            fcode = frame.f_code
            fonction = fcode.co_name
            filename = Path(fcode.co_filename)#.name
            toclick = clickableLink(filename, frame.f_lineno)#, fcode.co_name)
            print>>output,u'    [%-20s] '%fonction + toclick
            k += 1
#             filename = Path(frame.f_code.co_filename).name
#             print>>output, u" ['%s':%d, (%s)] "%(filename, frame.f_lineno, frame.f_code.co_name)
        except :
            break
    print>>output,10*u"="+u"      Fin pile      "+10*u"="
    print>>output,40*u"="
#     print>>output,u"======  Fin pile  ======"
#     print>>output,u"========================"
    return

def clickableLink(filename, lineno):#, name):
#     filename = filename.replace(SOURCES_DIR,'').strip('/\\')
    return unicode(u'File "%s", line %d'%(filename, lineno))#, name))

# clickableLink = _clickableLink if not config.TEST_MODE else lambda f,l : f.name

def _strace(*args,**kargs) :
    u"""_trace simplifié, sans le nom de la classe"""
    frame = sys._getframe(2)
    fcode = frame.f_code
    fonction = fcode.co_name
    filename = Path(fcode.co_filename)#.name
    if spleenconfig.TEST_MODE :
        toclick = u"(module %s)"%filename.name
    else :
        toclick = clickableLink(filename, frame.f_lineno)
    output = args[0]
    args = args[1:]
    try :
        titre = sursouligne(kargs.pop(u'titre'),u'=')
    except KeyError :
        titre = None
    try :
        paragraphe = souligne(kargs.pop(u'paragraphe'),u'-',shift=4)
    except KeyError :
        paragraphe = None

    lmsg = [unicode(arg) for arg in args]+\
           [unicode(key)+u" = "+unicode(value) for key,value in kargs.iteritems()]
    msg = u'[%s] '%fonction + toclick
    b, l, s = u'    ', u'\n', u' ; '#blanc, saut de ligne, separateur
    if titre and paragraphe :
        msg = msg + 2*l + titre + l + paragraphe
    elif paragraphe :#pas de titre
        msg = msg + 2*l + paragraphe
    elif titre :
        msg = msg + l + titre + l
    msg = msg + l + s.join(lmsg)
    if 0:#output in (logger.info, logger.warning):
        output(msg)
    else :
        print>>output, msg

def debug(*args,**kargs):
    u"""trace simplifié, stdout, identique à strace"""
    _strace(sys.stdout, *args, **kargs)
#     _strace(logger.info, *args, **kargs)

def rdebug(*args, **kargs):
    u"""trace simplifié, logger.warning, identique à salert"""
    _strace(sys.stderr, *args, **kargs)
#     _strace(logger.warning, *args, **kargs)

# def strace(*args,**kargs) :
#     u"""trace simplifié, sur logger.info"""
#     _strace(logger.info, *args, **kargs)
#
#
# def salert(*args,**kargs) :
#     u"""alert simplifié, sur logger.warning"""
#     _strace(logger.warning, *args, **kargs)

def whoami(objet=None):
    if objet is None:
        msg=''
    else:
        try : msg=objet.__class__.__name__+'.'
        except AttributeError : msg=''

    return '%s::'%(msg+sys._getframe(1).f_code.co_name)#, id(object)

def toDict(cles,valeurs):
    """
    - cles = "cle1 cle2 cle3 ..."
    - valeurs = "val1 val2 val3...", les valeurs sont des entiers ou des reels
    retourne un dictionnaire cle,valeurs
    """
    d={}
    for key,value in zip(cles.split(),valeurs.split()) :
        try: w=int(value)
        except ValueError : w=float(value)
        d[key.lower()]=w
    return d

def findAll(tag,lines,first=0,last=-1):
    """
    Retourne une liste triée des numeros de ligne de
    TOUTES les occurences de tag dans lines[first:last]
    """
    #
    n=first-1
    n0=findRel (tag,lines,n+1,last)
    N=[]
    while n0>=0 :
        n=n0+n+1
        N.append(n)
        n0=findRel(tag,lines,n+1,last)
    return tuple(N)

def findAllLines(tag,lines,first=0,last=-1):
    '''
    Comme findAll(), mais il faut que la ligne complete (nettoyée) soit égale à tag, pas seulement une partie de la ligne.
    Par exemple : line = 'TOTO_EST_CONTENT' ne match pas avec tag='TOTO'
    '''
    liste=findAll(tag,lines,first,last)
    newlist=[]
    for n in liste :
        if lines[n].strip()==tag :
            newlist.append(n)
    return tuple(newlist)

def findRel(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i-first de la premiere occurence trouvee,
          c'est à dire le numéro de ligne RELATIF : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    found=find(tag,lines,first,last)
    if found is None : return None
    else : return found-first

def find0(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i-first de la premiere occurence trouvee,
          c'est à dire le numéro de ligne RELATIF : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    found=findRel(tag,lines,first=first,last=last)
    if found is not None : return  found+first
    else : return None

def find(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i de la premiere occurence trouvee,
          c'est à dire le numéro de ligne ABSOLU : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    if last<0 : last=len(lines)+last
    if first is None : first=0
    elif first<0 : first=len(lines)+first
    if not tag : return None
    i=first-1
    while i<last:
        i+=1
        try :
            if lines[i].find(tag)>=0 :
                return i
        except IndexError :
            return None
    return None

def rreplace(chaine,sub,replace):
    '''
    Remplace récursivement dans chaine, toutes les occurences de sub  par replace
    >>> rreplace("aaaaab","aa","a")
    ... "ab"
    '''
    while 1:
        oc=chaine
        chaine=oc.replace(sub,replace)
        if oc==chaine : break
    return chaine


def find_names(obj):
    """Renvoit les noms de variable de l'objet 'obj' passé en paramètre"""
    frame=sys._getframe()
    for frame in iter(lambda: frame.f_back,None):
        frame.f_locals
    result=[]
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer,dict):
            for k,v in referrer.iteritems():
                if v is obj:
                    result.append(k)
    return result

def structuredGrid(nx,ny,nz,delta=-1):
    """
    retourne un tableau ndarray de shape(nx,ny,nz,3)

    produit tensoriel de X=(x_1,x_2,...x_nx), Y=(y_1,y_2,...y_ny)
    et Z=(z_1,z_2,...z_nz), aléatoirement espacés, mais triés.
    ou bien régulièrement espacés
    Contrainte : nx>1, ny>1, nz>1
    """
#    nx,ny,nz = shape
    if delta==-1 :
        X=random(nx)
        Y=random(ny)
        Z=random(nz)
        X.sort()
        Y.sort()
        Z.sort()
    else :#TODO: delta=dx,dy,dz
        dx,dy=1.0/(nx-1),1.0/(ny-1)
        X=arange(0.0,1.0+dx,dx)
        Y=arange(0.0,1.0+dy,dy)
    if nz==1 :
        Z=asarray([0.0])
    else :
        dz=1.0/(nz-1)
        Z=arange(0.0,1.0+dz,dz)
    G=ndarray((nx,ny,nz,3),dtype=float)
    for i,x in enumerate(X) :
        for j,y in enumerate(Y) :
            for k,z in enumerate(Z) :
                G[i,j,k,:]=[x,y,z]
    return G


def my2dPlot(XYs,legends=[],equal=False,cosmetic=[], title='No title'):
    u'''
    Tracer des points en 2d
    XYs est une liste de (6 maximum) tableaux de points chaque tableau est de shape (n,2) ou (n,3),
    comporte n points en 2d ou 3d. seules les deux premieres coord de chaque point sont prises en compte.
    '''

    from matplotlib import pyplot
#    matplotlib.use('MacOSX')
    nbcourbes=len(XYs)
#     print 'nb courbes', nbcourbes
    if not cosmetic :
        cosmetic=(
                  'r-o','g-o','b-o',
                  'r-^','g-^','b-^',
                  'r:.','g:.','b:.',
                  'r*','g*','b*',
                  'r^','g^','b^',
                  'r.','g.','b.',
                  )
        cosmetic=cosmetic[:nbcourbes]
    if legends in ([],) :
        legends=[str(k) for k in range(nbcourbes)]
#     if legends is None:
#         legends = []
    xmin=ymin=10^10
    xmax=ymax=-10^10
    colors=cosmetic
    for k,xy in enumerate(XYs) :
        try: color=colors[k]
        except : color='b-'
        pyplot.plot(xy[:,0],xy[:,1],color)
        xmin=min(xmin,min(xy[:,0]))
        ymin=min(ymin,min(xy[:,1]))
        xmax=max(xmax,max(xy[:,0]))
        ymax=max(ymax,max(xy[:,1]))
    w,h=xmax-xmin,ymax-ymin
    # dizaine=int(math.log10(w))#10^dizaine <= w < 10^(1+dizaine)
    ax=pyplot.axes()
    if equal : ax.set_aspect('equal')
    ax.set_xlim(xmin-w/10,xmin+w+w/10)
    ax.set_ylim(ymin-h/10,ymin+h+h/10)
    ax.grid(which='major',axis='x',linewidth=0.75,linestyle='-',color='0.75')
    ax.grid(which='minor',axis='x',linewidth=0.25,linestyle='-',color='0.75')
    ax.grid(which='major',axis='y',linewidth=0.75,linestyle='-',color='0.75')
    ax.grid(which='minor',axis='y',linewidth=0.25,linestyle='-',color='0.75')
#    ax.xaxis.set_major_locator(pyplot.MultipleLocator(10**dizaine))
#    ax.xaxis.set_minor_locator(pyplot.MultipleLocator(10**(dizaine-1)))
#    ax.yaxis.set_major_locator(pyplot.MultipleLocator(10**dizaine))
#    ax.yaxis.set_minor_locator(pyplot.MultipleLocator(10**(dizaine-1)))
    if legends is not None:
        pylab.legend(legends,shadow=True)#, loc = (0.01, 0.55))
        ltext=pylab.gca().get_legend().get_texts()
        for k, _ in enumerate(legends) :
            pylab.setp(ltext[k],fontsize=10)#, color = 'b')

#     print ltext
    pylab.title(title)
    pyplot.show()

def load(filinname):
    """pickle => python"""
    try :
        filin=open(filinname,'r')
        return pickle.load(filin)
    except IOError :
        print>>sys.stderr,"Impossible d'ouvrir le fichier dump %s, pas de lecture."%filinname
#    raise PyglideIOError("Impossible d'ouvrir le fichier dump %s, pas de lecture."%filinname)
# def doc(object):
#     for key, value in object.__dict__.iteritems() :
#         if type(value).__name__ in ('function',) :
#             print "%s:\n"%(key)+len(key)*'='+ "\n%s"%value.__doc__



def goodName(name):
    allowed=string.ascii_letters+string.digits+'_-+=#'
    for letter in name :
        if letter not in allowed :
            return False
    return True

def diff(A):
    u"""
    J'en ai marre de taper tout le temps 'dA = A[1:]-A[:-1]'
    :param A: un tableau de points ndarray de shape (n, dim)
    :return d:les différences premières de A:
        dA[i] = A[i+1]-A[i]
    :rtype dA: ndarray((n-1,dim))
    """
    return A[1:]-A[:-1]

def XY(xy):
    u"""
    J'en ai marre de taper tout le temps (pour le graphique) 'X,Y = xy[:,0],xy[:,1]'
    :param A: un tableau de points ndarray de shape (n, 2)
    :return (X, Y) : (ndarray((n,1)), ndarray((n,1))) les deux colonnes de A
    """
    return xy[:,0],xy[:,1]

def rayCross(polygon, point):
    u""" détermine si point est à l'intérieur ou a l'extérieur du polygône
    par la technique de "ray-crossing"
        The essence of the ray-crossing method is as follows.
        Think of standing inside a field with a fence representing the polygon.
        Then walk north. If you have to jump the fence you know you are now
        outside the poly. If you have to cross again you know you are now
        inside again; i.e., if you were inside the field to start with, the total
        number of fence jumps you would make will be odd, whereas if you were
        ouside the jumps will be even.

        The code below is from Wm. Randolph Franklin <wrf@ecse.rpi.edu>
        (see URL below) with some minor modifications for speed.  It returns
        1 for strictly interior points, 0 for strictly exterior, and 0 or 1
        for points on the boundary.  The boundary behavior is complex but
        determined; in particular, for a partition of a region into polygons,
        each point is "in" exactly one polygon.
        (See p.243 of [O'Rourke (C)] for a discussion of boundary behavior.)

        int pnpoly(int npol, float *xp, float *yp, float x, float y)
        {
          int i, j, c = 0;
          for (i = 0, j = npol-1; i < npol; j = i++) {
            if ((((yp[i]<=y) && (y<yp[j])) ||
                 ((yp[j]<=y) && (y<yp[i]))) &&
                (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))

              c = !c;
          }
          return c;
        }
    """
    print 'todo...'

def isInsideOrFrontier(p, Q, eps=1.0e-6):
    u"""
    Utilise le package shapely (Polygon et Point)
    :param p : (x,y) = ndarray((1,2), float) un point en 2d
    :param Q : est un polygone = ndarray((n,2), float) supposé fermé i.e.
        le premier point == le dernier point,
    :param eps : float = precision
    :return :
        True si p est à l'intérieur du polygone Q ou sur la frontiere
            (bande de largeur eps) et
        False s'il est à l'extérieur.
    """
    P = Polygon(Q)
    p = Point(p)
    return P.distance(p)<eps

def pointsDoubles(points, eps=1.0e-10):
    u"""retourne la liste des numéros de points doubles CONSECUTIFS à eps pres.
    points est de shape(n,2), n>=2"""
    if not isinstance(points, (ndarray,list,tuple)) :
        raise NotImplementedError('points doit etre de type numpy.ndarray')
    if len(points)<2 : return []
    ac = absCurv(points)
    dac = ac[1:]-ac[:-1]
    avirer = []
    for k, d in enumerate(dac) :
        if abs(d)<=eps :
            avirer.append(k)

    return avirer#, agarder


def eliminerPointsDoublesConsecutifs(points, eps=0.0, vires=False):
    avirer = pointsDoubles(points, eps)
#     debug(len_avirer=len(avirer))
    if avirer : points = delete(points, avirer, 0)
    return points, avirer if vires else points
#         msg1 = 'Les points de controle %s sont des doublons consecutifs.'%avirer
#         msg1 += 'Suppression des doublons'
#         rdebug(msg1,'')
#                 rdebug(avirer=avirer, delta_abscurv=dac.tolist(), controle=sum(dac))
#         for k in reversed(avirer) :
#             self.removePoint(k, update=False)

def projSurProfil(prof,iba,pos,cote):
    ''' Calcul des caractéristiques d'un point projeté sur un profil (à l'échelle ou non) à partir :
            de sa position (% corde)
            de son côté (intrados/extrados)
        Entrées :
            prof : profil en entrée
            iba : indice du point de BA
            pos : position du point en % de corde à partir du BA
            cote : 0:intrados, 1:extrados
        Sorties :
            i,t : indice point profil et coordonnée barycentrique tels que
                  pt se trouve à t entre les points profil i et i+1
            pt : coordonnées du point projeté

        On applique la méthode de updateAncrages (cloison).

        Le profil est supposé horizontal et à l'échelle : BA=(0,0) et BF=(corde,0)
        On recherche le segment [P(i), P(i+1)] dans lequel se trouve le point projeté :
            A = P(i) + t.(P(i+1) - P(i)) ce qui donne deux équations :
            a = x(i) + t.(x(i+1) - x(i))
                d'ou on tire t = (a-x(i))/(x(i+1)-x(i)) que l'on reporte dans la deuxième équation :
            b = y(i) + t.(y(i+1) - y(i))
    '''
    corde = prof[0][0]-prof[iba][0]
    ''' abscisse du point '''
    a = prof[iba][0]+(pos/100.0)*corde
    ptout = [a, 0.0]

    if cote==0:#intrados
        i=iba
        while i<len(prof)-1 and a>prof[i,0] :
            i+=1
        i=i-1
        t=(a-prof[i,0])/(prof[i+1,0]-prof[i,0])
        ptout[1]=prof[i,1]+t*(prof[i+1,1]-prof[i,1])
    elif cote==1:#extrados
        i=0
        while i<iba-1 and a<prof[i,0] :
            i+=1
        i=i-1
        t=(a-prof[i,0])/(prof[i+1,0]-prof[i,0])
        ptout[1]=prof[i,1]+t*(prof[i+1,1]-prof[i,1])
    elif 0<cote<1:
        ''' Point intrados '''
        ptoutInt = [a, 0.0]
        i=iba
        while i<len(prof)-1 and a>prof[i,0] :
            i+=1
        i=i-1
        t=(a-prof[i,0])/(prof[i+1,0]-prof[i,0])
        ptoutInt[1]=prof[i,1]+t*(prof[i+1,1]-prof[i,1])
        ''' Point extrados '''
        ptoutExt = [a, 0.0]
        i=0
        while i<iba-1 and a<prof[i,0] :
            i+=1
        i=i-1
        t=(a-prof[i,0])/(prof[i+1,0]-prof[i,0])
        ptoutExt[1]=prof[i,1]+t*(prof[i+1,1]-prof[i,1])
        ''' Point milieu '''
        ptout = [ptoutInt[0]+cote*(ptoutExt[0]-ptoutInt[0]), ptoutInt[1]+cote*(ptoutExt[1]-ptoutInt[1])]

    return i,t,ptout

def vect2d(u, v):
    u"""
    Retourne un réel, produit vectoriel de deux vecteurs (2d), u et v
    C'est aussi le déterminant
    """
    return u[0]*v[1] - u[1]*v[0]

det = vect2d

def rcercle(A, B, C, eps=0.0):
    u"""retourne le rayon dun cercle passant par les 3 points A,B,C, distincts et non alignes.
    Si les 3 points sont presque alignés retourne inf
    si deux points sont presque confondus, retourne nan"""
    A, B, C = asarray(A), asarray(B), asarray(C)
    AB, BC, CA = B-A, C-B, A-C
#     print AB, CA
    c, a, b = sqrt(AB[0]**2 + AB[1]**2), sqrt(BC[0]**2 + BC[1]**2), sqrt(CA[0]**2 + CA[1]**2),
    abc = a*b*c
    d = det(AB,CA)
    s = abs(d)
#     print s, abc
    if abc <= eps:#2 points presque confondus : infinité de cercles passant par deux points
        return nan
    elif s <= eps:# trois points presque alignés
        return inf
    else :
        return 0.5*abc/d

def rayonDeCourbure(P):
    u"""
    Retourne la liste des rayons de courbure du polyGONE P.
    S'il n'est pas fermé, on le ferme.
    Il y a exactement n=len(P) rayons
    """
    rayons = zeros(len(P))
    if norm(P[0]-P[-1])>1.0e-9 :#debut != fin
        rayons[ 0] = rcercle(P[-1], P[0], P[1])
        rayons[-1] = rcercle(P[-2], P[-1], P[0])
    else :
        rayons[0] = rayons[-1] = rcercle(P[-2], P[0], P[1])
    ABC = zip(P[0:-2], P[1:-1], P[2:])
    for k, (A, B, C) in enumerate(ABC) :
        rayons[1+k] = rcercle(A, B, C)
#     print 'dernier', k, A, B, C, 1.0/rcercle(A, B, C)
    return rayons

def courbure(P):
    u"""
    Retourne la liste des courbures du polyGONE P.
    S'il n'est pas fermé, on le ferme.
    Il y a exactement n=len(P) valeurs de la courbure
    """
    courbature = zeros(len(P))
    if norm(P[0]-P[-1])>0 :#debut != fin
        courbature[ 0] = 1.0/rcercle(P[-1], P[0], P[1])
        courbature[-1] = 1.0/rcercle(P[-2], P[-1], P[0])
    else :
        courbature[0] = courbature[-1] = 1.0/rcercle(P[-2], P[0], P[1])
    ABC = zip(P[:-2], P[1:-1], P[2:])
    for k, (A, B, C) in enumerate(ABC) :
        courbature[1+k] = 1.0/rcercle(A, B, C)
    return courbature

def scourbure(S, T):
    u"""
    ATTENTION, donne des resultats fantaisistes sur un cercle.????? a vérifier ????
    utiliser plutôt la version discrete courbure() ??????
    Parametres:
    ----------
    - S = (sx, sy) : 0<=t<=1 --> S(t)=(sx(t), sy(t) est une spline numpy,
    - T = [t0, t1, ...tn] : les valeurs du parametre t pour lesquels on calcule la courbure.
    retourne:
    --------
        un ndarray((n,)) avec les valeurs de la courbure c(ti) aux pointx sx(ti), sy(ti), ti dans T
        La courbure est l'inverse du rayon de courbure.
        Pour un arc x(t), y(t), le rayon de courbure est r(t)=((x'^2+y'^2)^3/2)/(x'y"-y'x")
        x'=x'(t), y'=y'(t) x"=x"(t), y"=y"(t).
        si ||(x",y")|| = 0, la courbure est nulle.
        cf https://fr.wikipedia.org/wiki/Rayon_de_courbure
    """
    sx, sy = S
    dx,  dy  = sx(T, 1), sy(T, 1)
    d2x, d2y = sx(T, 2), sy(T, 2)
    norm3_d2 = sqrt(dx**2+dy**2)**3
    # si norm_d2=0, x"(t)=y"(t)=0, c'est une droite, courbure nulle
    sc = (dx*d2y-dy*d2x)/(norm3_d2)
    sc[where(norm3_d2 < 1.0e-12)] = 0.0
    return sc

def simpson(f, a, b, n=10):#n doit être pair, integration precise ordre 3
    u"""
    Integrale de f sur [a,b], méthode de Simpson composite. (ordre 3)
    n DOIT être pair, on découpe le segment d'integration en n sous segments
    On applique la méthode de Simpson sur chaque sous segment"""
    h = float(b-a)/n
    T = linspace(a, b, n+1)
    C = f(T)
    A1 = C[0] + C[-1]
    A2 = 2*sum(C[i] for i in range(2,n) if i%2==0)
    A4 = 4*sum(C[i] for i in range(1,n) if i%2==1)
#         debug (h, A1, A2, A4, (h/3)*(A1 + A2 + A4))
    return (h/3)*(A1 + A2 + A4)

def dictsAreNotEqual(d1,d2):
    d1 = OrderedDict(d1)
    d2 = OrderedDict(d2)
    if not d1.keys() == d2.keys() :
        return 'diff :\n    keys1 = %s,\n    keys2 = %s'%(d1.keys(),d2.keys())
    for (k1,v1), (k2,v2) in zip(d1.items(),d2.items()) :
        if not type(v1) == type(v2) :
            return 'diff :\n    %s : %s\n    %s : %s'%(k1,className(v1),k2,className(v2))
        if isinstance(v1, (dict,OrderedDict)) :
            return dictsAreNotEqual(v1, v2)
        elif not v1==v2 :
            return 'diff :\n    %s : %s\n    %s : %s'%(k1,v1,k2,v2)
#     debug('equal')
    return False
if __name__=="__main__":
    dd1 = dict(aa=10, bb=1000)
    dd2 = dict(aa=10, bb=100)
    d1 = dict(a=1,b=2,c=3,d=dd1)
    d2 = dict(a=1,b=2,c=3,d=dd2)
    print dictsAreNotEqual(d1, d2)
#     P = asarray([[0.,0.],[1.,2.],[2.,3.],[3.,4.]])
#     print P
#     RP = rotate(P,radians(90), centre=[0,0], keep=[2,3])
#     print RP
#     from matplotlib import pyplot as plt
#     plt.plot(P[:,0],P[:,1],'r-o')
#     plt.plot(RP[:,0],RP[:,1],'b-o')
#     plt.show()
