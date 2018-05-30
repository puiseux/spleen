#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Created on 11 mai 2012

@author: puiseux
'''
import datetime
import gc  # garbage colector
# from PyQt4 import QtGui,QtCore
# from PyQt4.QtCore import (Qt,QSize,QVariant,QPointF,QPoint,)
# from PyQt4.QtGui import (QTableWidget,QTableWidgetItem,QPushButton,QLayout,
#                          QVBoxLayout,QHBoxLayout,QGridLayout,QSpacerItem,QSizePolicy,
#                          QWidget,QScrollArea,QApplication,QKeySequence,
#                          QShortcut,QFont,QFontMetrics,QHeaderView, QPolygonF, QPolygon)
# from PyQt4.Qt import SIGNAL,SLOT
# from inout.format import formData
#from gui.graphicsbase.graphicscommon import pointsFromPolygon
import math
import pickle
import string
import sys
# from exceptions import *
from math import acos, atan, cos, sin, tan
from pprint import pprint
from random import choice

import numpy as np
import scipy as sp
from numpy import _distributor_init, asarray, ones, zeros
from numpy.matlib import rand
from path import Path
from scipy import interpolate, ones, prod, sum
from scipy.interpolate import (CubicSpline, InterpolatedUnivariateSpline,
                               LSQUnivariateSpline, Rbf, UnivariateSpline,
                               interp1d, splev, splrep)
# from config import SOURCES_DIR
# from gui.graphicsbase.graphicscommon import qpolygonFrom
# import logging
# DEBOG=True
# VERBOSITY_LEVEL=2
# DEBOG_OUTPUT=sys.stderr
# DEBOG_OUTPUT=sys.stdout
# output=DEBOG_OUTPUT
from scipy.optimize import minimize_scalar
# from model.basicobjects.splineInterpolation import splineInterpolation
# from model.decoration.contrainte import Contrainte
from shapely.geometry import LinearRing, Point, Polygon

# import logging
# logger = logging.getLogger('')
# logger.setLevel(logging.DEBUG)

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
            points = np.vstack((points, points[0]))
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
#     trace(None, methode=methode, tension=tension, degre=degre)
    if methode in ('ius','interpolatedunivariatespline') :
        try :
            sx = InterpolatedUnivariateSpline(T, X, k=degre)#s=la précision de l'ajustement s=0 <=> interpolation
            sy = InterpolatedUnivariateSpline(T, Y, k=degre)
        except Exception as msg:
            trace(None)
            print unicode (msg)
#             trace(None,u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
            sx = sy = None
    elif methode in ('us','univariatespline') :
        try :
            weights = np.ones(N)
            W = 1000.0
            # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
            # le choix de s suivant implique
            # abs(xi-f(ti))<eps et
            # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
#             eps = 10.0**(-tension)
            weights[0] = weights[-1] = W
            weights /= np.sum(weights)
            s = eps/(N*W)
            sx = UnivariateSpline(T, X, w=weights, k=degre, s=s)#s=la précision de l'ajustement s=0 <=> interpolation
            sy = UnivariateSpline(T, Y, w=weights, k=degre, s=s)
        except Exception as msg:
            trace(None)
            print unicode(msg)
#             trace(None,u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
            sx = sy = None
    elif methode in ('interp1d',) :
        try :
            sx = interp1d(T, X, kind=degre)
            sy = interp1d(T, Y, kind=degre)
        except ValueError as msg:
            trace(None)
            print unicode(msg)
            sx = sy = None
#     elif methode in ('periodic',) :
#         try :
#             sx = PeriodicSpline(T, X, k=degre, s=eps)
#             sy = PeriodicSpline(T, Y, k=degre, s=eps)
#         except ValueError as msg:
#             trace(None)
#             print unicode(msg)
#             sx = sy = None
    elif 'cubic' in methode :#or isinstance(methode, (tuple, list, np.ndarray)):
        if methode == 'p cubic' : bc_type='periodic'
        elif methode == 'c cubic' : bc_type='clamped'
        elif methode == 'n cubic' : bc_type='natural'
        else : bc_type = 'not-a-knot'

        try :
#             trace(None, T)
            sx = CubicSpline(T, X, bc_type=bc_type)
            sy = CubicSpline(T, Y, bc_type=bc_type)
        except ValueError as msg:
            print unicode(msg)
            sx = sy = None

    elif isinstance(methode, (tuple, list, np.ndarray)):
        bc_type = methode
        try :
#             trace(None, T)
            sx = CubicSpline(T, X, bc_type=bc_type)
            sy = CubicSpline(T, Y, bc_type=bc_type)
        except ValueError as msg:
            print unicode(msg)
            sx = sy = None
    return T, sx, sy

def rotate(points, alfa, centre):
    u'''
    alfa en radians
     Retourne une COPIE de points rotationnée(!).
    points est supposé stocké par ligne (shape=(n,2)), chaque point est de shape (1,2),
    il les faudrait en colonne (shape=(2,1)) pour faire le produit matriciel.
    Donc on transpose tout et on ecrit Xi' = C' + (Xi'-C')*A' au lieu de
    Xi = C + A*(Xi-C), pour i= 0, 1,...
    '''
    Ct = np.asarray(centre).reshape((1,2))
    cosa, sina = np.cos(alfa), np.sin(alfa)
    At = np.matrix([[cosa,-sina], [sina,cosa]]).transpose()
    Xt = points - Ct
    Xt = Xt*At + Ct
#     trace('', Xt.shape)
    return np.asarray(Xt)

def symetriser(points, sym):
    u"""
    Symetrie in situ de 'points' par rapport à un des plan Oxy, Oyz ou Oxz
    Parametres:
    ----------
        - points = np.ndarray de shape (m,n,...,3) ou (m,n,...,2)
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
#     centre=np.asarray(centre).reshape((1,2))
    points *= np.asarray([echelle[0],echelle[1]])
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
    points=np.copy(points)#sinon effet de bord difficile à identifier !!!
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
    molecule=np.asarray(molecule)#np.asarray([1., 2., 1.])/4.
    molecule/=np.sum(np.absolute(molecule))
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
        if surface : return np.asarray([np.nan, np.nan]), 0.0
        else : return np.asarray([np.nan, np.nan])
    elif len(points) == 1 :
        if surface : return np.asarray(points[0]), 0.0
        else : return np.asarray(points[0])

    points=np.asarray(points)
    Sa = 0.0
    xG, yG = 0.0, 0.0
    G = np.asarray([xG, yG])
    A = points[0]
    T = zip(points[:-1], points[1:])+[(points[-1],points[0])]
    for b,c in T :
#    for b, c in zip(P[1:-2], P[2:-1]) :?????
        Gt = (A + b + c)/3.0
        Sat = (b[0] - A[0])*(c[1] - A[1]) - (b[1] - A[1])*(c[0] - A[0])
        G += Sat*Gt
#         trace(self, "pouet, a faire")
        Sa += Sat
    if Sa == 0.0 :
        if surface :  return np.asarray((np.nan, np.nan)), 0
        else : return np.asarray((np.nan, np.nan))
    else :
        if surface : return G / Sa, Sa*0.5
        else : return G / Sa

def aire(points):
    '''
    Calcul de l'aire algébrique d'un polygone.
    Si le polygone ne se recoupe pas, (intérieur connexe), alors l'aire donne le sens de rotation :
    si elle est positive, ses points tournent dans le sens trigonométrique, sinon, le sens des aiguilles d'une montre.
    La méthode utilisée est
    - fermer le polygone : P(n)=P(0)
    - de fixer un point 'O' quelconque, ici on prend O=P[0]
    - sommer les aires des triangles A = som(A(ti)) 0<=i<n) où
    A(ti) = (1/2).OP(i)^OP(i+1) est l'aire algébrique du triangle ti=OP(i)P(i+1)
    '''
    if len(points) <= 2 : return 0.0
    points=np.asarray(points)
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
    points=np.asarray(points)
#     trace('', points.shape, len(points))
    N=len(points)
    if masses is None :
        X,Y = points[:,0], points[:,1]
#         debug (X=X, Y=Y)
    else :
        masses=masses/np.sum(masses)
        X,Y=np.prod([masses,points[:,0]],axis=0),np.prod([masses,points[:,1]],axis=0)
    bc = [[np.sum(X)/N, np.sum(Y)/N]]
    return np.asarray(bc)


def dist2(p1,p2,n=2):
    '''retourne le carré de la distance de p1 à p2 en norme n=2'''
    try :
        x1,y1=p1[0],p1[1]
        x2,y2=p2[0],p2[1]
        return (x2-x1)**2+(y2-y1)**2
    except TypeError : # Si p1 ou p2=None
        return np.nan

def dist(p1,p2,n=2):
    '''retourne la distance de p1 à p2 en norme n=2'''
    return math.sqrt(dist2(p1,p2))

'''vectors 3D '''
def length3(v):
    x,y,z = v
    return math.sqrt(x*x + y*y + z*z)
def vector3(b,e):
    x,y,z = b
    X,Y,Z = e
    return (X-x, Y-y, Z-z)
def unit3(v):
    x,y,z = v
    mag = length3(v)
    if mag==0:
        return v
    else:
        return (x/mag, y/mag, z/mag)
def distance3(p0,p1):
    return length3(vector3(p0,p1))
def scale3(v,sc):
    x,y,z = v
    return (x * sc, y * sc, z*sc)
def add3(v,w):
    x,y,z = v
    X,Y,Z = w
    return (x+X, y+Y, z+Z)

def rot3(pt,u,theta):
    ''' 3D rotation of point pt around axis u of angle theta '''
    matR=zeros((3,3),dtype=float)
    ux,uy,uz=(u[0],u[1],u[2])
    matP=[[ux**2.0,ux*uy ,ux*uz],[ux*uy,uy**2.0,uy*uz],[ux*uz,uy*uz,uz**2.0]]
    matI=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    matQ=[[0.0,-uz,uy],[uz,0.0,-ux],[-uy,ux,0.0]]
    for i in range(3):
        for j in range(3):
            matR[i,j]=matP[i][j]+(matI[i][j]-matP[i][j])*cos(theta)+matQ[i][j]*sin(theta)

    rp=[float() for i in range(3)]
    for i in range(3):
        s=0.0
        for j in range(3):
            s=s+matR[i,j]*pt[j]
        rp[i]=s

    return rp

def vect2d(u, v):
    u"""
    Retourne un réel, produit vectoriel de deux vecteurs (2d), u et v
    C'est aussi le déterminant
    """
    return u[0]*v[1] - u[1]*v[0]
det = vect2d
# def sensDeRotaion(points):
#     u"""détermine le sens de rotation d'un polygone qui de se croise pas. C'est le signe de la surface algebrique """
#     pass


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
    return 0.5*np.sqrt(r*r+s*s+t*t)

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
        return np.nan,np.nan,np.asarray((np.nan,np.nan))
    else  :
        r=(w[0]*v[1]-w[1]*v[0])/denominateur
        s=(w[0]*u[1]-w[1]*u[0])/denominateur
        return r,s,A+r*u

def segmentPlusProche(points,P):
        u'''
        Paramètres:
        ----------
            - points : np.ndarray((n,2)) est un polyligne dans lequel on cherche le segment
            - P : np.ndarray((1,2)) est un point quelconque
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
        points = np.asarray(points)
        u = -(points-P)[:-1]
        v = points[1:]-points[:-1]
#         debug(v=v)
#         distances = [dist(point,(X,Y)) for point in points]
#         debug(points-P)
#         debug(distances=distances)
#         debug(P, len(points), distances)
#         distances = np.linalg.norm(u, axis=1)
        longueurs = np.linalg.norm(v, axis=1)
#         debug(distances_P_points=distances)
#         debug(longueurs_segments=longueurs)
        ps = [np.dot(ui, vi) for (ui, vi) in zip(u,v)]
        r = ps/(longueurs*longueurs)
#         debug(r_entre_0_et_1=r)
#         psn = [np.dot(ui, vi)/np.dot(vi,vi) for (ui, vi) in zip(u,v)]
#         debug(psn=psn)
#         distances = r*longueurs
#         debug(distances_P_Segment=distances)
        candidats, = np.where(np.logical_and(0.0<=r,r<=1.0))
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
            winner = np.argmin(distances)
#             debug(winner=winner)
#             debug(candidats_winner=candidats[winner])
            i = candidats[winner]
            projete = projetes[winner]
#         debug(i=i, projete=projete)
            return i, projete

def pointLePlusProche(points,P,return_distances_array=False):
    '''
    Parametres:
    ----------

    :param points: np.ndarray, tableau de n points 2d de shape (n,2)
    :param P: un point 2d de type QPointF, QPoint ou tuple ou liste (x,y)
    :param return_distances_array: bool, retourne ou non le tableau des distances.

    :return: (i, dist) avec

        - i : l'indice dans 'points' du point qui réalise le min
        - dist : la distance de P à points si return_distances_array est False
        - dist : le tableau des distances de P à points si return_distances_array est True
    '''
    X,Y=P[0],P[1]
    distances=[dist(point,(X,Y)) for point in points]
    index=np.argmin(distances)
    if return_distances_array :
        return index,distances
    else :
        return index,distances[index]

def pointLePlusProche3d(points,P,return_distances_array=False):
    '''
    Parametres:
    ----------
    - points : np.ndarray, tableau de n points 3d de shape (n,3)
    - P : un point 3d de type QPointF, QPoint ou tuple ou liste (x,y)
    retourne l'indice dans points, du point le plus proche de P
    Retourne (i, distances)
    --------
    - i : l'indice dans 'points' du point qui réalise le min
    - distance : la distances de P à points
    '''
    X,Y,Z=P[0],P[1],P[2]
    distances=[distance3(point,(X,Y,Z)) for point in points]
    index=np.argmin(distances)
    if return_distances_array :
        return index,distances
    else :
        return index,distances[index]

def pointsLesPlusProches3d(pt,tab):
    ''' Recherche de l'indice i tel que tab(i) et tab(i+1) soient les deux points les plus proches de pt
        - on cherche le point le plus proche, d'indice i
        - on détermine lequel des points i+1 ou i-1 est le plus proche de pt
    '''
    i1,dpt1 = pointLePlusProche3d(tab,pt,True)
    if (0):
        print 'pointsLesPlusProches3d:pt',pt
        print 'pointsLesPlusProches3d:tab',tab
        print 'pointsLesPlusProches3d:i1,dpt1',i1,dpt1
    if (0<i1<(len(tab)-1)):
        dpt1p1 = distance3(pt,tab[i1+1])
        dpt1m1 = distance3(pt,tab[i1-1])
        if (dpt1p1 <= dpt1m1):
            i = i1
        else:
            i = i1-1
    elif  i1==0:
        i = i1
    elif i1==len(tab)-1:
        i = i1-1

    return i

def longueur2d(polyline):
    P = np.asarray(polyline)
    return sum([dist(p1,p2) for (p1, p2) in zip(P[1:], P[:-1])])

def longueur2d_array(tab):
    return sum([dist(p1,p2) for (p1, p2) in zip(tab[1:], tab[:-1])])

def longueur3d(tab):
    return sum([distance3(p1,p2) for (p1, p2) in zip(tab[1:], tab[:-1])])

def encombrement(piece, dim=3):
    u"""
    piece doit être un np.ndarray de shape (N, dim) quelconque,
    retourne le parallelepipede d'encombrement du nuage de points
    """
    if isinstance(piece, (list, tuple, )) :
        return encombrement(np.asarray(piece),dim)#recursif
    elif isinstance(piece,(np.ndarray,)):
#         trace("", piece_shape=piece.shape)
        points = piece
        Max, Min = np.max, np.min #ca reste local
        if dim == 1 :
            return Min(points), Max(points)
        points = piece.view()
        points.shape = -1, dim
        if dim == 3 :
            M = [Max(points[:,0]),Max(points[:,1]),Max(points[:,2])]
            m = [Min(points[:,0]),Min(points[:,1]),Min(points[:,2])]
#             M = np.asarray([Max(points[:,0]),Max(points[:,1]),Max(points[:,2])])
#             m = np.asarray([Min(points[:,0]),Min(points[:,1]),Min(points[:,2])])
            return(m,M)
        elif dim == 2 :
            M = [Max(points[:,0]),Max(points[:,1])]
            m = [Min(points[:,0]),Min(points[:,1])]
            return(m,M)
    else :
        raise NotImplementedError

def encombrement1(piece, dim=3):
    u"""retourne l'encombrement de la piece sous forme dx[,dy[,dz]]"""
    m, M = encombrement(piece, dim)
    return np.asarray(M) - np.asarray(m)


def maintenant():
    u'''La date et l'heure formatées "human readable"'''
    return unicode(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))

def doublons(points,voisinages=[],eps=1.0e-8,sym=True):
    u"""
    Verifier qu'il n'y a pas de point double.
    Paramètres :
    ----------
    - points : ndarray(nbp,3) nbp points 3d. nbp peut être une
        shape (nl, nc, 3) et dans ce cas nbp = nl*nc
    - voisinages ndarray(shape=(np,1+nv), dtype=int) np numéros de points,
        faisant référence au tableau 'points', chacun des np points à
        vérifier possède nv voisins. voisinage est donc de la forme
        [
         [k_0,  v_0_1,  v_0_2, ...,  v_0_nv],
         [k_1,  v_1_1,  v_1_2, ...,  v_1_nv],
         ...
         [k_np, v_np_1, v_np_2, ..., v_np_nv]
        ]
        les valeurs k_i, v_i_j sont des entiers entre 0 et nbp-1
    - eps : réel. deux points p,q sont considérée comme confondus si dist(p,q)<eps
    - sym : booleen. True si i voisin de j <==> j voisin de i. Divise par deux
        le travail.

    Retourne :
    --------
    une liste de doublons (d_0,dd_0), ... (d_n, dd_n)
    """
    def cond(v,k):
        if sym : return v>k
        else : return v!=k

    qoints=points.view()
    qoints.shape=(qoints.size/3,3)
    nbp=qoints.shape[0]
    if voisinages!=[] :
        voisinages=np.asarray(voisinages)
#        debog(whoami(), voisinages.shape)
        nbp=voisinages.shape[-1]
        voisinages.shape=(-1,nbp)
    else :
        if sym :
            voisinages=np.asarray([[i]+range(i,nbp) for i in range(nbp)])
        else :
            voisinages=np.asarray([[i]+range(nbp) for i in range(nbp)])

    doublons=[]
    for voisinage in voisinages:
        k,voisins=voisinage[0],voisinage[1:]
        point=qoints[k]
        for voisin in voisins :
            if cond(voisin,k) :
                if np.linalg.norm(point-qoints[voisin])<eps :
                    doublons.append((k,voisin))
    return doublons

u"""
Quelques fonctions de débogage qui ont évolué au cours du temps.

Ecriture dans le fichier log 'AXile.log'
----------------------------------------

1- trace et alert

    >>> trace(x, y=yzt) écrit dans le fichier log niveau INFO
    >>> alert(x, y=yzt) écrit dans le fichier log niveau WARNING

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

    >>> strace(x, y='toto') écrit dans le fichier log niveau INFO
    >>> salert(x, y='toto') écrit dans le fichier log niveau WARNING

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
    # if output is sys.stdout :
    #     logger.info(msg0 + u' ; '.join(lmsg))
    # elif output is sys.stderr :
    #     logger.warning(msg0+u' ; '.join(lmsg))
#     try : print>>output, msg0, u' ; '.join(lmsg)
#     except UnicodeEncodeError: print msg0, lmsg

def rstack(commentaire=''):
    _stack(sys.stderr, commentaire)
def stack(commentaire='') :
    _stack(sys.stdout, commentaire)

def _stackOld(output, commentaire=''):
    u"""Impression de la pile des appels sur output"""
    print>>output,u"========================"
    print>>output,u"=== Pile des appels ===="
#     k = 2
#     frame = sys._getframe(2)
#     fcode = frame.f_code
#     fonction = fcode.co_name
#     filename = Path(fcode.co_filename)#.name
#     toclick = clickableLink(filename, frame.f_lineno)#, fcode.co_name)
#
#     print>>output,u'[%s] '%fonction + toclick
    k = 2
    if commentaire :
        print>>output, commentaire
    while 1 :
        try :
            frame = sys._getframe(k)
            filename = Path(frame.f_code.co_filename).name
            k += 1
            print>>output, u" ['%s':%d, (%s)] "%(filename, frame.f_lineno, frame.f_code.co_name)
        except :
            break
    print>>output,u"======  Fin pile  ======"
    print>>output,u"========================"
    return

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

def _strace(*args,**kargs) :
    u"""_trace simplifié, sans le nom de la classe"""
    frame = sys._getframe(2)
    fcode = frame.f_code
    fonction = fcode.co_name
    filename = Path(fcode.co_filename)#.name
#     filename = u'./source/gui/'+filename.name
#     msg0 = u"\n'+unicode(u" [%s:%d, (%s)] '%(filename, frame.f_lineno, fcode.co_name))#+u" : \n    '
#     msg0 = unicode(u'File "%s", line %d, (%s)]'%(filename, frame.f_lineno, fcode.co_name))#+u" : \n    '
#     msg0 = unicode(u'File "%s", line %d\n[%s]'%(filename, frame.f_lineno, fcode.co_name))#+u" : \n    '
    toclick = clickableLink(filename, frame.f_lineno)#, fcode.co_name)
    output = args[0]
    args = args[1:]
    lmsg = [unicode(arg) for arg in args]+\
           [unicode(key)+u" = "+unicode(value) for key,value in kargs.iteritems()]
    msg = u'[%s] '%fonction + toclick + u'\n    ' + u" ; ".join(lmsg)
    # if output in (logger.info, logger.warning):
    #     output(msg)
    # else :
    print>>output, msg
#     try : print>>output, msg0, u" ; '.join(lmsg)
#     except UnicodeEncodeError: print msg0, lmsg
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
    n0=_findRel (tag,lines,n+1,last)
    N=[]
    while n0>=0 :
        n=n0+n+1
        N.append(n)
        n0=_findRel(tag,lines,n+1,last)
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

def _findRel(tag,lines,first=0,last=-1):
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
    found=_findRel(tag,lines,first=first,last=last)
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
    retourne un tableau np.ndarray de shape(nx,ny,nz,3)

    produit tensoriel de X=(x_1,x_2,...x_nx), Y=(y_1,y_2,...y_ny)
    et Z=(z_1,z_2,...z_nz), aléatoirement espacés, mais triés.
    ou bien régulièrement espacés
    Contrainte : nx>1, ny>1, nz>1
    """
#    nx,ny,nz = shape
    if delta==-1 :
        X=np.random.random(nx)
        Y=np.random.random(ny)
        Z=np.random.random(nz)
        X.sort()
        Y.sort()
        Z.sort()
    else :#TODO: delta=dx,dy,dz
        dx,dy=1.0/(nx-1),1.0/(ny-1)
        X=np.arange(0.0,1.0+dx,dx)
        Y=np.arange(0.0,1.0+dy,dy)
    if nz==1 :
        Z=np.asarray([0.0])
    else :
        dz=1.0/(nz-1)
        Z=np.arange(0.0,1.0+dz,dz)
    G=np.ndarray((nx,ny,nz,3),dtype=float)
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

    import math,pylab
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
#        print xy.shape
        try: color=colors[k]
        except : color='b-'
        pyplot.plot(xy[:,0],xy[:,1],color)
        xmin=min(xmin,np.min(xy[:,0]))
        ymin=min(ymin,np.min(xy[:,1]))
        xmax=max(xmax,np.max(xy[:,0]))
        ymax=max(ymax,np.max(xy[:,1]))
    w,h=xmax-xmin,ymax-ymin
    dizaine=int(math.log10(w))#10^dizaine <= w < 10^(1+dizaine)
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
        for k, legend in enumerate(legends) :
            pylab.setp(ltext[k],fontsize=10)#, color = 'b')

#     print ltext
    pylab.title(title)
    pyplot.show()

def load(filinname):
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

# def testContraindre():
#     import matplotlib
#     P=np.asarray([(3.0,-1.0),(2,-1),(0,0),(2,2),(3,2),(5,2),(10,0)])
#     P=P[::-1]
#     S=[((0.5,0.25),(0.5,-1)),((1.5,3),(1.5,-2)),((2.5,3),(2.5,-2)),((4,3),(4,-2)),((4.5,3),(4.5,-2)),((5,3),(5.01,-1))]
# #     S.reverse()
# #     print PC[0].shape
# #     print PC
#     Sp=[]
#     for A,B in S : Sp+=[A,B]
#     Sp=np.asarray(Sp)
#     S1=dupliquerContraintes(S,3)
#     S1p=[]
#     for A,B in S1 : S1p+=[A,B]
#     S1p=np.asarray(S1p)
# #     my2dPlot((Sp,S1p, ), equal=False, cosmetic=('r-o','b-o', 'g-o', 'g*', 'r*','r^'))
# #     return
#     PC=contraindre(P,S)
#     PC2=contraindre(P,S,3)
#     PC1=contraindre(P,S1)
#     my2dPlot((P,PC[0],PC1[0],PC2[0]),equal=False,cosmetic=('r-o','b-o','g-*','r-^','r*','r^'))
# #     P = np.asarray([(0, -1.0),(-1,0),(1,0),(0,1)])

def diff(A):
    u"""
    J'en ai marre de taper tout le temps 'dA = A[1:]-A[:-1]'
    :param A: un tableau de points np.ndarray de shape (n, dim)
    :return d:les différences premières de A:
        dA[i] = A[i+1]-A[i]
    :rtype dA: np.ndarray((n-1,dim))
    """
    return A[1:]-A[:-1]


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


def isInside(p, Q):
    """
    Paramètres:
    ==========
        - p = (x,y) = np.ndarray((1,2), float) un point en 2d
        - Q est un polygone = np.ndarray((n,2), float) supposé fermé i.e. le premier point == le dernier point,
    Retourne :
    ========
        - vrai si p est à l'intérieur du polygone Q et
        - faux s'il est à l'extérieur ou frontiere
    Fonctionnement :
    ==============
    p est intérieur, si et seulement si les produits vectoriels AB^Ap ont tous le même signe.
    (on note A=Q[i], B=Q[i+1], 0<= i < n )
    """
#    trace('', p,Q)
    (x, y) = p
#    trace('', x,y)
    w0 = vect2d(Q[1]-Q[0], p-Q[0])
    for (A, B) in zip(Q[:-1], Q[1:]) :
        #[A,B] parcourt les segments de Q
        w =  vect2d(B-A, p-A)
#        trace('', w0, w)
        if w*w0 <= 0 : return False

    return True

def isInside0(p, Q, frontiere=False):
    """
    Paramètres:
    ==========
        - p = (x,y) = np.ndarray((1,2), float) un point en 2d
        - Q est un quadrilatère = np.ndarray((5,2), float) avec deux cotés parallèles à l'axe des y
            Q est supposé fermé i.e. le premier point == le dernier point,
            on n'utilise que les points 0 à 3 notés A, B, C, D
        - ATTENTION à l'ordre des points : cf dessin ci-dessous.
          A x¨ ¨ - - _ _
            |            ¨ ¨ -x B
            |                 |
            |      . p        |
            |                 |
            |         _ _ --  x C
          D x - - ¨ ¨
        - TODO : frontiere=False pour indiquer qu'un point frontière est considéré comme extérieur.
    Retourne :
    ========
        - vrai si p est à l'intérieur du quadrilatere Q = ABCDA et
        - faux s'il est à l'extérieur.
    Fonctionnement :
    ==============
    pour être intérieur, p doit vérifier les conditions suivantes :
        - xA < x <xB
        - les deux produits vectoriels AB^Ap et CD^Cp  ont même signe.
    """
#    trace('', p,Q)
    (x, y) = p
#    trace('', x,y)
    A, B, C, D = Q[0], Q[1], Q[2], Q[3]
#    trace('', A,B,C,D)
    if x <= A[0] or x >= B[0] :
#        print '****x'
        return False

    wa = vect2d(B-A, p-A)
    wf = vect2d(D-C, p-C)
#    print wa, wf
    if wa*wf > 0 : return True# meme signe  => point interieur
    return False

def isOn(P, segment, eps=1.0e-10, debug=False):
    u'''
    Retourne True si P est sur le segment, à eps près
    Paramètres :
    ----------
    - segment : [A,B] ou A et B de type np.ndarray((1,2),dtype=float)
    - P : de type np.ndarray((1,2),dtype=float)
    - eps : de type float >= 0.0
    Subject 1.02: How do I find the distance from a point to a line?
    ===============================================================
    Projection orthogonale d'un point sur un segment

    Let the point be C (Cx,Cy) and the line be AB (Ax,Ay) to (Bx,By).
    Let P be the point of perpendicular projection of C on AB.  The parameter
    r, which indicates P's position along AB, is computed by the dot product
    of AC and AB divided by the square of the length of AB:

    (1)     AC dot AB
        r = ---------
            ||AB||^2

    r has the following meaning:

        r=0      P = A
        r=1      P = B
        r<0      P is on the backward extension of AB
        r>1      P is on the forward extension of AB
        0<r<1    P is interior to AB

    The length of a line segment AB is computed by:
        L = sqrt( (Bx-Ax)^2 + (By-Ay)^2 )
    and the dot product of two vectors , U dot V is computed:
        D = (Ux * Vx) + (Uy * Vy)
    So (1) expands to:
            (Cx-Ax)(Bx-Ax) + (Cy-Ay)(By-Ay)
        r = -------------------------------
                          L^2

    The point P can then be found:

        Px = Ax + r(Bx-Ax)
        Py = Ay + r(By-Ay)
    '''
    A, B = segment[0], segment[1]
    x, y = P - A
    u, v = B - A
    det_, dot_ = abs(u*y - v*x), (x*u + y*v)/(u*u + v*v)
    val = det_ < eps and -eps <= dot_ <= 1.0 + eps
#    if debug : return val, det_, dot_
#    else :
    return val

def isInsideOrFrontier(p, Q, eps=1.0e-6):
    u"""
    Paramètres:
    ==========
        - p = (x,y) = np.ndarray((1,2), float) un point en 2d
        - Q est un polygone = np.ndarray((n,2), float) supposé fermé i.e. le premier point == le dernier point,
        - eps : precision
    Retourne :
    ========
        - vrai si p est à l'intérieur du polygone Q ou sur la frontiere (bande de largeur eps) et
        - faux s'il est à l'extérieur.
    """
    P = Polygon(Q)
    p = Point(p)
    #trace('', p=Point(p), lr=lr)
#     trace('', distance=P.distance(p) - p.distance(P))
    return P.distance(Point(p))<eps


def isInsideOrFrontierOld(p, Q, eps=1.0e-10):
    u"""
    Paramètres:
    ==========
        - p = (x,y) = np.ndarray((1,2), float) un point en 2d
        - Q est un polygone = np.ndarray((n,2), float) supposé fermé i.e. le premier point == le dernier point,
        - eps : precision
    Retourne :
    ========
        - vrai si p est à l'intérieur du polygone Q ou sur la frontiere (bande de largeur eps) et
        - faux s'il est à l'extérieur.
    Fonctionnement :
    ==============
    p est intérieur, si et seulement si les produits vectoriels AB^Ap ont tous le même signe.
    (on note A=Q[i], B=Q[i+1], 0 <= i < n )
    TODO : optimiser avec raytracing car très consommateur de cpu.
    """
    A, B = Q[0], Q[1]
    AB, Ap = B-A, p-A
    w0 = vect2d(AB, Ap)
#    if w0 == 0.0 :
    if abs(w0) < eps :
#        trace('', AB, Ap)
        return True if 0.0 <= det(AB, Ap) <= det(AB, AB)  else False
#         return True if 0.0 <= np.det(AB, Ap) <= np.det(AB, AB)  else False
#    if number is 0 :
#    trace('', '(x,y)=%s, w0=%.2g'%(str(p), w0))
    for (A, B) in zip(Q[1:-1], Q[2:]) :
        #[A,B] parcourt les segments de Q
        AB, Ap = B-A, p-A
        w =  vect2d(AB, Ap)
#        if w == 0.0 :
        if abs(w) < eps :
            return True if 0.0 <= det(AB, Ap) <= det(AB, AB)  else False
#             return True if 0.0 <= np.det(AB, Ap) <= np.det(AB, AB)  else False

#        if number == 0 :
#        trace('', w0, w)
        if w*w0 < 0 : return False

    return True

def testEmpilementSegments():
    """
    _____1_____________________________________________
    |          |       4    |            |            |
    0          |            |   6        |            8
    |          2    3       5            |            |
    |          |            |            7            |
    |          |            |            |            |
    |          |            |            |            |
    |          |            |            |            |
    |          |            |            |            |
    |          |            |            |            |
    |          |           12            |    10      9
    |          13           |            |            |
    |          |            |           11            |
    14         |            |            |            |
    |          |            |            |            |
    ---------------------------------------------------
    [(0,1,2),(13,14)], [(2,3,4,5),(12,13)], [(),()], [(),()], [(),()]


    """
#def eliminerPointsDoubles0(points, eps=0.0):
#    '''
#    Marche bien pour un tableau de points python (liste de points)
#    Ne marche pas pour ndarray
#    '''
#    i = 0
#    while i<len(points)-1 :
#        if dist(points[i], points[i+1])<=eps :
#            points.remove(points[i+1])
#        else : i+=1
#    return points

def pointsDoubles(points, eps=1.0e-10):
    u"""retourne la liste des numéros de points doubles CONSECUTIFS à eps pres.
    points est de shape(n,2), n>=2"""
    if not isinstance(points, (np.ndarray,list,tuple)) :
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
    if avirer : points = np.delete(points, avirer, 0)
    return points, avirer if vires else points
#         msg1 = 'Les points de controle %s sont des doublons consecutifs.'%avirer
#         msg1 += 'Suppression des doublons'
#         rdebug(msg1,'')
#                 rdebug(avirer=avirer, delta_abscurv=dac.tolist(), controle=sum(dac))
#         for k in reversed(avirer) :
#             self.removePoint(k, update=False)

def testPointsDoubles():
#     from gui.graphicsbase.graphicscommon import qpolygonFrom
    points=np.asarray([[0,0],[0,0],[0,0],])
    debug('points_sales=%s'%points.tolist(), avirer=pointsDoubles(points, 0.0))
    debug(points_propres=eliminerPointsDoublesConsecutifs(points, 0, True))
    X = [1,2,3,4,4,5,5,5,6,7,7]
    Y = [1,2,3,4,4,5,5,5,6,7,7]
    points = zip(X,Y)
    print 'points initiaux :\n',points
    print 'a virer, version Python', pointsDoubles(points)
    points = np.asarray(points)
    print 'a virer, version  numpy', pointsDoubles(points)
    print 'points nettoyé :\n', eliminerPointsDoublesConsecutifs(points, 0)
#     points = qpolygonFrom(points)
#     print 'avirer, version Qt', pointsDoubles(points)

def testIsInside():
    """
    str::isInsideOrFrontier ; (x,y)=[-5.8630197   0.11092514], w0=0.77
    str::isInsideOrFrontier ; 0.766641878852 ; 3.56936145336
    str::isInsideOrFrontier ; 0.766641878852 ; 0.846764334388
    str::isInsideOrFrontier ; 0.766641878852 ; 0.0
    Caisson::maille ; 0 ; maille ;
    [[-5.8630197  -0.37438482]
     [-5.8630197  -0.37264027]
     [-4.2833244  -0.82439736]
     [-4.2833244  -0.82825686]
     [-5.8630197  -0.37438482]]
    Caisson::maille ; inside
    """
    print 10*'#'
    p = [-5.79201521, 0.61218603]
    Q = np.asarray([[-5.92606363, -0.47661251],
                    [-5.78006363, -0.51037401],
                    [-5.78006363,  0.62300842],
                    [-5.92606363,  0.49080254],
                    [-5.92606363, -0.47661251]])
    print isInsideOrFrontier(p, Q, 1.0e-9)
    return
    Q = np.asarray([(0.0,0.0), (1.0,0.5), (1.0,2.0), (0.0, 0.5), (0.,0.)])
    Q = np.asarray([
                    [-5.8630197,  -0.37438482],
                    [-5.8630197,  -0.37264027],
                    [-4.2833244,  -0.82439736],
                    [-4.2833244,  -0.82825686],
                    [-5.8630197,  -0.37438482]
                    ])
    p = np.asarray([-5.8630197,   0.11092514])
    print isInsideOrFrontier(p, Q)
    Q = np.asarray([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0, 1.0), (0.,0.)])
    p = np.asarray([1.0, 1.0])
    print isInsideOrFrontier(p, Q)
#    alert('', 'bizarre, a verifier...')
    return
    p = np.asarray((0.5, 1.24999))
    print isInsideOrFrontier(p, Q)
    print isInside0(p, Q)
    p = np.asarray((0.0, 2.0))
    print isInsideOrFrontier(p, Q)
    print isInside(p, Q)
    print isInside0(p, Q)

def addActions(menu, actions):
    '''
    Pour inserer les actions de la liste 'actions' dans 'menu', avec None qui sert de séparateur dans 'actions',
    et eviter un message d'alerte 'QWidget::insertAction: Attempt to insert null action'
    menu peut etre un QMenu ou bien une QToolBar
    '''
    for action in actions :
        if action is None :
            menu.addSeparator()
        else :
            menu.addAction(action)

def clearLayout(layout):
    ''' Purge d'un layout '''
    if layout is not None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
            else:
                clearLayout(item.layout())

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

# def rcercle(A, B, C, eps=1.0e-8):
#     u"""retourne le rayon dun cercle passant par les 3 points A,B,C, distincts et non alignes.
#     Si les 3 points sont presque alignés retourne np.inf
#     si deux points sont presque confondus, retourne np.nan"""
#     A, B, C = np.asarray(A), np.asarray(B), np.asarray(C)
#     AB, BC, CA = B-A, C-B, A-C
# #     print AB, CA
#     c, a, b = np.sqrt(AB[0]**2 + AB[1]**2), np.sqrt(BC[0]**2 + BC[1]**2), np.sqrt(CA[0]**2 + CA[1]**2),
#     abc = a*b*c
#     s = abs(det(AB,CA))
# #     print s, abc
#     if abc < eps:#2 points presque confondus : infinité de cercles passant par deux points
#         return np.nan
#     elif s<eps:# trois points presque alignés
#         return np.inf
#     else :
#         return 0.5*abc/(det(AB, CA))
#         return 0.5*abc/abs(det(AB, CA))


def pointsFromFile(filename):
    filename = unicode(filename)
    filename = Path(filename)
    ext = filename.ext.lower()
#     debug(filename=filename, ext=ext)
    if ext in ('.dxf',) :
        raise NotImplementedError
        # lecteur = LecteurDXFNervures(filename)
        # polygon = lecteur.lire()
        # return polygon
    elif ext in ('.gnu',) :
        raise NotImplementedError
        # lecteur = LecteurGnuplot(filename,patch=False)
        # lecteur.dim = 2
        # lecteur.patch = False
        # polygon = lecteur.lire()
        # return polygon
    elif ext in ('.pts',) :
        raise NotImplementedError
        # lecteur = LecteurNervures(filename)
        # polygon = lecteur.lire()
        # if lecteur.genre.lower()=='voile' :
        #     raise IOError(whoami()+u"%s n'est pas un polygon2d mais 3d (voile?, profil3d?)"%filename.name)
        # elif lecteur.dimension == 3 :
        #     polygon = lecteur.points[:,:-1]
        # return polygon
    elif ext in ('.spl',) :
        """Une spline"""
        with open(filename,'r') as f :
            dump = eval(f.read())
#         if len(lines)>1 :
#             line = ' '.join([l.strip() for l in lines])
#         else :
#             line = lines[0]
#         dump = eval(line)
#         debug(dump=dump)
        try :
            return np.asarray(dump['cpoints'])
        except :
            raise IOError(u"je ne sais pas extraire les points de ce fichier : %s"%filename.name)
    else :
#         rdebug()
        raise IOError(u'Je ne sais pas lire ce fichier "%s"'%filename)
        return np.zeros((0,2))


def pointsFrom(x):#,close_=False):
    u'''
    Paramètre x:
        - un np.ndarray((n,2))
        - un chaine de caractères evaluable par eval(x)
        - un nom de fichier
    Retourne :
    --------
    - points : np.ndarray, de shape (N,2) si close_==False
    '''
    debug(x)
    if x is None :
        return np.zeros((0,2))
    elif isinstance(x, (str, unicode)) :
        try :
            debug(x=x)
            return pointsFromFile(x)#Fichier
        except Exception as msg:
            try :
                return pointsFrom(eval(x))#chaine de caracteres p.ex.'[1,2,5,12]'
            except Exception as msg : #pass
                debug('Exception:', str(msg))
                raise 
    elif isinstance(x, (list, tuple)) :
        npa = np.asarray(x)
#         debug(shape=npa.shape)
        sh = npa.shape
        if len(sh) != 2 or sh[1] != 2 :
            npa.shape=(len(npa)/2,2)
        return npa
    elif isinstance(x, np.ndarray) :
#         debug(u'x est deja un np.ndarray. Je ne modifie pas sa shape...%s'%str(x.shape))
        return x
    else :
        raise TypeError (u"pointsFrom() : cas non pris en charge")

def rcercle(A, B, C, eps=0.0):
    u"""retourne le rayon dun cercle passant par les 3 points A,B,C, distincts et non alignes.
    Si les 3 points sont presque alignés retourne np.inf
    si deux points sont presque confondus, retourne np.nan"""
    A, B, C = np.asarray(A), np.asarray(B), np.asarray(C)
    AB, BC, CA = B-A, C-B, A-C
#     print AB, CA
    c, a, b = np.sqrt(AB[0]**2 + AB[1]**2), np.sqrt(BC[0]**2 + BC[1]**2), np.sqrt(CA[0]**2 + CA[1]**2),
    abc = a*b*c
    d = det(AB,CA)
    s = abs(d)
#     print s, abc
    if abc <= eps:#2 points presque confondus : infinité de cercles passant par deux points
        return np.nan
    elif s <= eps:# trois points presque alignés
        return np.inf
    else :
        return 0.5*abc/d
        return 0.5*abc/abs(det(AB, CA))

def rayonDeCourbure(P):
    u"""
    Retourne la liste des rayons de courbure du polyGONE P.
    S'il n'est pas fermé, on le ferme.
    Il y a exactement n=len(P) rayons
    """
    rayons = zeros(len(P))
    if np.linalg.norm(P[0]-P[-1])>1.0e-9 :#debut != fin
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
    if np.linalg.norm(P[0]-P[-1])>0 :#debut != fin
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
        un np.ndarray((n,)) avec les valeurs de la courbure c(ti) aux pointx sx(ti), sy(ti), ti dans T
        La courbure est l'inverse du rayon de courbure.
        Pour un arc x(t), y(t), le rayon de courbure est r(t)=((x'^2+y'^2)^3/2)/(x'y"-y'x")
        x'=x'(t), y'=y'(t) x"=x"(t), y"=y"(t).
        si ||(x",y")|| = 0, la courbure est nulle.
        cf https://fr.wikipedia.org/wiki/Rayon_de_courbure
    """
    sx, sy = S
    dx,  dy  = sx(T, 1), sy(T, 1)
    d2x, d2y = sx(T, 2), sy(T, 2)
    norm3_d2 = np.sqrt(dx**2+dy**2)**3
    # si norm_d2=0, x"(t)=y"(t)=0, c'est une droite, courbure nulle
    sc = (dx*d2y-dy*d2x)/(norm3_d2)
    sc[np.where(norm3_d2 < 1.0e-12)] = 0.0
    return sc
#     return np.abs(dx*d2y-dy*d2x)/(np.sqrt(dx**2+dy**2)**3)
def simpson(f, a, b, n=10):#n doit être pair, integration precise ordre 3
    u"""Integrale de f sur [a,b], méthode de Simpson composite. (ordre 3)
    n DOIT être pair"""
    h = float(b-a)/n
    T = np.linspace(a, b, n+1)
    C = f(T)
    A1 = C[0] + C[-1]
    A2 = 2*sum(C[i] for i in range(2,n) if i%2==0)
    A4 = 4*sum(C[i] for i in range(1,n) if i%2==1)
#         debug (h, A1, A2, A4, (h/3)*(A1 + A2 + A4))
    return (h/3)*(A1 + A2 + A4)

def testSimpson():
    def f0(T): return ones(len(T))
    def f1(T): return T
    def f2(T): return T*T
    def f3(T): return T*T*T
    def f4(T): return T*T*T*T
    debug(simpson(f0, 0, 1, 1000))
    debug(simpson(f1, 0, 1, 1000))
    debug(simpson(f2, 0, 1))
    debug(simpson(f3, 0, 1))
    debug(simpson(f4, 0, 1))

def testSegmentPlusProche():

    points = [[1.00000000e+02, 0.00000000e+00],
            [2.38750710e+01, 1.12195000e+01],
            [2.06107410e+01, 1.11545910e+01],#Le plus proche
            [1.75275990e+01, 1.09294100e+01],
            [1.20666720e+01, -5.76084900e+00],
            [1.46258840e+01, -5.92738400e+00],
            [1.77233340e+01, -6.10722700e+00],#trouvé...
            [2.07627560e+01, -6.22899600e+00],
            [2.41177090e+01, -6.36394300e+00],
            [1.00000000e+02, 0.00000000e+00]]
    pt = (20.,20.0)
    points = [[ 0., 0. ], [ 1., -1. ],[ 3., -0.5],[ 5., 0. ],[ 2.0, 0.5 ]]#, [ 0., 0. ]]
    pt = (0.0,0.0)
    pt = [ 1., -1. ]
    pt = [ 1.8, 0.5 ]
    pt = (5.20,-0.20)
    pt = (2.0,-1.0)
#     pt1 = 21,11
    idx, pj = segmentPlusProche(points, pt)
    A, B = points[idx],points[idx+1]
    debug(indice_segment=idx, projection=pj, point_a_inserer=pt)
    from matplotlib import pyplot as plt
    points = np.asarray(points)
    fig = plt.figure()
    a = fig.add_subplot(1,1,1)
    plt.plot(points[:,0], points[:,1],'b-o',label=u'polyligne')
    plt.plot([A[0],B[0]], [A[1],B[1]],'r-',label=u'segment %d : winner'%idx)
    if pj is not None :
        plt.plot([pt[0],pj[0]],[pt[1],pj[1]],'g-o',label=u'projection')
    else :
        plt.plot([pt[0]],[pt[1]],'r-o',label='??')
    plt.plot()
    a.set_aspect(1)
#     plt.plot([pt[0],pt1[0]],[pt[1],pt1[1]],'r-o',label='winner')
    plt.legend()
    plt.show()
    return
    p7 = 7,0#None, None
    p2 = 2,0#(1, array([ 2.17647059, -0.70588235]))
    for x in range(7) :
        S = segmentPlusProche(points, (x,0))
        print 'Point = ', (x,0) , ' ; segmentPlusProche =',S

def testCourbure():
    p = np.asarray([[0.0, 0.0],
                    [0.0020691201283767205, -0.0042941912136635636], [0.007625888799913062, -0.011329446395562288], [0.015694418444218296, -0.019910545432815785],
                    [0.025298821490901697, -0.02884226821254368], [0.035463210369572534, -0.03692939462186559], [0.04525122291119395, -0.0430513056100728],
                    [0.05446593471151137, -0.04748301229776649], [0.06356139613433973, -0.05171627435197123], [0.07281050074621176, -0.05640916969890315],
                    [0.08221422031141906, -0.06105196236599214], [0.0917576221148263, -0.06524123911924189], [0.10143185809958596, -0.06895424232522496],
                    [0.1112288171535668, -0.07221431753663529], [0.12114038816463749, -0.07504481030616675], [0.1311584600206667, -0.07746906618651322],
                    [0.1412749216095232, -0.0795104307303686], [0.15148166181907566, -0.0811922494904268], [0.1617705695371927, -0.0825378680193817],
                    [0.17213353365174314, -0.0835706318699272], [0.18256244305059566, -0.08431388659475722], [0.19304918662161893, -0.0847909777465656],
                    [0.20358565325268163, -0.08502525087804626], [0.2141637318316525, -0.08504005154189312], [0.2247753112464002, -0.08485872529080002],
                    [0.2354122803847935, -0.08450461767746092], [0.24606652813470106, -0.08400107425456962], [0.2567299433839916, -0.08337144057482014],
                    [0.26739441502053374, -0.0826390621909063], [0.2780518319321962, -0.08182728465552197], [0.28869408300684785, -0.08095945352136108],
                    [0.2993130571323571, -0.08005891434111753], [0.30990192473004063, -0.07914582201669142], [0.3204600849636786, -0.07822482366739393],
                    [0.33098882667127716, -0.0772958616655942], [0.34148943911682783, -0.0763588773230801], [0.35196321156432236, -0.07541381195163939],
                    [0.3624114332777522, -0.07446060686305987], [0.372835393521109, -0.07349920336912937], [0.38323638155838424, -0.07252954278163569],
                    [0.39361568665356955, -0.07155156641236664], [0.4039745980706565, -0.07056521557311002], [0.4143144050736366, -0.06957043157565367],
                    [0.42463639692650157, -0.06856715573178536], [0.43494186289324277, -0.06755532935329289], [0.445232092237852, -0.0665348937519641],
                    [0.45550837422432056, -0.06550579023958679], [0.4657719981166402, -0.06446796012794875], [0.4760242531788025, -0.06342134472883781],
                    [0.48626642867479897, -0.06236588535404178], [0.49649981386862113, -0.061301523315348445], [0.5067256980242607, -0.06022819992454561],
                    [0.5169453704057091, -0.059145856493421106], [0.5271601202769579, -0.058054434333762735], [0.537371236901999, -0.056953874757358296],
                    [0.5475800095448234, -0.0558441190759956], [0.5577877274694232, -0.05472510860146246], [0.5679956799397896, -0.05359678464554667],
                    [0.5782051562199142, -0.05245908852003607], [0.588417445573789, -0.051311961536718416], [0.5986338372654051, -0.050155345007381565],
                    [0.6088554653393123, -0.048989213246300733], [0.6190825327819455, -0.04781373852770645], [0.6293148970074011, -0.046629166600818006],
                    [0.6395524148593865, -0.04543574333612987], [0.6497949431816094, -0.04423371460413647], [0.660042338817777, -0.04302332627533221],
                    [0.670294458611597, -0.04180482422021157], [0.6805511594067771, -0.040578454309268965], [0.6908122980470244, -0.03934446241299885],
                    [0.7010777313760468, -0.03810309440189562], [0.7113473162375514, -0.03685459614645376], [0.7216209094752458, -0.035599213517167674],
                    [0.7318983679328378, -0.034337192384531805], [0.7421795484540348, -0.03306877861904061], [0.7524643078825441, -0.031794218091188486],
                    [0.7627525030620732, -0.030513756671469904], [0.7730439908363299, -0.029227640230379295], [0.7833386280490214, -0.02793611463841107],
                    [0.7936362715438553, -0.026639425766059696], [0.8039367781645393, -0.025337819483819583], [0.8142400047547806, -0.024031541662185196],
                    [0.824545808158287, -0.02272083817165094], [0.8348540452187656, -0.02140595488271128], [0.8451645727799243, -0.020087137665860646],
                    [0.8554772476854705, -0.018764632391593452], [0.8657919267791115, -0.017438684930404163], [0.8761084669045551, -0.016109541152787184],
                    [0.8864267249055084, -0.014777446929236986], [0.8967465576256795, -0.01344264813024797], [0.9070678219087754, -0.012105390626314615],
                    [0.9173903745985037, -0.010765920287931328], [0.9277140725385721, -0.00942448298559254], [0.938038772572688, -0.008081324589792714],
                    [0.9483643315445588, -0.006736690971026247], [0.958690606297892, -0.005390827999787622], [0.9690174536763952, -0.004043981546571245],
                    [0.9793447305237761, -0.0026963974818715497], [0.9896722936837418, -0.001348321676182991], [0.9999999999999999, -3.903127820947816e-18]])

    c=np.abs(courbure(p))
    r=rayonDeCourbure(p)

    print c
    print r

def testSCourbure():
    C = np.asarray([(1,0),(0,1),(-1,0),(0,-1),(1,0)])
    ac = absCurv(C, normalise=True)
    k = 3
    sx4 = UnivariateSpline(ac, C[:,0], k=k, s=0.0)
    sy4 = UnivariateSpline(ac, C[:,1], k=k, s=0.0)
    T = np.linspace(0,1, 10000)
    X4, Y4 = sx4(T), sy4(T)
    curv4 = scourbure((sx4,sy4), T)-1

    pprint(curv4)
    from matplotlib import pyplot as plt
    plt.plot(T, sx4(T), T, sy4(T))
    plt.show()
    plt.plot(T, curv4)
    plt.plot(X4,Y4)
    plt.show()

def testRayonCercle():
    T = np.linspace(0, 2*math.pi, 100)
#     T.shape = 10,1
#     print T,
    R = np.random.rand()
    print R
    X, Y = R*np.cos(T), R*np.sin(T)
    P = np.asarray([(x,y) for x,y in zip(X,Y)])
    print P
    print rcercle(P[-1], P[0], P[1])


if __name__=="__main__":
    testPointsDoubles()
#     exit()
    testSegmentPlusProche()

    p = [[2.97370904e+00, 6.44820360e-19],
        [2.62290948e+00, 6.02959508e-02],
        [2.24308167e+00, 1.23744268e-01],
        [1.87452472e+00, 1.83613921e-01],
        [1.50115462e+00, 2.42394363e-01],
        [1.22747086e+00, 2.84350243e-01],
        [1.12399937e+00, 2.98561198e-01],
        [1.00689675e+00, 3.12871734e-01],
        [9.34115493e-01, 3.20342559e-01],
        [8.67134243e-01, 3.25454296e-01],
        [7.98757537e-01, 3.28922015e-01],
        [7.32512484e-01, 3.30579575e-01],
        [6.69915006e-01, 3.30470431e-01],
        [6.12428610e-01, 3.28699942e-01],
        [5.58750286e-01, 3.25359730e-01],
        [5.07959735e-01, 3.20500775e-01],
        [4.60353384e-01, 3.14257943e-01],
        [4.15956345e-01, 3.06735929e-01],
        [3.73947713e-01, 2.97892075e-01],
        [3.33827496e-01, 2.87700041e-01],
        [2.96006551e-01, 2.76350898e-01],
        [2.61205671e-01, 2.64151008e-01],
        [2.28902585e-01, 2.51027363e-01],
        [1.98298707e-01, 2.36743322e-01],
        [1.69361638e-01, 2.21350097e-01],
        [1.43241296e-01, 2.05567685e-01],
        [1.19840851e-01, 1.89463177e-01],
        [9.85093007e-02, 1.72704901e-01],
        [7.88006468e-02, 1.55037597e-01],
        [6.12129425e-02, 1.37020921e-01],
        [4.64766811e-02, 1.19548694e-01],
        [3.41003061e-02, 1.02253781e-01],
        [2.37764244e-02, 8.49053918e-02],
        [1.55363360e-02, 6.77601586e-02],
        [8.82760258e-03, 4.97397618e-02],
        [3.45630886e-03, 3.02145976e-02],
        [7.77361884e-04, 1.51843334e-02],
        [7.55648859e-21, 0.00000000e+00],
        [9.68226599e-04, -1.28975276e-02],
        [3.25212510e-03, -2.40215827e-02],
        [7.46226482e-03, -3.44277207e-02],
        [1.25640953e-02, -4.28441260e-02],
        [1.91184423e-02, -4.98146745e-02],
        [2.97370904e-02, -5.77199739e-02],
        [4.47227984e-02, -6.73131096e-02],
        [5.61194425e-02, -7.43857830e-02],
        [6.90344014e-02, -8.20925698e-02],
        [8.52837290e-02, -9.14662261e-02],
        [1.16833061e-01, -1.09259597e-01],
        [1.35057261e-01, -1.19683057e-01],
        [1.48685452e-01, -1.27730648e-01],
        [1.65211057e-01, -1.37333205e-01],
        [1.79501219e-01, -1.43984114e-01],
        [1.95710381e-01, -1.50032156e-01],
        [2.14137505e-01, -1.55449459e-01],
        [2.35700604e-01, -1.60332396e-01],
        [2.61869977e-01, -1.64782729e-01],
        [2.95823959e-01, -1.69070402e-01],
        [3.40834223e-01, -1.73219635e-01],
        [4.17064411e-01, -1.78399785e-01],
        [5.26630638e-01, -1.84648210e-01],
        [6.52240905e-01, -1.89484615e-01],
        [7.74966197e-01, -1.94548433e-01],
        [8.45738690e-01, -1.96600587e-01],
        [9.22197728e-01, -1.97437555e-01],
        [1.03117425e+00, -1.97001292e-01],
        [1.12509145e+00, -1.95546293e-01],
        [1.20976462e+00, -1.92949955e-01],
        [1.30714611e+00, -1.88476358e-01],
        [1.42330517e+00, -1.81902516e-01],
        [1.50908393e+00, -1.75841562e-01],
        [1.59606381e+00, -1.68363291e-01],
        [1.70815592e+00, -1.57186925e-01],
        [1.79449328e+00, -1.47547612e-01],
        [1.91366126e+00, -1.32171863e-01],
        [2.01694955e+00, -1.18217320e-01],
        [2.08412098e+00, -1.07736892e-01],
        [2.19909217e+00, -8.79476287e-02],
        [2.27284129e+00, -7.59140432e-02],
        [2.37543690e+00, -6.09027024e-02],
        [2.47881269e+00, -4.69689912e-02],
        [2.57710869e+00, -3.50125893e-02],
        [2.67456383e+00, -2.44616999e-02],
        [2.77305513e+00, -1.51044063e-02],
        [2.87236864e+00, -6.97302358e-03],
        [2.97370904e+00, 6.22623175e-19]]
    p=[[  2.97370904e+00,   6.44820360e-19],
 [  2.62290948e+00,   6.02959508e-02],
 [  2.24308167e+00,   1.23744268e-01],
 [  1.87452472e+00,   1.83613921e-01],
 [  1.50115462e+00,   2.42394363e-01],
 [  1.22747086e+00,   2.84350243e-01],
 [  1.12399937e+00,   2.98561198e-01],
 [  1.00689675e+00,   3.12871734e-01],
 [  9.34115493e-01,   3.20342559e-01],
 [  8.67134243e-01,   3.25454296e-01],
 [  7.98757537e-01,   3.28922015e-01],
 [  7.32512484e-01,   3.30579575e-01],
 [  6.69915006e-01,   3.30470431e-01],
 [  6.12428610e-01,   3.28699942e-01],
 [  5.58750286e-01,   3.25359730e-01],
 [  5.07959735e-01,   3.20500775e-01],
 [  4.60353384e-01,   3.14257943e-01],
 [  4.15956345e-01,   3.06735929e-01],
 [  3.73947713e-01,   2.97892075e-01],
 [  3.33827496e-01,   2.87700041e-01],
 [  2.96006551e-01,   2.76350898e-01],
 [  2.61205671e-01,   2.64151008e-01],
 [  2.28902585e-01,   2.51027363e-01],
 [  1.98298707e-01,   2.36743322e-01],
 [  1.69361638e-01,   2.21350097e-01],
 [  1.43241296e-01,   2.05567685e-01],
 [  1.19840851e-01,   1.89463177e-01],
 [  9.85093007e-02,   1.72704901e-01],
 [  7.88006468e-02,   1.55037597e-01],
 [  6.12129425e-02,   1.37020921e-01],
 [  4.64766811e-02,   1.19548694e-01],
 [  3.41003061e-02,   1.02253781e-01],
 [  2.37764244e-02,   8.49053918e-02],
 [  1.55363360e-02,   6.77601586e-02],
 [  8.82760258e-03,   4.97397618e-02],
 [  3.45630886e-03,   3.02145976e-02],
 [  7.77361884e-04,   1.51843334e-02],
 [  8.81590336e-21,   6.44820360e-19],
 [  9.68226599e-04,  -1.28975276e-02],
 [  3.25212510e-03,  -2.40215827e-02],
 [  7.46226482e-03,  -3.44277207e-02],
 [  1.25640953e-02,  -4.28441260e-02],
 [  1.91184423e-02,  -4.98146745e-02],
 [  2.97370904e-02,  -5.77199739e-02],
 [  4.47227984e-02,  -6.73131096e-02],
 [  5.61194425e-02,  -7.43857830e-02],
 [  6.90344014e-02,  -8.20925698e-02],
 [  8.52837290e-02,  -9.14662261e-02],
 [  1.16833061e-01,  -1.09259597e-01],
 [  1.35057261e-01,  -1.19683057e-01],
 [  1.48685452e-01,  -1.27730648e-01],
 [  1.65211057e-01,  -1.37333205e-01],
 [  1.79501219e-01,  -1.43984114e-01],
 [  1.95710381e-01,  -1.50032156e-01],
 [  2.14137505e-01,  -1.55449459e-01],
 [  2.35700604e-01,  -1.60332396e-01],
 [  2.61869977e-01,  -1.64782729e-01],
 [  2.95823959e-01,  -1.69070402e-01],
 [  3.40834223e-01,  -1.73219635e-01],
 [  4.17064411e-01,  -1.78399785e-01],
 [  5.26630638e-01,  -1.84648210e-01],
 [  6.52240905e-01,  -1.89484615e-01],
 [  7.74966197e-01,  -1.94548433e-01],
 [  8.45738690e-01,  -1.96600587e-01],
 [  9.22197728e-01,  -1.97437555e-01],
 [  1.03117425e+00,  -1.97001292e-01],
 [  1.12509145e+00,  -1.95546293e-01],
 [  1.20976462e+00,  -1.92949955e-01],
 [  1.30714611e+00,  -1.88476358e-01],
 [  1.42330517e+00,  -1.81902516e-01],
 [  1.50908393e+00,  -1.75841562e-01],
 [  1.59606381e+00,  -1.68363291e-01],
 [  1.70815592e+00,  -1.57186925e-01],
 [  1.79449328e+00,  -1.47547612e-01],
 [  1.91366126e+00,  -1.32171863e-01],
 [  2.01694955e+00,  -1.18217320e-01],
 [  2.08412098e+00,  -1.07736892e-01],
 [  2.19909217e+00,  -8.79476287e-02],
 [  2.27284129e+00,  -7.59140432e-02],
 [  2.37543690e+00,  -6.09027024e-02],
 [  2.47881269e+00,  -4.69689912e-02],
 [  2.57710869e+00,  -3.50125893e-02],
 [  2.67456383e+00,  -2.44616999e-02],
 [  2.77305513e+00,  -1.51044063e-02],
 [  2.87236864e+00,  -6.97302358e-03],
 [  2.97370904e+00,   6.45922348e-19]]

    my2dPlot([np.asarray(p),])
#     exit()
#     app=QtGui.QApplication(sys.argv)
#     testTableDataDiag()
# #     tde = TableDataEditor()
#     sys.exit(app.exec_())

#     testRayonCercle()
#     exit()
#     testSCourbure()
#     testCourbure()
#     exit()
#     testSegmentPlusProche()
#     exit()
    testIsInside()
#     exit()
    A, B, C = [0,     0], [1,     0], [0,      1]
    print rcercle(A, B, C)
    A, B, C = [0.0, 0.0], [-1.0, 0.0], [1.01, 1.0]
    print rcercle(A, B, C)
    print '***'
    P = [[1.00000000e+00,   0.00000000e+00],
        [9.37800919e-01,  -4.53983927e-03],
        [8.35970247e-01,  -1.52885157e-02],
        [6.92659143e-01,  -3.72445481e-02],
        [5.04164573e-01,  -5.87701128e-02],
        [3.13878806e-01,  -6.54986341e-02],
        [1.46258259e-01,  -5.92691821e-02],
        [6.34134460e-02,  -4.84279108e-02],
        [1.19604073e-02,  -1.94216706e-02],
        [2.22044605e-16,   0.00000000e+00],
#         [2.22044605e-16,   0.00000000e+00],
        [3.80610473e-02,   6.31931730e-02],
        [1.46449137e-01,   1.05380236e-01],
        [3.08658678e-01,   1.09106172e-01],
        [5.00560799e-01,   8.28268707e-02],
        [6.91521129e-01,   5.22806851e-02],
        [8.53153837e-01,   2.53405643e-02],
        [9.39723441e-01,   1.05280598e-02],
        [1.00000000e+00,   0.00000000e+00],
        ]
#     print courbure(np.asarray(P[:-1]))
    ac = absCurv(P)
    ac /= ac[-1]
#     print ac
    c = courbure(np.asarray(P))
    T, sx,sy = splineInterpolation(np.asarray(P), 'c cubic')
#     print T-ac
#     sc = scourbure((sx,sy), T)
#     print np.linalg.norm(c-sc)/ac[-1]
#     print sc
    print "***cercle rayon 1"
    T = np.linspace(0, 2*np.pi,100)#[:-1]
#     P = np.asarray(zip(T, T))
    P = np.asarray(zip(np.sin(T), np.cos(T)))
#     print P
    Ts, sx,sy = splineInterpolation(np.asarray(P), 'c cubic')
    c = courbure(P)
    print c
    T1 = np.linspace(0, 1, 20)
    P1 = np.asarray(zip(sx(T1),sy(T1)))
    c1 = courbure(P1)
    print c1
    from matplotlib import pyplot
    print len(Ts), len(c), Ts[-1]
    print len(T1), len(c1)
    pyplot.plot(Ts, c,'r.',T1, c1,'b.')
    pyplot.show()
#     pyplot.plot(sx(Ts,2),'r.', sy(Ts,2),'b.')
#     pyplot.show()
#     pyplot.plot(sx(Ts,2), sy(Ts,2),'.', )
#     pyplot.show()
#     pyplot.plot(sc ,'.')
#     pyplot.show()
#     pyplot.axes().set_aspect('equal')
#     pyplot.plot(P[:,0], P[:,1] ,'r.', sx(Ts) ,sy(Ts),'-', )
#     pyplot.show()

#     exit()
    points  = np.arange(1,25)
    points.shape = 2,4,3
    points.shape = 3,4,2
    print '\n==> x\n', symetriser(points, 'x')
    print '\n==> y\n', symetriser(points, 'y')
    print '\n==> z\n', symetriser(points, 'z')
#     exit()

#     help()
#     testPointsDoubles()
#     testIsInside()
#     exit()
#     testPretty()
    import matplotlib
#     testContraindre()
#     testCurv()
#     exit()
#     print goodName('name')
#     print goodName('nâme')
# #     exit()
    points=np.asarray((np.arange(5),np.arange(5,10))).reshape((-1,2))
    print points
    print hardScale(points, (2,0), centre=None, translate=False)
# #     exit()
#     dico=dict(locals())
#     for key,value in  dico.iteritems():
#         if type(value).__name__ in ('function',) :
#             print "%s:\n"%(key)+len(key)*'='+"\n%s"%value.__doc__
#     print phrase()
#     from config import VALIDATION_DIR
#     filename=Path(VALIDATION_DIR,'dxf','cot_d03.DXF')
#     lines=open(filename).readlines()
#     print findAll('LINE',lines)
#     p1=QPointF(1,0)
#     p2=QPointF(1,1)
#     print dist2(p1,p2)
#     p1=[1,0.]
#     p2=[1,1]
#     print dist2(p1,p2)
#     #### segmentPlusProche(points, P):####
    P=np.asarray([(3.0,-1.0),(2,-1),(0,0),(2,2),(3,2),(5,2),(10,0)])
#     print P
#     P=np.asarray([(0,-1.0),(-1,0),(1,0),(0,1)])
#     print 'barycentre',baryCentre(P)
#     exit()
    a=np.asarray([2.0,0.0])
    i=segmentPlusProche(P,a)
    a.shape=(1,2)
    my2dPlot((P,a),equal=True,cosmetic=('b-','bo','r-','ro','g*','r*','r^'))
#     print P[i],P[i+1]
#    exit()
#    print rreplace("kze,$$$$$l$$$izefh","$$","$")
#    datadir = '/Users/puiseux/Documents/workspace/pyglide/test'
#    filename = Path(datadir, 'simple.OUT')
#    lines =  open(filename).readlines()
#    print find('yyy', lines)
#    print findAll('xxxxxx', lines)
#    tag = 'TIME STEP'
#    deb = _findRel(tag, lines)
#    print "relatif", _findRel(tag, lines, 100)
#    print "absolu",find(tag, lines, 4000)
#    print findAll(tag, lines)
#    points = np.arange(21, dtype='float')
#    ppp = points
#    points.shape = (7,3)
#    print points
#    print find_names(points)
#    print structuredGrid(2,3,1)
#    print 'doublons : '
#    nbp = 10
#    a = np.random.randint(0, 3, 3*nbp).reshape((-1,3))
##    voisinages = [range(i,nbp) for i in range(nbp)]#[[0,...n-1],[1,...n-1],[2,...n-1],...[n-1]]
#    voisinages = [[i]+range(nbp) for i in range(nbp)]#[[0,0,...n-1],[1,0,...n-1],[2,0,...n-1],...[n-1, 0,...n-1]]
##    print a
#    for k0,k1 in doublons(a, voisinages=voisinages, sym=False) :
#        print k0, k1
#        print a[k0],a[k1]
#

#    print encombrement(points)

#    a = 100.0
#    points = [(-a,-a),(a,-a),(a,a),(-a,a)]
#    print encombrement(points)
#    tags = ['zozo','TIME STEP','COMPONENT COORDINATE SYSTEMS','AERODYNAMIC DATA FOR PATCH','UNIT NORMAL VECTORS AND PANEL AREAS FOR PATCH']
#
#    print mark(tags, lines)????
#     sys.exit(app.exec_())
