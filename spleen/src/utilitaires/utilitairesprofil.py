#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Created on 28 octobre 2016

@author: puiseux
'''
import sys,string
import numpy as np
from path import Path
import pickle
# from exceptions import *
# import datetime
# import scipy as sp
# from scipy import interpolate, sum, prod, ones
# from random import choice
# from PyQt4 import QtGui,QtCore
# from PyQt4.QtCore import Qt,QSize,QVariant
# from PyQt4.QtGui import (QTableWidget,QTableWidgetItem,QPushButton,QLayout,QVBoxLayout,QHBoxLayout,QGridLayout,QSpacerItem,QSizePolicy,
#                          QWidget,QScrollArea,QApplication,QKeySequence,QShortcut,QFont,QFontMetrics,QHeaderView)
# from PyQt4.Qt import SIGNAL,SLOT
# from inout.format import formData
#from PyQt4.QtCore import QPointF
#from gui.graphicsbase.graphicscommon import pointsFromPolygon
# import math
# from scipy.interpolate import (Rbf, InterpolatedUnivariateSpline,
#                                LSQUnivariateSpline, UnivariateSpline,
#                                splrep, splev, interp1d)
from utilitaires.utilitaires import absCurv, courbure, splineInterpolation, dist2
from shapely.geometry import LineString, Polygon, Point
# from geos import LineString
# DEBOG=True
# VERBOSITY_LEVEL=2
# DEBOG_OUTPUT=sys.stderr
# DEBOG_OUTPUT=sys.stdout
# output=DEBOG_OUTPUT

def computeCordeAndNBA(points):
    u"""
    retourne la distance et le numéro du point le plus éloigné de BF=points[0]
    Dés que la distance décroit, on y est.
    """
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :
            corde, nba = d, k
    return np.sqrt(corde), nba

def ajustePoints(n, s, dt=[0,1], prec=500):
    u"""positionne n points (= un polyligne) sur la spline s[dt], de sorte que
    la distance entre s et le polyligne soit minimale"""
    sx, sy = s
    t0, t1 = dt
    Ts = np.linspace(t0, t1, prec)
    S = np.asarray(zip(sx(Ts), sy(Ts)))
    T = np.linspace(t0, t1, n)
    P = np.asarray(zip(sx(T),sy(T)))
#     print P
    P = LineString(P)
#     .....TODO
    return S, P,
if __name__=="__main__":
    print '***'
    P0 = [[  1.00000000e+00,   0.00000000e+00],
        [  9.37800919e-01,  -4.53983927e-03],
        [  8.35970247e-01,  -1.52885157e-02],
        [  6.92659143e-01,  -3.72445481e-02],
        [  5.04164573e-01,  -5.87701128e-02],
        [  3.13878806e-01,  -6.54986341e-02],
        [  1.46258259e-01,  -5.92691821e-02],
        [  6.34134460e-02,  -4.84279108e-02],
        [  1.19604073e-02,  -1.94216706e-02],
        [  2.22044605e-16,   0.00000000e+00],
#         [  2.22044605e-16,   0.00000000e+00],
        [  3.80610473e-02,   6.31931730e-02],
        [  1.46449137e-01,   1.05380236e-01],
        [  3.08658678e-01,   1.09106172e-01],
        [  5.00560799e-01,   8.28268707e-02],
        [  6.91521129e-01,   5.22806851e-02],
        [  8.53153837e-01,   2.53405643e-02],
        [  9.39723441e-01,   1.05280598e-02],
        [  1.00000000e+00,   0.00000000e+00],
        ]
    P0 = list(reversed(P0))
    print P0
#     print courbure(np.asarray(P[:-1]))
    P0 = np.asarray(P0)
    ac = absCurv(P0)
    ac /= ac[-1]
#     c = courbure(np.asarray(P))
#     print c
    Ts, sx,sy = splineInterpolation(np.asarray(P0), 'c cubic')
    dt = 0.4, 0.6
    S, P = ajustePoints(10, (sx,sy), dt)
    P = np.asarray(P.coords)
    print list(P)
#     print T-ac
#     sc = scourbure((sx,sy), T)
#     print np.linalg.norm(c-sc)/ac[-1]
#     print sc
#     print "***cercle rayon 1"
#     T = np.linspace(0, 2*np.pi,100)#[:-1]
# #     P = np.asarray(zip(T, T))
#     P = np.asarray(zip(np.sin(T), np.cos(T)))
# #     print P
#     Ts, sx,sy = splineInterpolation(np.asarray(P), 'c cubic')
#     c = courbure(P)
#     print c
    T1 = np.linspace(dt[0], dt[1], 200)
    P1 = np.asarray(zip(sx(T1),sy(T1)))
#     c1 = courbure(P1[1:-1])
#     print c1
    from matplotlib import pyplot
#     print len(Ts), len(c), Ts[-1]
#     print len(T1), len(c1)
#     pyplot.plot(Ts, c,'r-',T1, c1,'b.')
    pyplot.plot(S[:,0], S[:,1], 'r-', P[:,0], P[:,1], 'o')
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

    exit()
