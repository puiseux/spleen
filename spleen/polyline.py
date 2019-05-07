#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures

    Classe NPolyline

@author:     Michel Le Berre, Pierre Puiseux
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 31 Jan 2013
'''
from path import Path
from utilitaires import (debug, absCurv, my2dPlot,)
import numpy as np
from numpy import asarray
from utilitaires.lecteurs import pointsFrom, pointsFromFile
from splinesimple import NSplineSimple
from matplotlib import pyplot as plt

class NPolyLine(NSplineSimple):
    u"""
    Un NPolyLine est une NSpline d'interpolation, de degré 1 : methode=('ius', 1)
    """
    SUPPORTED_FILES = ('*.gnu', '*.dxf', '*.pts', '*.pkl')
#     def __init__(self, parent=None, points=None, name='NoName', role='polyline', **dump):
    def __init__(self, **dump):
#         dump['methode'] = ('ius', 1)#c'estt fait dans load
#         dump['mode']  = 'cpoints'
        super(NPolyLine, self).__init__(**dump)
        self.precision = None
#         self._update()#deja fait dans NSplineSimple a supprimer ?
#         debug(u' _update() de trop ?')
#         debug(mode=self.mode, methode=self.methode)
    def load(self,dump):
        dump['methode'] = 'ius',1
        dump['mode']  = 'cpoints'
        super(NPolyLine, self).load(dump)
        self.precision = None

    @property
    def precision(self):
        u"""precision=Le nb de points de discrétisation pour affichage.
        Les points de discrétisation sont les points de contrôle"""
        return len(self.cpoints)

    @precision.setter
    def precision(self, prec):
        u"""prec peut prendre n'importe quelle valeur, ça n'a aucune influence.
        La precision est toujours le nb de points de self.cpoints,
        et self.dpoints==self.cpoints"""
        try :
            self._precision = len(self._cpoints)
            self._dpoints = self.cpoints.copy()
        except (AttributeError, TypeError) :
            pass

    def close_(self) :
        u'''Si le qpolygon est EXACTEMENT fermé (eps=0.0), on ne rajoute pas de point.'''
#        self.isclosed = True
        if len(self)<=1 : return False
        if self.isClosed(0.0) :
            return False
        else :
            self.appendPoint(self[0])
            return True
    def echantillonner(self, nbp=0, mode='cpoints', ta=0, tb=1, onlyT=False):
        u"""On doit préciser le mode d'échantillonnage à 'cpoints',
        car dans NSplineSimple, il est à 'segment'
        """
        return super(NPolyLine, self).echantillonner(nbp, mode, ta, tb, onlyT)
############################
    @property
    def points(self):
        return self.cpoints
    @points.setter
    def points(self, points):
        u"""Appelle la property setter.cpoints de NSplineSimple"""
        self.cpoints(points)
    @property
    def dpoints(self):
        if hasattr(self, '_dpoints') :
            del self._dpoints
        return self.cpoints.copy()
#     @points.setter
#     def points(self, points):
#         self.cpoints(points)
#     @property
#     def qpolygon(self):
#         return self.qcpolygon
#     @qpolygon.setter
#     def qpolygon(self, points):
#         self.qcpolygon=points
############################

class NPolygone(NPolyLine):
    """Polyligne fermé"""
    def __init__(self, **dump):
        super(NPolygone, self).__init__(**dump)
        self.close_()

def test1NPolyLine():
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    points = pointsFromFile(filename)
    dump = {'classename': 'NGraphicsShapeItem',
            'color': (100, 100, 100, 0),
            'precision': 1000,
            'methode': ('cubic', 'not-a-knot'),
            'nbpe': 30,
            'verrou': False,
            'name': u'polyline-1',
            'cpoints': points,
#             'cpoints': [[0.0, 0.0], [-1.1313737630844116, 0.11398803442716599], [-0.8625661730766296, -0.015311826020479202]],
            'role': u'polyline',
            'mode': 'segment',
            'position': '(1, 0)'}
    p = NPolyLine(**dump)
    debug(dump)
    debug(p)
    p.cpoints = points
    p.plot(plt, titre='original')
    p = NPolyLine()
    p.load(dump)
    debug(dump)
    debug(p)
    p.plot(plt, titre='load')

def testNPolyLine():
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    p=NPolyLine()
    p=NPolyLine(points=pointsFrom(filename))
    p.appendPoint((5,10))
    print " => Constructeur filename :",p
    points=p.points
    p=NPolyLine(points=points)
    print " => constructeur np.ndarray  :",p
    # p=NPolyLine(points=pointsFrom(p.qcpolygon))
    # print " => constructeur QPolygon  :",p
    print p.scaled(asarray((2,2)))
    #print "constructeur QPolygon :",p
    #print '########## NPolyLine divers (BA, absCurv, extrados, ligne moyenne...) ##########'
    p.removePoint(0)
    curv=absCurv(p)
    dcurv=curv[1:]-curv[:-1]
    print dcurv
    print 'p.absCurv()',absCurv(p)
    print '########## NPolyLine geometrie ##########'
    p.hardScale((2,2))
    print p.points

    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    p = NPolyLine(points=pointsFrom(filename))
    mp = NPolyLine(points=pointsFrom(filename))
    centre=np.asarray((7,0))
    mp.hardRotate(30,centre)#en degrés
    mp.hardScale((0.5,0.5),centre)
#     p = NPolyLine(points=([0,0],[1,0],[1,1]), parent=None)
    from matplotlib import pyplot as plt
    mp.plot(plt).show()

    for k in range(1) :
        k
        XY0=p.points
        XY=mp.points
#         trace('',type(XY))
        my2dPlot((XY0,XY0,XY,XY,np.asarray(centre).reshape((1,2))),equal=True,cosmetic=('b-','bo','r-','ro','g^','r*','r^'))

if __name__=="__main__":
    from spleenconfig import VALIDATION_DIR

    testNPolyLine()
    test1NPolyLine()
