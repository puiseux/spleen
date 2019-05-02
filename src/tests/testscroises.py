#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on Feb 1, 2019

@author: puiseux
'''
from numpy.linalg import norm
from splinesimple import NSplineSimple
from matplotlib import pyplot as plt
from utilitaires.utilitairesprofil import computeCordeAndNBA
from pprint import pprint
from testsplinesimple import testDivers
plt.rcParams["figure.figsize"] = (20,10)
from utilitaires.lecteurs import pointsFrom, Path, pointsFromFile
from config import VALIDATION_DIR
from splinecomposee import NSplineComposee
from utilitaires import (debug, XY, rdebug, dist2, className, diff)
# from numpy import (asarray, min, max, argmin, linspace,sqrt,ndarray,zeros,ones,
#           vstack, hstack, nan, NaN, empty, where, isnan, zeros_like, isfinite)
import config

from profils import Profil, ProfilNormalise
from testsprofil import testProfil as tp
from testsprofil import testSaveAndOpen as tso

from testsplinecomposee import testConstructeurs as tc
from testsplinecomposee import testMethodesGlobales as tmg
from testsplinecomposee import testMethodesLocales as tml

filenames = [
             Path(VALIDATION_DIR,'points-86pts.gnu'),
             Path(VALIDATION_DIR,'splinesimple-86pts.pkl'),
             Path(VALIDATION_DIR,'shark-profilnormalise-26pts.spl'),
             Path(VALIDATION_DIR,'diamirE-profilnormalise-24pts.spl'),
             
            ]
def testConstructeursVides():
    debug(titre='P = NSplineSimple()')
    P = NSplineSimple()
    print P
    
    debug(titre='P = NSplineComposee()')
    P = NSplineComposee()
    print P
    
    debug(titre='P = Profil()')
    P = Profil()
    pprint(P.Default().dump)
    print P.Default().dump
    print P
    
    debug(titre='P = ProfilNormalise()')
    P = ProfilNormalise()
    pprint(P.Default().dump)
    print P

def testConversions(filename):
    name = filename.name
    debug(titre='testConversions(%s)'%name)

    debug(paragraphe='NSplineSimple.open(%s)'%name)
    P = NSplineSimple()
    P.open(filename)
#     dump = P.toDump()
#     dump.pop('cpoints')
    print P

    debug(paragraphe='NSplineComposee.open(%s)'%name)
    P = NSplineComposee()
    P.open(filename)
    print P

    debug(paragraphe='Profil.open(%s)'%name)
    P = Profil()
    P.open(filename)
    print P
    
    debug(paragraphe='ProfilNormalise.open(%s)'%name)
    P = ProfilNormalise()
    P.open(filename)
    print P
    
def testsDivers(filename, show=True):
    name = filename.name
    debug(titre='testsDivers(%s)'%name)
    S = NSplineComposee()
    S.open(filename)
#     debug('Avant split',S.name,[s.name for s in S.splines])
    a, b = len(S)/3, 2*len(S)/3
    S.split([a,b])
#     debug('Apres split',S.name,[s.name for s in S.splines])
    debug(S)
    debug(epoints=S.epoints.tolist())
    P = ProfilNormalise()
    P.open(filename)
    if show:P.plot(numbers=['2c'])
    P[12] = (0.024,-0.03)
    if show:P.plot(numbers=['2c'])
    debug('verif')
    for v in P.verification(): print v
    debug(paragraphe='_getT()')
    pourcent = -2.5
    t0, r0 = P._getT(pourcent, nbit=True)
    t1, r1 = P._getT(pourcent, nbit=True, t0=t0*2)
#     t1, r1 = P._getT(pourcent, nbit=True, t0=t0*10)
    print '==> t0 = %.3g'%t0
    print type(r0)
    print '==> t1 = %.3g'%t1
    print r1
    ispl = 0 if pourcent>0 else 1
    p0 = P.splines[ispl](t0) 
    p1 = P.splines[ispl](t1)
    print 'p0 =', p0
    print 'p1 =', p1
    if show : 
        more = [(p0[0],p0[1],'rv','p0'),(p1[0],p1[1],'bv','p1')]
        P.plot(more=more)
    
config.TEST_MODE = False
if 1 : testConstructeursVides()
if 1 : testsDivers(filenames[3], True)
for filename in filenames[:] :
    if 0 : tc(filename, show=False)
    if 0 : tmg(filename, show=False)
    if 0 : tml(filename, show=False)
    if 0 : testConversions(filename)
    if 0 : tso(filename, show=False)