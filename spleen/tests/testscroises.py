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

from profil import Profil
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
    print P
# exit()

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
#     debug('  >> dump sauf cpoints')
#     pprint(dump)
    print P

def testsDivers(filename, show=True):
    name = filename.name
    debug(filename.name)
    S = NSplineComposee()
    S.open(filename)
#     debug('Avant split',S.name,[s.name for s in S.splines])
    a, b = len(S)/3, 2*len(S)/3
    S.split([a,b])
#     debug('Apres split',S.name,[s.name for s in S.splines])
    debug(S)
    debug(epoints=S.epoints.tolist())
    
config.TEST_MODE = False
if 1 : testsDivers(filename=filenames[3])
if 1 : testConstructeursVides()
for filename in filenames[:] :
    if 1 : tc(filename, show=False)
    if 1 : tmg(filename, show=False)
    if 1 : tml(filename, show=False)
    if 1 : testConversions(filename)
