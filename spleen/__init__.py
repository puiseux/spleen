#!/usr/local/bin/python2.7
# encoding: utf-8
__updated__='2019-02-15'
__version__ = "0.1-"+__updated__
__all__ = ['utilitaires', 'profils', 'tests']
from splinesimple import NSplineSimple
from splinecomposee import NSplineComposee
from splineabstraite import NSplineAbstract
from polyline import NPolygone, NPolyLine
from utilitaires import *
from profils import Profil, ProfilNormalise, naca4, naca5
import spleenconfig
