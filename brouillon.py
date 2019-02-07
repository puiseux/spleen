#!/usr/local/bin/python2.7
#coding: utf-8
from setuptools import setup
import spleen
from spleen import NSplineSimple
from spleen.tests import config
from utilitaires import debug
s = NSplineSimple(cpoints=[[0,1],[0,0],[1,1]])
debug(s)