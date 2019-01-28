#!/usr/local/bin/python2.7
# coding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016-2017-2018 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
__updated__="2018-07-01"
from .. import config
import os, sys
from path import Path
from utilitaires import debug, rdebug
import numpy as np
# print os.getcwd()
if __name__ == '__main__' :
#     print os.getcwd()

    cmd = ' '.join([u'python', Path('..','splinesimple.py')])
    rdebug(cmd)
    os.system(' '.join([u'python', Path('..','splinesimple.py')]))
