#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe xxx
Description :
@module : programminghelpers.drawing
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 19 janv. 2013
'''
__updated__="2019-01-22"
VERSION = 'v1-'+__updated__
NAME = 'AXile'
TEMPLATE = 'temp'
import math
from utilitaires import trace
#from config import LOG_FILES
#from config import NB_BACKUP
class SplinePrefs(object) :
    precision = 1000#nb points affichage
    tension = 5.0
    poids = 1.0#le poids des noeuds par defaut
    degre = 3
    methode = ('cubic','not-a-knot')
    modes = ['linear', 'cosinus', 'courbure','segment','cpoints']#comment echantillonner
    nbpe = 30 #nb points echantillonnage
    mode = modes[4]
    elagage = 1.0#precision elagage, en mm/m

class ProfilPrefs(object) :
    pouverture = (
                  -1.0, # extrados : position de l'ouverture en %de corde, si <0, ouvext sur l'intrados
                  -5.0  #idem intrados
                  )
    precisionaffichage = 1000 #nb de points a afficher
    profparamnew = (86, 43, 53, 37)#nptprof, iouvext, iouvint, iba
    EPS = 1.0e-5

if __name__ == "__main__" :
    print u'Rien a faire'


