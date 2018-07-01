#!/usr/local/bin/python2.7
# encoding: utf-8
from numpy import asarray
u'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe ProfilNormalise

@author:     Michel Le Berre, Pierre Puiseux
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 31 Jan 2013
'''

from pprint import pprint
from utilitaires import Path
from profil import Profil
from utilitaires import (debug, rdebug,)
from lecteurs import pointsFrom#, qpolygonFrom
from profilnormalise import ProfilNormalise
from config import DATA_DIR,VALIDATION_DIR, RUNS_DIR, WORK_DIR

def testProfilNormalise():
    from matplotlib import pyplot as plt
    dump = {'classename': 'ProfilNormalise',
           'cpoints': [[1.0, 1.3552527156068805e-20],
                       [0.7084010949941603, 0.04443863894727284],
                       [0.4284347004138754, 0.08826190999847804],
                       [0.23674398594773688, 0.10760076664620725],
                       [0.09207463882331422, 0.08905579820615639],
                       [0.01637063370876783, 0.03853089493038632],
                       [0.0, 0.0],
                       [0.002371876527552363,-0.007975320245593812],
                       [0.013417354346383206,-0.014834590263683093],
                       [0.03378250361398215, -0.02084196973495454],
                       [0.05401765843464571, -0.038677498922310316],
                       [0.07141124918455544, -0.04659119673690394],
                       [0.09262143929059875, -0.05112526462431936],
                       [0.25257255205423007, -0.06300340478020656],
                       [0.4221571793400746, -0.057482467002892425],
                       [0.7519294850346585, -0.02637922651126774],
                       [1.0, 1.3552527156068805e-20]],
           'graphicsstate': {},
           'methode': [('cubic', ((2, 0.0, 0.0), (1, 0, -1.0))),
                       ('cubic', ((1, 0, -1.25), (2, 0.0, 0.0)))],
           'mode': ['courbure', 'courbure'],
           'name': u'E-P01',
           'nbpe': [41, 46],
           'pouverture': (-3.0, -3.0),
           'precision': [1000, 1000],
           'role': 'gprofil',
           'ruptures': [0, 6, 16]}
    dump2 = {'classename': 'ProfilNormalise',
           'cpoints': [[1.0, 0.0],
                       [0.7071902169056651, 0.04483465725889014],
                       [0.42811066069700837, 0.0887659973792878],
                       [0.236294834409799, 0.10816259984228063],
                       [0.1313314565503818, 0.10073384297862178],
                       [0.027883090329988702, 0.051699110582146095],
                       [0.005415505333442232, 0.021367022007105285],
                       [0.0, 0.0],
                       [0.023078710408487593,-0.016725499215111874],
                       [0.028786127653164142,-0.018032064869529244],
                       [0.04310599499863427, -0.029182612249222292],
                       [0.06478092132361202, -0.042961243975068356],
                       [0.092577693242155, -0.0501852851518773],
                       [0.14921811143141142, -0.057938824515836304],
                       [0.25429096518199634, -0.06245277820071929],
                       [0.35087672191875136, -0.061065851086376605],
                       [0.4226884306798086, -0.057002810226840214],
                       [0.6483460163197147, -0.037385199262451686],
                       [0.7531716291754094, -0.02605979826089128],
                       [1.0, 0.0]],
           'graphicsstate': {},
           'methode': [('cubic', ((2, 0.0, 0.0), (1, 0, -1.0))),
                       ('cubic', ((1, 0, -1.15), (2, 0.0, 0.0)))],
           'mode': ['courbure', 'courbure'],
           'name': u'E-P06-5',
           'nbpe': [41, 46],
           'pouverture': (-2.88, -4.06),
           'precision': [1000, 1000],
           'role': 'gprofil',
           'ruptures': [0, 7, 19]}
    p = Profil(**dump)
    p.hardRotate(30)
    p.plot(plt, titre=u'%s : profil original, à normaliser'%p.name)
    dump = p.toDump()
    p = ProfilNormalise(**dump)
    p.iba = 40
    p.iouverture = (50, 60)
    print p.name, p.epoints.shape
    pprint (p.epoints.tolist())
#     print p
#     plt.figure(1)
    p.plot(plt, show=False, titre=u'%s normalisé'%p.name)
    plt.draw()#

    p = ProfilNormalise(**dump2)
    p.iba = 40
    p.iouverture = (50, 60)
    print p.name, p.epoints.shape
    pprint (p.epoints.tolist())
#     print p
#     plt.figure(1)
    p.plot(plt, show=False, titre=u'%s normalisé'%p.name)
    plt.draw()#


    plt.show()
    return
    filename=Path(VALIDATION_DIR,'1.gnu')
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename=Path(VALIDATION_DIR,'faial2.dxf')
    filename=Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    print " => Constructeur presque vide, puis append(5,10):\n    +++++++++++++++++++++++++++++++++++++"
    p = ProfilNormalise()
    print p
    print " => Constructeur filename :\n    +++++++++++++++++++++"
    p = ProfilNormalise(points=pointsFrom(filename))
    print p
    cpoints = p.cpoints
#     print " => constructeur QPolygon  :\n    +++++++++++++++++++++"
#     try : p = ProfilNormalise(points=p.qpolygon)
#     except TypeError as msg :
#         print"    => " , msg
    print " => constructeur np.ndarray  :\n    +++++++++++++++++++++++"
#     print(cpoints.shape)
    p = ProfilNormalise(points=cpoints)
    print p.normalite()
#     print " => constructeur recopie  :\n    +++++++++++++++++++++"
#     p.iouverture = p.nba+4,p.nba+5
#     q = ProfilNormalise(profil=p)
#     print q
    print '########## ProfilNormalise geometrie ##########'
    print p.verification()
    print p.normalite()
    p.normalise()
    for msg in p.verification() :
        print msg
    p.iouverture=p.nba+2,p.nba+3
    print p.normalite()
#     exit()
    print p.cpoints
    print p.normalite()
    filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
    p = ProfilNormalise(points=pointsFrom(filename))
    p.verification()
    mp = ProfilNormalise(points=pointsFrom(filename))
    dump = p.toDump()
#     debug()
#     pprint(dump)
    p = ProfilNormalise(**dump)
    p.verification()

    debug(p.normalite())
    p.plot(plt, titre='p = Profil(**p.toDump())')
    p.plotCourbure()#(plt, titre='p = Profil(**p.toDump())')
    plt.show()

    print p.normalite()
    p.normalise()
    print p.normalite()

    print '########## ProfilNormalise opérations interdites ##########'
    print '    #hardScale'
    try : p.hardScale((2,2))
    except RuntimeError as msg : print msg
    print p.normalite()
#     return
    centre=asarray((7,0))
    print '    #hardRotate'
    try : mp.hardRotate(30,centre)
    except RuntimeError as msg : print msg
    print p.normalite()
    print '    #appendPoint'
    try : p.appendPoint((5,10))
    except RuntimeError as msg : print msg
    print p.normalite()
    print '    #p[0] = (1,2)'
    try : p[0] = (1,2)
    except ValueError as msg : print msg
    print p.normalite()
    print '    #insertPoint'
    try : p.insertPoint((1.1, 0.2))
    except Exception as msg : print msg
    print p.normalite()
    print '################## FIN testProfilNormalise'
#     exit()
#     p = ProfilNormalise(points=([0,0],[1,0],[1,1]), parent=None)
def testDivers():
#
    from numpy import asarray as array
    PN = ProfilNormalise()
    debug(PN)
    filename = Path(RUNS_DIR,'D2V1-Pb.txt')
    filename = Path(VALIDATION_DIR,'P0.spl')
    with open(filename) as f :
        dump = eval(f.read())
    print dump
    p = ProfilNormalise(**dump)
    print p

def testMain():
    p=ProfilNormalise()
    debug(p)
    testDivers()
    testProfilNormalise()
    print '################## FIN main #################'

if __name__=="__main__":
    testMain()
#     sys.exit(app.exec_())
