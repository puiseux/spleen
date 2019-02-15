#!/usr/local/bin/python2.7
# encoding: utf-8
u'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe ProfilNormalise

@author:     Michel Le Berre, Pierre Puiseux
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 31 Jan 2013
'''

from numpy import asarray, abs, vstack
from pprint import pprint
from utilitaires import Path
from profil import Profil, ProfilNormalise, ProfsParam1
from utilitaires import (debug, rdebug,)
from utilitaires import pointsFrom#, qpolygonFrom
from config import VALIDATION_DIR, RUNS_DIR
from matplotlib import pyplot as plt

def testProfilNormalise(filename, show):
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
    p = Profil()
    p.open(filename)
    debug(p)
#     p = Profil(**dump)
#     p.hardRotate(30)
    if show : p.plot(titre=u'%s : profil original, à normaliser'%p.name)
    dump = p.toDump()
    p = ProfilNormalise(**dump)
    p.iba = 40
    p.iouverture = (50, 60)

    print p.name, p.epoints.shape
    debug(p.epoints.tolist())
    if show : p.plot(titre=u'%s normalisé'%p.name)

    p = ProfilNormalise(**dump2)
    p.iba = 40
    p.iouverture = (50, 60)
    print p.name, p.epoints.shape
    debug(p.epoints.tolist())
#     print p
#     plt.figure(1)
    if show : p.plot(titre=u'%s normalisé'%p.name)


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
def testDivers(filename, show):
    name = filename.name
    debug(titre=name)
#     with open(filename) as f :
#         dump = f.read()
#         debug(dump)
#         dump = eval(f.read())
#     print dump
    p = ProfilNormalise()
    p.open(filename)
    print p
    debug(p.echantillonner().tolist())


def testPinces(filename, show):
    """
    Les couples de pourcentage concernent les % gauche et droite de
    BF_xt, BA_ext, BA_int, BF_int
    On bouge 4 points de l'échantillonnage du profil pour les positionner à
        - pbai% du BA (intrados)
        - pbfi% du BF(intrados)
        - pbae% du BA (extrados)
        - pbfe% du BF(extrados)
    """
    pourcentages=((10,10.1), (5,5), (7.3,7.3), (5.5,5.5))
    P = ProfilNormalise(cpoints=pointsFrom(filename), name=filename.name)
    pp = ProfsParam1(nptprof=150, iouvext=70, iouvint=76, iba=60)
    P.profparam = pp
    debug(pp=pp)
    debug(P)
#     debug(P.pouverture[1])
    dpoints = P.dpoints
    Xd, Yd = dpoints[:,0], dpoints[:,1]
    epoints = P.epoints
#     techext0, techint0 = P.techext.copy(), P.techint.copy()#Pour debug
    Xe, Ye = epoints[:,0], epoints[:,1]
    plt.plot(Xd,Yd,'k-', label='dpoints')
    plt.plot(Xe,Ye,'go', label='epoints')
    def moveOnePoint(tech, t0, forbiden=[]):
        """
        Dans le tableau tech (les t d'échantillonage de la spline),
        on repère le point le plus proche de t0, son indice est k c'est l'indice ou
        A = abs(tech[k]-t0) est minimal
        SI on a le droit, (k n'appartient pas à forbiden) on modifie : tech[k] <- t0
        SI on n'a le droit, (k appartient à forbiden) on prend k <- k-1 ou k <- k+1 et on modifie tech[k]"""
#         tech = P.techext #les t échantillonnage, extrados théorique
#         t = P._getT(pc)         # le t du point de pince sur son trados
        A = abs(tech-t0)
        k = A.argmin()# le t echantillon actuel le plus proche de t
        if k in forbiden :
            A[k] = 2*max(A)
        k = A.argmin()
        tech[k] = t0             #On modifie le t echantillonnage dans P
        return k
    numeros = []
    for i, (pg, pd) in enumerate(pourcentages) :#pg=%gauche, pd=%droite
        #BF si i = 0 ou 3
        #BA si i = 1 ou 2
        #ext si i = 0 ou 1
        #int si i = 2 ou 3
        if i in (0,1) : #extrados
            tech = P.techext #les t échantillonnage, extrados théorique
            sgn = 1
        else : #intrados
            tech = P.techint
            sgn = -1
        if i in(0,3) :#BF
            pg = (100 - pg)*sgn
            pd = (100 - pd)*sgn
        else :#BA
            pg =  (pg - P.pouverture[1])*sgn
            pd =  (pd - P.pouverture[1])*sgn
        tg = P._getT(pg)         # le t du point de pince sur son trados
        k1 = moveOnePoint(tech, tg, forbiden=[])
        if pd == pg :
            k2 = k1
        else :
            td = P._getT(pd)     # le t du point de pince sur son trados
            k2 = moveOnePoint(tech, td, forbiden=[k1])
        numeros.append((k1,k2))

    P._epoints = vstack((P.splines[0](P.techext), P.splines[1](P.techint)))
    epoints = P.epoints
    Xe, Ye = epoints[:,0], epoints[:,1]
    plt.plot(Xe,Ye,'rx', label='epoints-new')
#     print 'n'.join((
#     print 'pourcentages=%s'%str(pourcentages),
#     print 'numeros points=%s'%str(numeros)
#     ))
    plt.title('pourcentages=%s'%str(pourcentages)+'\n'+'numeros points=%s'%str(numeros))
#     debug(deltaTechext=P.techext-techext0)
#     debug(deltaTechint=P.techint-techint0)
    print(numeros)
    plt.legend()
    plt.axis('equal')
#     plt.title(u"4 points de l'échantillonnage déplacés \n(BF Ext=%g%%, BA Ext=%g%%, BA Int=%g%%, BF Int%g%%)"%(pbfe,pbae,pbai,pbfi))
    plt.show()
    return

def testMain():
    files = [
            Path(VALIDATION_DIR,'splinesimple-86pts.spl'),
            Path(VALIDATION_DIR,'profilnormalise-86pts.spl'),
            Path(VALIDATION_DIR,'blocjonc-splinesimple.spl'),
            Path(VALIDATION_DIR,'shark-profilnormalise-86pts.spl'),
            Path(VALIDATION_DIR,'diamirE-profilnormalise-86pts.spl'),
            Path(VALIDATION_DIR,'profilnormalise-21pts.spl'),
            Path(VALIDATION_DIR,'splinesimple-21pts.spl'),
            Path(VALIDATION_DIR,'shark-profilnormalise-26pts.spl'),
            Path(VALIDATION_DIR,'diamirE-profilnormalise-24pts.spl'),
            Path(VALIDATION_DIR,'NACA2415-100pts.spl'),
            Path(VALIDATION_DIR,'points-86pts.gnu'),
            Path(VALIDATION_DIR,'splinesimple-86pts.pkl'),
            ][::-1]
    
    p = ProfilNormalise()
    debug('constructeur vide\n',p)
    for filename in files[:2] :
        if 1 : testPinces(filename, show=True)
        if 1 : testDivers(filename, show=True)
        if 1 : testProfilNormalise(filename, show=True)
    print '################## FIN main #################'

if __name__=="__main__":
    testMain()
#     sys.exit(app.exec_())
