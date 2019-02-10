#!/usr/local/bin/python2.7
# encoding: utf-8
u'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSplineComposee
Description :
@author:      puiseux
@copyright:   2016 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
import numpy as np
from path import Path
from config import VALIDATION_DIR#, SOURCES_DIR
from splinesimple import NSplineSimple
from config import WORK_DIR, TEST_MODE
# from splinesimple import NSplineSimple
from splinecomposee import NSplineComposee
# from utilitaires.utilitairesdivers import className
from utilitaires import (debug, rdebug)
# debug(config.TEST_MODE)
from utilitaires.utilitairesprofil import computeCordeAndNBA
from utilitaires.lecteurs import pointsFrom
from pprint import pprint
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (20,10)
import config
config.TEST_MODE = True

def testBlocJonc(filename, show=True):
#     show=True
    debug(titre="testBlocJonc()")
    filename = Path(VALIDATION_DIR,'blocjonc.spl')
    name = filename.name
    debug(paragraphe="Construction et dump '%s'"%name)
    with open(filename,'r') as f :
        dump = eval(f.read())['dmodel']
        pprint (dump)
    S = NSplineSimple(**dump)
#     debug(S)
    if show : S.plot(numbers=['1c'])
    filename = Path(VALIDATION_DIR,'blocjonc-compose.spl')
    Sc = NSplineComposee()
#     with open(filename,'r') as f :
#         dump = eval(f.read())
#     debug('dump')
#     pprint (dump)

    Sc.open(filename)
    dump['methode'] = ('ius',3)
#     Sc = NSplineComposee(**dump)
# #     for s in Sc.splines :
# #         s.methode = 'ius',3
# #         s._update()
# #     Sc._update()
    debug(Sc)
#     exit()
#     Sc.plot(numbers=['1c'])
    debug(paragraphe="Split() et dump de '%s'"%name)

#     Sc.split([29,38,39,40,41,50,52,54,55,63,66,74,94,105])
#     debug(Sc)
    if show : Sc.plot()#, numbers=['1c'])
    debug("dump de '%s' apres split()"%name)
    pprint(Sc.toDump())
    filename = Path(WORK_DIR,'blocjonc-compose-work.spl')
    Sc.save(filename)
    name = filename.name
    debug(paragraphe="Construction NSplineComposee a partir de '%s' (=le dump apres split)"%name)
    with open(filename,'r') as f :
        dump = eval(f.read())
    debug('>>>yeap, lecture de dump')
    pprint(dump)
    Sc = NSplineComposee(**dump)
    for i,S in enumerate(Sc.splines) :
        debug("%d-eme composante de '%s'"%(i,name))
        print S
        if show : S.plot(titre='spline-%d'%i, numbers='c')
    if show :
        Sc.plot(numbers=['3c'])
        Sc.plotCourbure()
    debug(titre="Fin testBlocJonc()")

def testParticular(filename, show):
    name = filename.name
    debug(titre='testParticular : (%s)'%name)
    methode = [#demi-cercle
               ('cubic',((1, 0, 3.14), (1, 0, -3.14))),
               ('cubic',((1, 0, 0), (1, 0, 0)))
               ]
    S = NSplineComposee(points=pointsFrom(filename),
                        ruptures=[0,40,-1],
                        methode=methode,
                        precision=[10000,10000],
                        mode=['courbure','courbure'],
                        nbpe=[100,100],
                        name='%s'%(name)
                        )
    debug(S)
    if show : S.plot()
    S.join(1)
    if show : S.plot()
    debug(titre='Fin testParticular : (%s)'%name)

def testNSplineComposee(filename, show=False):
#     show=True
    name = filename.name
    debug(titre='testNSplineComposee : (%s)'%name)
    debug(paragraphe='1. Constructeur vide')
    p = NSplineComposee(name='vide')
    debug(p)
    k=0
    for methode in (
                    [#moche
                        ('cubic','natural'),
                        ('cubic','natural')
                    ],
                    [#demi-cercle
                        ('cubic',((1, 0, 3.14), (1, 0, -3.14))),
                        ('cubic',((1, 0, 0), (1, 0, 0)))
                    ],
                    [#demi-cercle
                        ('cubic','not-a-knot'),
                        ('cubic','not-a-knot')
                    ],
                    [
                        ['cubic',((1, 0, 2), (1, 0, -2))], #extrados
                        ['cubic',((1, 1, 0), (1, 1, 0))]#intrados
                    ],
                    [
                        ('cubic',((2, 0, 0), (1, 0, -5))), #extrados
                        ('cubic',((1, 0, -5), (2, 0, 0))) #intrados
                    ],
                    ):
        debug(paragraphe='methode-%d : %s (%s)'%(k,str(methode), filename.name))
        k+=1
        points=pointsFrom(filename)
        corde, nba = computeCordeAndNBA(points)
        S1 = NSplineComposee(points=points,
                             ruptures=[0,nba,-1],
                             methode=[methode[0],methode[1]],
                             precision=[1000,1000],
                             mode=['courbure','courbure'],
                             nbpe=[100,100],
                             name='%s-%d'%(name,k)
                             )
        debug(S1)
        if show : S1.plot()
#         S2 = NSplineComposee(points=filename, methode=methode, ruptures=[0,2,-1])
    debug(paragraphe='hardScale : %s'%(name))
    S1.hardScale((0.2, 0.2))
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='translate((-1.5,0.05)) : %s'%(name))
    S1.translate((-1.5,0.05))
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='hardRotate(30) : %s'%(name))
    S1.hardRotate(30)
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='symetriser(1) : %s'%(name))
    S1.symetriser(1)
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='symetriser(1) : %s'%(name))
    points = pointsFrom(filename)
    _, nba = computeCordeAndNBA(points)
    S = NSplineComposee(points=points,
                         ruptures=[0,nba,-1],#simple
#                          ruptures=[0,40,-1],#reference
                         methode = [#demi-cercle
#                                     ('cubic','not-a-knot'),
#                                     ('cubic','not-a-knot')
                                ['cubic',((2, 0, 0), (1, 0, -20))], #extrados
                                ['cubic',((1, 0, -20), (2, 0, 0))]#intrados
                                    ],
                         precision=[1000,1000],
                         mode=['courbure','courbure'],
                         nbpe=[20,10],
                         name=name
                         )
    debug(S)

    if show : S.plot( titre='Original '+str(S.methode))

    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    if show : S.plot( titre='S.hardScale((0.2,0.2))')

    S.hardRotate(10)
    debug(S)
    if show : S.plot(titre='S.hardRotate(10)')

    dump = S.toDump()
    S = NSplineComposee(**dump)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    if show : S.plot(titre='S = NSplineComposee(**S.toDump())')

    dump = S.toDump()
    S.load(dump)
    R = S.cpoints[S.ruptures]
    debug(S=S,R=R)
    if show : S.plot(titre='S.load(S.toDump())')

    S.join(1)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    if show : S.plot(titre='join(1)'+str(S.methode))

    S.split(
            nba,
#             ['cubic',((1, 0, 20), (1, 0, -2))], #extrados
#             ['cubic',((1, 20, 0), (1, 1, 0))]#intrados
            )
    R = S.cpoints[S.ruptures]
    debug(split_nba=S,R=R)
    if show :
        S.plot(titre='S.split(40)'+str(S.methode))
        S.plotCourbure()
#     numfig += 1
#     plt.figure(numfig)

#     S.translate((100,0))
#     S.plot(titre='translate((100,0))')
    npt = S.ruptures[1]
    try :
        S.removePoint(npt)
    except NotImplementedError as msg :
        rdebug('normal', msg)
#         raise
    if show : S.plot(titre='removePoint(%d)'%npt)
#     S[4] = (0,0.8)
    try :
        pt = (20,20)
        S.insertPoint(pt)
        if show : S.plot(titre='insertPoint(%s)'%str(pt))
    except ValueError as msg :
        print msg
        raise
    debug(points=S.cpoints.tolist())
    i = len(S)/2
    S.insertPoint(pt,i)
    if show : S.plot(titre='insertPoint(%s,%d)'%(pt,i))
    if show : S.plot(titre='S[4] = (0,0.8)')
    S[0] = (2,0.5)
    if show : S.plot(titre='S[0] = (2,0.5)')
    S[len(S)-1] = S[0]
    if show : S.plot(titre='S[len(S)-1] = S[0]')
#     return
    S.appendPoint((2,1.5))
#     print S
#     print S.cpoints
    if show : S.plot(titre='S.appendPoint((2,1.5))')
    S.insertPoint((2,-2))
    if show : S.plot(titre='S.insertPoint((2,-2))')
    S.insertPoint((6,0))
    print S
    if show : S.plot(titre='S.insertPoint((6,0))')
    S.removePoint(3)
    if show : S.plot(titre='S.removePoint(3)')
    debug(titre='Fin testNSplineComposee : (%s)'%name)


def debugNSplineComposee(filename, show):
    name = filename.name
    debug(titre='debugNSplineComposee : (%s)'%name)
    debug(paragraphe='Constructeur vide')
    S0 = NSplineComposee()
    print S0
    k=0
    for methode in (
#                     [
#                         ['cubic',((1, 0, 2), (1, 0, -2))], #extrados
#                         ['cubic',((1, 1, 0), (1, 1, 0))]#intrados
#                     ],
                    [
                        ('cubic',((2, 0, 0), (1, 0, -5))), #extrados
                        ('cubic',((1, 0, -5), (2, 0, 0))) #intrados
                    ],
                    ):
        k+=1
        points=pointsFrom(filename)
        corde, nba = computeCordeAndNBA(points=points)
        S1 = NSplineComposee(points=points,
                             ruptures=[0,nba,-1],
                             methode=[methode[0],methode[1]],
                             precision=[1000,1000],
                             mode=['courbure','courbure'],
                             nbpe=[50,40],
                             name=name
                             )
        debug(S1)
#         S1.plot()
        debug(paragraphe='elaguer(eps=1, replace=True) : (%s)'%name)
        S1.elaguer(eps=1, replace=True)
        debug(S1)
        if show : S1.plot(titre=u'élaguage')
    debug(paragraphe='echantillonner() (%s)'%name)
    pts = S1.echantillonner()
    debug(echantillon=pts.tolist())
    debug(paragraphe='hardScale((0.2, 0.2)) : (%s)'%name)
    S1.hardScale((0.2, 0.2))
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='translate((-1.5,0.05)) : (%s)'%name)
    S1.translate((-1.5,0.05))
    debug(S1)
    if show : S1.plot()
    debug(paragraphe='hardRotate(30) (%s)'%name)
    debug(S1)
    S1.hardRotate(30)
    if show : S1.plot()
    debug(paragraphe='symetriser(1) (%s)'%name)
    debug(S1)
    S1.symetriser(1)
    if show : S1.plot()
    points = pointsFrom(filename)
    corde, nba = computeCordeAndNBA(points)
    debug(paragraphe='reconstruction (%s)'%name)
    S = NSplineComposee(points=points,
                        ruptures=[0,nba,-1],#reference
                        methode = [['cubic',((2, 0, 0), (1, 0, -20))], #extrados
                                   ['cubic',((1, 0, -20), (2, 0, 0))]],#intrados
                        precision=[1000,1000],
                        mode=['courbure','courbure'],
                        nbpe=[50,50],
                        name=name
                        )
    debug(S)
#     S[4] = (95,5)
#     R = S.cpoints[S.ruptures]
#     debug(original=S,R=R)

    if show : S.plot(titre='Original '+str(S.methode))
    debug(S)

    debug(paragraphe='hardScale((0.2,0.2), centre=np.asarray((0,0))) (%s)'%name)
    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    if show : S.plot(titre='S.hardScale((0.2,0.2))')

    debug(paragraphe='S.hardRotate(10) (%s)'%name)
    S.hardRotate(10)
    debug(S)
    if show : S.plot(titre='S.hardRotate(10)')


    debug(paragraphe='S.toDump() (%s)'%name)
    dump0 = S.toDump()

    pprint(dump0)
    debug(paragraphe='NSplineComposee(**dump0) (%s)'%name)
    S = NSplineComposee(**dump0)
    debug(join_1=S)
    debug(Rupture=S.cpoints[S.ruptures].tolist())
    if show : S.plot(titre='S = NSplineComposee(**S.toDump())')

    debug(paragraphe='dump1 = S.toDump() (%s)'%name)
    dump1 = S.toDump()
    debug(dump_equal_dump1=dump0==dump1)
    debug(paragraphe='S.load(dump1) (%s)'%name)
    S.load(dump1)
    R = S.cpoints[S.ruptures]
    debug(S)
    debug(R=R)
    dump2 = S.toDump()
    debug(dump1_equal_dump2=dump1==dump2)

    if show : S.plot(titre='S.load(**S.toDump())')

    debug(paragraphe='S.join(1) (%s)'%name)
    S.join(1)
    debug(S)
    debug(Ruptures=S.cpoints[S.ruptures])
    if show : S.plot(titre='join(1)'+str(S.methode))

    debug(paragraphe='S.split(len(S)/2) (%s)'%name)
    S.split(len(S)/2)
    debug(S)
    debug(ruptures=S.cpoints[S.ruptures])
    if show : S.plot(titre='S.split(40)'+str(S.methode))

    debug(paragraphe='removePoint(rupture[1]) (%s)'%name)
    npt = S.ruptures[1]
    try : S.removePoint(npt)
    except NotImplementedError as msg :
        debug('normal', msg)
#         raise
    debug(S)
    if show : S.plot(titre='removePoint(%d)'%npt)

    pt = (20,20)
    debug(paragraphe='S.insertPoint(%s) %s'%(str(pt),name))
    try :
        S.insertPoint(pt)
        if show : S.plot(titre='insertPoint(%s)'%str(pt))
    except ValueError as msg :
        print msg
        raise
    debug(S)
    debug(points=S.cpoints.tolist())

    i = len(S)/2
    debug(paragraphe='S.insertPoint(pt,%d) %s'%(i,name))
    S.insertPoint(pt,i)
    debug(S)
    if show : S.plot(titre='insertPoint(%s,%d)'%(pt,i))

    debug(paragraphe='S[0] = (2,0.5) (%s)'%name)
    S[0] = (2,0.5)
    debug(S)
    if show : S.plot(titre='S[0] = (2,0.5)')
    debug(paragraphe='S[len(S)-1] = S[0] (%s)'%name)
    S[len(S)-1] = S[0]
    debug(S)
    if show : S.plot(titre='S[len(S)-1] = S[0]')
#     return
    debug(paragraphe='S.appendPoint((2,1.5)) (%s)'%name)
    S.appendPoint((2,1.5))
    debug(S)
    if show : S.plot(titre='S.appendPoint((2,1.5))')

    debug(paragraphe='S.insertPoint((2,-2)) (%s)'%name)
    S.insertPoint((2,-2))
    debug(S)
    if show : S.plot(titre='S.insertPoint((2,-2))')

    debug(paragraphe='S.insertPoint((6,0)) (%s)'%name)
    S.insertPoint((6,0))
    debug(S)
    if show : S.plot(titre='S.insertPoint((6,0))')

    debug(paragraphe='S.removePoint(3) (%s)'%name)
    S.removePoint(3)
    debug(S)
    if show : S.plot(titre='S.removePoint(3)')
    debug(paragraphe='S._update() (%s)'%name)
    S._update()
    debug(S)
    if show : S.plot(titre='S._update()')
#     rdebug('**********************\ndouteux, a verifier\n**********************')
    debug(titre='Fin debugNSplineComposee : %s'%filename.name)


def testMain(show=False):
    files = []
    files.append(Path(VALIDATION_DIR,'unenervure2d.gnu'))
    files.append(Path(VALIDATION_DIR,'diamirE.spl'))
    files.append(Path(VALIDATION_DIR,'E-diamirE.spl'))
    files.append(Path(VALIDATION_DIR,'E-shark.spl'))
    files.append(Path(VALIDATION_DIR,'shark.spl'))
    files.append(Path(VALIDATION_DIR,'naca2415.spl'))
    files.append(Path(VALIDATION_DIR,'spline-0#.pkl'))
    files.append(Path(VALIDATION_DIR,'spline-0#.spl'))
    files.append(Path(VALIDATION_DIR,'spline-0.spl'))
    if 1 : testBlocJonc(Path(VALIDATION_DIR,'blocjonc.spl'),show=show)
#     return
    if 1 : testParticular(files[6], show=show)
    for k,filename in enumerate(files[:]) :
        if 1:debugNSplineComposee(filename, show=show)#k in [6,7])
        if 1:testNSplineComposee(filename, show=show)#k in [6,7])
if __name__=="__main__":
    testMain()
