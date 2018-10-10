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
# import sys,os,math
import numpy as np
# from numpy import asarray as array
from path import Path
from config import VALIDATION_DIR#, SOURCES_DIR
# from splineabstraite import NSplineAbstract
# from splinesimple import NSplineSimple
from lecteurs import pointsFrom
from utilitaires import (debug, rdebug)
from splinecomposee import NSplineComposee

def testNSplineComposee():
#     def testSplit(s):
#         s.split(40,
#                 ['cubic',((1, 0, 20), (1, 0, -2))], #extrados
#                 ['cubic',((1, 20, 0), (1, 1, 0))]#intrados
#                 )
# #         debug(s)
#     def testJoin(s):
#         s.join(1)
    from matplotlib import pyplot as plt

    p = NSplineComposee(name='vide')
    debug(p=p)
    # for point in p.qcpolygon :
    #     print point
    # for point in p.qdpolygon :
    #     print point
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename=Path(VALIDATION_DIR,'1.gnu')
    filename=Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    filename=Path(VALIDATION_DIR,'faial2.dxf')
    filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename=Path(VALIDATION_DIR,'simple','intrados.gnu')
    filename=Path(VALIDATION_DIR,'simple','anguleux.gnu')
    filename=Path(VALIDATION_DIR,'simple','extrados.gnu')
    filename=Path(VALIDATION_DIR,'simple','demi-cercle.gnu')
    filename=Path(VALIDATION_DIR,'simple','profil.gnu')
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename=Path(VALIDATION_DIR,'reference.pts')

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
        k+=1
        S1 = NSplineComposee(points=pointsFrom(filename),
                             ruptures=[0,4,-1],
                             methode=[methode[0],methode[1]],
                             precision=[1000,1000],
                             mode=['courbure','courbure'],
                             nbpe=[20,10],
                             name='Profil-%d'%k
                             )
        print S1
        S1.plot(plt)
#         S2 = NSplineComposee(points=filename, methode=methode, ruptures=[0,2,-1])
    S1.hardScale((0.2, 0.2))
    S1.plot(plt)
    S1.translate((-1.5,0.05))
    S1.plot(plt)
    S1.hardRotate(30)
    S1.plot(plt)
    S1.symetriser(1)
    S1.plot(plt)
    S = S1
    S = NSplineComposee(points=pointsFrom(filename),
                         ruptures=[0,4,-1],#simple
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
                         name='Profil-%d'%k
                         )
    debug(S)
#     S[4] = (95,5)
#     R = S.cpoints[S.ruptures]
#     debug(original=S,R=R)

    S.plot(plt, titre='Original '+str(S.methode))

    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    S.plot(plt, titre='S.hardScale((0.2,0.2))')

    S.hardRotate(10)
    S.plot(plt, titre='S.hardRotate(10)')


    dump = S.toDump()
    S = NSplineComposee(**dump)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='S = NSplineComposee(**S.toDump())')

    dump = S.toDump()
    S.load(dump)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='S.load(S.toDump())')

    S.join(1)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='join(1)'+str(S.methode))

    S.split(
            len(S)/2,
#             ['cubic',((1, 0, 20), (1, 0, -2))], #extrados
#             ['cubic',((1, 20, 0), (1, 1, 0))]#intrados
            )
    R = S.cpoints[S.ruptures]
    debug(split_40=S,R=R)
    S.plot(plt, titre='S.split(40)'+str(S.methode))
#     numfig += 1
#     plt.figure(numfig)

#     S.translate((100,0))
#     S.plot(plt, titre='translate((100,0))')
    npt = S.ruptures[1]
    try : S.removePoint(npt)
    except NotImplementedError as msg :
        rdebug(msg)
    S.plot(plt, titre='removePoint(%d)'%npt)
#     S[4] = (0,0.8)
    try :
        pt = (20,20)
        S.insertPoint(pt)
        S.plot(plt, titre='insertPoint(%s)'%str(pt))
    except ValueError as msg :
        print msg
    debug(points=S.cpoints)
    i = len(S)/2
    S.insertPoint(pt,i)
    S.plot(plt, titre='insertPoint(%s,%d)'%(pt,i))
    S.plot(plt, titre='S[4] = (0,0.8)')
    S[0] = (2,0.5)
    S.plot(plt, titre='S[0] = (2,0.5)')
    S[len(S)-1] = S[0]
    S.plot(plt, titre='S[len(S)-1] = S[0]')
#     return
    S.appendPoint((2,1.5))
#     print S
#     print S.cpoints
    S.plot(plt, titre='S.appendPoint((2,1.5))')
    S.insertPoint((2,-2))
    S.plot(plt, titre='S.insertPoint((2,-2))')
    S.insertPoint((6,0))
    print S
    S.plot(plt, titre='S.insertPoint((6,0))')
    S.removePoint(3)
    S.plot(plt, titre='S.removePoint(3)')
    return plt.show()
    # plt.show()
    # return


def debugNSplineComposee():
    from lecteurs import pointsFrom
    S0 = NSplineComposee()
    print S0
#     return
#     def testSplit(s):
#         s.split(40,
#                 ['cubic',((1, 0, 20), (1, 0, -2))], #extrados
#                 ['cubic',((1, 20, 0), (1, 1, 0))]#intrados
#                 )
# #         debug(s)
#     def testJoin(s):
#         s.join(1)
    from matplotlib import pyplot as plt

    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename=Path(VALIDATION_DIR,'1.gnu')
    filename=Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    filename=Path(VALIDATION_DIR,'faial2.dxf')
    filename=Path(VALIDATION_DIR,'simple','intrados.gnu')
    filename=Path(VALIDATION_DIR,'simple','anguleux.gnu')
    filename=Path(VALIDATION_DIR,'simple','extrados.gnu')
    filename=Path(VALIDATION_DIR,'simple','demi-cercle.gnu')
    filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename=Path(VALIDATION_DIR,'simple','profil.gnu')
    filename=Path(VALIDATION_DIR,'reference0.pts')
    filename=Path(VALIDATION_DIR,'unenervure2d.gnu')

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
        S1 = NSplineComposee(points=pointsFrom(filename),
#                              ruptures=[0,4,-1],#simple
                             ruptures=[0,40,-1],#reference
                             methode=[methode[0],methode[1]],
                             precision=[1000,1000],
                             mode=['courbure','courbure'],
                             nbpe=[50,40],
                             name='Profil-%d'%k
                             )
        print S1
#         S1.plot(plt)
        S1.elaguer(eps=1, replace=True)
        S1.plot(plt,titre=u'élaguage')
#         S2 = NSplineComposee(points=filename, methode=methode, ruptures=[0,2,-1])
    pts = S1.echantillonner()
    debug(echantillon=pts.tolist())
    S1.hardScale((0.2, 0.2))
    S1.plot(plt)
    S1.translate((-1.5,0.05))
    S1.plot(plt)
    S1.hardRotate(30)
    S1.plot(plt)
    S1.symetriser(1)
    S1.plot(plt)
    S = S1
    S = NSplineComposee(points=pointsFrom(filename),
#                          ruptures=[0,4,-1],#simple
                        ruptures=[0,40,-1],#reference
                         methode = [#demi-cercle
#                                     ('cubic','not-a-knot'),
#                                     ('cubic','not-a-knot')
                                ['cubic',((2, 0, 0), (1, 0, -20))], #extrados
                                ['cubic',((1, 0, -20), (2, 0, 0))]#intrados
                                    ],
                         precision=[1000,1000],
                         mode=['courbure','courbure'],
                         nbpe=[50,50],
                         name='Profil-%d'%k
                         )
    debug(S)
#     S[4] = (95,5)
#     R = S.cpoints[S.ruptures]
#     debug(original=S,R=R)

    S.plot(plt, titre='Original '+str(S.methode))

    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    S.plot(plt, titre='S.hardScale((0.2,0.2))')

    S.hardRotate(10)
    S.plot(plt, titre='S.hardRotate(10)')


    dump = S.toDump()
    S = NSplineComposee(**dump)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='S = NSplineComposee(**S.toDump())')

    dump = S.toDump()
    S.load(dump)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='S.load(**S.toDump())')

    S.join(1)
    R = S.cpoints[S.ruptures]
    debug(join_1=S,R=R)
    S.plot(plt, titre='join(1)'+str(S.methode))

    S.split(
            len(S)/2,
#             ['cubic',((1, 0, 20), (1, 0, -2))], #extrados
#             ['cubic',((1, 20, 0), (1, 1, 0))]#intrados
            )
    R = S.cpoints[S.ruptures]
    debug(split_40=S,R=R)
    S.plot(plt, titre='S.split(40)'+str(S.methode))
#     numfig += 1
#     plt.figure(numfig)

#     S.translate((100,0))
#     S.plot(plt, titre='translate((100,0))')
    npt = S.ruptures[1]
    try : S.removePoint(npt)
    except NotImplementedError as msg :
        rdebug(msg)
    S.plot(plt, titre='removePoint(%d)'%npt)
#     S[4] = (0,0.8)
    try :
        pt = (20,20)
        S.insertPoint(pt)
        S.plot(plt, titre='insertPoint(%s)'%str(pt))
    except ValueError as msg :
        print msg
    debug(points=S.cpoints.totxt())
    i = len(S)/2
    S.insertPoint(pt,i)
    S.plot(plt, titre='insertPoint(%s,%d)'%(pt,i))
    S.plot(plt, titre='S[4] = (0,0.8)')
    S[0] = (2,0.5)
    S.plot(plt, titre='S[0] = (2,0.5)')
    S[len(S)-1] = S[0]
    S.plot(plt, titre='S[len(S)-1] = S[0]')
#     return
    S.appendPoint((2,1.5))
#     print S
#     print S.cpoints
    S.plot(plt, titre='S.appendPoint((2,1.5))')
    S.insertPoint((2,-2))
    S.plot(plt, titre='S.insertPoint((2,-2))')
    S.insertPoint((6,0))
    print S
    S.plot(plt, titre='S.insertPoint((6,0))')
    S.removePoint(3)
    S._update()
    print 'removepoint(3)', S
    S.plot(plt, titre='S.removePoint(3)')
    rdebug('**********************\ndouteux, a verifier\n**********************')
#     return plt.show()
    plt.show()
    return

def testMain():
    debugNSplineComposee()
    testNSplineComposee()

if __name__=="__main__":
    testMain()
