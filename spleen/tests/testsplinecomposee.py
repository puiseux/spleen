#!/usr/local/bin/python2.7
# encoding: utf-8
from filenames import filenames
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
from utilitaires import (debug, rdebug,dictsAreNotEqual)
# debug(config.TEST_MODE)
from utilitaires.utilitairesprofil import computeCordeAndNBA
from utilitaires.lecteurs import pointsFrom
from pprint import pprint
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (20,10)
import config
config.TEST_MODE = True

def testBlocJonc(filename, show):
#     show=True
    debug(titre="testBlocJonc()")
#     filename = Path(VALIDATION_DIR,'blocjonc.spl')
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

def testParticulier(filename, show):
    name = filename.name
    debug(titre='testParticulier : (%s)'%name)
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
                        name=name
                        )
    debug(S)
    if show : S.plot()
    S.join(1)
    if show : S.plot()
    debug(titre='Fin testParticulier : (%s)'%name)

def testConstructeurs(filename, show):
#     show=True
    name = filename.name
    debug(titre='testConstructeurs : (%s)'%name)
    debug(paragraphe='Constructeur vide')
    S0 = NSplineComposee()
    print S0
    pprint(S0.toDump())
    methodes = [('cubic','not-a-knot'),('ius',3), ('cubic','natural'), ('ius',1)]
    points=pointsFrom(filename)
    _, nba = computeCordeAndNBA(points=points)
    S = NSplineComposee(points=points,
                        ruptures=[0,nba/2,nba, nba+nba/2,-1],
                        methode=methodes,
                        precision=4*[500],
                        mode=4*['courbure'],
                        nbpe=[50,40,40,30],
                        name=name
                         )
    debug(S)
    if show : S.plot()

def testMethodesGlobales(filename,show):
    name = filename.name
    debug(titre='testMethodesGlobales : (%s)'%name)
    S = NSplineComposee()
    S.open(filename)
    print S
    debug(paragraphe='echantillonner() (%s)'%name)
    pts = S.echantillonner()
    debug(echantillon=pts.tolist())
    debug(epoints=S.epoints.tolist())

    debug(paragraphe='hardScale((0.2,0.2), centre=np.asarray((0,0))) (%s)'%name)
    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S.hardScale((0.2,0.2),centre=(0,0))')

    debug(paragraphe='translate((-1.5,0.05)) : (%s)'%name)
    S.translate((-1.5,0.05))
    debug(S)
    if show : S.plot(titre='translate((-1.5,0.05)) : (%s)'%name)

    debug(paragraphe='hardRotate(30) centre=(0,0)(%s)'%name)
    debug(S)
    S.hardRotate(30, centre=(0,0))
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='hardRotate(30) centre=(0,0)(%s)'%name)

    debug(paragraphe='symetriser(1) (%s)'%name)
    debug(S)
    S.symetriser(1)
    if show : S.plot(titre='symetriser(1) (%s)'%name)


def testNSplineComposee(filename, show):
    name = filename.name
    debug(titre='testNSplineComposee : (%s)'%name)
    debug(paragraphe='Constructeur vide')
    S0 = NSplineComposee()
    print S0
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
        _, nba = computeCordeAndNBA(points=points)
        S = NSplineComposee(points=points,
                             ruptures=[0,nba,-1],
                             methode=[methode[0],methode[1]],
                             precision=[1000,1000],
                             mode=['courbure','courbure'],
                             nbpe=[50,40],
                             name='%s-%d'%(name,k)
                             )
        debug(S)

        debug(paragraphe='elaguer(eps=1, replace=True) : (%s)'%name)
        S.elaguer(eps=1, replace=True)
        debug(S)
        debug(cpoints=S.cpoints.tolist())
        debug(epoints=S.epoints.tolist())
        if show : S.plot(titre=u'élaguage')

    debug(paragraphe='echantillonner() (%s)'%name)
    pts = S.echantillonner()
    debug(echantillon=pts.tolist())
    debug(epoints=S.epoints.tolist())

    debug(paragraphe='hardScale((0.2,0.2), centre=np.asarray((0,0))) (%s)'%name)
    S.hardScale((0.2,0.2), centre=np.asarray((0,0)))
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S.hardScale((0.2,0.2),centre=(0,0))')

    debug(paragraphe='translate((-1.5,0.05)) : (%s)'%name)
    S.translate((-1.5,0.05))
    debug(S)
    if show : S.plot(titre='translate((-1.5,0.05)) : (%s)'%name)

    debug(paragraphe='hardRotate(30) centre=(0,0)(%s)'%name)
    debug(S)
    S.hardRotate(30, centre=(0,0))
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='hardRotate(30) centre=(0,0)(%s)'%name)

    debug(paragraphe='symetriser(1) (%s)'%name)
    debug(S)
    S.symetriser(1)
    if show : S.plot(titre='symetriser(1) (%s)'%name)

    titre='reconstruction (%s)'%name
    debug(paragraphe=titre)
    points = pointsFrom(filename)
    corde, nba = computeCordeAndNBA(points)
    S = NSplineComposee(points=points,
                        ruptures=[0,nba,-1],#reference
                        methode = [['cubic',((2, 0, 0), (1, 0, -10))], #extrados
                                   ['cubic',((1, 0, -10), (2, 0, 0))]],#intrados
                        precision=[1000,1000],
                        mode=['courbure','courbure'],
                        nbpe=[20,40],
                        name=name
                        )
    debug(S)
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='Original '+str(S.methode))

    debug(paragraphe='S.toDump() (%s)'%name)
    dump0 = S.toDump()
    pprint(dump0)

    debug(paragraphe='NSplineComposee(**dump0) (%s)'%name)
    S = NSplineComposee(**dump0)
    debug(join_1=S)
    debug(Rupture=S.cpoints[S.ruptures].tolist())
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S = NSplineComposee(**S.toDump())')

    debug(paragraphe='dump1 = S.toDump() (%s)'%name)
    dump1 = S.toDump()
    debug(dump0_not_equal_dump1=dictsAreNotEqual(dump0,dump1))

    debug(paragraphe='S.load(dump1) (%s)'%name)
    S.load(S.toDump())
    R = S.cpoints[S.ruptures]
    debug(S)
    debug(R=R)
    dump2 = S.toDump()
    debug(dump1_not_equal_dump2=dictsAreNotEqual(dump1,dump2))
    if show : S.plot(titre='S.load(**S.toDump())')

    debug(paragraphe='S.join(1) (%s)'%name)
#     debug(paragraphe='S.join(1) (%s) ici je perd le nom de la spline'%name)
    S.join(1)
    debug(S)
    debug(Ruptures=S.cpoints[S.ruptures])
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='join(1)'+str(S.methode))

    debug(paragraphe='S.split(len(S)/2) (%s)'%name)
    S.split(len(S)/2)
    debug(S)
    debug(ruptures=S.cpoints[S.ruptures])
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S.split(len(S)/2)'+str(S.methode))

    debug(paragraphe='removePoint(rupture[1]) (%s)'%name)
    npt = S.ruptures[1]
    try : S.removePoint(npt)
    except NotImplementedError as msg :
        debug('normal', msg)
#         raise
    debug(S)

    pt = (0.33, 0.1)
    titre='S.insertPoint(%s) %s'%(str(pt),name)
    debug(paragraphe=titre)
    try :
        S.insertPoint(pt)
        if show : S.plot(titre=titre)
    except ValueError as msg :
        print msg
        raise
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())

    i = len(S)/2
    titre = 'S.insertPoint(pt,%d) pt=%s, %s'%(i,str(pt),name)
    debug(paragraphe=titre)
    S.insertPoint(pt,i)
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre=titre)

    pt = (1.5,0.2)
    titre = 'S[0] = %s (%s)'%(str(pt),name)
    debug(paragraphe=titre)
    S[0] = (1.5,0.2)
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre=titre)

    titre = 'S[len(S)-1] = S[0] (%s)'%name
    debug(paragraphe=titre)
    S[len(S)-1] = S[0]
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S[len(S)-1] = S[0]')
    for pt in((1.5,0.1),(2,2),(6,0)) :
        titre = 'S.appendPoint(%s) (%s)'%(pt,name)
        debug(paragraphe=titre)
        S.appendPoint(pt)
        debug(S)
        debug(cpoints=S.cpoints.tolist())
        debug(epoints=S.epoints.tolist())
        if show : S.plot(titre=titre)

    debug(paragraphe='S.removePoint(3) (%s)'%name)
    S.removePoint(3)
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    debug(epoints=S.epoints.tolist())
    if show : S.plot(titre='S.removePoint(3)')

    debug(paragraphe='S._update() (%s)'%name)
    S._update()
    debug(S)
    debug(cpoints=S.cpoints.tolist())
    k = len(S.splines)
    S.nbpd = k*[100]
    debug(dpoints=S.dpoints.tolist())
    debug(epoints=S.epoints.tolist())

    if show : S.plot(titre='S._update()')
#     rdebug('**********************\ndouteux, a verifier\n**********************')
    debug(titre='Fin testNSplineComposee : %s'%filename.name)

def testMain(show=False):
    files = []
    files.append(Path(VALIDATION_DIR,'points-86pts.gnu'))
    files.append(Path(VALIDATION_DIR,'diamirE-profilnormalise-86pts.spl'))
    files.append(Path(VALIDATION_DIR,'diamirE-profilnormalise-24pts.spl'))
    files.append(Path(VALIDATION_DIR,'shark-profilnormalise-26pts.spl'))
    files.append(Path(VALIDATION_DIR,'shark-profilnormalise-86pts.spl'))
    files.append(Path(VALIDATION_DIR,'NACA2415-100pts.spl'))
    files.append(Path(VALIDATION_DIR,'splinesimple-86pts.pkl'))
    files.append(Path(VALIDATION_DIR,'splinesimple-86pts.spl'))
    files.append(Path(VALIDATION_DIR,'splinesimple-21pts.spl'))

    if 1 : testBlocJonc(Path(VALIDATION_DIR,'blocjonc-splinesimple.spl'),show=show)
    if 1 : testParticulier(files[6], show=show)
    for filename in files[:] :
        if 1:testNSplineComposee(filename, show=show)
    for filename in files[:] :
        if 1 : testConstructeurs(filename, show=show)
    debug(titre='Fin testMain : %s'%filename.name)
if __name__=="__main__":
    testMain()
