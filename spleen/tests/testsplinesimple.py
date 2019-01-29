#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016-2017-2018 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
'''
from matplotlib import pyplot as plt
from difflib import context_diff
from collections import OrderedDict
from scipy.interpolate import dfitpack
__updated__="2019-01-29"
import os, sys
from path import Path
# print os.getcwd()
from splinesimple import (NSplineSimple, placementReperesMontage,
                          elaguer, correctionReperesMontage)
from utilitaires.utilitaires import (diff,XY)
from splineabstraite import (distance2PointSpline,)
from utilitaires.lecteurs import pointsFrom
import sys,os,math
# from array import array

import numpy as np
from numpy import log, linspace, asarray, sqrt, zeros
from numpy.linalg import  norm
# import scipy as sp
from config import VALIDATION_DIR, RUNS_DIR
# from scipy.optimize import newton, minimize
from pprint import pprint
from utilitaires.utilitaires import (debug, debug,absCurv,dist2)
import cPickle
def mesures(S0) :
    def fmt(r):return int(1000*r)/1000.0
    x, y = S0.barycentre[0]
    L = S0.longueur
    l = S0.largeur
    h = S0.hauteur
    c = S0.integraleCourbure()
    last = -2 if S0.isClosed() else -1
    s = norm(S0[0]-S0[last])/L
    d = {'barycentre':(fmt(x), fmt(y)),
         'largeur':fmt(l),
         'hauteur':fmt(h),
         'longueur':fmt(L),
         'courbure':fmt(c),
         'sinuosite':fmt(s)}
    return OrderedDict(sorted(d.items(), key=lambda t: t[0]))

def testSequenceBasique(filename, show=True):
    """Une séquence basique d'utilisation de spline"""
    debug(titre='testSequenceBasique %s'%filename.name)
    debug(paragraphe='0. NSplineSimple() ')
    S = NSplineSimple()
    debug(S=S)

    c, d, knots = S.cpoints, S.dpoints, S.knots
    print 'cpoints =', c, type(c), c.shape
    print 'dpoints =', d, type(d), d.shape
    print '  knots =', knots, type(knots), knots.shape
    debug(titre='differentes methodes %s'%filename.name)
    for k, methode in enumerate((
                                ('us',1),#'linear'),
#                                 ('us',2),#'linear'),
                                ('us',3),#'courbure'),
#                                 ('us',4),#'courbure'),
#                                 ('us',7),#'courbure'),
                                ('ius',1),#'courbure'),
                                ('ius',3),#'courbure'),
#                                 ('lsqus',1),#'linear'),
#                                 ('lsqus',2),#'linear'),
#                                 ('lsqus',3),#'courbure'),
#                                 ('lsqus',4),#'courbure'),
                                ('cubic','not-a-knot'),
                                ('cubic','periodic'),
                                ('cubic','natural'),
                                ('cubic',((1, 0, -5), (1, 0, -5))),
                                ('cubic',((2, 0, 0), (1, 0, -5))), #extrados
                                ('cubic',((1, 0, -5), (2, 0, 0))), #intrados
#                                 ('ius',7),#'courbure'),
                                )):
        debug(paragraphe='1. methode-%d=%s (fichier %s)'%(k, str(methode),filename.name))
        S = NSplineSimple(cpoints=pointsFrom(filename), methode=methode, name=str(methode))
        debug(S)
        try : debug(S._dpoints.shape)
        except AttributeError as msg:debug(u'    Normal pas de S._dpoints : %s'%str(msg))
        except Exception as msg : debug(str(msg))
        try : debug(u'    Normal, S.dpoint est calcule',lendpoints=len(S.dpoints))
        except Exception as msg : debug(str(msg))
        debug(paragraphe='2. S._update() methode-%d=%s  (fichier %s)'%(k, str(methode),filename.name))
        S._update()
        try :debug(S._dpoints)
        except AttributeError as msg:debug(u'    Normal pas de S._dpoints : %s'%str(msg))
        try : debug(u'    Normal, S.dpoint est calcule',lendpoints=len(S.dpoints))
        except Exception as msg : debug(str(msg))
        debug(S)
        if show :
            try : S.plot(plt, show=True)
            except Exception as msg :
                debug(str(msg))
                plt.plot(title=str(msg))
                plt.title(str(msg)+'\n'+str(methode))
                plt.show()
        debug('    Longueurs (c,d)=(%.5g,%.5g) ; nspline=%d'%(S.longueur('c'),S.longueur('d'), S.nbspline))
#         debug('Longueurs (c,d,e)=(%.5g,%.5g,%.5g)'%(S.longueur('c'),S.longueur('d'),S.longueur('e')))

        i, p = len(S)/2, S.centregravite
        debug(paragraphe='3. S[%d]=[%.3g, %.3g], methode-%d=%s  (fichier %s)'%(i, p[0], p[1], k, str(methode),filename.name))
        oldp = S[i].copy()
        S[i] = p
        debug(S[i],oldp=oldp)
        debug('    Longueurs (c,d)=(%.5g,%.5g) ; nspline=%d'%(S.longueur('c'),S.longueur('d'), S.nbspline))
        debug('    S(1)=%s'%S(1.0))#ne donne pas tout a fait cpoints[-1]
        ac = absCurv(S._cpoints, normalise=True)
        print '    abscurv=',ac.shape, 'knots=',S.knots.shape, 'cpoints=',S._cpoints.shape
        print '    norm(S(ac)-S._cpoints)=',norm(S(ac)-S._cpoints)
        print '    norm(ac-S.knots)',norm(ac-S.knots)
        if show : S.plot(plt,show=True, titre='%s : modif S[%d]=[%.3g, %.3g]'%(filename.name, i, p[0], p[1]))
        debug(S[i],oldp=oldp)
        S[i] = oldp
        debug(paragraphe='4. echantillonnage methode-%d=%s  (fichier %s)'%(k, str(methode),filename.name))
        S.nbpe = 120
        print 'mode=',S.mode, ', epoints=',S.epoints.shape
        e0, mode0 = S.epoints, S.mode
        print 'mode=',S.mode, ', epoints=',e0.shape
        print 'T = %s'%S.tech.tolist()
        S.mode = 'cos'#efface tech et epoints
        e1, mode1 = S.epoints, S.mode
        print 'mode=',S.mode, ', epoints=',e0.shape
        print 'T = %s'%S.tech.tolist()
        #on regarde si ca se recalcule bien apres suppression epoints et tech
        del S._epoints, S._tech
        e2, mode2 = S.epoints, S.mode
        print 'mode=',S.mode, ', epoints=',e2.shape
        print 'T = %s'%S.tech.tolist()
        print 'epoints norme(e1-e2)=%.3g'%norm(e1-e2)
        S.mode = 'x3'#efface tech et epoints
        e2, mode2 = S.epoints, S.mode
        print 'mode=',S.mode, ', epoints=',e2.shape
        print 'T = %s'%S.tech.tolist()
        S.mode = 'courbure'#efface tech et epoints
        mode3 = S.mode
        try :
            e3, mode3 = S.epoints, S.mode
            print 'mode=',S.mode, ', epoints=',e3.shape
            print 'T = %s'%S.tech.tolist()
            X3, Y3 = XY(e3)
        except ValueError as msg :
            print 'methode %s, mode echantillonnage = %s : %s'%(str(S.methode), str(S.mode), str(msg))
            X3,Y3 = zeros((0,)), zeros((0,))
        X0, Y0 = XY(e0)
        X1, Y1 = XY(e1)
        X2, Y2 = XY(e2)
#         if show : S.plot(plt,show=True)
        if show :
#             plt.plot(X0,Y0,'r.-',label='%s : mode=%s'%(S.name,mode0))
            plt.plot(X1,Y1,'b.-',label='%s : mode=%s'%(S.name,mode1))
#             plt.plot(X2,Y2,'k.-',label='%s : mode=%s'%(S.name,mode2))
            plt.plot(X3,Y3,'m.-',label='%s : mode=%s'%(S.name,mode3))
            plt.legend()
            plt.axis('equal')
            plt.show()






def testConstructeurs(filename, show=False) :
    debug(titre='testConstructeurs %s'%filename.name)
    from matplotlib import pyplot as plt
    p = NSplineSimple()
    debug(paragraphe='1. constructeur vide')
    debug(show=show)
    debug(p=p)
    p = NSplineSimple(cpoints=np.zeros((0,2)),
                      parent=None,
                      methode=('cubic','periodic'),
                      mode='rayon',
                      name='SPLINE1')
    debug(paragraphe='2. constructeur presque vide')
    print p
    # for point in p.qcpolygon :
    #     print point
    # for point in p.qdpolygon :
    #     print point
    debug(paragraphe='3. p vide => p.cpoints=points de %s'%filename.name)
    points = pointsFrom(filename)
    p.cpoints = points
    if show : p.plot(plt, titre='cpoints=pointsFrom()')

    S0 = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic','periodic'),
#                        mode='courbure',
                       name=filename.name)
    debug(S0)
    if show :S0.plot(plt, titre=filename.name)

    debug(paragraphe='4. S = NSplineSimple(**dump00) (points:%s)'%filename.name)
    dump00 = S0.toDump()
    debug("dump00")
    pprint(dump00)
    S = NSplineSimple(**dump00)
    debug("dump00==S.toDump() : %s"%dump00==S.toDump())

    debug(paragraphe='4. S1 = S0.copy() (points:%s)'%filename.name)
    S = S0.copy()
    debug("dump00==S1.toDump() : %s"%dump00==S.toDump())







    filename = Path(RUNS_DIR, 'spline.npkl')
    S0.save(filename=filename)
    S0.open(filename)
#     debug('S0.toDump()')
#     pprint(S0.toDump())
    debug(S0)
    debug('S0.toDump()==dump00',S0.toDump()==dump00)
    if show : S0.plot(plt, titre='open pkl')
    if show : plt.show()
    S1 = S0.copy()
    d0,d1 = dump00, S1.toDump()
    print sorted(d0.keys())
    print sorted(d1.keys())
    for key in sorted(set(d0.keys()).union(set(d1.keys()))):
        v0, v1 = d0[key], d1[key]
        print '%10s : %5s : %20s : %-20s :'%(key,v0==v1,v0,v1)
    debug('S1.toDump()==dump00',d0==d1)
    if not d0==d1 :
        s0,s1 = str(d0).split(),str(d1).split()
        while s0 and s1 :
            w0,w1 = s0.pop(),s1.pop()
            print w0==w1,w0,w1
#             if not w0==w1:break
#         print w0, w1
    if show :S1.plot(plt, titre='copy')
#     exit()

def testMethodes(filename, show=True):
    debug(titre='testMethodes(%s)'%filename.name)
    from matplotlib import pyplot as plt
    k = 0
    for methode in (
                    ('cubic','not-a-knot'),
                    ('cubic','periodic'),
                    ('cubic','natural'),
                    ('cubic',((1, 0, -5), (1, 0, -5))),
                    ('cubic',((2, 0, 0), (1, 0, -5))), #extrados
                    ('cubic',((1, 0, -5), (2, 0, 0))), #intrados
                    ):
        k +=1
        S = NSplineSimple(points=pointsFrom(filename),
                          methode=methode,
                          mode='courbure',
                          name='SPLINE%d'%(k+1))
        debug(S)
        if show : S.plot(plt, titre=filename.name+str(methode))
        if show : S.plotCourbure()

#     filename=Path(VALIDATION_DIR,'simple','simple2d2.gnu')
    for methode, mode in (
#                     (('ius',1),'linear'),
#                     (('ius',2),'linear'),
#                     (('ius',3),'courbure'),
                    (('ius',4),'courbure'),
                    (('ius',7),'courbure'),
                    (('us',1),'linear'),
                    (('us',2),'linear'),
                    (('us',3),'courbure'),
                    (('us',4),'courbure'),
                    (('us',7),'courbure'),
                    ):
        try :
            S = NSplineSimple(points=pointsFrom(filename), methode=methode, mode=mode)
            debug(S)
#             debug(u'\n'+70*u'*'+u'\n'+ u' TODO : methode=("us",k) ne marche pas. Va y avoir plein d\'exceptions'+u'\n'+70*u'*'+u'\n')

            if show : S.plot(plt, titre=filename.name+str(methode))
        except Exception as msg:
            debug('(methode=%s; mode=%s) => pas de spline'%(str(methode),str(mode)),
                   '%s : %s'%(type(msg).__name__,unicode (msg)))
    if show : plt.show()
#     exit()
def testModifGlobales(filename, show=True) :
    debug(titre="testModifGlobales : %s"%filename.name)
    from matplotlib import pyplot as plt
#     numfig = 0
#     plt.figure(numfig)
    S0 = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                       precision=1000,
                       mode='courbure',
                       name='SPLINE01')
    if show : S0.plot(plt, titre='NSplineSimple(points=pointsFrom(filename),...')
#     numfig += 1
#     plt.figure(numfig)
    S0.cpoints = pointsFrom(filename)
    debug(initial=mesures(S0))
    if show : S0.plot(plt, titre=filename.name+'cpoints = pointsFrom()')

#     numfig += 1
#     plt.figure(numfig)
    S0.symetriser(0)
    debug(symetriser_0=mesures(S0))
    if show : S0.plot(plt, titre='symetriser(0)')

#     numfig += 1
#     plt.figure(numfig)
    S0.symetriser(1)
    debug(symetriser_1=mesures(S0))
    if show : S0.plot(plt, titre='symetriser(1)')

#     numfig += 1
#     plt.figure(numfig)
    S0.hardRotate(30)#en degres par defaut
    debug(hardRotate_30=mesures(S0))
    if show : S0.plot(plt, titre='hardRotate(30)')

#     numfig += 1
#     plt.figure(numfig)
    S0.hardScale((2,0.5))
    debug(hardScale_2_0point5=mesures(S0))
    if show : S0.plot(plt, titre='hardScale((2,0.5))')

#     numfig += 1
#     plt.figure(numfig)
    S0.translate((2,3))
    debug(translate_2_3=mesures(S0))
    if show : S0.plot(plt, titre='translate((2,3))')
#     if show : plt.show()

def testModifLocales(filename, show=True)  :
    debug(titre="testModifLocales : %s"%filename.name)
#     exit()
    from matplotlib import pyplot as plt
    S = NSplineSimple(points=pointsFrom(filename),
                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                       mode='courbure',
                       name='SPLINE01')
    mesure = mesures(S)
    keys = mesure.keys()
    for key, value in mesure.items() :
        exec('%s=%s'%(key,value))
    for key in keys :
        exec('print %s'%key)
    if show : S.plot(plt, titre=filename.name)

    try : S.insertPoint((0,0))
    except ValueError as msg : debug(msg)

    S.appendPoint((2,1.5))
    if show : S.plot(plt, titre='appendPoint((2,1.5))')

    S.insertPoint((2,-2))
    if show : S.plot(plt, titre='insertPoint((2,-2))')

    try :
        S.insertPoint((6,0))
        if show : S.plot(plt, titre='insertPoint((6,0))')
    except ValueError as msg :
        debug(msg)
    S.removePoint(3)
    if show : S.plot(plt, titre='removePoint(3)')

    S[1] = (4,-0.1)
    if show : S.plot(plt, titre='S[1] = (4,-0.1)')
    print S[1]
#     if show : S.plot(plt, )
    try : S._update()
    except AttributeError as msg : debug(msg)
    pprint(S.__dict__)
    dump = S.toDump()
    debug()
    pprint(dump)
    dump['methode'] = ('cubic','natural')
    S.load(dump)
    S.name='SPLINE3'
    if show : S.plot(plt, titre='load(toDump())')
    print S
    if show : plt.show()

def testSaveRead(filename, show=True):
    u"""Diverses manières de lire, construire, sauvegarder une NSplineSimple"""
    # from matplotlib import pyplot as plt
    debug(titre="testSaveRead : %s"%filename.name)
    print "    ==> NSplineSimple::__init__()"
    S = NSplineSimple(points=pointsFrom(filename),
                      methode=('cubic',((2, 0, 0), (1, 0, -5))),
                      mode='courbure',
                      name='SPLINE01')
    S4 = NSplineSimple()
    S4.load(S.toDump())
    dump4 = S4.toDump()
    fname = Path(filename.dirname, filename.namebase+'#.pkl')
    print fname
    print "    ==> cPickle::dump()/load()"
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = NSplineSimple(**dump)
#     S1.plot(plt)
#     if show : plt.show()
    dump  = S.toDump()
    dump1 = S1.toDump()
    print 'dump==dump1 ?', dump==dump1

    print "    ==> file::read()/write() "
    fname = Path(filename.dirname, filename.namebase+'#.spl')
    S2 = NSplineSimple(**dump)
    dump2 = S2.toDump()
    with open(fname,'w') as f :
        f.write(str(S2.toDump()))
    with open(fname,'r') as f :
        dump = eval(f.read())
    print "    ==> NSplineSimple::load()"
    S3 = NSplineSimple(**dump)
    dump3 = S3.toDump()
#     dump3_1 = S3.toDump()
    print 'dump==dump1=dump2=dump3=dump4 ?', dump==dump1==dump2==dump3==dump4
    debug(dump=dump)
#     debug(dump1=dump1)

#     debug(dump2=dump2)
#     debug(dump3=dump3)
#     debug(dump4=dump4)
#     exit()

def testEchantillonnage(filename, trados, show=True):
    u"""echantillonnage entre ta et tb. [ta,tb]=[0,1] => echantillonnage complet"""
#     s = NSplineSimple(points=pointsFrom(filename))
    u"""Extraction intrados et extrados, normalisation"""
    debug(titre="testEchantillonnage : %s"%filename.name)
    from matplotlib import pyplot as plt

    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points*=(1000/corde) #en mm
    corde = 1000
    if trados == 'e' : #Extrados
        methode = ('cubic',((2, 0, 0), (1, 0, -corde/2)))#extrados
        c0 = points[:1+nba]
    elif trados == 'i': #Intrados
        methode = ('cubic',((1, 0, -corde/4), (2, 0, 0)))#intrados
        c0 = points[nba:]
    u"""Construction et élagage de la spline d'interpolation du trados"""
    precision = 1000
#     T0 = linspace(0,1,precision)
    s0 = NSplineSimple(cpoints=c0, precision=precision, methode=methode, mode='courbure')
#     s0.elaguer(1, replace=True)
    ne, ta, tb = 30, 0.2, 0.8
    E = s0.echantillonner(nbp=ne, ta=ta, tb=tb, mode='courbure')
    if show :
        s0.plot(plt, more=[(E[:,0],E[:,1],'y^','echantillon \nn=%d, %.2f=>%.2f'%(ne,ta,tb)),])
#     plt.plot(E[:,0],E[:,1],'y^','echantillon')
        plt.show()

def testDivers(filename, show=True):
    debug(titre="testDivers : %s"%filename.name)

    msg=u"""Dans cet exemple, on a des points doubles p4=p5=p6 et p8=p9,
    donc la spline x(t) ne peut pas être construite.
    Que faire ?
    - supprimer p[6] <=== Oui c'est la solution adoptée ici
    - déplacer légèrement le point
    - faire une spline composée, splittée en 6
    """
    debug(msg)
    from matplotlib import pyplot as plt
#     cpoints = [
#                [ 1.00827138,  0.08545793],
#                [ 1.0063704 ,  0.06720821],
#                [ 1.00541204,  0.04849635],
#                [ 1.00541204,  0.02962959],
#                [ 1.0063704 ,  0.01091773],
#                [ 1.0063704 ,  0.01091773],
#                [ 1.0063704 ,  0.01091773],
#                [ 1.00827138, -0.00733198],
#                [ 1.12226291,  0.03906297],
#                [ 1.12226291,  0.03906297],
#                [ 1.12178274,  0.02023498],
#                [ 1.00827138,  0.08545793],
#                [ 1.00827138,  0.08545793],
#                ]
    cpoints = pointsFrom(filename)
    S0 = NSplineSimple(cpoints=cpoints, name=filename.name)
    debug(S0)
    if show : S0.plot(plt, show=True)
    ########################################
    msg=u"""S est une NSplineSimple. Teste les méthodes :
    - S.knots() qui retourne les knots de la spline.
        Normalement on doit avoir knot==abscurv()==sx.x==sy.x sauf à l'initialisation
    - S.abscurv(T) qui retourne les abscisses curvilignes VRAIES des points S(T)
    - S.abscurv(None) qui retourne les abscisses curvilignes NORMALISÉES des points S(T)
    - S.longueur par intégration numérique (voir la méthode S.abscurv(T))
    - S.distance2To(p) et affiche le résultat pour 3 points
        Calcule la distance du point 'p' à la spline 'S' et le point 'h' de la spline
        le plus proche de 'p'. Retourne également la precision et le nb d'appels à fonction
        Utilise :
        "https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html
    """
    dump = {'classename': 'NSplineSimple',
            'cpoints': cpoints,
            'methode': ('cubic', 'not-a-knot'),
#             'mode': 'segment',
#             'methode': ('ius', 1),
            'mode': 'rayon',
#             'mode': 'cpoints',
            'name': filename.name,
            'nbpe': 30,
            'precision': 1000,
            'role': 'NSplineSimple'}
    dump = {'cpoints': [[0.5279636979103088, 0.07332829385995865], [0.7137287259101868, 0.3275330364704132], [1.268811468596966, 0.3365727279966774], [1.1390328407287598, 0.07332829385995865], [1.2571670716747245, 0.2148051133421408], [1.2206038430660453, -0.0507648238639925], [1.5545598268508911, -0.0048885527066886425]],
#             'methode': ('cubic', 'not-a-knot'),
            'methode': ('ius', 3),
            'name': u'test',
            'nbpe': 30,
            'precision': 1000,
            'mode': 'lineaire'}
    debug(msg)
    S = NSplineSimple(**dump)
    if show : S.plot(plt)
    knots = S.knots
    acn = S.abscurv(knots)#T Normalisé == knots
    acr = S.abscurv(acn)
    debug('S.abscurv(None)',acn.tolist())
    debug('S.knots',knots.tolist())
    l = S.longueur
    debug('(acr-acn)*longueur',(acr-acn*l).tolist())
    debug('longueur=%.10g'%l)
    lr = acr[-1]
    debug('vraie abscurv[-1]=%.10g'%lr)
    ln = acn[-1]
    debug('abscurv[-1]=%.10g'%ln)
    S.precision=1000
    ld3 = S.dlongueur
    debug('S.dlongueur=%g'%ld3+u"(calculé sur %d points)"%S.precision)
    S.precision=10000
    ld4 = S.dlongueur
    debug('S.dlongueur=%g'%ld4+u"(calculé sur %d points)"%S.precision)
    S.precision=100000
    ld5 = S.dlongueur
    debug('S.dlongueur=%g'%ld5+u"(calculé sur %d points)"%S.precision)
    le = S.elongueur
    debug('S.elongueur=%g'%le)
    lc = S.clongueur
    debug('S.clongueur=%g'%lc)
    diff = [l-lr,l-ln,l-ld3,l-ld4,l-ld5,l-le,l-lc]
    plt.bar(range(len(diff)),diff)
    plt.xticks(range(len(diff)),'lr,ln,ld3,ld4,ld5,le,lc'.split(','))
    plt.title('differentes longueurs calculees')

    if show : plt.show()
    p1, p2, p3 = [0.2,0.01], [1.1,0.2], [-0.1,0.3]
    res1, res2, res3 = S.distance2To(p1), S.distance2To(p2), S.distance2To(p3)
    pj1, pj2, pj3 = S(res1.x), S(res2.x), S(res3.x)
#     debug(str(res))
    debug("****  DISTANCE POINT-SPLINE  ****")
    debug('message="%s"'%res1.message,'d(S,p1)=%g'%res1.fun, 't=%g'%res1.x, 'nb appels a fct=%d'%res1.nfev)
    debug('message="%s"'%res2.message,'d(S,p2)=%g'%res2.fun, 't=%g'%res2.x, 'nb appels a fct=%d'%res2.nfev)
    debug('message="%s"'%res3.message,'d(S,p3)=%g'%res3.fun, 't=%g'%res3.x, 'nb appels a fct=%d'%res3.nfev)
    more = [([p1[0],pj1[0]],[p1[1],pj1[1]],'g-o', 'd(S,p1) : %.2g'%sqrt(res1.fun)),
            ([p2[0],pj2[0]],[p2[1],pj2[1]],'g-o', 'd(S,p2) : %.2g'%sqrt(res2.fun)),
            ([p3[0],pj3[0]],[p3[1],pj3[1]],'g-o', 'd(S,p3) : %.2g'%sqrt(res3.fun))]
    projections = [(pj1[0],pj1[1],'y.',''), (pj2[0],pj2[1],'y.', ''), (pj3[0],pj3[1],'y.', '')]
    more += projections
    if show : S.plot(plt, more=more, titre='distance point-spline', show=True)
#     plt.plot(p1[0],p1[1],'go', label='p1 : %g'%res1.fun)
#     plt.legend()
    if show : plt.show()
#     return
    debug('ORIGINAL', S=S)
    debug(rect=S.boundingRect())
#     debug('CLONE', clone=S.clone())
    debug('COPY', clone=S.copy())
    S.hardScale((0.5,0.5), S.centre)
    debug('APRES HARDSCALE',S=S)
    debug(rect=S.boundingRect())
    S.hardMove((10,10))
    debug('APRES HARDMOVE',S=S)
    debug(rect=S.boundingRect())
#     debug(top='%g'%rect.top(), left='%g'%rect.left(), width='%g'%rect.width(), height='%g'%rect.height(), )
#     debug(bottom='%g'%rect.bottom(), top='%g'%rect.top())
    if show : S.plot(plt, titre="**dump")
#     S.elaguer(1, replace=True)
#     if show : S.plot(plt, titre="elaguer")

def testPlacementRM(show=True):
    ##############################################
    ### debut construction exemple
    ##############################################
    debug()
    B1 = asarray([[0,0],[0.5,0.25],[0.5,0.8],[0,1.1],[0,1.7],[0.5,2]], float) #un S
    B2 = asarray([[2,0],[2.5,1],[2,2]], float)#un C à l'envers
#     delta = 0.0186#distance entre reperes de montage (en m)
    from matplotlib import pyplot as plt
    s1 = NSplineSimple(cpoints=B1, parent=None,
#                        methode=('ius',1),
                        methode=('cubic','natural'),
                       name='B1',
                       )
    s2 = NSplineSimple(cpoints=B2, parent=None,
#                        methode=('ius',1),
                        methode=('cubic','natural'),
                       name='B2',
                       )

#     ac1 = absCurv(s1.cpoints, normalise=False)
#     ac2 = absCurv(s2.cpoints, normalise=False)
    l1, l2 = s1.longueur, s2.longueur
    debug(u'avant hardScale',l1=l1, l2=l2)
    #Je les ramène à la même longueur
    s2.hardScale((l1/l2, l1/l2))
    l1, l2 = s1.longueur, s2.longueur
    lref = min(l1,l2)
    debug(u'après hardScale',l1=l1, l2=l2)
    B1, B2 = s1.epoints, s2.epoints#le tableaux echantillonnés que va me passer ML
    ##############################################
    ### fin construction exemple
    ##############################################
    delta = 0.2#distance entre reperes de montage (en m)
    s1, s2, T = placementReperesMontage(B1,B2,delta,)
    n=len(T)
    rm1, rm2 = s1(T), s2(T)
    ac1, ac2 = s1.abscurv(T), s2.abscurv(T)#les abs curv. vraies des reperes
    delta1 = ac1[1:]-ac1[:-1]
    delta2 = ac2[1:]-ac2[:-1]
    dd = np.max(np.abs(delta1-delta2))
    debug("Erreur : %.1e mm"%(1000*dd/lref))

    plt.plot(s1.dpoints[:,0], s1.dpoints[:,1], 'b-')
    plt.plot(s1.cpoints[:,0], s1.cpoints[:,1], 'bx', label=u'Points de contrôle des splines')
    plt.plot(s2.dpoints[:,0]-1, s2.dpoints[:,1], 'b-')
    plt.plot(s2.cpoints[:,0]-1, s2.cpoints[:,1], 'bx')
    plt.plot(rm1[:,0], rm1[:,1], 'ro', label=u'Repères de montage')
    plt.plot(rm2[:,0]-1, rm2[:,1], 'ro')
    plt.axis('equal')
    plt.legend()
    msg = u"""Deux bords de pièces à assembler, de longueur $l_1$=%.4g m $\simeq$ $l_2$=%.4g m
    Placements de %d repères de montage espacés de %.2g m.
    L'erreur est de %.2f mm soit %.2f‰"""%(l1,l2,n,delta,1000*dd,1000*dd/lref)
    plt.title(msg)
    if show : plt.show()

def testElaguer(filename, trados, fonction='elaguer', show=True):
    debug(titre="testElaguer : %s"%filename.name)
    from matplotlib.widgets import CheckButtons
    from matplotlib import pyplot as plt
    if fonction=='ajuster' :
        msg=u"""L'élagage par "ajustement" ne marche pas, et il est très long.
        Je le supprime des tests.
        :TODO: voir au besoin si on peut l'améliorer"""
        debug(msg)
        return

    def compareCourbures(s0,s1, corde, trados):
        u"""s0 et s1 sont des intrados (ou des extrados), s0:original, s1:elagué
        On trace les deux splines et les log des abs(courbures), mises à l'échelle pour y voir qque chose"""
#         from matplotlib import pyplot as plt
#         from matplotlib.widgets import CheckButtons
    #         nbpd = self.precision
        s1.nbpe = s0.nbpe
        titre = ' '.join(('courbures', s0.name, s1.name))
        plt.title(titre)
    #     E = s0.epoints
        _, ax = plt.subplots()
        T = linspace(0,1,100)

#         s0.echantillonner(nbpe)
        D0 = s0.dpoints
        C0 = s0.cpoints
        spline0, = ax.plot(D0[:,0], D0[:,1], 'r-', lw=1)
    #         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
        control0, = ax.plot(C0[:,0], C0[:,1], 'ro', lw=1, label=u'%d points contrôle'%len(C0))
        courbure0 = s0.courbure(T)
        courbure1 = s1.courbure(T)
        minc = min(min(courbure0), min(courbure1))
        maxc = max(max(abs(courbure0)), max(abs(courbure1)))
        courbure0 += (1.0 + abs(minc))
        courbure0 = log(courbure0)
        courbure0 /= maxc
        if trados == 'e' :
            courbure0 = courbure0[::-1]
        courbure0, = ax.plot(corde*T[1:], corde*courbure0[1:])
        D1 = s1.dpoints
        C1 = s1.cpoints
        spline1, = ax.plot(D1[:,0], D1[:,1], 'b-', lw=1)
        control1, = ax.plot(C1[:,0], C1[:,1], 'bo', lw=1, label=u'%d points contrôle'%len(C1))

        courbure1 += (1.0 + abs(minc))
        courbure1 = log(courbure1)
        courbure1 /= maxc
        if trados == 'e' :
            courbure1 = courbure1[::-1]
        courbure1, = ax.plot(corde*T[1:], corde*courbure1[1:])
        ax.legend(loc='upper right')
        buttons = ['control0', 'control1', 'courbure0', 'courbure1', 'spline0', 'spline1']
        values = [True, True, True, True, True, True]
        draws = [control0, control1, courbure0, courbure1, spline0, spline1]
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        rax = plt.axes([0.05, 0.4, 0.1, 0.15])
        check = CheckButtons(rax, buttons, values)

        def func(label):
            if label == 'spline0':
                spline0.set_visible(not spline0.get_visible())
            elif label == 'spline1':
                spline1.set_visible(not spline1.get_visible())
            elif label == 'control0':
                control0.set_visible(not control0.get_visible())
            elif label == 'control1':
                control1.set_visible(not control1.get_visible())
            elif label == 'courbure0':
                courbure0.set_visible(not courbure0.get_visible())
            elif label == 'courbure1':
                courbure1.set_visible(not courbure1.get_visible())
            else :
                draw = draws[buttons.index(label)]
                draw.set_visible(not draw.get_visible())
            plt.draw()
        check.on_clicked(func)
        if show : plt.show()
        return plt

    def comparePointsSpline(P,s):
        u"""comparaison graphique, P famille de points s spline
        On trace les projections de P sur la spline s"""
        projections = []
        for p in P :
            optres = distance2PointSpline(p, s)
            #t=abscisse scurviligne de la pj p' de p sur s, et d=distance d(p,p')
            t, d = optres.x, optres.fun
            projections.append(s(t))
        projections = asarray(projections)
        if show : s.plot(plt, more=([P[:,0], P[:,1],'g.', 'C0'], [projections[:,0], projections[:,1],'y.','projections']))
        if show : plt.show()

    u"""Extraction intrados et extrados, normalisation"""
#     if '.spl' in filename :
#         with open(filename,'r') as f :
#             lines = f.readlines()
#         for line in lines :
#             dump = eval(line)
#             points = dump['cpoints']
#         print S

    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points*=(1000/corde)#en mm
    corde = 1000
    if trados == 'e' :#Extrados
        methode = ('cubic',((2, 0, 0), (1, 0, -corde/2))) #extrados
        c0 = points[:1+nba]
    elif trados == 'i':#intrados
        methode = ('cubic',((1, 0, -corde/4), (2, 0, 0)))#intrados
        c0 = points[nba:]
    u"""Construction de la spline d'interpolation du trados"""
    precision = 1000
    T0 = linspace(0,1,precision)
    s0 = NSplineSimple(cpoints=c0, precision=precision, methode=methode, mode='courbure')
    u"""La spline elaguée"""
    if fonction == 'elaguer' :
        s1, d, (n0,n1) = s0.elaguer(1)
    elif fonction == 'ajuster' :
        s1, d, (n0,n1) = s0.ajuster(1)
    d0 = s0(T0)
    d1 = s1(T0)
    print u'  max dist(C0, s1)      = %.2e max des distances aux points de contrôle (‰) <  ################'%(d)
    print u'  max norme 2     D0-D1 = %.2e max des distances point à point des splines discrétisées'%max(norm(d0-d1, 2, axis=1))
    print u'  norme frobenius D0-D1 = %.2e somme des distances point à point des splines discrétisées'%norm(d0-d1, 'fro')
    print u'  norme frobenius D0-D1 = %.2e moyenne des distances point à point des splines discrétisées'%(norm(d0-d1, 'fro')/len(d0))
    print u'  nb pts contrôle       : %d => %d'%(len(s0),len(s1))
    compareCourbures(s0, s1, corde, trados)
#     comparePointsSpline(s0.cpoints, s1)
    return

def testElagage1(filename, show=True):
    u"""Elagage d'une spline simple=suppression de points de contrôle.
    filename doit contenir un profil.
    1-On détermine son extrados, son intrados, le ba, le bf, la corde...
    2-On construit la spline de l'intrados
    3-On l'élague
    3-On affiche"""
    debug(titre="testModifGlobales : %s"%filename.name)
    points = pointsFrom(filename)
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf,point)
        if d > corde :corde, nba = d, k
    corde = math.sqrt(corde)
    debug (corde=corde, nba=nba)
    points/=corde
    corde = 1.0
#     return
# #     #Extrados
#     methode = ('cubic',((2, 0, 0), (1, 0, -corde/2))) #extrados
#     s0 = NSplineSimple(points=points[:nba+1], precision=100, methode=methode)
    #intrados
    methode = ('cubic',((1, 0, -corde/2), (2, 0, 0)))#intrados
    s0 = NSplineSimple(points=points[nba:], precision=100, methode=methode, mode='courbure')
    print s0
    s1 = elaguer(s0, eps=1.0e-4)
#     s1 = elaguer(s1, eps=1.0e-4)
#     s1 = elaguer(s1, eps=1.0e-4)

    Td = linspace(0,1,1000)
    c0 = s0.cpoints.copy()*1000#en mm
    d0 = asarray(s0(Td))*1000
#     debug(d0_shape=d0.shape)
    c1 = s1.cpoints.copy()*1000
    d1 = asarray(s1(Td))*1000
#     debug(d1_shape=d1.shape)
#     d1 = s1.dpoints.copy()
#     plt = s0.plot()
    from matplotlib.widgets import CheckButtons
    from matplotlib import pyplot as plt
    _, ax = plt.subplots()
    spline0, = ax.plot(d0[:,0], d0[:,1], 'b-', lw=1)
#     echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
    control0, = ax.plot(c0[:,0], c0[:,1], 'b.', lw=1)
    spline1, = ax.plot(d1[:,0], d1[:,1], 'r-', lw=1)
#     echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
    control1, = ax.plot(c1[:,0], c1[:,1], 'r.', lw=1)

    plt.subplots_adjust(left=0.2)
    plt.axis('equal')
    rax = plt.axes([0.05, 0.4, 0.1, 0.15])
    check = CheckButtons(rax, ('control0', 'control1', 'spline0', 'spline1'), (True, True, True, True))

    def func(label):
        if label == 'spline0':
            spline0.set_visible(not spline0.get_visible())
        elif label == 'control0':
            control0.set_visible(not control0.get_visible())
        elif label == 'spline1':
            spline1.set_visible(not spline1.get_visible())
        elif label == 'control1':
            control1.set_visible(not control1.get_visible())
        plt.draw()
    check.on_clicked(func)

    if show : plt.show()
#     for precision in range(1000, 1001, 1000):
#         print 'precision=%d, distance(s0, s1)=%.2e'%(precision, distance(s0, s1, precision))
#     p1 = elaguer2(p0)

def testCorrectionRM(show=True):
    ##############################################
    ### debut construction exemple
    ##############################################
    from matplotlib import pyplot as plt
    debug()
    vc = 0.25
    B = asarray([[0,0],[0.5,0.25],[0.5,0.8],[0,1.1],[0,1.7],[0.5,2],
                 [0.6,1.8],[1.5,1.0]], float) #un S
    # construction valeurs couture
    SB1 = NSplineSimple(cpoints=B, parent=None,
#                         methode=('ius',1),
                        methode=('cubic','natural'),
                        name='B',
                        )
#     SB1.plot(plt)
    knots = SB1.knots#==B
    tangentes = SB1(knots,1)#dérivées premières
#     debug(tangentes=tangentes)
    normales = asarray(zip(-tangentes[:,1],tangentes[:,0]))
    for k, normale in enumerate(normales):
        normales[k] = normale/np.linalg.norm(normale)
#     u = normales[0]
#     debug (u[0]**2+u[1]**2)
#     debug(normales=normales)
    #points couture
    Bc = B + vc*normales
#     debug(Bc=Bc)
    #spline coutures/reperes
    C = NSplineSimple(cpoints=Bc,
                      methode=('cubic','natural'),
                      name='v.coutures')
#     debug(C)
#     C.plot(plt)
    DC = C.dpoints
    DB = SB1.dpoints
    ####################
    ## reperes mal placés
    ####################
    nrep = len(B)#nb reperes
    TR = absCurv(B, normalise=True)
    R = C(TR)#repere sur spline couture, mal placés
    ##############################################
    ### fin construction exemple
    ### les distances dans B ne sont pas les mêmes
    ###    que les distances dans R
    ##############################################
#     B,R = R,B#permutation B et R
#     DB,DC = DC, DB#permutation tracés
    ###################
    lR = absCurv(R, normalise=False)[-1]
    lB = absCurv(B, normalise=False)[-1]
    ##############################################
    #  correction
    ##############################################

    SCR, T = correctionReperesMontage(B, R, 'test')
    CR = SCR(T)

    ##############################################
    #tracé
    ##############################################
    plt.subplot(1,2,1)
    plt.plot(DC[:,0], DC[:,1],'g,')
#     plt.plot(CC[:,0], CC[:,1],'r.', label='spline coutures')
    plt.plot(DB[:,0], DB[:,1],'b,')
#     plt.plot(CB[:,0], CB[:,1],'b.', label='spline bord')
    plt.suptitle(u'Cas ou longueur(R) = %.5g >  %.5g = longueur(B)'%(lR, lB), fontsize=16, fontname='serif')
    plt.plot(B[:,0], B[:,1],'b-o',label=u'Points du bord (%d points)'%nrep)
    plt.plot(R[:,0], R[:,1],'g-o',label=u"Initialement %d repères mal placés"%nrep)
    plt.plot(CR[:,0], CR[:,1],'ro',label=u'Correction : %d repères '%nrep)
    plt.legend()
    plt.axis('equal')
#     if show : plt.show()
    def deltaLongueur(B, R, sCR, T):
        Bac, Rac, CRac = absCurv(B), absCurv(R), sCR.abscurv(T)
        lB, lR, lCR = diff(Bac), diff(Rac), diff(CRac)
        return lB-lR, lB-lCR
    lBR, lBCR = deltaLongueur(B, R, SCR, T)
    lB, lR, lCR = absCurv(B)[-1], absCurv(R)[-1], SCR.abscurv(T[-1])
    plt.subplot(1,2,2)
    plt.title(u"$\Delta$longueurs. (Longueur totale R Corrigé=%.5g)"%(lCR))
    plt.plot(lBR,'b.-', label=u'initial : max($\delta l)=$%.2g'%max(abs(lBR)))
    plt.plot(lBCR,'r.-', label=u'corrigé : max($\delta l)=$%.2g'%max(abs(lBCR)))
    plt.plot(lBR,'b.')
    plt.plot(lBCR,'r.')

    plt.legend()
    if show : plt.show()
#     P7 = asarray(C(T[7]))
#     C.precision*=10
#     pts = C.dpoints-P7
#     debug([(i, norm(p)) for i, p in enumerate(pts[5261:5565])])
#     debug(C)


#     msg = u"""Deux bords de pièces à assembler, de longueur $l_1$=%.4g m $\simeq$ $l_2$=%.4g m
#     Placements de %d repères de montage espacés de %.2g m.
#     L'erreur est de %.2f mm soit %.2f‰"""%(l1,l2,n,delta,1000*dd,1000*dd/lref)
#     plt.title(msg)
#     if show : plt.show()

def mainTest(show=False):
    files = [Path(VALIDATION_DIR,'unenervure2d.gnu'),
            Path(VALIDATION_DIR,'1.gnu'),
            Path(VALIDATION_DIR,'simple','simple2d1.gnu'),
            Path(VALIDATION_DIR,'unenervure2d.gnu'),
            Path(VALIDATION_DIR,'simple','anguleux.gnu'),
            Path(VALIDATION_DIR,'simple','intrados.gnu'),
            Path(VALIDATION_DIR,'simple','extrados.gnu'),
            Path(VALIDATION_DIR,'simple','profil.gnu'),
            Path(VALIDATION_DIR,'simple','demi-cercle.gnu'),
            Path(VALIDATION_DIR,'reference.pts'),
            Path(VALIDATION_DIR,'ARBIZON.PTS'),
            Path(VALIDATION_DIR,'shark.pts'),
            Path(RUNS_DIR,'profils','D2v5p2.pts'),
            Path(VALIDATION_DIR,'diamirM.pts'),
            Path(VALIDATION_DIR,'blocjonc.spl'),
            Path(VALIDATION_DIR,'spline-0.spl'),
            Path(VALIDATION_DIR,'P0.spl'),
            Path(VALIDATION_DIR,'P0#.spl'),
            ][::-1]

    debug(titre='TODO : methode=("us",k) ne marche pas')
    for filename in files[0:3] :
        if 1:
            testSequenceBasique(filename, show=show)
            debug(titre='fin testSequenceBasique %s'%filename.name)
        if 1:
            testConstructeurs(filename, show=show)
            debug(titre='fin testConstructeurs %s'%filename.name)
            exit()
        if 1:
            testSaveRead(filename, show=show)
            debug('fin testSaveRead %s'%filename.name)
        if 1:
            testMethodes(filename, show=show)
            debug('fin testMethodes %s'%filename.name)
        if 1:
            testModifGlobales(filename, show=show)
            debug('fin testModifGlobales %s'%filename.name)
        if 1:
            testModifLocales(filename, show=show)
            debug('fin testModifLocales %s'%filename.name)
        if 1:
            testEchantillonnage(filename, trados='e',show=show)
            testEchantillonnage(filename, trados='i',show=show)
            debug('fin testEchantillonnage %s'%filename.name)
        if 1:
            testElaguer(filename, trados='e', fonction='ajuster',show=show)
            debug('fin testElaguer extrados %s'%filename.name)
        if 1:
            testElaguer(filename, trados='i', fonction='ajuster',show=show)
            debug('fin testElaguer intrados %s'%filename.name)
        if 1:
            testElaguer(filename, trados='e', fonction='elaguer',show=show)
            debug('fin testElaguer extrados %s'%filename.name)
        if 1:
            testElaguer(filename, trados='i', fonction='elaguer',show=show)
            debug('fin testElaguer intrados %s'%filename.name)
#         if 1:
#             testElagage1(filename, show=show)
#             debug('fin testElagage1 %s'%filename.name)
        if 0:
            testCorrectionRM(show=show)
            debug('fin testCorrectionRM')
        if 0:
            testPlacementRM(show=show)
            debug('fin testPlacementRM')
        if 1:
            try :
                testDivers(filename, show=show)
            except AttributeError as msg:
    #             stack()
                debug('AttributeError :'+str(msg))
        #         raise
            debug('fin testDivers %s'%filename.name)
    #     sys.exit(app.exec_())
        debug(titre='fin tests %s'%filename.name)
    debug(titre='fin tests')
if __name__=='__main__':
    mainTest(show=False)

    # sys.exit(app.exec_())
#     os.system(' '.join(['python', Path('..','splinesimple.py')]))