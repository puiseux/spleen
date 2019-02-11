#!/usr/local/bin/python2.7
# encoding: utf-8
'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe Profil

@author:     Pierre Puiseux
@copyright:  2018 Nervures. All rights reserved.
__updated__="2019-02-11"
'''
import sys, cPickle
from path import Path
from pprint import pprint
from config import VALIDATION_DIR,RUNS_DIR
from utilitaires import (stack, debug, rdebug)
from utilitaires.lecteurs import pointsFrom
from numpy import asarray as array
from profil import Profil
from matplotlib import pyplot as plt
import config
config.TEST_MODE = True

def testProfil(filename,show=False):
    name = filename.name
    debug(titre='testProfil : %s'%name)

    debug(paragraphe='p = Profil()')
    p = Profil()
    pprint(p.toDump())
    exit()
    debug(paragraphe='p.open(%s)'%name)
    p.open(filename)
    print p
    debug(paragraphe="p = Profil(points=None, parent=None, naca=['2415', 50], name=None)")
    p = Profil(naca=['2415', 50], name=None)
    print p.verification()
    print p
#     print " => Constructeur presque vide, puis append(5,10):\n    +++++++++++++++++++++++++++++++++++++"
#     debug( p)
#     try : p.appendPoint((5,10))
#     except RuntimeError as msg : rdebug(msg)
#     numfig = -1
#     print " => Constructeur filename :\n    +++++++++++++++++++++"
    debug(" => constructeur np.ndarray  :\n    +++++++++++++++++++++++")
    p = Profil(cpoints=pointsFrom(filename))
    p.normalise()
    debug(p)
    for msg in p.verification() :
        print msg
#     p.plot(plt,titre='Profil(points=pointsFrom(filename))')

    debug(" => test toDump/load  :\n    +++++++++++++++++++++++")
    dump = p.toDump()
    debug(p)
#     pprint(dump)
    p = Profil(**dump)
#     p.verification()
#     p.plot(plt, titre='p = Profil(**p.toDump())')

#     p[46] = (3.3, -1.9)
#     p.plot(plt, titre='p[46] = (3.3, -1.9)')

#     p = Profil(points=points)
#     debug( p)
#     print " => constructeur QPolygon  :\n    +++++++++++++++++++++"
#     p = Profil(points=pointsFrom(p.qpolygon))
#     print p
    debug(" => constructeur recopie  :\n    +++++++++++++++++++++")
#     rdebug(pp=p.profparam)
#     try :
#         p.iouverture = p.iba-10,p.iba+10
#         p.plot(plt, titre='PROBLEM:p.iouverture = p.iba-10,p.iba+10')
#     except UnboundLocalError as msg:
#         rdebug('nombre de points de retour nptret=%d<0 ? '%p.profparam.nptret,msg=msg)
    debug(pp=p.profparam)
    p.iouverture = p.iba+5,p.iba+10
    p.pouverture = -1.0,-5.0
    q = Profil(profil=p)
#     debug(q)
#     rdebug(pp=p.profparam)
    debug(qp=q.profparam)
#     exit()
#     return
    debug('scaled')
    debug(p.scaled(2))
    #print "constructeur QPolygon :",p
    #print '########## Profil divers (BA, absCurv, extrados, ligne moyenne...) ##########'
#     p.removePoint(0)
#     curv = absCurv(p)
#     dcurv = curv[1:]-curv[:-1]
#     print dcurv
#     print 'p.absCurv()',absCurv(p)
    debug('########## Profil geometrie ##########')
#     rdebug(p.profparam)
#     exit()
    p.iouverture = p.iba+3,p.iba+10#sinon il sait pas echantillonner
    debug(p)
    if show : p.plot(plt, titre='p.iouverture = p.iba+3,p.iba+10')

    p.hardScale((2,2))
    debug(p)
    pprint(p.toDump())
    if show : p.plot(plt, titre='p.hardScale((2,2))')
    p.hardRotate(30,centre=(2,2))
    if show : p.plot(plt, titre='p.hardRotate(30,centre=(2,2))')
    p.translate((100,100))
    if show : p.plot(plt, titre='p.translate((100,100)')
#     print p.verification()
    debug( p)
    p.normalise()
    debug(p)
    if show : p.plot(plt, titre='normalise()')
    p[1] = (0.6,0.12)
    if show : p.plot(plt, titre='p[1] = (0.6,0.12)')
    p.iouverture = p.iba+2,p.iba+5
    if show : p.plot(plt, titre='ouverture = %d, %d'%(p.iba+2,p.iba+5))
    p.insertPoint((0.23,0.16))
    if show : p.plot(plt, titre='insertPoint((0.23,0.16))')
    p.removePoint(1)
    if show : p.plot(plt, titre='p.removePoint(1)')
    if show : plt.show()
#     return

    p.iouverture=p.nba+2,p.nba+3
    debug( p)
#     exit()
    debug(p.points)
    debug( p)
#     filename=Path(VALIDATION_DIR,'ARBIZON.PTS')
#     p=Profil(points=pointsFromFile(filename))
    mp=Profil(points=pointsFrom(filename))
    centre=array((7,0))
    mp.hardRotate(30,centre)
    mp.hardScale((0.5,0.5),centre)
    p.dump(Path(RUNS_DIR,'profildump.pkl'))
    f=open(Path(RUNS_DIR,'profildump.pkl'),'r')
    d=cPickle.load(f)
    pprint(d)
    p.load(d)
    debug(p)
    p.normalise()
    debug(p)
    rdebug('################## Fin testprofil ##################')
#     exit()

def testElaguer(filename,show=False):
    name = filename.name
    debug(titre='testElaguer : %s'%name)
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba)
    debug("Extrados : Se'(1.0)=%s"%(p.splines[0](1.0,1)))
    debug("Intrados : Si'(0.0)=%s"%(p.splines[1](0.0,1)))
    sin0 = p.splines[0].integraleCourbure(0.01,1,1000)
    sin1 = p.splines[1].integraleCourbure(0.01,1,1000)
    debug('sinuosite totale avant elagage=',(sin0,sin1))
#     return
    p.elaguer(eps=1, replace=True)
    debug(rba=p.rba)
    sin10 = p.splines[0].integraleCourbure(0.01,1,1000)
    sin11 = p.splines[1].integraleCourbure(0.01,1,1000)
#     sin10 = p.splines[0].integraleCourbure(0.01,1,1000)
    debug('sinuosite totale apres elagage=',(sin10,sin11))
    debug(variation_relative_sinuosite = ((sin10-sin0)/sin0,(sin11-sin1)/sin1))
#     debug('sinuosite totale =',p.splines[1].integraleCourbure(0.01,1,1000))
    if show is False:
        p.plot(plt, nbpd=[1000,1000],titre='elagage')
        p.plotCourbure()
        plt.show()

def testDivers(filename,show=False):
    name = filename.name
    debug(titre='testDivers : %s'%name)
    if 'spl' in filename.ext :
        pass
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba, corde=p.corde)
    p.normalise()
    debug(rba=p.rba, corde=p.corde)

def testSaveAndOpen(filename,show=False):
    name = filename.name
    debug(titre='testSaveAndOpen : %s'%name)
    if '.spl' in filename :
        debug(paragraphe="Ouverture fichier (.spl) %s"%name)
        with open(filename,'r') as f :
            dump = eval(f.read())
#             lines = f.read()
#         dump = eval(lines)
        S = Profil(**dump)
        print S
    else :
        debug(paragraphe="Ouverture fichier %s"%name)
        S = Profil(points=pointsFrom(filename),
    #                methode=('cubic',((2, 0, 0), (1, 0, -5))),
                   mode=['courbure'],
                   precision=[3000])

    debug(paragraphe="S.normalise() %s"%name)
    debug(ouvertures=S.pouverture)
    S.normalise()
    debug(ouvertures=S.pouverture)
    debug(paragraphe="pickle.dump() et pickle.load() %s"%name)
    fname = Path(filename.dirname, filename.namebase+'Test.spl')
    debug(fname)
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = Profil(**dump)
    S1.elaguer(1, True)
    S1.rba = (-1,-1)
    if show :
        S1.plot(plt,titre='rba=%s'%str(S1.rba))
        plt.show()
    dump  = S.toDump()
    dump1 = S1.toDump()
    debug('dump==dump1 ?', dump==dump1)
    pprint(dump1)
    debug('S.rba   = %s\n'%str(S.rba)+\
          '    S1.rba  = %s\n'%str(S1.rba)+\
          '    S.erba  = %s\n'%str(S.erba)+\
          '    S1.erba = %s\n'%str(S1.erba)+\
          '    S.drba  = %s\n'%str(S.drba)+\
          '    S1.drba = %s'%str(S1.drba))


def testEchantillonner(filename,show=False):
    name = filename.name
    debug(titre='testEchantillonner : %s'%name)
    P = Profil(points=pointsFrom(filename),
#                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
#                       mode=['linear'],
                      precision=[3000])
    P.normalise()
    P.elaguer(2, True)
    P.pouverture = pam, pav = -10, -40#en % de corde
    debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)
    debug(ouverture=P.ouverture)

    debug('P._getT(%.f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
#     P.echantillonner()
    if show : P.plot(plt, titre='echantillonage : pouv=%s, touv=%s'%(str(P.pouverture),str(P.touverture)))
    P.hardScale(2, centre=array([0,0]))
    P.hardRotate(180, centre=(0,0))
#     debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)#,P_t=sint(t))
    debug(ouverture=P.ouverture)
    debug('P._getT(%.1f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
    if show : P.plot(plt,titre='echantillonage : rotation 180')#%(str(P.pouverture),str(P.touverture)))

def testOuverture(filename,show=False):
    name = filename.name
    debug(titre='testOuverture : %s'%name)
    P = Profil(points=pointsFrom(filename), name='Ouverture')
    if show : P.plot0(plt, titre='testOuverture:%s'%P.name)
    pouv = (-20, -50)
    P.pouverture = pouv
    if show : P.plot0(plt, titre='ouverture=%s'%(str(pouv)))

def testMain(show = False):
    p = Profil()#constructeur vide
    p.name = u'Vide (3 points)'
    debug(titre='Constructeur Vide')
    print p
    pprint(p.toDump())
    files = [
            Path(VALIDATION_DIR,'diamirE.spl'),
            Path(VALIDATION_DIR,'unenervure2d.gnu'),
            Path(VALIDATION_DIR,'E-diamirE.spl'),
            Path(VALIDATION_DIR,'E-shark.spl'),
            Path(VALIDATION_DIR,'shark.spl'),
            Path(VALIDATION_DIR,'naca2415.spl'),
            Path(VALIDATION_DIR,'spline-0#.pkl'),
            Path(VALIDATION_DIR,'spline-0#.spl'),
            Path(VALIDATION_DIR,'spline-0.spl'),

            ]

    for filename in files[:1] :
        if 1:testProfil(filename, show=show)
        if 1:testSaveAndOpen(filename, show=show)
        if 1:testOuverture(filename, show=show)
        if 1:testEchantillonner(filename, show=show)
        if 1:testElaguer(filename, show=show)
        if 1:testDivers(filename, show=show)

if __name__=="__main__":
    testMain()