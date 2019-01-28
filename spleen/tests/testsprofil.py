#!/usr/local/bin/python2.7
# encoding: utf-8
'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe Profil

@author:     Pierre Puiseux
@copyright:  2018 Nervures. All rights reserved.
'''
__updated__="2018-07-01"
import sys, cPickle
from path import Path
from pprint import pprint
from config import VALIDATION_DIR,RUNS_DIR
from utilitaires import (stack, debug, rdebug)
from lecteurs import pointsFrom, pointsFromFile
from numpy import asarray as array
from profil import Profil

def testProfil(filename):
#     p = Profil(points=None, parent=None, naca=['2415', 50], name=None)
#     print p.verification()
#     print p
    from matplotlib import pyplot as plt
#     print " => Constructeur presque vide, puis append(5,10):\n    +++++++++++++++++++++++++++++++++++++"
    p = Profil()
#     debug( p)
#     try : p.appendPoint((5,10))
#     except RuntimeError as msg : rdebug(msg)
#     numfig = -1
#     print " => Constructeur filename :\n    +++++++++++++++++++++"
    debug(" => constructeur np.ndarray  :\n    +++++++++++++++++++++++")
    p = Profil(points=pointsFrom(filename))
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
    p.plot(plt, titre='p.iouverture = p.iba+3,p.iba+10')

    p.hardScale((2,2))
    debug(p)
    pprint(p.toDump())
    p.plot(plt, titre='p.hardScale((2,2))')
    p.hardRotate(30,centre=(2,2))
    p.plot(plt, titre='p.hardRotate(30,centre=(2,2))')
    p.translate((100,100))
    p.plot(plt, titre='p.translate((100,100)')
#     print p.verification()
    debug( p)
    p.normalise()
    debug(p)
    p.plot(plt, titre='normalise()')
    p[1] = (0.6,0.12)
    p.plot(plt, titre='p[1] = (0.6,0.12)')
    p.iouverture = p.iba+2,p.iba+5
    p.plot(plt, titre='ouverture = %d, %d'%(p.iba+2,p.iba+5))
    p.insertPoint((0.23,0.16))
    p.plot(plt, titre='insertPoint((0.23,0.16))')
    p.removePoint(1)
    p.plot(plt, titre='p.removePoint(1)')
    plt.show()
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

def testElaguer(filename):
    from matplotlib import pyplot as plt
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba)
    debug('Extrados : Se\'(1.0)=%s'%(p.splines[0](1.0,1)))
    debug('Intrados : Si\'(0.0)=%s'%(p.splines[1](0.0,1)))
    debug('sinuosite=',p.splines[1].integraleCourbure(0.01,1,1000))
#     return
    p.elaguer(eps=1, replace=True)
    debug(rba=p.rba)
    debug('sinuosite=',p.splines[1].integraleCourbure(0.01,1,1000))
    p.plot(plt, nbpd=[1000,1000],titre='elagage')
    plt.show()

def testDivers(filename):
    if 'spl' in filename.ext :
        pass
    p = Profil(points=pointsFrom(filename), precision=[1000])
    debug(rba=p.rba, corde=p.corde)
    p.normalise()
    debug(rba=p.rba, corde=p.corde)

def testSaveAndOpen(filename):
    from matplotlib import pyplot as plt
    if '.spl' in filename :
        with open(filename,'r') as f :
            lines = f.read()
#         for line in lines :
        dump = eval(lines)
        S = Profil(**dump)
        print S
    else :
        S = Profil(points=pointsFrom(filename),
#                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
                      mode=['courbure'],
                      precision=[3000])
    S.normalise()
    fname = Path(filename.dirname, filename.namebase+'Test.spl')
    debug(fname)
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = Profil(**dump)
    S1.elaguer(1, True)
    S1.rba = (-1,-1)
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


def testEchantillonner(filename):
    from matplotlib import pyplot as plt
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
    P.plot(plt, titre='echantillonage : pouv=%s, touv=%s'%(str(P.pouverture),str(P.touverture)))
    P.hardScale(2, centre=array([0,0]))
    P.hardRotate(180, centre=(0,0))
#     debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)#,P_t=sint(t))
    debug(ouverture=P.ouverture)
    debug('P._getT(%.1f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
    P.plot(plt,titre='echantillonage : rotation 180')#%(str(P.pouverture),str(P.touverture)))

def testOuverture(filename):
    from matplotlib import pyplot as plt
    P = Profil(points=pointsFrom(filename), name='Ouverture par defaut')
    P.plot0(plt, titre='testOuverture:%s'%P.name)
    pouv = (-20, -50)
    P.pouverture = pouv
    P.plot0(plt, titre='ouverture=%s'%(str(pouv)))

def testMain():
#     sys.path.append(Path('..').abspath)
    p = Profil()#constructeur vide
    debug('Constructeur Vide', p)
    filename = Path(VALIDATION_DIR,'1.gnu')
    filename = Path(VALIDATION_DIR,'ARBIZON.PTS')
    filename = Path(VALIDATION_DIR,'faial2.dxf')
    filename = Path(VALIDATION_DIR,'simple','simple2d1.gnu')
    filename = Path(VALIDATION_DIR,'unenervure2d.gnu')
    filename = Path(VALIDATION_DIR,'simple','profil.gnu')
    filename = Path(RUNS_DIR,'profils','D2v5p2.pts')
    filename = Path(VALIDATION_DIR,'reference.pts')
    filename = Path(RUNS_DIR,'P0.spl')
    testOuverture(filename)
#     exit()
    testProfil(filename)
    testEchantillonner(filename)
    testElaguer(filename)
    testSaveAndOpen(filename)
    testDivers(filename)

if __name__=="__main__":
    testMain()