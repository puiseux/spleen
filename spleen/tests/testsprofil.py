#!/usr/local/bin/python2.7
# encoding: utf-8
'''
AXile -- Outil de conception/simulation de parapentes Nervures

    Classe Profil

@author:     Pierre Puiseux
@copyright:  2018 Nervures. All rights reserved.
__updated__="2019-02-17"
'''
import cPickle
from pprint import pprint
from spleen.path import Path
from spleen.utilitaires import (debug, rdebug, dictsAreNotEqual)
from spleen.utilitaires.lecteurs import pointsFromFile
from spleen.profils import Profil
from numpy import asarray as array, where
from matplotlib import pyplot as plt

def testProfil(filename,show=True):

    def expose(par) :
        print p
        print 'cpoints = asarray(%s)'%p.cpoints.tolist()
        ep = p.epoints.tolist()
        print 'epoints = asarray(%s)'%ep
        print 'techext = asarray(%s)'%p.techext.tolist()
        print 'techint = asarray(%s)'%p.techint.tolist()
        if show : p.plot(titre=par, show=True)
        return ep
    name = filename.name
    debug(titre='testProfil : %s'%name)

    par='p = Profil()'
    debug(paragraphe=par)
    p = Profil()
    expose(par)

    par = "p = Profil(naca=['2415', 50], name=None)"
    debug(paragraphe=par)
    p = Profil(naca=['2415', 20], name=None)
    expose(par + 'verification')
    for v in p.verification() : print v

    par = 'p.open("%s")'%name
    debug(paragraphe=par)
    try :
        p.open(filename)
        expose(par)
    except ValueError as msg :
        rdebug("TODO : %s, %s"%(par, msg))
    except IndexError as msg :
        rdebug("TODO : %s, %s"%(par, msg))
#     exit()
    par = 'Profil(points) (%s)'%name
    debug(paragraphe=par)
    p = Profil(cpoints=pointsFromFile(filename), name=name)
    expose(par)

    par = "toDump/load (%s)"%name
    debug(paragraphe=par)
    dump = p.toDump()
    p = Profil(**dump)
    expose(par)
    dump1 = p.toDump()
    debug(dump_NOT_equal_dump1=dictsAreNotEqual(dump,dump1))

    iouv = p.iba+5, p.iba+10
    ouv = -1.0, -10.0
    par = 'p.pouverture = %s; p.iouverture = %s'%(ouv,iouv)
    debug(paragraphe=par)
    debug(pp=p.profparam)
    p.iouverture = iouv
    p.pouverture = ouv
    ep0 = expose(par)

    par = "constructeur recopie (%s)"%name
    p = Profil(profil=p)
    ep = expose(par)
    debug(where(ep!=ep0))

    pt, i = (0.55, 0.1), 10
    par = "%s : insert(%s,%d)"%(name,pt,i)
    debug(paragraphe=par)
    p.insertPoint(pt, i)
    expose(par)

    par = '%s.scaled(2.0)'%name
    debug(paragraphe=par)
    p = p.scaled(2.0)
    expose(par)

    par = '%s.hardScale((2,2))'%name
    debug(paragraphe=par)
    p.hardScale((2,2))
    expose(par)

    par = '%s.hardRotate(30,centre=(2,2))'%name
    debug(paragraphe=par)
    p.hardRotate(30,centre=(2,2))
    expose(par)

    par = 'p.translate((100,100)'
    debug(paragraphe=par)
    p.translate((100,100))
    expose(par)
#     print p.verification()
    par = 'normalise()'
    debug(paragraphe=par)
    p.normalise()
    expose(par)

    par = 'p[1] = (0.6,0.12)'
    debug(paragraphe=par)
    p[1] = (0.6,0.12)
    expose(par)
    par = 'p.removePoint(1)'
    debug(paragraphe=par)
    p.removePoint(1)
    expose(par)
#     return

    par = "dump('profildump.pkl') puis cPickle.load()"
    debug(paragraphe=par)
    p.dump(Path(RUNS_DIR,'profildump.pkl'))
    f=open(Path(RUNS_DIR,'profildump.pkl'),'r')
    d=cPickle.load(f)
    pprint(d)
    p.load(d)
    expose(par)
    debug(p)
    p.normalise()
    expose('p.normalise')
    debug(p)
    debug(titre='Fin testprofil')
#     exit()

def testElaguer(filename,show=False):
    name = filename.name
    debug(titre='testElaguer : %s'%name)
    p = Profil(points=pointsFromFile(filename), precision=[1000,1000])
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
    if show :
        p.plot(titre='elagage')
        p.plotCourbure()
        plt.show()
    debug(titre='Fin testElaguer')

def testDivers(filename,show=False):
    name = filename.name
    debug(titre='testDivers : %s'%name)
    if 'spl' in filename.ext :
        pass
    p = Profil(points=pointsFromFile(filename), precision=[1000, 1000])
    debug(rba=p.rba, corde=p.corde)
    p.normalise()
    debug(rba=p.rba, corde=p.corde)
    debug(titre='Fin testDivers')

def testSaveAndOpen(filename,show=False):
    name = filename.name
    titre ='testSaveAndOpen : %s'%name
    debug(titre=titre)
    if '.spl' in filename :
        par = "Ouverture fichier (.spl) %s"%name
        debug(paragraphe=par)
#         with open(filename,'r') as f :
#             dump = eval(f.read())
#             lines = f.read()
#         dump = eval(lines)
        S = Profil()
        try :
            S.open(filename)
            print S
        except IndexError as msg:
            rdebug('TODO:%s, %s'%(par, msg))
            return#On peut pas continuer S est verolé
    else :
        debug(paragraphe="Ouverture fichier %s"%name)
        S = Profil(points=pointsFromFile(filename),
    #                methode=('cubic',((2, 0, 0), (1, 0, -5))),
                   mode=['courbure', 'courbure'],
                   precision=[3000,3000])

    debug(paragraphe="S.normalise() %s"%name)
    debug(ouvertures=S.pouverture)
    S.normalise()
    debug(ouvertures=S.pouverture)
    debug(paragraphe="pickle.dump() et pickle.load() %s"%name)
    fname = Path(WORK_DIR, filename.namebase+'Test.spl')
    debug(fname)
    with open(fname,'w') as f :
        cPickle.dump(S.toDump(),f)
    with open(fname,'r') as f :
        dump = cPickle.load(f)
        S1 = Profil(**dump)
    S1.elaguer(1, True)
    S1.rba = (-1,-1)
    if show :
        S1.plot(titre='rba=%s'%str(S1.rba))
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
    debug(titre='Fin testSaveAndOpen')


def testEchantillonner(filename,show=False):
    name = filename.name
    debug(titre='testEchantillonner : %s'%name)
    P = Profil(points=pointsFromFile(filename),
#                       methode=('cubic',((2, 0, 0), (1, 0, -5))),
#                       mode=['linear'],
                      nbpd=[1000,1000])
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
    if show : P.plot(titre='echantillonage : pouv=%s, touv=%s'%(str(P.pouverture),str(P.touverture)))
    P.hardScale(2, centre=array([0,0]))
    P.hardRotate(180, centre=(0,0))
#     debug(P)
    touv  = P.touverture #= (kam, tam), (kav,tav)
    debug(touv=touv)#,P_t=sint(t))
    debug(ouverture=P.ouverture)
    debug('P._getT(%.1f%%)=%f'%(pam,P._getT(pam)))
    debug('P._getT(%.1f%%)=%f'%(pav,P._getT(pav)))
    if show : P.plot(titre='echantillonage : rotation 180')#%(str(P.pouverture),str(P.touverture)))
    debug(titre='Fin testEchantillonner')

def testOuverture(filename,show=False):
    name = filename.name
    debug(titre='testOuverture : %s'%name)
    P = Profil(points=pointsFromFile(filename), name='Ouverture')
    if show : P.plot(titre='testOuverture:%s'%P.name)
    pouv = (-20, -50)
    P.pouverture = pouv
    if show : P.plot(titre='ouverture=%s'%(str(pouv)))
    debug(titre='Fin testOuverture')

def testMain(show = False):
    p = Profil()#constructeur vide
    p.name = u'Vide (3 points)'
    debug(titre='Constructeur Vide')
    print p
    pprint(p.toDump())
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
            ]

    for filename in files[:1] :
        if 1:testProfil(filename, show=show)
        if 1:testSaveAndOpen(filename, show=show)
        if 1:testOuverture(filename, show=show)
        if 1:testEchantillonner(filename, show=show)
        if 1:testElaguer(filename, show=show)
        if 1:testDivers(filename, show=show)
    debug(titre='Fin testMain')

if __name__=="__main__":
    from spleen.spleenconfig import VALIDATION_DIR, RUNS_DIR, WORK_DIR
    from spleen import spleenconfig

    spleenconfig.TEST_MODE = True
    testMain()