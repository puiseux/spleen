#!/usr/local/bin/python2.7
# encoding: utf-8
from utilitaires.utilitairesdivers import dictsAreNotEqual
from testsplinecomposee import pprintSaufCpoints
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
from profils import ProfilNormalise, ProfsParam1
from utilitaires import (debug, rdebug,)
from utilitaires import pointsFrom#, qpolygonFrom
from config import VALIDATION_DIR
from matplotlib import pyplot as plt
import config
config.TEST_MODE = False
def testProfilNormalise(filename, show):
    name = filename.name
    debug(titre="testProfilNormalise(%s)"%name)
    debug(paragraphe="p=ProfilNormalise()  (%s)"%name)
    p = ProfilNormalise()
    debug(paragraphe="(%s).open(filename)"%name)
    p.open(filename)
    debug(p)
    if show : p.plot(titre=u'%s : profil original, à normaliser'%p.name)
    debug(paragraphe="p=ProfilNormalise(p.toDump())  (%s)"%name)
    dump = p.toDump()
    p = ProfilNormalise(**dump)
    debug(p)
    iba = 40
    iouverture = (50, 60)
    par = u"p.iba=%d; p.iouverture=%s  (%s)"%(iba, iouverture, name)
    debug(paragraphe=par)
    p.iba = iba
    p.iouverture = iouverture
    if show : p.plot(titre=par)

    debug(p.name)
    debug(epoints=p.epoints.tolist())
    if show : p.plot(titre=u'%s normalisé'%p.name)

    debug(normalite=p.normalite())

    par = u"q = ProfilNormalise(**p.toDump())  (%s)"%(name)
    debug(paragraphe=par)
    dumpp = p.toDump()
    q = ProfilNormalise(**dumpp)
#     print p
#     pprintSaufCpoints(dumpp)
#     print q
    debug(p_info_egale_q_info=(str(p)==str(q)))
    dumpq = q.toDump()
    debug(dump_p_DIFFERENT_DE_dump_q=dictsAreNotEqual(dumpp, dumpq))
#     pprintSaufCpoints(dumpq)
    par = u"q = ProfilNormalise(**p.toDump())  (%s)"%(name)
    debug(paragraphe=par)
    p = ProfilNormalise(**p.toDump())
    dumpp1 = p.toDump()
    debug(dumpp1_DIFFERENT_DE_dumpp=dictsAreNotEqual(dumpp, dumpp1))
    debug(paragraphe=par)
    par = u"p.verification  (%s)"%(name)
    for val in p.verification():
        print val

    print p.normalite()
    p.normalise()
    p.iouverture=p.nba+2,p.nba+3
    print p.normalite()
#     exit()
    print p.cpoints.tolist()
    mp = ProfilNormalise(points=pointsFrom(filename), name=name)
    if show : 
        p.plot(show=show, titre='p = Profil(**p.toDump())')
        p.plotCourbure()#(plt, titre='p = Profil(**p.toDump())')
        
    titre = u'ProfilNormalise opérations interdites (%s)'%name
    debug(titre=titre)
    
    par = '    #hardScale'
    debug(paragraphe=par)
    try : p.hardScale((2,2))
    except RuntimeError as msg : print msg
    print p.normalite()
#     return
    centre=asarray((7,0))
    par = '#hardRotate'
    debug(paragraphe=par)
    try : mp.hardRotate(30,centre)
    except RuntimeError as msg : print msg
    print p.normalite()
    par = '#appendPoint'
    debug(paragraphe=par)
    try : p.appendPoint((5,10))
    except RuntimeError as msg : print msg
    print p.normalite()
    par = '#p[0] = (1,2)'
    debug(paragraphe=par)
    try : p[0] = (1,2)
    except ValueError as msg : print msg
    print p.normalite()
    par = '#insertPoint hors limites'
    debug(paragraphe=par)
    try : p.insertPoint((1.1, 0.2))
    except Exception as msg : print msg
    print p.normalite()
    par = '#insertPoint autorise'
    debug(paragraphe=par)
    try : p.insertPoint((0.5, 0.2))
    except Exception as msg : print msg
    print p.normalite()
    debug(titre="Fin testProfilNormalise(%s)"%name)
    
#     exit()
#     p = ProfilNormalise(points=([0,0],[1,0],[1,1]), parent=None)
def testDivers(filename, show):
    name = filename.name
    debug(titre="testDivers (%s)"%name)
    p = ProfilNormalise()
    p.open(filename)
    print p
    if show : p.plot()
    debug(p.echantillonner().tolist())


def testPinces(filename, show):
    pourcentages=((10,10.1), (5,5), (7.3,7.3), (5.5,5.5))
    u"""
    pourcentages = Les couples de pourcentage concernent les % gauche et droite
    de BF_xt, BA_ext, BA_int, BF_int
    On bouge 4 points de l'échantillonnage du profil pour les positionner à
        - pbai% du BA (intrados)
        - pbfi% du BF(intrados)
        - pbae% du BA (extrados)
        - pbfe% du BF(extrados)
    """
    name = filename.name
    debug(titre="testPinces(%s)"%name)
    P = ProfilNormalise(cpoints=pointsFrom(filename), name=name)
    pp = ProfsParam1(nptprof=150, iouvext=70, iouvint=76, iba=60)
    P.profparam = pp
    debug(pp=pp)
    debug(P)
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
        if 0 : testPinces(filename, show=True)
        if 0 : testDivers(filename, show=True)
        if 1 : testProfilNormalise(filename, show=False)
    print '################## FIN main #################'

if __name__=="__main__":
    testMain()
#     sys.exit(app.exec_())
