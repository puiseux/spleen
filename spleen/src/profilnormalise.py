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
from utilitaires import (debug, rdebug,dist2,dist)
import numpy as np
# import scipy as sp
from numpy import linspace, log
from scipy.optimize import newton
from lecteurs import pointsFrom#, qpolygonFrom
from preferences import ProfilPrefs
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

class ProfilNormalise(Profil):
    prefs = ProfilPrefs()
    u"""
    Un ProfilNormalise est un Profil tel que :
    - la corde vaut 1.0,
    - le BA est en (0.0,0.0) est INNAMOVIBLE.
    - le BF en (1.0, 0.0) est INNAMOVIBLE.
    on ne peut pas le normaliser, sauf une fois à l'initialisation.
    Tant qu'on n'ajoute ni ne retranche de point à l'extrados, son nba reste constant.
    """

    def __init__(self,**dump):
        u"""
        """
        super(ProfilNormalise, self).__init__(**dump)
        self.nb_normalisations = 0
        if len(self) <= 3 :
            self._nba = 1
#             self.qpolygon = QPolygonF([QPointF(), QPointF(1,0), QPointF()])
            self._cpoints = self._epoints = asarray([(0.,0.),(1.,0.),(0.,0.)])
            self.nb_normalisations = 1
        else :
            self.__normalise()
            self.normalite()

    def scaled(self,scale):
        '''
        Retourne une COPIE de self, mise à l'échelle scale.
        Pas de modif de self
        '''
        prof = ProfilNormalise(self)
        prof.hardScale((scale,scale))#on a le droit de le normaliser une seule fois
        prof.name = '%.2gx'%scale
        return prof

    def hardScale(self, scale, centre=None):
        if self.nb_normalisations == 0 :
#             debug("*****", scale=scale)
            Profil.hardScale(self, scale, centre)
        else :
            raise RuntimeError("%s.hardScale() : operation impossible"%self.classname)
            rdebug("%s.hardScale() : operation impossible"%self.classname)

    def hardRotate(self, angle,centre, unit='degres'):
        if self.nb_normalisations == 0 :
            Profil.hardRotate(self, angle, centre, unit)
        else :
            rdebug("%s.hardRotate() : operation impossible"%self.classname)

    def __setitem__(self, k, value):
        if k in (0, self.nba, -1, len(self)-1) and self.nb_normalisations>0:
            rdebug("%s.__setitem__() : operation impossible"%self.classname)
        else :
            Profil.__setitem__(self, k, value)

    def _getT(self, x, t0=None, nbit=False, dt=[-0.1, 1.1], tol=1.0e-10, maxiter=50):
        u"""
        Retourne la valeur du paramètre t correspondant à l'abscisse |x|/100.
        Si t est non trouvé, retourne np.nan
        Plus précisement la recherche se fait dans
            S = extrados si x>0
            S = intrados si x<0
            le t retourné est l'unique t tel que |x|/100 == S.sx(t)
        - t0 : la recherche se fait par iterations de Newton, elle démarre à t=t0
        - nbit : si nbit est True, retourne le couple (t, nb iterations)
                si nbit est False, retourne t
        - dt est l'intervalle de recherche.
        """
        if x == 0.0 : return (np.nan, 0) if nbit else np.nan
        ax = abs(x)/100.0#abscisse du point situé à ax% de corde
#         BA, BF = self[self.nba], self[0]
#         P = BA + ax*(BF-BA)
        if ax > 1.0 :
            raise ValueError("∣x∣=%g devrait etre dans [0,100]. C'est une abscisse en %% de corde"%abs(x))
        if x>=0 : S = self.splines[0]
        else    : S = self.splines[1]
        if t0 is None :
            if x>0 : t0 = 1.0 - ax
            elif x<0 : t0 = ax
        k = 0
        while 1 :
            t = newton(lambda t: S.sx(t)-ax, t0, lambda t:S.sx(t,1),
                       tol=tol, maxiter=50, fprime2=lambda t:S.sx(t,2))
            k += 1
#             #pour que le BA et le BF ne bougent pas, à 0.01 mm
#             if abs(t)<=1.0e-5 :
#                 t = 0.0#0.01 mm
#             elif abs(t-1)<=1.0e-5 :
#                 t = 1.0
            return (t,k) if nbit else t
            if not dt[0] <= t <= dt[1] : #dépasse les bornes
                return (np.nan, k) if nbit else np.nan
            t0 += 0.1
    def insertPoint(self, pos, k=None):
#         i = super(Profil, self).insertPoint(pos)
        if dist(pos, self[0]) >= 1 :
            rdebug("impossible d'inserer un point a une distance du BF >= 1")
        else :
            return Profil.insertPoint(self, pos, k)

    def pointValide(self, p):
        u""" retourne True si p peut faire partie du profil, i.e.
        - distance(p, BF)<1
        - xp compris entre 0 et 1
        """
        return dist2((1,0),p) < 1.0 and 0<p[0]<1

    def normalise(self):
        """on ne renormalise pas un profil normalisé mais
        normalise() est parfois appelé par Profil (dans le update en particulier)"""
#         debug("%s.normalise() : ne fait rien"%self.classname)
        return

    def normalite(self):
        u"""Vérifie que le profil est bien normalisé, EXACTEMENT normalisé,
        leve une exception si ca n'est pas le cas"""
        bf0, bf1, ba = self[0], self[len(self)-1], self[self.nba]
#         debug(bf0, bf1,ba)
        msg = []
        if (bf0[0],bf0[1]) != (1.0, 0.0) :
#             raise RuntimeError("%s non normalisé: self[0]=(%.5g,%.5g) != (0,1) "%(self.classname,bf0[0],bf0[1]))
            msg.append("%s non normalise: self[0]=(%.5g,%.5g) != (0,1) "%(self.classname,bf0[0],bf0[1]))
        elif (bf1[0],bf1[1])!= (1.0, 0.0) :
#             raise RuntimeError("%s non normalisé: self[-1]=(%.5g,%.5g) != (0,1) "%(self.classname,bf1[0],bf1[1]))
            msg.append("%s non normalise: self[-1]=(%.5g,%.5g) != (0,1) "%(self.classname,bf1[0],bf1[1]))
        elif (ba[0],ba[1]) != (0.0, 0.0):
#             raise RuntimeError("%s non normalisé: self[0]=(%.5g,%.5g) != (0,0) "%(self.classname,ba[0],ba[1]))
            msg.append("%s non normalise: self[0]=(%.5g,%.5g) != (0,0) "%(self.classname,ba[0],ba[1]))
        corde, nba = self.corde, self.nba#self.computeCordeAndNBA()
        if abs(corde-1.0)>=self.prefs.EPS:# or nba != self.nba :
#             raise RuntimeError("%s non normalisé: corde=%.5g != 1.0 ou nba=%d != self.nba=%d"%(self.classname,corde,nba, self.nba))
            msg.append("%s non normalise: corde=%.5g != 1.0 ou nba=%d != self.nba=%d"%(self.classname,corde,nba, self.nba))
        else :
            msg.append('%s : normalite OK'%self.classname)
        return '\n'.join(msg)

    def echantillonner(self):
        u"""l'échantillonnage de Profil peut déplacer légèrement le BA et le BF.
        Ici, on les repositionne à (0,0) et (1,0)"""
        Profil.echantillonner(self)
#         rdebug(iba=self.iba, len_epoints=len(self._epoints))
        self._epoints[self.iba] = [0,0]
        self._epoints[-1] = self._epoints[0]= [1,0]
        return self._epoints

    def plotCourbure(self):
        from matplotlib import pyplot as plt
        from matplotlib.widgets import CheckButtons
#         nbpd = self.precision
        nbpe = self.nbpe
        self.echantillonner()#nbpe)
        D = self.dpoints
        C = self.cpoints
        E = self.epoints
        _, ax = plt.subplots()
        titre = self.name+' courbure'
        plt.title(titre)
        T = linspace(0,1,100)
        spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
#         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
        control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        cext = self.splines[0].courbure(T)
        cext += (1.0 + abs(min(cext)))
        cext = log(cext)
        cext /= max(abs(cext))
        extcourbure, = ax.plot(T, cext[::-1])
        cint = self.splines[1].courbure(T)
        cint += (1.0 + abs(min(cint)))
        cint = log(cint)
        cint /= max(abs(cint))
#         cint /= max(cint)
        intcourbure, = ax.plot(T, cint)
        buttons = ['points controle','courbure (extrados)','courbure (intrados)', 'spline',]
        values = [True, True, True, True]
        draws = [control, extcourbure, intcourbure, spline]
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        rax = plt.axes([0.05, 0.4, 0.1, 0.15])
        check = CheckButtons(rax, buttons, values)

        def func(label):
            if label == 'spline':
                spline.set_visible(not spline.get_visible())
            elif label == 'points controle':
                control.set_visible(not control.get_visible())
            elif label == 'courbure (extrados)':
                extcourbure.set_visible(not extcourbure.get_visible())
            elif label == 'courbure (intrados)':
                intcourbure.set_visible(not intcourbure.get_visible())
            else :
                draw = draws[buttons.index(label)]
                draw.set_visible(not draw.get_visible())
            plt.draw()
        check.on_clicked(func)
        plt.show()
        return plt

    def update(self):
        u'''
        Mise à jour de nba et profparam.
        Doit être appelé à chaque modification (suppression, insertion, deplacement) d'un point du profil
        - suppression, insertion de point : on reconstruit profparam entier.
        - déplacement de point : seul nba peut changer.
        On appelle le _update de NSplineComposee (i.e. celui de NSplineAbstraite),
         ************************
         mais PAS celui de Profil (Qui n'existe pas à l'heure ou je parle)
         ************************
        '''
        return super(Profil,self)._update()
#         #Le nba a possiblement changé => non, il a ete modifié
# #         if hasattr(self, '_nba') :
# #             del self._nba
#         try :
#             ppold = self.profparam
# #            trace(self, 'Prof, ppold', ppold)
#         except AttributeError :
#             return
#         #mise à jour de iba
#         if ppold.nptprof == len(self) :#déplacement de noeud
#             self.normalite()
#         else : #suppression ou insertion de noeud : on ne peut pas maintenir les parametres
#             nbanew = self.nba#a faire en premier car recalcule self._nba
#             nptext = 1 + nbanew# + nptret
#             nptint = len(self) - nptext
#             self.profparam = ProfsParam(nptext, nptint, nbanew)
# #            trace(self, 'Prof:self.profparam', self.profparam)
#         try : del self.smoothpoints
#         except AttributeError : pass

    def __normalise(self):
        if self.nb_normalisations>0 :
            raise RuntimeError("%s.__normalise() : profil deja normalise %d fois"%(self.classname,self.nb_normalisations))
            rdebug("%s.__normalise() : profil deja normalise %d fois"%(self.classname,self.nb_normalisations))
        else :
            res = super(ProfilNormalise, self).normalise()
            self.nb_normalisations += 1
            return res

    def toDump(self, format_='new'):
        '''Lors d'un dump d'une base de profils, le profil self est sauvée sous cette forme'''
        dump = super(ProfilNormalise, self).toDump(format_)
#         dump['profparam'] = self.profparam.dump
        return dump

def testProfilNormalise():
    from config import DATA_DIR,VALIDATION_DIR, WORK_DIR
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
    filename = Path(WORK_DIR,'D2V1-Pb.txt')
    filename = Path(VALIDATION_DIR,'P0.spl')
    with open(filename) as f :
        dump = eval(f.read())
    print dump
    p = ProfilNormalise(**dump)
    print p
if __name__=="__main__":
    ##########
#     import matplotlib
#     app=QtGui.QApplication(sys.argv)
    p=ProfilNormalise()#constructeur vide
    debug(p)
    testDivers()
    testProfilNormalise()
    print '################## FIN main #################'
#     sys.exit(app.exec_())
