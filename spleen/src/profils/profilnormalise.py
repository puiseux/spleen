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

from profils import Profil
from utilitaires import (debug,rdebug,dist2,dist,className)
import numpy as np
from numpy import linspace, log, asarray
from scipy.optimize import newton
from preferences import ProfilPrefs
import config

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
    class Default(Profil.Default) :
        u"""Un dictionnaire avec les valeurs par défaut"""
        def __init__(self) :
#             prefs = ProfilPrefs
            Profil.Default.__init__(self)
            self.update(name='ProfilNormalise', role='ProfilNormalise')
            self.name      = 'ProfilNormalise'
            self.role      = 'ProfilNormalise'

    def __init__(self,**dump):
        u"""
        """
        super(ProfilNormalise, self).__init__(**dump)
#         self.nb_normalisations = 0
#         if len(self) <= 3 :
#             self._nba = 1
#             self._cpoints = self._epoints = asarray([(0.,0.),(1.,0.),(0.,0.)])
#             self.nb_normalisations = 1
#         else :
#             self.__normalise()
#             self.normalite()
    
    def load(self, dump):
        super(ProfilNormalise,self).load(dump)
        self.nb_normalisations = 0
        if len(self) <= 3 :
            self._nba = 1
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
            raise RuntimeError("%s.hardScale() : operation impossible"%className(self))
            rdebug("%s.hardScale() : operation impossible"%className(self))

    def hardRotate(self, angle,centre, unit='degres'):
        if self.nb_normalisations == 0 :
            Profil.hardRotate(self, angle, centre, unit)
        else :
            rdebug("%s.hardRotate() : operation impossible"%className(self))

    def __setitem__(self, k, value):
        if k in (0, self.nba, -1, len(self)-1) and self.nb_normalisations>0:
            rdebug("%s.__setitem__(%d) : operation impossible"%(className(self),k))
        else :
            Profil.__setitem__(self, k, value)

    def _getT(self, x, t0=None, nbit=False, dt=[-0.1, 1.1], tol=1.0e-10, 
              maxiter=50,full_output=False):
        u"""
        Retourne la valeur du paramètre t correspondant à l'abscisse |x|/100.
        Si t est non trouvé, retourne np.nan
        Plus précisement la recherche se fait dans
            S = extrados si x>0
            S = intrados si x<0
            le t retourné est l'unique t tel que |x|/100 == S.sx(t)
        :param t0 : float, en % de corde, la recherche se fait par iterations 
            de Newton, elle démarre à t=t0% de corde
        :param nbit : si nbit est True, retourne le couple (t, nb iterations)
                si nbit est False, retourne t
        :param dt: est l'intervalle de recherche.
        """
        if x == 0.0 : return (np.nan, 0) if nbit else np.nan
        ax = abs(x)/100.0#abscisse du point situé à ax% de corde
        if ax > 1.0 :
            raise ValueError(u"∣x∣=%g devrait etre dans [0,100]. C'est une abscisse en %% de corde"%abs(x))
        S = self.splines[0] if x >= 0 else self.splines[1]
        if t0 is None :
            t0 = 1.0 - ax if x>0 else ax
        k = 0
        while 1 :
            t, r = newton(lambda t: S.sx(t)-ax, t0, lambda t:S.sx(t,1),
                       tol=tol, maxiter=50, fprime2=lambda t:S.sx(t,2),
                       full_output=True)
            
            k += 1
#             debug(r)
            return (t,r) if nbit else t
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
        msg = u"""Aucune action. On ne renormalise pas un profil normalisé"""
#         debug("%s.normalise() : ne fait rien"%className(self))
        rdebug(msg)
        return

    def normalite(self):
        u"""Vérifie que le profil est bien normalisé, EXACTEMENT normalisé,
        leve une exception si ca n'est pas le cas"""
        bf0, bf1, ba = self[0], self[len(self)-1], self[self.nba]
#         debug(bf0, bf1,ba)
        msg = []
        if (bf0[0],bf0[1]) != (1.0, 0.0) :
#             raise RuntimeError("%s non normalisé: self[0]=(%.5g,%.5g) != (0,1) "%(className(self),bf0[0],bf0[1]))
            msg.append("%s non normalise: self[0]=(%.5g,%.5g) != (0,1) "%(className(self),bf0[0],bf0[1]))
        elif (bf1[0],bf1[1])!= (1.0, 0.0) :
#             raise RuntimeError("%s non normalisé: self[-1]=(%.5g,%.5g) != (0,1) "%(className(self),bf1[0],bf1[1]))
            msg.append("%s non normalise: self[-1]=(%.5g,%.5g) != (0,1) "%(className(self),bf1[0],bf1[1]))
        elif (ba[0],ba[1]) != (0.0, 0.0):
#             raise RuntimeError("%s non normalisé: self[0]=(%.5g,%.5g) != (0,0) "%(className(self),ba[0],ba[1]))
            msg.append("%s non normalise: self[0]=(%.5g,%.5g) != (0,0) "%(className(self),ba[0],ba[1]))
        corde, nba = self.corde, self.nba#self.computeCordeAndNBA()
        if abs(corde-1.0)>=self.prefs.EPS:# or nba != self.nba :
#             raise RuntimeError("%s non normalisé: corde=%.5g != 1.0 ou nba=%d != self.nba=%d"%(className(self),corde,nba, self.nba))
            msg.append("%s non normalise: corde=%.5g != 1.0 ou nba=%d != self.nba=%d"%(className(self),corde,nba, self.nba))
        else :
            msg.append('%s : normalite OK'%className(self))
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
        '''
        return super(ProfilNormalise,self)._update()

    def __normalise(self):
        if self.nb_normalisations>0 :
            raise RuntimeError(u"%s.__normalise() : profil deja normalise %d fois"%(className(self),self.nb_normalisations))
            rdebug(u"%s.__normalise() : profil deja normalise %d fois"%(className(self),self.nb_normalisations))
        else :
            res = super(ProfilNormalise, self).normalise()
            self.nb_normalisations += 1
            return res

    def toDump(self, format_='new'):
        u'''Lors d'un dump d'une base de profils, le profil self est sauvée sous cette forme'''
        dump = super(ProfilNormalise, self).toDump(format_)
#         dump['profparam'] = self.profparam.dump
        return dump

if __name__=="__main__":
    from testsprofilnormalise import testMain
    config.TEST_MODE = False
    testMain()
