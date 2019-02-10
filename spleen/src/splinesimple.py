#!/usr/local/bin/python2.7
#-*-coding: utf-8 -*-
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016-2017-2018 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
__updated__="2019-02-10"
'''
import sys, os, math
import cPickle
from array import array
from matplotlib import pyplot as plt
# from pprint import pprint
from numbers import Number
from matplotlib.widgets import CheckButtons
import config
# from brouillon4 import plot
# from splinecomposee import NSplineComposee
plt.rcParams["figure.figsize"] = (20,10)
from utilitaires import (rstack, eliminerPointsDoublesConsecutifs,
                         className, centreGravite, baryCentre, XY,
                         Path, segmentPlusProche, stack, debug, rdebug, dist,
                         hardScale, absCurv,dist2,rotate,courbure,symetrieAxe)
from utilitaires.lecteurs import pointsFrom, LecteurUniversel
from splineabstraite import computeSpline, NSplineAbstract
from numpy import (log, linspace, asarray, sqrt, arange, zeros, cos,
                   tan, pi, loadtxt, ndarray, argmin, where, savetxt, ones,
                   vstack, insert, delete, isnan, NaN, nan, NAN, empty,
                   empty_like, zeros_like, sign, inf)
from numpy.linalg import  norm
# from scipy.optimize import newton, minimize
from scipy.integrate.quadpack import quad
from scipy.optimize._minimize import minimize_scalar

class NSplineSimple(NSplineAbstract):
    u"""
    - Les attributs _cpoints, _dpoints, _epoints ne valent JAMAIS None,
        ils sont TOUJOURS :
        # soit inexistants
        # soit des ndarray de shape (N,2)
    """
    class Default(object):
        _nbpd = 1000
        _methode  = ('ius',1)
        _mode      = 'linear'
        _nbpe      = 30
        _cpoints  = zeros((0,2))
        gparent = None
        eps = 1.0e-5 # 10 microns si les distances sont en metres

    def __init__(self, **dump):# cpoints=None, name='spline', role='spline'):
        u"""
        * dump DOIT contenir les clés :
             - aucune
        * dump PEUT contenir les clés suivantes (sinon il y a des valeurs par défaut)
            - 'points' dont la valeur est un .ndarray, ou un QPolygon
                ou une liste de points, ou un filename ou encore None
            - 'precision' = nb de points pour le tracé de la spline
            - 'methode' : cf computeSpline() paramètre hyper important.
                détermine s'il s'agit d'une spline cubique, de degré autre,
                périodique, si les dérivées aux extrémités sont données, etc...
            - 'mode' : le mode de répartition des points en cas d'échantillonage
            - 'name' : str ou unicode pour affecter un nom
            - 'role' : str ou unicode pour affecter un role...
        """
        super(NSplineSimple,self).__init__(**dump)
#         self.setDefaultValues()
#         self.load(dump)
        if self.mode in ('rayon', 'courbure') :
            if self.methode in (('ius',1),('us',1)) :
                raise ValueError(u"Une spline lineaire ne peut pas etre echantillonnee avec le mode '%s'"%self.mode )
        self._update()

    def __len__(self):
        return len(self._cpoints)

    def __call__(self, T, diff=0):
        u"""
        # Paramètre :
            * T peut être
                - un réel entre 0 et 1
                - ou un .ndarray de n réels dans [0,1] (de shape (n,1)? ou (1,n)?)
                - ou une liste de n réels entre 0 et 1
            * diff est l'ordre de dérivation
        # Retourne :
        les 2 ndarrays de shape (n,1), zippés, contenant les points (sx(ti, diff), sy(ti, diff)) ou ti parcourt T
            * si diff=0, retourne les coordonnées du point x(t), y(t),
            * si diff=1, retourne les dérivées x'(t), y'(t)
            * si diff=2, retourne les dérivées secondes x"(t),y"(t)
        # Utilisation :
            >>> S = NSpline(....)
            >>> T = .linspace(0.0, 1.0, 11)
            >>> S(T)
            ... : les 11 points [S(0), S(0.1), ...,S(0.9), S(1)] de la spline
            N.B. : la fonction zip fait ceci:
            >>> z = zip(range(3),range(3))
            >>> z
            ... [(0, 0), (1, 1), (2, 2)]
        """
        if not hasattr(self,'_sx') :
            self.computeSpline()

        try :
            X, Y = self.sx(T, diff), self.sy(T, diff)
            try : return asarray(zip(X,Y))
            except TypeError : return asarray([X,Y])

        except TypeError as msg :#sx pas calculé
            debug('Spline inconstructible : %s'%str(msg))
            return zeros((0,2))

    def __setitem__(self,k,value):
        u"""
        mise à jour du polygone de contrôle et du polygone seulement
        Typiquement, on pourra ecrire
        >>> S = NSplineSimple(cpoints=[(1,1),(2,2),(3,3)])
        >>> S[1] = 20,30
        >>> print S[1]
        [20.0, 30.0]
        """
        if all(self._cpoints[k]==value) : #inutile de déclencher un _update()
            return
        self._cpoints[k] = value
        self._update()#indispensable pour mettre à jour self._xxx

    @property
    def info(self):
        infos=[
                u"<%s>"%className(self),
                u'%25s = '%u'name'+u'%s'%self.name,
                u'%25s = '%u'role'+u'%s'%self.role,
                u'%25s = '%u'nb pts controle'+u"%d"%len(self.cpoints),
                u'%25s = '%u'closed'+u'%s'%self.isClosed(),
                u'%25s = '%u'nb_splines, n_bech'+"%d, %d"%(self.nbspline, self.nbech),
                u'%25s = '%u'methode'+u"%s"%str(self.methode),
                u'%25s = '%u'degre'+u"%s"%str(self.degre),
                u'%25s = '%u'nb points discretisation'+"%s"%str(self.nbpd),
                u'%25s = '%u'mode echantillonage'+u"%s"%self.mode,
                u'%25s = '%u'nb pts echantillon'+u"%s"%self.nbpe,
#                 u'%25s = '%u'nb updates'+u"%s"%self.nbupdate,
                ]
        i = u'largeur, hauteur'
        try :
            i1 = u"%g, %g"%(self.width, self.height)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)
#
        i = u'position cg'
        try :
            i1 = u"%g, %g"%(self.centregravite[0], self.centregravite[1])
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)
#
        i = u'longueur'
        try :
            i1 = u"%g"%(self.longueur('r'))
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%25s = '%i+i1)

        return infos

    def setDefaultValues(self):
        """Valeurs par defaut:"""
        self._cpoints   = self.Default._cpoints
        self._methode   = self.Default._methode
        self._mode      = self.Default._mode
        self._nbpe      = self.Default._nbpe
        self._nbpd      = self.Default._nbpd
        self.role       = className(self)
        self.name       = className(self)
        self.nbspline   = 0# nb appels à computeSpline
        self.nbdpoint   = 0# nb calculs dpoints
        self.nbech      = 0# nb echantillonnages (epoints)
        self.nbupdate   = 0
        self.eps        = self.Default.eps#precision en mètres
        self.gparent    = None

    def load(self, dump):
        u"""Ce qui est commun à toutes les splines."""
#         super(NSplineSimple, self).load(dump)
        u"""
        Attention, la cle est 'points'  dans les anciens projets, et 'cpoints' dans les nouveaux.
        OUI=>On affecte _cpoints directement qui évite l'appel au setter de cpoints (qui appelle _update()).
        cpoints DOIT être liste, tuple, array ou ndarray, sans points doubles consécutifs
        """
        keys = dump.keys()
        if 'name'      in keys : self.name     = dump.pop('name')
        if 'gparent'   in keys : self.gparent  = dump.pop('gparent')
        if 'role'      in keys : self.role     = dump.pop('role')
        if 'methode'   in keys : self._methode = dump.pop('methode')#type spline
        if 'mode'      in keys : self._mode    = dump.pop('mode')#mode echantillonnage
        if 'nbpe'      in keys : self._nbpe    = dump.pop('nbpe')#nbp echantillonnage
        if 'nbpd'      in keys : self._nbpd    = dump.pop('nbpd')#nbp discretisation
        if 'precision' in keys : self._nbpd    = dump.pop('precision')#idem
#         debug(self._methode)
        for key in ('points', 'cpoints') :#la cle 'cpoints' est prioritaire, en cas
            if key in dump:
                cpoints = dump.pop(key)
                if isinstance(cpoints, (list, tuple, array)) :#types Python de base
                    #obligé de considerer ce cas, car pour le load d'un pkl cpoints est une liste
                    cpoints = asarray(cpoints)
                if not isinstance(cpoints, ndarray) :
                    msg = u'Les points de controle doivent etre de type ndarray((*,2)), et non pas %s.'%cpoints.__class__.__name__
                    msg = msg + u' \nUtilisez la fonction points=pointsFrom(points) pour transformer a peut pres n\'importe quoi en ndarray((n,2)).'
                    msg = msg + u' \nLes points de controle ne doivent pas avoir de points doubles consecutifs.'
                    msg = msg + u' \nUtilisez points = eliminerPointsDoublesConsecutifs(points) en cas de doute.'
                    raise TypeError(msg)
#                 debug(self.methode)
                if self.methode[1]=='periodic' and len(cpoints) > 0 :#and any(cpoints[0] != cpoints[-1]) :
                    self._cpoints = cpoints#pas de update pour le moment
                    self.close()#Ca appelle le _update si besoin
                self.cpoints = cpoints#Ca appelle le _update()

        try : del self._epoints
        except AttributeError : pass

    def toDump(self):
        u"""Ce qui est commun à toutes les splines"""
        return {
                'classename' : className(self),#besoin de savoir quel type de spline.
                'cpoints'    : self.cpoints.tolist(),
                'role'       : self.role,
                'name'       : self.name,
                'methode'    : self.methode,
                ############################
                'nbpd'       : self.nbpd,
                }

    def copy(self):
        u"""retourne une copie de self"""
        dump = self.toDump()
        return type(self)(**dump)

    @property
    def cpoints(self):
        u"""
        points de contrôle, sous forme ndarray((n,2)), None si non initialisé
        """
        return self._cpoints

    @cpoints.setter
    def cpoints(self, points):
        u"""points de contrôle, sous forme ndarray((n,2))
        - points est un ndarray((n,2)) et ne doit pas avoir de points doubles consécutifs."""
#         stack('points=%s'%points)
        if points is None or len(points)<1 :
            self._cpoints = zeros((0,2))
        else :
            self._cpoints = points
        self._update()

    @property
    def dpoints(self):
        u"""
        Les points discretises, sont recalculés
            - si _update() a été appelé (_dpoints a été supprimé) ou bien
            - si self.precision a changé
        Si on veut des dpoints aux abscisses T=(t1,...tn), on appelle directement
        X=self(T) plutôt que X=self.dpoints
        :return : ndarray((self.presision,2), les points de discrétisation
        """
        if not hasattr(self, '_dpoints') : #self._dpoints a été supprimé
            T = linspace(0.0, 1.0, self.nbpd)
            # si sx et sy n'existent pas, self(T) les recalcule
            self._dpoints = self(T)
            try : del self._dac
            except AttributeError : pass
        return self._dpoints

#     @property
#     def dac(self):
#         u"""
#         les abscisses curvilignes normalisees du POLYGONE dpoints
#         """
#         if not hasattr(self, '_dac') :
#             self._dac = absCurv(self.dpoints, True)
#         return self._dac

    @property
    def epoints(self):
        u"""points échantillonnés, sous forme ndarray((n,2))"""
        if not hasattr(self,'_epoints') :
            self._epoints = self(self.tech)
        return self._epoints

    @property
    def tech(self):
        u"""Les parametres T de l'echantillonnage"""
        if not hasattr(self, '_tech') :
            self.echantillonner()
        return self._tech

    @property
    def sx(self):
        u"""La spline en x, recalculée si _sx n'existe pas"""
        if not hasattr(self, '_sx') or self._sx is None:
            if len(self._cpoints)>=2 :
                _, self._sx, self._sy = computeSpline(self.cpoints, self.methode)
            else : self._sx, self._sy = None, None
        return self._sx

    @property
    def sy(self):
        u"""La spline en y, recalculée si _sy n'existe pas"""
        if not hasattr(self, '_sy') or self._sy is None:
            if len(self._cpoints)>=2 :
                _, self._sx, self._sy = computeSpline(self.cpoints, self.methode)
            else : self._sx, self._sy = None, None
        return self._sy

    @property
    def nbpe(self):return self._nbpe
    @nbpe.setter
    def nbpe(self,ne):
        if ne==self._nbpe :
            return
        else :
            self._nbpe = ne
            try : del self._tech
            except AttributeError : pass
            try : del self._epoints
            except AttributeError : pass

    @property
    def nbpd(self):
        """nb points discrétisation => self._dpoints"""
        return self._nbpd
    @nbpd.setter
    def nbpd(self,nd):
        if nd==self._nbpd :
            return
        else :
            self._nbpd = nd
            try : del self._dpoints
            except AttributeError : pass
            try : del self._dac
            except AttributeError : pass
    precision=nbpd

    @property
    def knots(self):
        u"""Ce sont les 'knots' au sens des splines sx et sy,
        i.e. les pseudo-abscisses curvilignes de self
        i.e. les T tels que self.sx(T[i]),self.sy(T[i]) = _cpoints[k]
        i.e. les abs. curv. des points du POLYGONE self._cpoints
        Ils sont stockés dans sx.XXX et sy.XXX"""
        try : return self.sx.x#si methode = (cubic','xx')
        except AttributeError as msg1:
#             debug(str(msg1))
            try :return self.sx._data[0]
            except AttributeError as msg2 :
#                 debug(str(msg2))
                try : return self.sx.get_knots()
                except AttributeError as msg3 :
#                     debug(str(msg3))
                    return zeros((0,))

    def _update(self):
        u"""
        Suppression de tous les attributs self._xxx volatiles, ils sont
            recalculés à la demande i.e. quand on appelle self.xxx
        """
        u'''
        Est appelé à chaque modification
        - (géométrique) d'un point de contrôle de la spline
        - ou bien du PolygonF de base
        - ou bien de methode spline (cubic, IUS, US,...), ou des dérivées aux extremites
        - ou bien de mode d'échantillonage
        ultra privée, ne pas l'appeler de l'extérieur.
                '''
#         try : del self._qcpolygon
#         except AttributeError : pass
#         try : del self._qepolygon
#         except AttributeError : pass
        try : del self._epoints
        except AttributeError : pass
        try : del self._dpoints
        except AttributeError : pass
        try : del self._height
        except AttributeError : pass
        try : del self._width
        except AttributeError : pass
        try : del self._longueur
        except AttributeError : pass
        try : del self._dac
        except AttributeError : pass
        try : del self._sx
        except AttributeError : pass
        try : del self._sy
        except AttributeError : pass

    def isClosed(self, eps=0.0):
        try : return dist2(self[0], self[-1]) <= eps*eps
        except IndexError : return None

    def close(self):
        """
        Fermeture de self._cpoints.
        Appelé automatiquement si self.methode[0] = 'periodic'
        :return : int,
            # 0 si rien touché,
            # 1 si ajustement mineur (<= self.eps) de self[-1]
            # 2 si ajout de point self[-1]=self[0]
        """
        if self.isClosed(0) :
            msg = u"    cpoints deja ferme"
        elif dist2(self[0],self[-1]) <= self.eps:
            #presque fermé, on n'ajoute pas de point
            p1, p2 = self[0], self[-1]
            self._cpoints[-1] = self._cpoints[0]
            msg = u"""    les points self[0]=%s et self[-1]=%s sont presque identiques.
            Ils sont à une distance d=%.3g
             => self[-1] est (legerement) modifie.
            """%(p1,p2, dist(p1,p2))
        else :
            #On ajoute un point
            p1, p2 = self[0], self[-1]
            self._cpoint.append[self[0]]
            msg = u"""    Les points self[0]=%s et self[-1]=%s sont distincts.
            Ils sont à une distance d=%.3g
             => le point self[0] a ete rajoute en fin de spline.
            """%(p1, p2, dist(p1,p2))
        debug(msg)

#     def echantillonner(self):
#         u"""
#         Méthode usuelle d'échantillonnage.
#         Utiliser cette méthode pour échantillonner avec les paramètres de self
#         et conserver le résultat dans self.epoints, self.tech.
#         :utilisation :
#             >>> S.mode='courbure'
#             >>> S.nbpe=27
#             >>> e = S.echantillonner()
#             >>>     #e et S.epoints contiennent les points echantillonnés
#             >>>     #les parametres des points d'echantillonnage sont dans S.tech
#             >>> T = S.echantillonner(True)
#             >>>     # T contient les parametres des points d'echantillonnage
#             >>>     # S.epoints contient les points echantillonnés
#         """
#         return self._echantillonner(self.nbpe, self.mode, 0, 1)

    def echantillonner(self, nbp=0, mode=None, ta=0, tb=1):
        u"""
        répartition (échantillonnage) de nbp points sur la spline self,
        entre les abscisses ta et tb, suivant le mode précisé par 'mode'
        *supprime self._epoints et modifie self._tech*
        :return : les abscisses T=ndarray(shape=(n,1)) (je crois!)
                des points échantillonnés

        :param mode : str ou unicode ou par defaut None
            - si mode==None => mode = self.mode
            - si mode is not None => self.mode = mode
            - si mode=='linear', les points d'échantillonnage sont régulièrement
                répartis tout au long de la spline
            - si mode=='rayon' ou 'courbure' la densité de points est
                approximativement proportionnelle à la courbure.
            - si mode=='cpoints' : retourne simplement les points de contrôle.
            - si mode=='telkel' : ta doit être un tableau des parametres t
                des points echantillonnés retourne self(ta).

        :param nbp : int, nombre de points d'échantillonnage.
            - Si mode='telkel' ou 'cpoints', ce paramètre est inutile.
            - Dans tous les autres cas il est indispensable, nbp>0

        :param ta : float dans [0,1] ou ndarray((n,1)) comme retourné par
            self.absCurv()
            - facultatif, ta=0 par defaut
            - l'abs. curv. du premier point d'échantillonnage

        :param tb : float dans [0,1] ou ndarray((n,1)) comme retourné par
            self.absCurv()
            - facultatif, tb=1 par defaut
            - l'abs. curv. du dernier point d'échantillonnage
        """

        if mode is None : mode = self._mode
        else : self._mode = mode

        if nbp==0 : nbp = self._nbpe
        else : self._nbpe = nbp
        #suppression _epoints qui sera recalculé à partir de _tech calculé ici
        try : del self._epoints
        except AttributeError : pass

        if mode == 'cos' :
            #Points serrés au debut et à la fin, lâches au milieu
            C = cos(linspace(pi, 2*pi, nbp))
            T = 0.5*(1+C)
            if (ta, tb) == (0,1) : self._tech = T
            else : self._tech =  ta + (tb-ta)*T

        elif mode == 'x3' :
            #Points serrés au milieu, lâches au debut et à la fin
            T = linspace(-1, 1, nbp)**3
            T = 0.5*(1+T)
            if (ta, tb) == (0,1) : self._tech =  T
            else : self._tech =  ta + (tb-ta)*T

        elif mode in ('linear', 'lineaire', 'lin') :
            self._tech =  linspace(ta, tb, nbp)#y compris ta et tb

        elif mode in ('rayon', 'courbure') :
            u"""On ne touche plus à RIEN !!!
            Calcul des points d'échantillonnage, de sorte que la densité de points
            soit localement proportionnelle à la courbure
            - la "COURBURE ABSOLUE" en un point t est ca(t) = sqrt(abs(self.courbure(t)))
            On la calcule en N points d'une discrétisation fine T, elle se trouve dans CA
            - La SINUOSITÉ sur un micro-intervalle [t, t+dt] est l'intégrale de t à t+dt de ca(s)ds
                On la calcule en chaque micro-intervalle de la discrétisation fine, elle est dans S
                La SINUOSITÉ TOTALE est la sinuosité sur [a,b], c'est la somme des sinuosités des micro-intervalles
            - On partage l'intervalle [ta,tb] en nbp-1 intervalles (de largeur variable)
                ta=te[0]<te[1]<te[2]<...<te[nbp-1]=tb;
            tels que sur chaque intervalle j, la sinuosité sj soit constante=s0=SINUOSITÉ TOTALE/(nbp-1)
            ou, ce qui revient au même, tq la sinuosité de ta à te[j] soit egale à Sj=j*s0
            """
#             rdebug(mode,self.mode)
            if self.methode[0] in ('us', 'ius') and self.methode[1] < 3 :
                raise ValueError,\
            u"""Impossible d'echantillonner, la courbure n'est pas definie pour une spline de degre %d< 3"""%self.methode[1]
            if nbp==1 :
                debug(nbp=nbp, ta=ta,tb=tb)
                stack()
                self._tech =  ta
#                 else : return self(ta)
            N = 100001
            T = linspace(ta,tb,N)#N-1 micro intervalles
            dt = T[1]-T[0]#micro intervalles de largeur dt
            CA = sqrt(abs(self.courbure(T)))
            S = (CA[1:]+CA[:-1])*(dt/2)#Les sinuosités des N micro-intervalles = dt*moyenne des courbures aux extrémités
            s0 = sum(S)/(nbp-1)#La sinuosité (constante) de chaque intervale Te[j],Te[j+1] défini ci dessous
            Te = ones(nbp)*ta# les te[j] cherchés
            Te[-1] = tb
#             debug(Te=Te)
            i, Sj = 0, 0.0#i=numéro du micro-intervalle
#             dbg = False
            for j in range(1, nbp) :#j=les numeros des points d'échantillonnage
                while Sj < j*s0 :#Sj=Sinuosité de ta à t[j]
                    try :
                        Sj += S[i]
                    except IndexError as msg :
                        msg = ['IndexError :', '[ta, tb]=[%.4g, %.4g]'%(ta,tb),
                               'nbp=%d, j=%d'%(nbp,j),
                               's0 = %.4g ; j*s0-Sj=%.4g > 0'%(s0,j*s0-Sj),
                               ]
                        if j*s0-Sj < s0/1000 :
                            msg = [u'\nPas grave %s'%self.name]
                        else :
                            msg = [u'\nAttention']
                            msg += [
                                   'i=%d'%i,
                                   'Te=%s'%(str(Te))]
                        debug('\n    '.join(msg))
                        break
                    i += 1#i=
#                 debug('Te[%d]=T[%d] : %.3g, %.3g'%(j,i,Te[j],T[i]))
                Te[j] = T[i]
            self._tech =  Te
        else :
            stack('mode echantillonnage inconnu ')
            rdebug(mode=mode, nbp=nbp, ta=ta, tb=tb)
            self._tech =  zeros((0,))
        return self._tech

    @property
    def methode(self):
        return self._methode
    @methode.setter
    def methode(self, newmethode):
        u"""
        Changement de spline (ius, cubic, periodic,...)
        """
        if len(self)==2 : #on fait une interpolation linéaire
#             debug('>>>>yeap')
            newmethode = 'ius', 1
        elif len(self) == 3 :#on fait une interpolation quadratique
            newmethode = 'ius', 2
        if newmethode[0] == self.methode[0] :
            self._methode = newmethode
            self._update()#tout est modifié
        else :
            dump = self.copy().toDump()
            dump['methode'] = newmethode
            self.load(dump)

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, newmode):
        u"""
        Changement mode d'echantillonnage
        """
        if newmode==self._mode :
            return
        else :
            self._mode = newmode
            try : del self._epoints
            except AttributeError : pass
            try : del self._tech
            except AttributeError : pass

    def computeSpline(self, methode=None):
        u"""
        Calcule la spline (sx, sy) considérée comme une courbe paramétrée sx(t), sy(t).
        sx et sy sont deux splines à une seule variable au sens scipy.
        """
        self.nbspline += 1
        if methode is None :
            methode=self.methode
        if methode[1] != self.methode[1] :
            debug(u'Attention, changement de parametres dans \'%s\': self.methode= %s => methode=%s'%(self.name,str(self.methode),str(methode)))
            self.methode = methode
        try :
#             debug('Avant computeSpline', self)
#             pprint(self.toDump())
            _, self._sx, self._sy = computeSpline(self.cpoints, methode)
#             rdebug('computeSpline:OK?')
#             debug('Apres computeSpline', self.name)
            try : #On regarde si sx et sy existent. Sinon, on remonte l'exception
                sx_05=self.sx(0.5)
                sy_05=self.sy(0.5)
#                 rdebug(sx_05=sx_05, sy_05=sy_05)
                if isnan(sx_05) or isnan(sy_05) :# in (nan,inf) or sy_05 in (nan,inf) :
                    raise ValueError('sx(0.5) ou sy(0.5) : not a number')
            except AttributeError as msg :
                rdebug(msg)
                raise msg
            except TypeError as msg :
                #sx=None, n'existe pas encore
                pass
        except ValueError as msg :
            # le message est `x` must be strictly increasing sequence.
            # il y a des points doubles certainement, on les vire
#             rdebug(self._knots)#existe pas
#             raise
#             rdebug(msg, u'"%s", Je tente d\'eliminer les points doubles consecutifs'%self.name)
#             debug(points_sales=self.cpoints.tolist())
            cpoints, avirer = eliminerPointsDoublesConsecutifs(self._cpoints, vires=True)
#             rdebug(avirer = avirer)
            if not avirer and len(cpoints)!=0:
                debug(BUG=u"Un update doit manquer..., len(cpoints)=%d"%len(cpoints))
#             debug(abscurv=absCurv(cpoints, False))
            if len(cpoints) < len(self._cpoints) :#sinon ca boucle
                try : role=self.role
                except : role='pas de role'
                rdebug(avirer=avirer, classename=className(self), name=self.name, role=role)
                self._cpoints = cpoints
                self.computeSpline(methode)

        except RuntimeWarning as msg :
#             rdebug(msg)
            raise(msg)
        except Exception as msg :
#             rdebug(type(msg).__name__, msg)
            raise(msg)

    def plot(self, figure=None, aspect={}, titre=None, more=[], texts=[], show=True,
             buttons=True, numbers=['3p']):
        """
        :param figure: une instance de matplotlib.figure.Figure
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
        :param texts : list, [(x,y, 'text', style=, fontsize=...), ...] texts en x,y
        :param numbers: list ou set ou tuple ['3c','12e','100d'] numéroter les points
            controle, echantillon, discretisation '3c' signifie que l'on numérote
            les cpoints par pas de 3, '100d' les dpoints par pas de 100 etc ...
        :param show: bool, True pour affichage de TOUS les plots en attente, False sinon
        :return figure: la figure passée en argument ou nouvelle figure si None

        """
    #     renderer = RendererBase()
        defaultaspect = {
                    'c':'ro',#control :        red
                    'd':'b-',#discretisation : blue
                    'e':'g.',#echantillon :    green
                    'r':'kv',#rupture:         black
                  }
        defaultaspect.update(aspect)
        if figure is None :#pas de figure => on la crée
            figure = plt.figure('plot(self)')
#         debug(figure=figure)

    #     if nbpd is None : nbpd = self.nbpd
    #     if nbpe is None : nbpe = self.nbpe
    #     if mode is None : mode = self.mode
        D = self.dpoints
        C = self.cpoints
        E = self.epoints
    #     debug(self)
    #     exit()
        axes = figure.get_axes()
        if axes : ax = axes[0]
        else : ax = figure.subplots()
    #    debug(axes=axes, ax=ax,gca=plt.gca())

    #     ax = figure.subplots()
#         debug(ax=className(ax),figure=className(figure))
        if titre : ax.set_title(titre)
#         ax.set_title('prout')
    #         titre = self.name+str(self.methode)

        fmtc, fmtd, fmte, fmtr = (defaultaspect['c'], defaultaspect['d'],
                                  defaultaspect['e'], defaultaspect['r'])
    #     debug('fmtc=%s, fmtd=%s, fmte=%s'%(fmtc, fmtd, fmte))
        if fmtc : control,     = ax.plot(C[:,0], C[:,1], fmtc, lw=1, label=u'Contrôle')
        if fmtd : spline,      = ax.plot(D[:,0], D[:,1], fmtd, lw=1, label=u'Spline')
        if fmte : echantillon, = ax.plot(E[:,0], E[:,1], fmte, lw=1, label=u'Échantillon')
#         if isinstance(self, NSplineComposee) :
#             R = asarray(self[self.ruptures])#des points
#     #         Rx, Ry = XY(R)
#             fmtr = defaultaspect['r']
#             if fmtr : ruptures,     = ax.plot(R[:,0], R[:,1], fmtr, lw=1, markersize=12, label=u'Rupture')

        for x, y, color, name in more:
            _, = ax.plot(x, y, color, label=name)
        figure.legend()

        for chars in numbers :
            if len(chars)>1 :
                step = eval(chars[:-1])
            else :
                step = 1#numerotation de un point sur step
            char = chars[-1]
            if char=='c' :
                for k, p in enumerate(C[::step]) :
                    ax.text(p[0], p[1], '%d'%(step*k))
            elif char=='e' :
                for k, p in enumerate(E[::step]) :
                    ax.text(p[0], p[1], '%d'%(step*k))

        for txt in texts :
            ax.text(*txt)
        ax.axis('equal')

        if buttons :
    #         figure.add_subplot('112')
    #         rax = plt.axes([0.05, 0.4, 0.1, 0.15])
            rax = figure.add_axes([0.05, 0.4, 0.1, 0.15])
#             debug(rax=className(rax))#,fsubplot=className(fsubplot))
            labels = [u'control', u'spline',u'echantillon']
            values = [True, True, True]
            draws = [control, spline, echantillon]

#             if isinstance(self, NSplineComposee) :
#                 labels.append(u'ruptures')
#                 values.append(True)
#                 draws.append(ruptures)

        #         for x, y, color, name in more:
        #             temp, = ax.plot(x, y, color)
        #             if not name : continue
        #             if buttons :
        #                 draws.append(temp)
        #                 labels.append(name)
        #                 values.append(True)
            figure.subplots_adjust(left=0.2)
            for item in draws :
                item.set_visible(True)
    #     if buttons :
            check = CheckButtons(rax, labels, values)

            def func(label):
                try :
                    draw = draws[labels.index(label)]
                    draw.set_visible(not draw.get_visible())
                except Exception as msg:
                    rdebug(u"label=%s marche po : %s"%(label,str(msg)))
                    pass

    #             if label == u'spline':
    #                 spline.set_visible(not spline.get_visible())
    #             elif label == u'control':
    #                 control.set_visible(not control.get_visible())
    #             elif label == u'echantillon':
    #                 echantillon.set_visible(not echantillon.get_visible())
    #             elif label == u'ruptures':
    #                 ruptures.set_visible(not ruptures.get_visible())
    #             else :
    #                 draw = draws[labels.index(label)]
    #                 draw.set_visible(not draw.get_visible())
                plt.draw()
            check.on_clicked(func)
    #         figure.subplots_adjust(left=0.2)
        axes = figure.get_axes()
#         debug(axes=axes)
        if show : plt.show()
        return figure

    def plotCourbure(self):
#         from matplotlib import pyplot as plt
        from matplotlib.widgets import CheckButtons
#         nbpd = self.precision
#         nbpe = self.nbpe
#         self.echantillonner()
        D = self.dpoints
        C = self.cpoints
        # E = self.epoints
        _, ax = plt.subplots()
        titre = self.name+' courbure'
        plt.title(titre)
        T = linspace(0,1,100)
        spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
#         echantillon, = ax.plot(E[:,0], E[:,1], 'bv', lw=1)
        control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        courbure = self.courbure(T)
        courbure += (1.0 + abs(min(courbure)))
        courbure = log(courbure)
        courbure /= max(abs(courbure))
        courbure, = ax.plot(T, courbure)
        buttons = ['points controle','courbure','spline',]
        values = [True, True, True, True]
        draws = [control, courbure, spline]
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        rax = plt.axes([0.05, 0.4, 0.1, 0.15])
        check = CheckButtons(rax, buttons, values)

        def func(label):
            if label == 'spline':
                spline.set_visible(not spline.get_visible())
            elif label == 'points controle':
                control.set_visible(not control.get_visible())
            elif label == 'courbure':
                courbure.set_visible(not courbure.get_visible())
            else :
                draw = draws[buttons.index(label)]
                draw.set_visible(not draw.get_visible())
            plt.draw()
        check.on_clicked(func)
        plt.show()
        return plt

    def longueur(self, p='r'):
        if p=='r':#longueur vraie
            if hasattr(self, '_longueur') :
                return self._longueur
            else :
                self._longueur = self.absCurv(1) if len(self)>1 else 0
                return self._longueur
        else :
            return super(NSplineSimple, self).longueur(p)

    @property
    def height(self):
        if not hasattr(self, '_height') :
            self._height = 0 if len(self)<=1 else\
                max(self.dpoints[:,1]) - min(self.dpoints[:,1])
        return self._height

    hauteur = height

    @property
    def width(self):
        if not hasattr(self, '_width') :
            self._width = 0 if len(self)<=1 else\
                max(self.dpoints[:,0]) - min(self.dpoints[:,0])
        return self._width

    largeur = width

    def boundingRect(self):
        dpts = self.dpoints
        xM, xm, yM, ym = dpts[:,0].max(), dpts[:,0].min(), dpts[:,1].max(), dpts[:,1].min()
        return xm,ym,xM,yM

    @property
    def gravitycenter(self):
        u"""Le vrai centre de gravité de la surface (plaque plane)
        délimitée par le polygone fermé."""
        return self[0] if len(self)==1 else centreGravite(self.dpoints)
    centregravite=gravitycenter

    @property
    def barycentre(self):
        u"""Le centre de gravité du nuage de points matériels cpoints."""
        return baryCentre(self.cpoints)
    @property
    def degre(self):
        m0,m1 = self.methode
        if m0 in ('cubic',) : return 3
        elif m0 in ('us','ius') : return m1
        else : raise NotImplementedError,u"TODO:degre de la spline (methode=%s)"%self.methode

    def courbure(self, T) :
        """
        Calcule et retourne la courbure en T.
        :param T: ndarray(n,).
        :return sc: ndarray(n,) la courbure aux points de T.

        :TODO :
            raffiner en autorisant le calcul pour les splines de degré 1,
            la dérivée première n'existe pas,
                -> retourner inf si point anguleux et 0 sinon.
            Analogue pour degré 2, la dérivée seconde n'existe pas,
                -> retourner ce qu'il faut
        """
        try :
            dx,  dy  = self.sx(T, 1), self.sy(T, 1)
        except ValueError as msg :
            rdebug(str(msg))
            sc = empty_like(T)
            sc[:] = nan
            return sc
        try :
            d2x, d2y = self.sx(T, 2), self.sy(T, 2)
        except ValueError as msg :
            rdebug(str(msg))
            sc = empty_like(T)
            sc[:] = nan
            return sc
        norm3_d2 = sqrt(dx**2+dy**2)**3
        v = dx*d2y-dy*d2x
        sc = v/norm3_d2
        iz = where(norm3_d2 < 1.0e-12)
        sc[iz] = sign(v)[iz]*inf
#         debug(sc=sc)
        return sc
        # si norm_d2=0, x"(t)=y"(t)=0, c'est une droite, courbure nulle
#         sc[where(norm3_d2 < 1.0e-12)] = 0.0

    def absCurv(self, T, witherr=False):
        u"""
        Calcule et retourne l'abscisse curviligne réelle des points self(T) sur
            la spline self.
        L'abscisse curviligne d'un point de paramètre t dans [0,1] est peu
            différente de t lorsqu'il y a beaucoup de points de contrôle.
        Un temps 't' de T est une abscisse curviligne le long du polygone cpoints
        le temps absCurv(t) est l'abscisse curviligne réelle du point self(t) le
            long de self.
        C'est l'intégrale de 0 à t de phi(s) = sqrt(self.sx(s)**2+self.sy(s)**2)
        L'intégration est assurée par scipy.integrate.quad()
        Si la spline self a trop de points de contrôle, ca rame et l'erreur est
        importante
        :param self: une NSplineAbstract
        :param T: les n paramètres t des points dont on veut l'abscisse curviligne.
            Ces n paramètres doivent être dans l'intervalle [0,1]
        :type T: au choix
            - un ndarray de shape (n,1) à valeurs réelles dans [0,1],
            - une liste de n valeurs réelles dans [0,1],
            - un tuple de n valeurs réelles dans [0,1],
            - un réel unique t dans [0,1]
        :return ac ou (ac, err):
            - si T est réel, ac et err sont réels, err est l'erreur estimée
            - si T est un ndarray((n,1)) ou ndarray((n,)) de n réels
                alors ac et err sont de même type
            (voir scipy.integrate)
        """
        if isinstance(T, Number) :
            #On integre sur les sous intervalles de self délimités par les knots
            phi = lambda s : sqrt(self.sx(s,1)**2+self.sy(s,1)**2) #la fonction à integrer
            bornes = [tk for tk in self.knots if tk<T]+[T]#bornes de sous intervalles
            intervalles = zip(bornes[:-1], bornes[1:])#les intervalles
            ac, err = 0, 0
            for (t1, t2) in intervalles :
                int_t1_t2, err12 = quad(phi, t1, t2)#integration intervalle [t1, t2]
                ac += int_t1_t2
                err = max(err,err12)
            return (ac, err) if witherr else ac

    #             return ac1+ac2, max(err1, err2)
        else :
            res = asarray([absCurv(self, t) for t in T])
    #         res = asarray([quad(lambda s : sqrt(self.sx(s,1)**2+self.sy(s,1)**2), 0.0, t) for t in T])
            return (res[:,0], res[:,1]) if witherr else res[:,0]

    def integraleCourbure(self, a=0.0, b=1.0, n=100):
        u"""Integrale de la valeur absolue de la courbure, caractérise la
        régularité de la courbe"""
        #     n DOIT être pair"""
        h = float(b-a)/n
        T = linspace(a, b, n+1)
        C = abs(self.courbure(T))
        A1 = C[0] + C[-1]
        A2 = 2*sum(C[i] for i in range(2,n) if i%2==0)
        A4 = 4*sum(C[i] for i in range(1,n) if i%2==1)
    #         debug (h, A1, A2, A4, (h/3)*(A1 + A2 + A4))
        return (h/3)*(A1 + A2 + A4)

    def symetriser(self, axe, centre=None):
        u'''
        Modif de la spline elle meme. self est modifié.
        C'est à dire
        >>> self[i,axe] = 2*centre -self[i,axe]
        si centre == None on prend l'isobarycentre[axe]
        '''
        if len(self)==0 : return
        if centre is None :
            centre = self.barycentre[0][axe]
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On fait tourner les dérivées du point p0, puis du dernier point.
            #(dérivées premiere ou seconde suivant methode)
            newderivees = []
            for k in range(2) :
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être symétrisée
                d = derivees[0]#ordre de derivation
                u = asarray(derivees[1:], dtype=float).reshape((1,2))
#                 debug('Avant sym:',u=u, axe=axe, centre=centre)
                u = symetrieAxe(u,axe,0.0)[0]
#                 debug('Apres sym:',u=u, axe=axe, centre=centre)
                newderivees.append((d, u[0], u[1]))

            self._methode = (self.methode[0], tuple(newderivees))
        self.cpoints = symetrieAxe(self.cpoints,axe,centre)
#         self._update()c'est fait par le setter de cpoints
#         self.qcpolygon=points#qpolygonFrom(points)
#         self.update()C'est fait dans le setter qcpolygon()

    def hardScale(self, scale=(1.0,1.0), centre=None):
        u'''
        Modif de la spline elle même. self est modifié.
        Mise à l'échelle d'un scale=(fx,fy), centrée sur centre. (en fait, une homothétie)
        C'est à dire

            - self[i][0] = centre[0] + scale[0]*(self[i]-centre)[0]
            - self[i][1] = centre[1] + scale[1]*(self[i]-centre)[1]

        Le cas échéant, les valeurs des dérivées aux extrémités de la spline sont
        mises à l'échelle, de sorte que self transformée soit réellement homothétique de self.
        En modifiant self.cpoints, le _update est appelé, tous les attributs de self
        (_epoints, _absCurv, _dpoints, etc...) sont mis à jour.
        :param centre: (float, float) les coordonées (x0, y0) du centre d'homothétie.
            si centre==None on prend l'isobarycentre des points de contrôle.
        :param scale: (float, float)
        :return: None

        '''
        if len(self) == 0 : return
        if centre is None :
            centre = self.barycentre
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On scale les dérivées du point p0, puis du dernier point.
            #(dérivées premiere ou seconde suivant methode)
            newderivees = []
            for k in range(2) :
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être mise a l'echelle
                d = derivees[0]#ordre de derivation
                #u = (x'(t), y'(t)) ou (x"(t),y"(t)) en t=0 ou 1
                u = asarray(derivees[1:], dtype=float)
                u[0] *= scale[0]
                u[1] *= scale[1]
                newderivees.append((d, u[0], u[1]))

            self._methode = (self.methode[0], tuple(newderivees))
#         debug(Avant=self.cpoints)
        self.cpoints = hardScale(self.cpoints, scale, centre)#ca fait le _update
#         debug(APRES=self.cpoints)
        return

    def translate(self,vecteur):
        self.cpoints = self.cpoints+vecteur
#         self._update()#C'est fait dans la property cpoints
    hardMove = translate

    def appendPoint(self,pos):
        i = self.cpoints.shape[0]#l'indice du dernier point
        if len(self) and dist2(pos, self.cpoints[-1]) == 0 :
            raise ValueError('Point double, Impossible d\'ajouter ce point:%s'%(pos))
            # return #eviter points doubles, les abs curv douvent être strictement croissantes.
        else :
            self.cpoints = vstack((self.cpoints, pos))
            self._update()
            return i

    def index(self, t):
        """Localise t dans self.abscurv:retourne les deux entiers
        k1,k2  tels que k1 < t <= k2"""
        if t<0 or t>1 : return nan, nan
        T = self.knots
        k = 0
        while T[k]<t : k+=1
        if T[k] == t : return k,k
        else : return k-1, k

    def insertPoint(self, pos, i=None):
        u"""
        segment_i = [self[i-1], self[i]]
        Si i = None :
            Calcule dans quel segment il est raisonable d'insérer le point pos
        si i != None,
            insertion du point pos entre i-1 et i (de sorte qu'il deviennet le point i)
        """
        # if isinstance(pos, QPointF) :
        # pos = pos.x(), pos.y()
        cpoints = self.cpoints
        if i is None :
            im1, _ = segmentPlusProche(cpoints, pos)
            i = im1 + 1
            #Les segments sont numerotes de 0 à len(cpoints)-2
            #le segment i est délimité par les points i-1 et i
#             debug('Calcul segment le plus proche de %s'%pos)
#             if im1 is None :#Ca n'arrive jamais, cf segmentPlusProche()
#                 debug('Je ne sais pas ou (dans quel segment) inserer ce point %s'%str(pos))
#                 return
        elif i == 0 :
            #si i est le premier point, on insere dans le premier segment
            im1, i = 0, 1
        else :
            im1 = i-1
        if dist2(pos, cpoints[i]) == 0 or dist2(pos, cpoints[im1]) == 0 :
#             debug('segment le plus proche',i=i,H=H,pos=pos)
            raise ValueError('Point double, impossible d\'inserer ce point:%s en position %d'%(str(pos),i))
            # return #eviter points doubles, car les abs curv doivent être strictement croissantes.
        else:
            cpoints = insert(cpoints, i, (pos,), axis=0)
#             debug(cpoints=cpoints)
            self.cpoints = cpoints
#             self._update()#deja fait par le cpoint.setter
            return i#La position du point inséré

    def removePoint(self, pnum, update=True):
        u"""
        - Suppression du point pnum de qcpolygon
        """
        point = self[pnum]
        if update : #le setter de cpoints fait un update
            self.cpoints = delete(self.cpoints, pnum, 0)
        else :
            #le setter de cpoints fait un update, parfois c'est indésirable
            #dans ce cas, on intervient directement sur _cpoints
            self._cpoints = delete(self._cpoints, pnum, 0)
#         self._update()#c'est fait par le setter de cpoints
        return point

    def hardRotate(self, alfa, centre=None, unit='degres'):
        u'''
        Rotation de self, alfa est en degres par défaut
        si les dérivées aux extremités sont précisées, il faut les tourner elles aussi
        methode est de la forme : ('cubic',(D0, Dn))
        avec D0 = (d, a, b) ou (a,b) est le vecteur dérivée d-ieme au point 0
        et Dn idem pour le dernier point.
        '''
#         debug('avant:', alfa=alfa)
        if unit == 'degres' :
            alfa = math.radians(alfa)
#         debug('apres:',alfa=alfa)
        if self.methode[0] == 'cubic' and isinstance(self.methode[1],(list, tuple)) :
            #On fait tourner les dérivées premiere et seconde du point 0, puis du dernier point.
            #Si les data sont x' et y'' par exemple, il faut faire tourner x',y' puis x'', y''
            newderivees = []
            for k in range(2) : #on tourne les deux vecteurs dérivés aux extrémités
                derivees = self.methode[1][k]# Point 0 ou n
                #la derivée doit être tournee de alfa
                d = derivees[0]#ordre de derivation
                u = derivees[1:], #u=(x', y') ou bien (x'', y'') c'est lui qui tourne
                #debug('Avant rotate:x',d=d, u=u)
                u = rotate(u, alfa, (0,0))[0]
                #idem pour la derivée k-ième
#                 #on ne garde que u[0] et v[1]
                newderivees.append((d, u[0], u[1]))

            self._methode = [self.methode[0], newderivees]

#         debug('Avant rotate', cpoints=self.cpoints, barycentre=self.barycentre, gravity=self.gravitycenter)
        if centre is None :
            centre = self.barycentre
        points = rotate(self.cpoints, alfa, centre)
#         debug('apres rotate', cpoints=self.cpoints, barycentre=self.barycentre, gravity=self.gravitycenter)
        self.cpoints = points
#         self._update()c'est fait dans le setter cpoints()

#     def elaguerNeMarcheQuePourUnTrados(self, eps=0.5, replace=False,debog=False):
#         u"""
#         On cherche une spline s1 avec un minimum de points de contrôle
#         et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)
#
#         La distance(self.cpoints, s1) est le plus grand écart entre
#         la spline calculée et les points de contrôle de la spline self.
#         autrement dit le max des distances d(self.cpoints[k], s1), k=0,1,...
#         où d(P, s1) est le min de la fonction t -> norme(P-s1(t)).
#         Voir la methode self.distanceTo().
#         On discrétise finement self (une seule fois) -> tableau de points D0
#         On initialise s1 partant des quelques points de contrôle dont
#         self.cpoints[0] et self.cpoints[-1].
#         On discrétise finement s1 (à chaque ajout de point) -> tableau de points D1
#         puis on rajoute des points de contrôle à s1 jusqu'a obtention de la précision désirée (eps)
#         Pour rajouter un point, on repère la distance point à point de D0 et D1,
#         disons qu'elle est au i-eme point et on insère dans s1 un point de contrôle
#         que l'on positionne exactement en D0[i].
#
#         [ne marche pas =>] Quand la précision désirée est atteinte, on met
#         un petit coup d'optimisation pour améliorer la position des points
#         de contrôle de s1.
#
#         :param eps: float, la spline resultante est à une distance eps au maximum de self.
#         :param replace: bool, si True, la spline élaguée remplace self.
#         :return: s1, pm, (n0,n1)
#
#             - s1 est la spline NSplineSimple élaguée
#             - pm = self.pourMille(d) est la précisison en valeur relative :
#                 la spline s1 ne s'éloigne pas à plus de pm ‰ de self.
#             - n0, n1 = nb de points de contrôle avant et après élagage
#
#         """
#         if len(self) < 10 :
#             rdebug(u"Élagage inutile, %d points de contrôle seulement"%len(self))
#             return self, self.pourMille(0),(len(self),len(self))
#         s0 = self
#         nd = 5000
# #         Td = linspace(0,1,nd)
# #         T0 = self.knots
#         c0 = s0.cpoints.copy()
# #         d0 = asarray(s0(Td))
#         t, m = self.methode
#         n = len(self)
#         if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
#             init = (c0[0], self(0.33), self(0.66), c0[0])
#             init = (c0[0], self(0.25), self(0.5), self(0.75), c0[0])
# #             debug('polyligne (4points) ferme')
#             seen = [0, n/4, n/2, 3*n/4, n-1]
#             init = self[seen]
#             c1 = pointsFrom(init)
#         elif dist2(c0[0], c0[-1]) < 1.0e-6 :
# #             init = (c0[0], self(0.5), c0[0])
#             seen = [0, n/3, 2*n/3, n-1]
#             init = self[seen]
# #             init = (c0[0], self(0.33), self(0.66), c0[0])
#             c1 = pointsFrom(init)
# #             debug('Initial : polyligne ferme =%s'%str(c1))
#         else:
# #             debug('polyligne non ferme')
#             seen = [0,n-1]
#             init = self[seen]
#             c1 = pointsFrom(init)
# #         debug(initial=c1)
#         s1 = NSplineSimple(cpoints=c1,
#                            methode=s0.methode,
#                            nbpd=s0.nbpd,
#                            name='%s-elaguee'%s0.name,
#                            mode=s0.mode)
# #         d1 = asarray(s1(Td))
#         if debog :
#             X0, Y0 = XY(c0)
#             more = [(X0,Y0, 'r.-','c0'),]
#         T, D, _ = s1.distanceTo(c0)
# #         D[0] = D[-1] = 0.0
# #         for k in range(len(c0)) :
#         while len(seen)<len(c0):
#             #indices du-des point-s ou la distance entre les deux splines est max
#             debug(seen=seen)
#             D[seen] = 0.0#On ne veut pas revisiter 2 fois le même point (=>point double)
#             d = max(D)
#             imaxd = (D == d).nonzero()
#             dm = self.pourMille(sqrt(d))
#             idx = imaxd[0][0]
#             seen = sorted(list(set(seen+[idx])))
#             pos0 = c0[idx]#position(x,y) du point à inserer dans s1
#             #On cherche maintenant en quelle position (k) il faut l'insérer
#             #
#             K1 = s1.knots.tolist()#Les t des noeuds de s1
#             t = T[idx]#le param t (sur s1) de la distance max
#             # La position de t dans K1
#             i = 0
#             while K1[i] <= t : i += 1
#             #i est l'indice d'insertion dans s1
#             if len(seen)<7 :
#                 debug('    t=%.3g'%t, knots=s1.knots.tolist(), )
#                 print '    idx=',idx, ' ; insertion dans s1=',i#, ' ; dist=',D[idx],' ; maxD=',d
#                 print '    dist=maxD ?',D[idx]==d
#                 print '    *** distances de s1 (%d points) à c0'%len(s1), D.tolist()
#             if debog :
#                 s1.plot(plt,
#                         titre=u'Spline élaguée : \nseen=%s, \ndist=%.1g ‰'%(str(seen),d),
#                         more=more,show=True)
#             try :
#                 s1.insertPoint(pos0, i)
#             except ValueError as msg :#Impossible d'inserer le point (point double ?)
#                 debug(u'Ca devrait pas arriver, je tente autre chose', msg, pos0=pos0)
#                 Td = linspace(0,1,nd)
#                 d0 = asarray(s0(Td))
#                 d1 = asarray(s1(Td))
#                 ad = norm(d0-d1, 2, axis=1)#[1:-1]#ecart point à point des deux splines discrétisées
#                 mad = (ad == max(ad)).nonzero()#indice des points ou la distance entre les deux splines est max
#                 idx = mad[0]
#                 t = Td[idx][0]
#                 pos0 = d0[idx][0]
#                 try :
#                     s1.insertPoint(pos0)
#                 except ValueError as msg :
#                     rdebug(msg, pos0=pos0)
#                     rdebug(u'Precision non atteinte, iteration %d : dist = %.2e '%(len(seen),d))
#                     break
# #                     c1 = s1.cpoints.copy()
#
# #             d = [distancePointSpline(c0[i], s1, t0=t).fun for i, t in enumerate(T0)]
#             D = s1.distanceTo(c0)[1]
#             D[0] = D[-1] = 0.0
#             dm = self.pourMille(sqrt(max(D)))
#             debug(u'dist-%d = %.2g‰ '%(len(seen),dm))
# #             d1 = asarray(s1(Td))
#             c1 = s1.cpoints.copy()
#             if dm<eps : break
#         if len(s1) == len(self) :#même nb de points : on garde la spline initiale.
#             s1.cpoints = self.cpoints.copy()
#             n0 = len(self)
#             n1 = len(s1)
#             return s1, d,(n0,n1)
#         #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de s1
#         #ca ne marche pas du tout !!!
#         n0 = len(self)
#         n1 = len(s1)
#         debug('Apres ELAGAGE : dist = %.2g mm/m ; %d => %d '%(d,n0,n1))
#
#         if replace :
#             self.cpoints = s1.cpoints
#         return s1, d,(n0,n1)
#
#     def elaguerBeurk(self, eps=0.5, replace=False, debog=False):
#         u"""
#         On cherche une spline s1 avec un minimum de points de contrôle
#         et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)
#
#         La distance(self.cpoints, s1) est le plus grand écart entre
#         la spline calculée et les points de contrôle de la spline self.
#         autrement dit le max des distances d(self.cpoints[k], s1), k=0,1,...
#         où d(P, s1) est le min de la fonction t -> norme(P-s1(t)).
#         Voir la methode self.distanceTo().
#         On discrétise finement self (une seule fois) -> tableau de points D0
#         On initialise s1 partant des quelques points de contrôle dont
#         self.cpoints[0] et self.cpoints[-1].
#         On discrétise finement s1 (à chaque ajout de point) -> tableau de points D1
#         puis on rajoute des points de contrôle à s1 jusqu'a obtention de la précision désirée (eps)
#         Pour rajouter un point, on repère la distance point à point de D0 et D1,
#         disons qu'elle est au i-eme point et on insère dans s1 un point de contrôle
#         que l'on positionne exactement en D0[i].
#
#         [ne marche pas =>] Quand la précision désirée est atteinte, on met
#         un petit coup d'optimisation pour améliorer la position des points
#         de contrôle de s1.
#
#         :param eps: float, la spline resultante est à une distance eps au maximum de self.
#         :param replace: bool, si True, la spline élaguée remplace self.
#         :return: s1, pm, (n0,n1)
#
#             - s1 est la spline NSplineSimple élaguée
#             - pm = self.pourMille(d) est la précisison en valeur relative :
#                 la spline s1 ne s'éloigne pas à plus de pm ‰ de self.
#             - n0, n1 = nb de points de contrôle avant et après élagage
#
#         """
#         if len(self) < 10 :
#             rdebug(u"Élagage inutile, %d points de contrôle seulement"%len(self))
#             return self, self.pourMille(0),(len(self),len(self))
#         self = self
#         nd = 100
# #         Td = linspace(0,1,nd)
# #         T0 = self.knots
#         c0 = self.cpoints.copy()
# #         d0 = asarray(self(Td))
#         t, m = self.methode
#         n = len(self)
#         if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
#             seen = [0, n/4, n/2, 3*n/4, n-1]
#             init = self[seen]
#             c1 = pointsFrom(init)
#         elif dist2(c0[0], c0[-1]) < 1.0e-6 :
#             seen = [0, n/3, 2*n/3, n-1]
#             init = self[seen]
#             c1 = pointsFrom(init)
#         else:
#             seen = [0,n-1]
#             init = self[seen]
#             c1 = pointsFrom(init)
#         s1 = NSplineSimple(cpoints=c1,
#                            methode=self.methode,
#                            nbpd=self.nbpd,
#                            name='%s-elaguee'%self.name,
#                            mode=self.mode)
# #         d1 = asarray(s1(Td))
#         if debog :
#             X0, Y0 = XY(c0)
#             more = [(X0, Y0, 'k.-','c0'),]
#             texts = [(X0[0],Y0[0]-0.05,u'%d'%0),(X0[-1],Y0[-1]-0.05,u'%d'%(len(X0)-1))]
#             texts+= [(c1[0,0],c1[0,1],u's1:%d'%0),(c1[-1,0],c1[-1,1],u's1:%d'%(len(c1)-1))]
#         T1, D01, _ = s1.distanceTo(c0,discret=nd)
#         while len(seen)<len(c0):
#             #indices du-des point-s ou la distance entre les deux splines est max
#             debug(seen=seen)
#             D01[seen] = 0.0#On ne veut pas remettre 2 fois le même point (=>point double)
#             d = max(D01)#le point de c0 le plus eloigné de s1
#             imaxd = (D01 == d).nonzero()
#             dm = self.pourMille(sqrt(d))
#             idx0 = imaxd[0][0]
#             seen = sorted(list(set(seen+[idx0])))
#             pos0 = c0[idx0]#position(x,y) du point à inserer dans s1
#             pj01 = s1(T1[idx0])#pj de pos0 sur s1
# #             print 'c0 = asarray(%s)'%c0.tolist()
#             print u's1 = NSplineSimple(**%s)'%s1.toDump()
#             print u'seen = %s       #deja visités'%seen
#             print u'idx0 = %d       #numero dans c0'%idx0
#             print u't1 = %.3g       #temps dans s1'%T1[idx0]
#             print u'pos0 = x, y = (%.3g,%.3g) #position c0[idx0]'%(pos0[0],pos0[1])
#             print u'pj0 = xj, yj = (%.3g,%.3g) #projeté sur s1'%(pj01[0],pj01[1])
#             #On cherche maintenant en quelle position (k) il faut l'insérer
#             #
#             K1 = s1.knots.tolist()#Les t des noeuds de s1
#             t = T1[idx0]#le param t (sur s1) de la distance max
#             # La position de t dans K1
#             print 'knots = %s'%K1
#             print 't = %g'%t
#             i = 0
#             while K1[i] < t : i += 1
# #             i -= 1
#             #i est l'indice d'insertion dans s1
#             if len(seen)<7 :
#                 debug('    t=%.3g'%t, knots=s1.knots.tolist(), )
#                 print '    idx=',idx0, ' ; insertion dans s1=',i#, ' ; dist=',D[idx],' ; maxD=',d
#                 print '    dist=maxD ?',D01[idx0]==d
#                 print '    *** les distances [d(c0[i],s1) i=0...len(c0) ](len(s1) = %d points) à c0'%len(s1), D01.tolist()
# #                 plt.plot(D)
# #                 plt.title(u"distances de c0 à s1 ; min=%.3g point num%d"%(d, idx))
# #                 plt.show()
#             if debog :
#                 pj = self(t)
#                 more += [((pos0[0],pj[0]),(pos0[1],pj[1]),'go-',u'à ajouter en position %d'%(i))]
#                 texts += [(pos0[0],pos0[1],u'%d'%(i))]
# #                 debug(more=more)
#                 s1.plot(plt,
#                         titre=u'Spline élaguée : \nseen=%s, \ndist=%.1g ‰'%(str(seen),d),
#                         more=more,texts=texts, show=True)
#                 more = more[:-1]
# #                 texts=texts[:-2]
#                 c1 = s1.cpoints
#                 texts = [(c1[0,0],c1[0,1],u's1:%d'%0),(c1[-1,0],c1[-1,1],u's1:%d'%(len(s1)))]
#
#             try :
#                 s1.insertPoint(pos0, i)
#             except ValueError as msg :#Impossible d'inserer le point (point double ?)
#                 debug(u'Ca devrait pas arriver, je tente autre chose', msg, pos0=pos0)
#                 raise msg
#                 Td = linspace(0,1,nd)
#                 d0 = asarray(self(Td))
#                 d1 = asarray(s1(Td))
#                 ad = norm(d0-d1, 2, axis=1)#[1:-1]#ecart point à point des deux splines discrétisées
#                 mad = (ad == max(ad)).nonzero()#indice des points ou la distance entre les deux splines est max
#                 idx = mad[0]
#                 t = Td[idx][0]
#                 pos0 = d0[idx][0]
#                 try :
#                     s1.insertPoint(pos0)
#                 except ValueError as msg :
#                     rdebug(msg, pos0=pos0)
#                     rdebug(u'Precision non atteinte, iteration %d : dist = %.2e '%(len(seen),d))
#                     break
# #                     c1 = s1.cpoints.copy()
#
# #             d = [distancePointSpline(c0[i], s1, t0=t).fun for i, t in enumerate(T0)]
#             D01 = s1.distanceTo(c0, discret=nd)[1]
# #             D01[0] = D01[-1] = 0.0
#             dm = self.pourMille(sqrt(max(D01)))
#             debug(u'dist-%d = %.2g‰ '%(len(seen),dm))
# #             d1 = asarray(s1(Td))
#             c1 = s1.cpoints.copy()
#             if dm<eps : break
#         if len(s1) == len(self) :#même nb de points : on garde la spline initiale.
#             s1.cpoints = self.cpoints.copy()
#             n0 = len(self)
#             n1 = len(s1)
#             return s1, d,(n0,n1)
#         #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de s1
#         #ca ne marche pas du tout !!!
#         n0 = len(self)
#         n1 = len(s1)
#         debug('Apres ELAGAGE : dist = %.2g mm/m ; %d => %d '%(d,n0,n1))
#
#         if replace :
#             self.cpoints = s1.cpoints
#         return s1, d,(n0,n1)

    def pourMille(self, longueur):
        u"""longueur convertie en ‰ de la longueur de self"""
        return 1000*longueur/self.longueur()

    def projection(self, p, discret=0, t0=0.0, t1=1.0):
        u"""
        Calcule le projeté du point p sur la portion S01 de la spline self
            S01 = {self(t), t dans I=[t0,t1]},
            i.e. le min de la fonction phi : s-> dist(p,self(s)), t dans I.
        Pour cela,
        - si discret=0, on appelle scipy.minimize_scalar() sur l'intervalle I.
            Ce faisant, on risque fort de ne pas tomber sur le bon min,
            lorsqu'il y a plusieurs minima locaux.
            Typiquement, si self est un profil (fermé) et p=centre de gravité de
            self, la courbe t -> dist(p,self(t)) présente 2 minima locaux.
            On doit donc localiser le min absolu en discrétisant la spline.
            C'est ce qui est fait si discret>0
        - si discret>0, les valeurs de t0 et t1 sont ignorées.
            On localise le min absolu en discrétisant la spline
            => self.dpoints = self(T) ou T contient 1+discret pas de temps
            linéairement répartis, i.e. discret intervalles de largeur dt.
            # on cherche tm=le min (discret) des distances dist(p, self.dpoints)
            # on raffine ensuite sur l'intervalle [tm-dt, tm+dt]
        :param p : tuple, list ou ndarray((2,)), les coord. du point p.
        :param discret : int, le nombre de points de discrétisation
        :param t0,t1: float, float, l'intervalle de recherche du min.
            on doit avoir 0.0 <= t0 < t1 <= 1.0
        :return tt, dd, nev, msg : float, float, int, str
            - tt dans [t0,t1] est le temps sur self du projeté p' = self(tt)
            - dd est la distance euclidienne dist(p,p')
            - nev le nombre d'évaluations de la fonction phi(s)
            - msg un message fournit par minimize_scalar()
        """
        debog=False
        if discret>0 :
            self.nbpd = discret#_dpoints est effacé et recalculé à la demande
            eps = 1.0/(discret-1)
            #les distances de p aux points de dpoints
            D = norm(self.dpoints-p,axis=1)
            idx = argmin(D)
            d = D[idx]
    #         debug(d=d,idx=idx)
            T = linspace(0,1,discret)
            #t est le temps du point dpoints[idx] sur le polygone self.dpoints
            t = T[idx]
            #raffinement
            t0, t1 = max([0,t-eps]), min([1,t+eps])
            try :
                tt, dd, nev, msg = self.projection(p, discret=0, t0=t0, t1=t1)
            except ValueError as msg :
                tt, dd, nev, msg = t, d, nan, u'Echec optimisation'
            if debog :
                debug(titre=u'==> p=c0[%d]=, discret = %-6d'%(debog, discret))
                print u' D%d = asarray(%s)'%(debog, D.tolist())
                print u' T%d = asarray(%s)'%(debog, T.tolist())
                print u' idx = %d  indice du point le plus éloigné de se dans c0#'%idx
                print u" Avant raffinement : t=%-6.3g, d=%-6.3g, (t0,t1)=(%.3g,%.3g)"%(t, d, t0, t1)
                print u" Après raffinement : t=%-6.3g, d=%-6.3g, ev=%2d, msg='%s'"%(tt, dd, nev, msg)
    #         return t, d, 0, 'pas de raffinement'
            return tt, dd, nev, msg
        else :#discret=0
            a, b = p[0], p[1]
            res = minimize_scalar(lambda t: (a-self.sx(t))**2 + (b-self.sy(t))**2,
                              bounds=(t0, t1),
                              method='bounded',
                              options={'xatol':1.0e-9})
            return res.x, sqrt(res.fun), res.nfev, res.message

    def projectionObj(self, obj, discret=0):
        u"""
        Calcule la projection de l'objet obj sur la spline self par appels répétés à
            projection(self,obj[k],...)
        :param obj: NSplineSimple ou ndarray((n,2),dtype=float) ou bien
                    un seul point (float, float) a projeter sur self.
            # si obj et une NSplineSimple on projete les points de obj.cpoints sur self
            # si obj est un ndarray((n,2)) on projete les points de obj sur self
        :param discret: int, un nombre de points de discrétisation pour le calcul
            des projections.
        :return tprojs, dists, nevs, msgs: cf doc scipy.optimize.OptimizeResult
            - tprojs est la liste des parametres t des (c)points de obj sur self
            - dists est la liste des distances à self des projetés
            - nevs est la liste des nb d'évaluation de phi
            - msgs sont des messages de convergence donnés par scipy.minimize_scalar
        La fonction phi(t) (distance de p à self, lorsque p est un point (x,y))
        peut admettre plusieurs minima locaux.
        Il faut les trouver tous et les comparer entre eux pour obtenir le vrai
        minimum.
        C'est ce que fait l'appel à projete()
        """
        if isinstance(obj, (tuple, list, ndarray)) and len(obj)==2 :
            return self.projection(p=obj, discret=discret)
        else :
            tprojs = zeros((len(obj),))
            dists  = zeros((len(obj),))
            nevs   = zeros((len(obj),),dtype=int)
            msgs   = zeros((len(obj),), dtype=str)
            if isinstance(obj, NSplineSimple) :
                P = obj.cpoints
            elif isinstance(obj, ndarray) and len(obj[0])==2 :
                P = obj
            for k, p in enumerate(P) :
    #             debog = False if k != NUM else k
                tprojs[k], dists[k], nevs[k], msgs = self.projection(p, discret)
            return tprojs, dists, nevs, msgs

    def elaguer(self, eps=0.5, replace=False, debog=False):
        u"""
        On cherche une spline se avec un minimum de points de contrôle
        et qui soit à une distance de self.cpoints < eps (en ‰ de la longueur de self)

        La distance(self.cpoints, se) est le plus grand écart entre
        la spline calculée et les points de contrôle de la spline self.
        autrement dit le max des distances d(self.cpoints[k], se), k=0,1,...
        où d(P, se) est le min de la fonction t -> norme(P-se(t)).
        Voir la methode self.distanceSplineToObj().
        On discrétise finement self (une seule fois) -> tableau de points D0
        On initialise se partant des quelques points de contrôle dont
        self.cpoints[0] et self.cpoints[-1].
        On discrétise finement se (à chaque ajout de point) -> tableau de points D1
        puis on rajoute des points de contrôle à se jusqu'a obtention de la précision désirée (eps)
        Pour rajouter un point, on repère la distance point à point de D0 et D1,
        disons qu'elle est au i-eme point et on insère dans se un point de contrôle
        que l'on positionne exactement en D0[i].

        [ne marche pas =>] Quand la précision désirée est atteinte, on met
        un petit coup d'optimisation pour améliorer la position des points
        de contrôle de se.

        :param eps: float, la spline resultante est à une distance eps au maximum de self.
        :param replace: bool, si True, la spline élaguée remplace self.
        :return: se, pm, (n0,n1)

            - se est la spline NSplineSimple élaguée
            - pm = self.pourMille(d) est la précisison en valeur relative :
                la spline se ne s'éloigne pas à plus de pm ‰ de self.
            - n0, n1 = nb de points de contrôle avant et après élagage

        """
        histoire = []
        if len(self) < 10 :
            rdebug(u"Élagage inutile, %d points de contrôle seulement"%len(self))
            return self, self.pourMille(0),(len(self),len(self))
        nd = 1000
        c0 = self.cpoints.copy()
        t, m = self.methode
        n = len(self)
        if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
            seen = [0, n/4, n/2, 3*n/4, n-1]
#             debug('cubic, periodic, 5 points',seen=seen)
            c1 = pointsFrom(self[seen])
        elif dist2(c0[0], c0[-1]) < 1.0e-6 :
            seen = [0, n/3, 2*n/3, n-1]
            c1 = pointsFrom(self[seen])
#             debug('ferme, 4 points',seen=seen)
        else:
            seen = [0,n-1]
            c1 = pointsFrom(self[seen])
#             debug('ouvert, 2 points',seen=seen)
        se = NSplineSimple(cpoints=c1,
                           methode=self.methode,
                           nbpd=self.nbpd,
                           name='%s-elaguee'%self.name,
                           mode=self.mode)
        if debog :
            X0, Y0 = XY(c0)
            more = [(X0, Y0, 'k.-','c0'),]
            texts = [(X0[0],Y0[0]-0.05,u'%d'%0),(X0[-1],Y0[-1]-0.05,u'%d'%(len(X0)-1))]
            texts+= [(c1[0,0],c1[0,1],u'se:%d'%0),(c1[-1,0],c1[-1,1],u'se:%d'%(len(c1)-1))]
        T1, D01, _, _ = se.projectionObj(c0, discret=nd)#T1, D01 de longueur len(c0)
        while len(seen)<len(c0):
#             debug(seen=seen)
            D01[seen] = 0.0#On ne veut pas remettre 2 fois le même point (=>point double)
            d = max(D01)#le point de c0 le plus eloigné de se
            imaxd = (D01 == d).nonzero()
            idx0 = imaxd[0][0]
            histoire.append((idx0, self.pourMille(d)))
            seen = sorted(list(set(seen+[idx0])))
            pos0 = c0[idx0]#position(x,y) du point à inserer dans se
            #On cherche maintenant en quelle position (k) il faut l'insérer
            K1 = se.knots.tolist()#Les t des noeuds de se
            t1 = T1[idx0]#le param t (sur se) de la distance max
            if debog :
                print u'    idx0 = %-6d       #numero dans c0'%idx0
                print u'    t1   = %-6.3f       #temps du projeté de c0[%d] sur se'%(t1,idx0)
                print u'    d01  = %-6.3f       #distance de c0[%d] à se'%(D01[idx0],idx0)
            # La position de t dans K1
            i = 0
            while K1[i] < t1 : i += 1
            #i est l'indice d'insertion dans se
            if debog :
                debug()
                print u'    knots = %s #les knots de se'%se.knots.tolist()
                print u'    insertion dans se en position i =',i#, ' ; dist=',D[idx],' ; maxD=',d
                pj = se(t1)
                more += [((pos0[0],pj[0]),(pos0[1],pj[1]),'go-',u'à ajouter en position %d'%(i))]
                texts += [(pos0[0],pos0[1],u'%d'%(i))]
                se.plot(titre=u'Spline élaguée : \nseen=%s, \ndist=%.1g ‰'%(str(seen),d),
                        more=more,texts=texts, show=True)
                more = more[:-1]
                c1 = se.cpoints
                texts = [(c1[0,0],c1[0,1],u'se:%d'%0),(c1[-1,0],c1[-1,1],u'se:%d'%(len(se)))]

            se.insertPoint(pos0, i)#ca devrait fonctionner sans exception...
            T1, D01, _, _ = se.projectionObj(c0, discret=nd)
            dm = self.pourMille(d)
#             debug(u'dist-%d = %.2g‰ '%(len(seen),dm))
            c1 = se.cpoints.copy()
            if dm<eps : break
#         histoire.append((idx0, self.pourMille(d)))
        if len(se) == len(self) :#même nb de points : on garde la spline initiale.
            se.cpoints = self.cpoints.copy()
            n0 = len(self)
            n1 = len(se)
            return se, d,(n0,n1)
        #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de se
        #ca ne marche pas du tout !!!
        n0 = len(self)
        n1 = len(se)
#         debug(historique_distances=["dm(%d)=%5.3g (pt %d)"%(k, dm, i) for k,(i,dm) in  enumerate(histoire)])
        debug(historique_distances=[float("%5.2f"%(dm)) for k,(_,dm) in  enumerate(histoire)])
        debug('Apres ELAGAGE : dist = %.2g mm/m ; %d points => %d points '%(dm,n0,n1))

        if replace :
            self.cpoints = se.cpoints
        return se, d, (n0,n1)


    ################################################################################
    #############                    Fonctions                     #################
    ################################################################################


def adet(u,v):
    u"""Retourne 2*surface du triangle défini par les vecteurs u et v"""
    return abs(u[0]*v[1]-u[1]*v[0])

def dist2s(P0, P):
    u"""
        - distance entre les polygones et P0 et P au sens suivant :
        la surface en valeur absolue entre les deux polygones P0 et P,
        - P et P0 ont le meme nombre de points.
        - La surface est calculée comme la somme des surfaces (en val. abs.)
        des n-1 quadrilatères P0[i-1], P0[i], P[i], P[i-1]"""
    if len(P0) != len(P) :
        raise ValueError("len(P0)=%d doit etre egal a len(P)=%d"%(len(P0), len(P)))
    AB = zip(P0, P)[:-1]
    surf = 0.0
    for (a0,b0), (a1,b1) in zip(AB[1:], AB[:-1]) :
#         a0a1, a1b0, b0a1, a1b1 = a1-a0, b0-a1, a1-b0, b1-a1
        surf += adet(a1-a0, b0-a1) + adet(a1-b0, b1-a1)
    return surf

def distance(s1, s2, precision):
    u"""
    distance separant deux splines, au sens :
    norme de s1(T)-s2(T)
    """
#     s2.precision = s1.precision = precision
    T = linspace(0, 1, precision)
    d = s1(T) - s2(T)
    return norm(d)
# #############################

def placementReperesMontage(B1,B2,delta,):
    u"""Besoin:
    ----------
    Calculer une correction sur les repères de montage des cloisons, intrados, extrados...
    pour éviter qu'ils se trouvent légèrement décalés à cause du décalage
    de montage et de la courbure des pièces
    Actuellement les repères sont décalés perpendiculairement au bords de pièce,
    du coup ils s'écartent si le bord est convexe et se rapprochent si le bord est concave.
    L'idée est de recalculer les repères de montage avec respect des distances
    des points initiaux des bords de pièces,
    ainsi les repères de toutes les pièces vont s'aligner parfaitement.

    Algo correction:
    -------------------
    Soit B un tableau de points initiaux 2D correspondant à un bord de pièce à plat.
    Soit RM1 le tableau des repères de montage calculés avec une fréquence d'échantillonage de 1
    (donc pas sous échatillonés, à un point de B correspond un point de RM1,
    je précise car après je les réchantillone).
    1. Calculer le tableau D des distances successives des points de B
    2. Calculer une spline SRM1 qui passe par les points de RM1
    3. Calculer un tableau de points RM2 sur la spline SRM1
       tel que les distances successives D soient respectées

    => le plus simple pour moi serait une fonction qui prenne en entrée B et RM1
    (donc deux tableaux de points avec le même nombre de points 2D qui
    représentent des courbes ayant à peu près la même longueur géométrique...
    mais la 2e peut être plus ou moins longue que la 1ère,
    je ne sais pas comment tu gères si les points dépassent au bout ou pas...)
    et qui sorte RM2

    Réponse : la fonction placementReperesMontage(B1, B2, delta)
    -------
    :param B1, B2: les tableaux de points des deux bords à assembler
        (quelconques mais de même longueur, à peu près)
    :type B1, B2 : ndarray((n,2))
    :param delta: (float) le  pas des repères
        (distance entre deux repères consécutifs, constant pour le moment mais
        ça peut se changer)
    :return S1, S2, T: deux splines NSplinesSimple, et un tableau T de paramètres t_i,
        de type ndarray((n,1)).
    :utilisation:
        - n=len(T) est le nombre de repères calculé, les reperes eux mêmes sont obtenus par
        - S1(T): un ndarray((n,2)) avec les coordonnées des repères sur B1
        - S2(T): un ndarray((n,2)) avec les coordonnées des repères sur B2
    :remarque: sur le cas test, l'erreur est de 0,4 mm pour des bords de longueur 2 m environ.
        Soit 2 pour 10 000.
    """
    ##############################################
    ### construction des repères de montage
    ##############################################
    s1 = NSplineSimple(cpoints=B1, parent=None,
#                         methode=('cubic','natural'),
                        methode=('ius',1),
                       name='s1')
    s2 = NSplineSimple(cpoints=B2, parent=None,
#                         methode=('cubic','natural'),
                        methode=('ius',1),
                       name='s2')
    l1, l2 = s1.longueur(), s2.longueur()
    lref = min(l1,l2)
    debug(longueur1=l1, longueur2=l2)
    n, r = divmod(lref, delta)#nb de reperes de montage
    debug(u'nb reperes=%d, reste=%.2g'%(n,r))
    n = int(n)
    T = arange(1,n+1)*(delta/lref)
    return s1, s2, T

def correctionReperesMontage(B,R, mode='production'):
    u"""Besoin:
    ----------
    Calculer une correction sur les repères de montage des cloisons, intrados, extrados...
    pour éviter qu'ils se trouvent légèrement décalés à cause du décalage
    de montage et de la courbure des pièces
    Actuellement les repères sont décalés perpendiculairement au bords de pièce,
    du coup ils s'écartent si le bord est convexe et se rapprochent si le bord est concave.
    L'idée est de recalculer les repères de montage avec respect des distances
    des points initiaux des bords de pièces,
    ainsi les repères de toutes les pièces vont s'aligner parfaitement.

    Algo correction:
    -------------------
    Soit B un tableau de points initiaux 2D correspondant à un bord de pièce à plat.
    Soit R le tableau des repères de montage calculés avec une fréquence d'échantillonage de 1
    (donc pas sous échatillonés, à un point de B correspond un point de R,
    je précise car après je les réchantillone).
    1. Calculer le tableau dB des distances successives des points de B
    2. Calculer une spline sR qui passe par les points de R
    3. Calculer un tableau de points R2 sur la spline sR
       tel que les distances successives dB soient respectées

    => le plus simple pour moi serait une fonction qui prenne en entrée B et R
    (donc deux tableaux de points avec le même nombre de points 2D qui
    représentent des courbes ayant à peu près la même longueur géométrique...
    mais la 2e peut être plus ou moins longue que la 1ère,
    je ne sais pas comment tu gères si les points dépassent au bout ou pas...)
    et qui sorte sR

    Réponse : la fonction correctionReperesMontage(B, R)
    -------
    qui fournit un tableau RC de n=len(B) repères de montage, répartis sur
    le **chemin polygonal** défini par R, de sorte que les distances entre points de RC
    soient les mêmes que les distances entre les points de B
    **Les distances entre points de RC sont mesurées le long du chemin R**
    (et non pas de point à point)

    :param B, R: points du bord et repères de montage d'origine (chemin R)
    :type B : ndarray((n,2))
    :type R : ndarray((m,2)) (n et m n'ont pas besoin d'être égaux)
    :return:
        - si mode = 'test' : (sR, T) = la spline polygonale passant par les
            repères d'origine, que j'appelle chemin polygonal ci-dessus,
            et les valeurs du paramètre t de sorte que les n repères de montage
            corrigés sont RC=sR(T).
        - si mode = 'production' :les repères de montage corrigés,
            sous forme d'un ndarray((n,2))
    :remarque: sur le cas test, construit par la fonction testCorrectionRM()
    la longueur du bord B est 4.303, la longueur du chemin R est 4.3365
    l'erreur max initiale sur les longueurs est 0.033 soit 8‰
    l'erreur max sur les longueurs corrigées est 4.5e-16 soit
        moins de 1.0e-15 (en valeur relative)
    """
    #spline repères
    sR = NSplineSimple(cpoints=R, methode=('ius',1))
    # distances à respecter (distances à B[0])
    distB = absCurv(B)
    # Comme sR est de degré 1, sur sR le paramètre t EST EGAL à
    # dist(sR(0.),sR(t))/sR.longueur == t
    #rdebug(u'verif dist(sR(0.),sR(t)) == t', longueur=sR.longueur)
    #l=sR.longueur
    #for t in linspace(0, 1., 20) :
    #    print (sR.abscurv(t)/(t*l), dist(sR(0),sR(t)))
    # OK : Sur sR : dist(sR(0.),sR(t)) == t*sR.longueur
    # paramètres T correspondant aux distances à respecter
    T = distB/sR.longueur()
    if mode == 'test' :
        return sR, T#pour tests
    elif mode == 'production' :
        return sR(T)#en production

if __name__=="__main__":
    from testsplinesimple import mainTest
    config.TEST_MODE = False
    mainTest()
#     placementReperesMontage(T=asarray([[2,-2],[1,-1],[0,0],[1,1],[2,2]]),
#                              TR=asarray([[1.5,-2.5],[0.5,-1.5],[-math.sqrt(2)*0.5,0],[0.5,1.5],[1.5,2.5]]))
#     exit()
    # app = QApplication(sys.argv)

