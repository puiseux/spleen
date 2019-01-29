#!/usr/local/bin/python2.7
#-*-coding: utf-8 -*-
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe NSpline
Description :
@author:      puiseux
@copyright:   2016-2017-2018 Nervures. All rights reserved.
@contact:    pierre@puiseux.name
__updated__="2019-01-29"
'''
from utilitaires.utilitaires import (rstack, eliminerPointsDoublesConsecutifs, diff,
    className, centreGravite, baryCentre)
# from splineabstraite import absCurvReal
from utilitaires.lecteurs import pointsFrom, LecteurUniversel
import sys,os,math
from array import array
#
import numpy as np
from numpy import log, linspace, asarray, sqrt, arange, zeros, cos, tan, pi,\
    loadtxt
from numpy.linalg import  norm
# import scipy as sp
# from config import VALIDATION_DIR, RUNS_DIR
from scipy.optimize import newton, minimize
from pprint import pprint
from utilitaires.utilitaires import (Path, segmentPlusProche, stack, debug, rdebug, dist,
                        hardScale, absCurv,dist2,rotate,courbure,symetrieAxe)
from splineabstraite import NSplineAbstract, computeSpline, distance2PointSpline
from numbers import Number
from scipy.integrate.quadpack import quad
import cPickle

class NSplineSimple(object):
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
        eps = 1.0e-5 # 10 microns si les distances sont en metres

    def __init__(self, **dump):# cpoints=None, name='spline', role='spline'):
        u"""
        * dump DOIT contenir les clés :
             - aucune
        * dump PEUT contenir les clés suivantes (sinon il y a des valeurs par défaut)
            - 'points' dont la valeur est un np.ndarray, ou un QPolygon
                ou une liste de points, ou un filename ou encore None
            - 'precision' = nb de points pour le tracé de la spline
            - 'methode' : cf computeSpline() paramètre hyper important.
                détermine s'il s'agit d'une spline cubique, de degré autre,
                périodique, si les dérivées aux extrémités sont données, etc...
            - 'mode' : le mode de répartition des points en cas d'échantillonage
            - 'name' : str ou unicode pour affecter un nom
            - 'role' : str ou unicode pour affecter un role...
        """
        super(NSplineSimple,self).__init__()
        self.setDefaultValues()
        self.load(dump)
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
                - ou un np.ndarray de n réels dans [0,1] (de shape (n,1)? ou (1,n)?)
                - ou une liste de n réels entre 0 et 1
            * diff est l'ordre de dérivation
        # Retourne :
        les 2 ndarrays de shape (n,1), zippés, contenant les points (sx(ti, diff), sy(ti, diff)) ou ti parcourt T
            * si diff=0, retourne les coordonnées du point x(t), y(t),
            * si diff=1, retourne les dérivées x'(t), y'(t)
            * si diff=2, retourne les dérivées secondes x"(t),y"(t)
        # Utilisation :
            >>> S = NSpline(....)
            >>> T = np.linspace(0.0, 1.0, 11)
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

    def __getitem__(self,k):
        u"""
        Pour accéder aux points du polygone de controle (en lecture)
        ne traite pas les slices
        Retourne un tuple (x,y)
        >>> S = NSplineSimple(points=....)
        >>> S[3]
        [10.0,3.14] le point points[3]
        ne gère pas les slices
        >>> S[1:3]
        ... AttributeError: 'QPolygonF' object has no attribute 'x'
        """
        return self._cpoints[k]

    def __setitem__(self,k,value):
        u"""
        mise à jour du polygone de contrôle et du polygone seulement
        Typiquement, on pourra ecrire
        >>> S = NSplineSimple(cpoints=[(1,1),(2,2),(3,3)])
        >>> S[1] = 20,30
        >>> print S[1]
        [20.0, 30.0]
        """
        if np.all(self._cpoints[k]==value) : #inutile de déclencher un _update()
            return
        self._cpoints[k] = value
        self._update()#indispensable pour mettre à jour self._xxx

    def __str__(self):
        return u'\n'.join(self.info)

    @property
    def info(self):
        infos=[
                u"<%s>"%className(self),
                u'%20s = '%u'name'+u'%s'%self.name,
                u'%20s = '%u'role'+u'%s'%self.role,
                u'%20s = '%u'nb pts controle'+u"%d"%len(self.cpoints),
                u'%20s = '%u'closed'+u'%s'%self.isClosed(),
                u'%20s = '%u'nb_splines, n_bech'+"%d, %d"%(self.nbspline, self.nbech),
                u'%20s = '%u'methode'+u"%s"%str(self.methode),
                u'%20s = '%u'nb points discretisation'+"%s"%str(self.nbpd),
                u'%20s = '%u'mode echantillonage'+u"%s"%self.mode,
                u'%20s = '%u'nb pts echantillon'+u"%s"%self.nbpe,
#                 u'%20s = '%u'nb updates'+u"%s"%self.nbupdate,
                ]
        i = u'largeur, hauteur'
        try :
            i1 = u"%g, %g"%(self.width, self.height)
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%20s = '%i+i1)
#
        i = u'position cg'
        try :
            i1 = u"%g, %g"%(self.centregravite[0], self.centregravite[1])
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%20s = '%i+i1)
#
        i = u'longueur'
        try :
            i1 = u"%g"%(self.longueur('r'))
        except Exception as msg :
            i1 = u"? (%s, %s)"%(className(self), msg)
        infos.append(u'%20s = '%i+i1)

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

    def load(self, dump):
        u"""Ce qui est commun à toutes les splines."""
#         super(NSplineSimple, self).load(dump)
        u"""
        Attention, la cle est 'points'  dans les anciens projets, et 'cpoints' dans les nouveaux.
        OUI=>On affecte _cpoints directement qui évite l'appel au setter de cpoints (qui appelle _update()).
        cpoints DOIT être liste, tuple, array ou np.ndarray, sans points doubles consécutifs
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

        for key in ('points', 'cpoints') :#la cle 'cpoints' est prioritaire, en cas
            if key in dump:
                cpoints = dump.pop(key)
                if isinstance(cpoints, (list, tuple, array)) :#types Python de base
                    #obligé de considerer ce cas, car pour le load d'un pkl cpoints est une liste
                    cpoints = asarray(cpoints)
                if not isinstance(cpoints, np.ndarray) :
                    msg = u'Les points de controle doivent etre de type np.ndarray((*,2)), et non pas %s.'%cpoints.__class__.__name__
                    msg = msg + u' \nUtilisez la fonction points=pointsFrom(points) pour transformer a peut pres n\'importe quoi en np.ndarray((n,2)).'
                    msg = msg + u' \nLes points de controle ne doivent pas avoir de points doubles consecutifs.'
                    msg = msg + u' \nUtilisez points = eliminerPointsDoublesConsecutifs(points) en cas de doute.'
                    raise TypeError(msg)
                if self.methode[1]=='periodic' and len(cpoints) > 0 :#and np.any(cpoints[0] != cpoints[-1]) :
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
    def save(self, filename):
        filename = Path(filename)
        ext = filename.ext
#         debug(filename=filename, ext=ext)
        try :
            if ext in (".gnu", '.txt'):
                #seulement les points échantillonnés
                np.savetxt(filename, self.epoints)
            # elif ext==".pts":
            #     #seulement les points échantillonnés
            #     writeProfil(filename, self.epoints)
            # elif ext=='.npkl':
                #Spline complète, pickle
                cPickle.dump(self.toDump('new'),open(filename,'w'))
            elif ext=='.pkl':
                #seulement les points échantillonnés
                cPickle.dump(self.toDump('old'),open(filename,'w'))
            elif ext=='.npkl':
                cPickle.dump(self.toDump('new'),open(filename,'w'))
            elif ext=='.spl':
                #Spline complète, dict
                with open(filename,'w') as f :
                    f.write(str(self.toDump()))
            elif ext=='.dxf':
                raise NotImplementedError('Format dxf')
            debug('Sauvegarde %s : OK'%filename)
        except Exception as msg:
            rdebug(msg, 'Sauvegarde %s : impossible'%filename)

    def open(self, filename):
        filename = Path(filename)
        ext = filename.ext
        #debug(filename=filename)
        if ext in (".gnu", '.txt'):
            dump = self.toDump()
            dump['cpoints'] = loadtxt(filename)
        elif ext==".pts":
            dump = self.toDump()
            dump['cpoints'] = LecteurUniversel(filename).points
        elif ext in ('.pkl', '.npkl'):
            dump = cPickle.load(open(filename,'r'))
            for key in ('points', 'cpoints') :
                if dump.has_key(key) :
                    dump[key] = pointsFrom(dump.pop(key))
#         #debug(dump=dump)
        self.load(dump)

    @property
    def cpoints(self):
        u"""
        points de contrôle, sous forme np.ndarray((n,2)), None si non initialisé
        """
        return self._cpoints

    @cpoints.setter
    def cpoints(self, points):
        u"""points de contrôle, sous forme np.ndarray((n,2))
        - points est un np.ndarray((n,2)) et ne doit pas avoir de points doubles consécutifs."""
#         stack('points=%s'%points)
        if points is None or len(points)<1 :
            self._cpoints = np.zeros((0,2))
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
        return self._dpoints

    @property
    def epoints(self):
        u"""points échantillonnés, sous forme np.ndarray((n,2))"""
        if not hasattr(self,'_epoints') :
            self._epoints = self(self.tech)
        return self._epoints

    @property
    def tech(self):
        u"""Les parametres T de l'echantillonnage"""
        if not hasattr(self, '_tech') :
            self._tech = self.echantillonner()
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
            del self._dpoints

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

    def echantillonner(self):
        u"""
        Méthode usuelle d'échantillonnage.
        Utiliser cette méthode pour échantillonner avec les paramètres de self
        et conserver le résultat dans self.epoints, self.tech.
        :utilisation :
            >>> S.mode='courbure'
            >>> S.nbpe=27
            >>> e = S.echantillonner()
            >>>     #e et S.epoints contiennent les points echantillonnés
            >>>     #les parametres des points d'echantillonnage sont dans S.tech
            >>> T = S.echantillonner(True)
            >>>     # T contient les parametres des points d'echantillonnage
            >>>     # S.epoints contient les points echantillonnés
        """
        return self._echantillonner(self.nbpe, self.mode, 0, 1)

    def _echantillonner(self, nbp, mode, ta=0, tb=1):
        u"""
        répartition (échantillonnage) de nbp points sur la spline self,
        entre les abscisses ta et tb, suivant le mode précisé par 'mode'
        *supprime self.epoints et modifie self._tech*
        :return : les abscisses T=np.ndarray(shape=(n,1)) (je crois!)
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

        :param ta : float dans [0,1] ou np.ndarray((n,1)) comme retourné par
            self.absCurv()
            - facultatif, ta=0 par defaut
            - l'abs. curv. du premier point d'échantillonnage

        :param tb : float dans [0,1] ou np.ndarray((n,1)) comme retourné par
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
            if (ta, tb) == (0,1) : return T
            else : return ta + (tb-ta)*T

        elif mode == 'x3' :
            #Points serrés au milieu, lâches au debut et à la fin
            T = linspace(-1, 1, nbp)**3
            T = 0.5*(1+T)
            if (ta, tb) == (0,1) : return T
            else : return ta + (tb-ta)*T

        elif mode in ('linear', 'lineaire', 'lin') :
            return linspace(ta, tb, nbp)#y compris ta et tb

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
            rdebug(mode,self.mode)
            if nbp==1 :
                debug(nbp=nbp, ta=ta,tb=tb)
                stack()
                return ta
#                 else : return self(ta)
            N = 100001
            T = linspace(ta,tb,N)#N-1 micro intervalles
            dt = T[1]-T[0]#micro intervalles de largeur dt
            CA = np.sqrt(np.abs(self.courbure(T)))
            S = (CA[1:]+CA[:-1])*(dt/2)#Les sinuosités des N micro-intervalles = dt*moyenne des courbures aux extrémités
            s0 = sum(S)/(nbp-1)#La sinuosité (constante) de chaque intervale Te[j],Te[j+1] défini ci dessous
            Te = np.ones(nbp)*ta# les te[j] cherchés
            Te[-1] = tb
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
                            deb, msg = debug, [u'\nPas grave %s'%self.name]
                        else :
                            deb, msg = rdebug, [u'\nAttention']
                            msg += [
                                   'i=%d'%i,
                                   'Te=%s'%(str(Te))]
                        deb('\n    '.join(msg))
                        break
                    i += 1#i=
                Te[j] = T[i]
            return Te
        else :
            rdebug('mode echantillonnage inconnu : ',mode=mode, nbp=nbp, ta=ta, tb=tb)
            return zeros((0,))

    @property
    def methode(self):
        return self._methode
    @methode.setter
    def methode(self, newmethode):
        u"""
        Changement de spline (ius, cubic, periodic,...)
        """
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
                if np.isnan(sx_05) or np.isnan(sy_05) :# in (np.nan,np.inf) or sy_05 in (np.nan,np.inf) :
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

    def plot(self, plt, control=True, nbpd=None, nbpe=None, mode=None,
             titre=None, more=[], show=True, buttons=True):
        """
        :param plt : une instance de pyplot, obtenue en amont par :
              >>> from matplotlib import pyplot as plt
        :param control : bool, affichage des points de controle True/False
        :param nbpd : int, nb points de discretisation
        :param nbpe : int, nb de points d'echantillonage
        :param titre : str ou unicode, titre
        :param more : list, [(X,Y, couleur,'nom'), ...] tracé supplementaire
          """
        from matplotlib.widgets import CheckButtons
#         plt.figure(numfig)
#         rdebug('***********')
        if nbpd is None : nbpd = self._nbpd
        if nbpe is None : nbpe = self.nbpe
        if mode is None : mode = self.mode
#         debug('appel echantillonnage', type(self))

        D = self.dpoints
        C = self.cpoints
        E = self.epoints
        _, ax = plt.subplots()
        if titre is None : titre = self.name+str(self.methode)
        plt.title(titre)
        spline, = ax.plot(D[:,0], D[:,1], 'b-', lw=1)
        echantillon, = ax.plot(E[:,0], E[:,1], 'g.', lw=1)
        control, = ax.plot(C[:,0], C[:,1], 'ro', lw=1)
        if buttons :
            butt = ['control', 'spline','echantillon',]
            values = [True, True, True]
            draws = [control, spline, echantillon]
        for x, y, color, name in more:
            temp, = ax.plot(x, y, color)
            if not name : continue
            if buttons :
                draws.append(temp)
                butt.append(name)
                values.append(True)
        plt.subplots_adjust(left=0.2)
        plt.axis('equal')

        if buttons :
            rax = plt.axes([0.05, 0.4, 0.1, 0.15])
            check = CheckButtons(rax, butt, values)

            def func(label):
                if label == 'spline':
                    spline.set_visible(not spline.get_visible())
                elif label == 'control':
                    control.set_visible(not control.get_visible())
                elif label == 'echantillon':
                    echantillon.set_visible(not echantillon.get_visible())
                else :
                    draw = draws[butt.index(label)]
                    draw.set_visible(not draw.get_visible())
                plt.draw()
            check.on_clicked(func)
        if show : plt.show()
        return plt

    def longueur(self, p='r'):
        if p=='c':
            return absCurv(self.cpoints, normalise=False)[-1]
        elif p=='d' :
            return absCurv(self.dpoints, normalise=False)[-1]
        elif p=='e' :
            return absCurv(self.epoints, normalise=False)[-1]
        else:#longueur vraie
            raise NotImplementedError('TODO : longueur(self, p='r') (longueur reelle)')

    @property
    def height(self):
        if not hasattr(self, '_height') :
            self._height = max(self.dpoints[:,1]) - min(self.dpoints[:,1])
        return self._height

    hauteur = height

    @property
    def width(self):
        if not hasattr(self, '_width') :
            self._width = max(self.dpoints[:,0]) - min(self.dpoints[:,0])
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
        return centreGravite(self.dpoints)
    centregravite=gravitycenter

    @property
    def barycentre(self):
        u"""Le centre de gravité du nuage de points matériels cpoints."""
        return baryCentre(self.cpoints)

    def courbure(self, T) :
        dx,  dy  = self.sx(T, 1), self.sy(T, 1)
        d2x, d2y = self.sx(T, 2), self.sy(T, 2)
        norm3_d2 = np.sqrt(dx**2+dy**2)**3
        sc = (dx*d2y-dy*d2x)/(norm3_d2)
        return sc
        # si norm_d2=0, x"(t)=y"(t)=0, c'est une droite, courbure nulle
#         sc[np.where(norm3_d2 < 1.0e-12)] = 0.0

    def absCurv(self, T):
        u"""
        Calcule et retourne l'abscisse curviligne réelle des points self(T) sur la spline self.
        L'abscisse curviligne d'un point de paramètre t dans [0,1] est
        l'intégrale de 0 à t de phi(t) = sqrt(self.sx(t)**2+self.sy(t)**2)
        L'intégration est assurée par scipy.integrate.quad()
        Si la spline self a trop de points de contrôle, ca rame et l'erreur est importante
        err
        :param self: une NSplineAbstract
        :param T: les n paramètres t des points dont on veut l'abscisse curviligne.
            Ces n paramètres doivent être dans l'intervalle [0,1]
        :type T: au choix
            - un ndarray de shape (n,1) à valeurs réelles dans [0,1],
            - une liste de n valeurs réelles dans [0,1],
            - un tuple de n valeurs réelles dans [0,1],
            - un réel unique t dans [0,1]
        :return ac, err:
            - deux ndarray de n réels avec ac = abscisses curvilignes, err erreur estimée
            - ou bien deux réels, si T est réel
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
            return ac, err

    #             return ac1+ac2, max(err1, err2)
        else :
            res = array([absCurv(self, t) for t in T])
    #         res = asarray([quad(lambda s : sqrt(S.sx(s,1)**2+S.sy(s,1)**2), 0.0, t) for t in T])
            return res[:,0],  res[:,1]

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
    valeur absolue de la surface de s1-s2
    """
#     s2.precision = s1.precision = precision
    T = linspace(0, 1, precision)
    d = s1(T) - s2(T)
    return np.linalg.norm(d)
# #############################

def elaguer(p0, eps=1.0e-3):
    Np = p0.precision#nb de points de discretisation
    D0 = p0.dpoints#les points discretisés de la spline de reference
    C0 = p0.cpoints#les points de contrôle de la spline de reference
#     L0 = p0.longueur()
    T = linspace(0, 1, Np)
    methode = p0.methode
#     D0courb = courbure(D0)
#     C0courb = courbure(C0)
    global nbiter
    nbiter = 0
    def f(C, T, methode):
        u"""Fonction à minimiser : distance de la spline S0 à la spline d'interpolation de C.
        C = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(C))"""
#         print T
        global nbiter
        nbiter += 1
        oldshape = C.shape#tableau a plat (np.ravel(C))
        C.shape = (len(C)/2,2)
#         Ts, Sx, Sy = splineInterpolation(C, 'c cubic')#la spline d'interpolation de C
        S = NSplineSimple(cpoints=C, methode=methode)
#         print Ts
        C.shape = oldshape
        D = S(T)#discretisation de la spline S
#         D = asarray(zip(Sx(T), Sy(T)))
#         return np.sqrt(np.linalg.norm(D0-D))
        E = D-D0
        E[ 0] *= 1.0e10#pénalisation des extrémités qui bouge (on ne veut pas qu'elles bougent)
        E[-1] *= 1.0e10

        norm = np.linalg.norm(E)
        if nbiter%100 == 0 : print 'iteration %d : %.2g'%(nbiter, norm)
        return norm

    def numPointsToKeep(C, lref, genre='shark'):
        curv = courbure(C)*lref
        dcurv = np.abs(curv[:-1]-curv[1:])
        tokeep = dcurv>0.1#pour sharknose : 0.12
#         tokeep = dcurv>0.15#sans sharknose : 0.15
        #Dans tous les cas on en supprime un peu moins que la moitié
        nptk = set([idx for idx, tk in enumerate(tokeep) if tk])
        nptkp = set([idx+1 for idx in tokeep if idx+1 < len(curv)]) - nptk
        nptkm = set([idx-1 for idx in tokeep if idx-1 >= 0 ]) - nptk
        nptk = nptk.union(nptkp.union(nptkm))
        return sorted(list(nptk))
#
    A = C0[ 0]
    B = C0[-1]
    corde = np.linalg.norm(B-A)
    numeros = list(set(numPointsToKeep(C0, corde) + [0,len(p0)-1,]))
    numeros.sort()
    fuera = list(set(range(len(p0))) - set(numeros))
    fuera.sort()
    print u'numeros des points conservé=%s'%str(numeros)
    print u'numéros des points supprimés=%s'%str(fuera)
    C0 = p0.cpoints[numeros].copy()
#     print 'C0=',C0.tolist()
    print 'f(C0,T) = %.2e'%f(np.ravel(C0),T,methode)
    res = minimize(f, np.ravel(C0), (T,methode), 'BFGS', options={'maxiter':0,'gtol':eps})#, 'disp':True})
    print 'res.fun=%.2e'%res.fun
    C1 = res.x.reshape((-1,2))
#     C1[ 0]   = A
#     C1[-1]   = B
    print 'x-C0', C1-C0
    p1 = NSplineSimple(cpoints=C1, precision=10*Np, methode=p0.methode, mode='courbure') #mode=p0.mode)#
    p0.precision = 10*Np
    p0._update()
    D0 = p0.dpoints
    D1 = p1.dpoints
    L0 = p0.dlongueur
    L1 = p1.dlongueur
    C1 = p1.cpoints
#     debug(len_C1=len(C1), len_T=len(T))

    print u'  longueur  L1          = %.2e '%abs(L1)
    print u'  δ-longueur            = %.2e '%abs(L1-L0)
    print u'  norme frobenius D0-D1 = %.2e '%norm(D0-D1, 'fro')
    print u'  max norme 2     D0-D1 = %.2e <=='%max(norm(D0-D1, 2, axis=1))
    print u'  distance(p0,p1)/L0    = %.2e '%(distance(p0, p1, Np)/L0)
    print u'  f(p1)                 = %.2e '%f(np.ravel(C1), linspace(0, 1, 10*Np),methode)
    print u'  nb de points éliminés = %d / %d'%(len(fuera), len(p0))

#     from matplotlib import pyplot
# #     C0 = p0.cpoints[numeros]
#     minx, maxx = np.min(C0[:,0])-0.1, np.max(C0[:,0])+0.1
#     miny, maxy = np.min(C0[:,1])-0.1, np.max(C0[:,1])+0.1
#     pyplot.plot(
#                 D0[:,0], D0[:,1],'r-',
#                 C0[:,0], C0[:,1],'ro',
#                 D1[:,0], D1[:,1] ,'b-',
#                 C1[:,0], C1[:,1] ,'b>',
#                 (minx,maxx,0.0,0.0), (0.0, 0.0,miny, maxy), 'w.',#pour cadre large
#                 )
#     pyplot.show()
    return p1
######################
global nbiter
nbfunc=0

def f1(c, s0, Td):
    u"""Fonction à minimiser :
    distance des points de contrôle s0.cpoints à la spline d'interpolation de c.
    c = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(c))
    s0 = spline de reference
    Td = les abscurv dans[0,1] des points de discrétisation pour la comparaison"""
#         print T
    global nbfunc
    nbfunc += 1
    oldshape = c.shape#tableau a plat (np.ravel(C))
#     debug (C.shape, len(C)/2,2)
    c.shape = (len(c)/2,2)
    s = NSplineSimple(cpoints=c, methode=s0.methode)
#     E = S0(T)-S(T)
    E = [distance2PointSpline(s0.cpoints[k], s, t0=t).fun for k, t in enumerate(s0.sx.x)]
    E[0]*=1.0e3
    E[-1]*=1.0e3
#     norme = np.sqrt(sum(E[1:-1]))
    norme = np.sqrt(sum(E))
    c.shape = oldshape
#     norm = np.linalg.norm(E)
#     if nbfunc%100 == 0 :
#         print 'nbfunc %d : %.2g'%(nbfunc, norme)
    return norme
def f(C, S0, T):
    u"""Fonction à minimiser : distance de la spline S0 à la spline d'interpolation de C.
    C = famille de points de R^2 = np.ndarray((2*n)), mis a plat (avec np.ravel(C))
    T = les abscurv dans[0,1] des points de discrétisation pour la comparaison"""
#         print T
    global nbiter
    nbiter += 1
    oldshape = C.shape#tableau a plat (np.ravel(C))
#     debug (C.shape, len(C)/2,2)
    C.shape = (len(C)/2,2)
    S = NSplineSimple(cpoints=C, methode=S0.methode)
    E = S0(T)-S(T)
    E[0]*=1.0e10
    E[-1]*=1.0e10
    C.shape = oldshape
    norm = np.linalg.norm(E)
    if nbiter%100 == 0 :
        print 'iteration %d : %.2g'%(nbiter, norm)
    return norm

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
    :type B1, B2 : np.ndarray((n,2))
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
    l1, l2 = s1.longueur, s2.longueur
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
    :type B : np.ndarray((n,2))
    :type R : np.ndarray((m,2)) (n et m n'ont pas besoin d'être égaux)
    :return:
        - si mode = 'test' : (sR, T) = la spline polygonale passant par les
            repères d'origine, que j'appelle chemin polygonal ci-dessus,
            et les valeurs du paramètre t de sorte que les n repères de montage
            corrigés sont RC=sR(T).
        - si mode = 'production' :les repères de montage corrigés,
            sous forme d'un np.ndarray((n,2))
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
    T = distB/sR.longueur
    if mode == 'test' :
        return sR, T#pour tests
    elif mode == 'production' :
        return sR(T)#en production

if __name__=="__main__":
    from testsplinesimple import mainTest
    mainTest()
#     placementReperesMontage(T=np.asarray([[2,-2],[1,-1],[0,0],[1,1],[2,2]]),
#                              TR=np.asarray([[1.5,-2.5],[0.5,-1.5],[-math.sqrt(2)*0.5,0],[0.5,1.5],[1.5,2.5]]))
#     exit()
    # app = QApplication(sys.argv)

