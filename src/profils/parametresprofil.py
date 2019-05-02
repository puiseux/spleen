#!/usr/bin/python
#-*-coding: utf-8 -*-
'''
Axile -- Outil de conception/simulation de parapentes Nervures

    Classe ParamGeneraux

@author:     Michel Le Berre
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    creation: 08 Jan 2013
__updated__ = "2019-02-06"
'''
import sys, os
from utilitaires import debug, rdebug, stack
from preferences import ProfilPrefs

class ProfsParamNew(object):
    u'''
    Attributs :
    ---------
    * Modifiables :
        self.nptprof    => nb de points profil
        self.iouvext, self.iouvint => index ouverture
        self.iba        => index BA = le point le plus éloigné du BF
        self.nbpbf      => nb de pts réservés au BF pour pinces
        self.pourcentbf => % de corde au BF occupé par les nbpbf points réservés

    * NON modifiables : les autres
        self.nptint  => nb points intrados
        self.nptext  => nb points extrados
        self.nptret  => nb points retour Y COMPRIS le point BA (BA -> iouvext)
        self.copy (property) fournit un clone de self.

    Méthodes & exemples :
    -------------------
    >>> pp1 = ProfsParam(10, 6, 7)
    >>> pp2 = ProfsParam(11, 5, 7)
    >>> pp1 == pp2 #tester si deux profparam sont identiques
    False
    >>> pp1 != pp2 #tester si deux profparam sont différents
    True
    >>> print pp1
    nb points [ext+int] : 16=[10+6]; iouverture=(9,10); iba=7, nb points retour=3
    '''
#     def __init__(self, npe=0, npi=0, iba=0):
    def __init__(self, nptprof=0, iouvext=0, iouvint=0, iba=0):
        """Suivant la date de création du projet, et l'endroit d'ou il est invoqué,
        l'instanciation peut être faite de différentes manières :
        - ProfsParamNew(pp) ou pp peut être
            * un dict {'nptprof':nptprof, 'iouvext':iouvext, 'iouvint':iouvint, 'iba':iba}
            * une tuple pp=(nptprof, iouvext, iouvint, iba)
            * une liste pp=[nptprof, iouvext, iouvint, iba]
            * un str ou unicode ou QString
        - ProfsParamNew(nptprof, [iouvext, [iouvint, [iba]]])
        """
        raise TypeError("Interdit d'instancier un ProfParamNew")
        if isinstance(nptprof, (tuple,list)):
            raise TypeError('Attention, il faut probablement ecrire pp=ProfsParamNew(*liste), au lieu de pp=ProfsParamNew(liste)')
        elif isinstance(nptprof, (ProfsParam, ProfsParamNew)) :
            raise TypeError('Attention, il faut probablement ecrire pp=x, au lieu de pp=ProfsParamNew(x)')
        elif isinstance(nptprof, dict) :
            raise TypeError('Attention, il faut probablement ecrire pp=ProfsParamNew(**x), au lieu de pp=ProfsParamNew(x)')
#         elif isinstance(nptprof, (str, QString, unicode)) :
#             raise TypeError('Attention, il faut probablement ecrire pp=ProfsParamNew(**x), au lieu de pp=ProfsParamNew(*x)')
        self._nptprof = nptprof
        self._iba = iba
        self.iouverture = ((iouvext,iouvint))
#         self._iouvext = iouvext
#         self._iouvint = iouvint

    def __eq__(self, other):
        return  self.nptprof    == other.nptprof and\
                self.iouvext    == other.iouvext and\
                self.iouvint    == other.iouvint and\
                self.iba==other.iba
    def __ne__(self, other):
        return not self.__eq__(other)

    def castToProfsParam1(self):
        """"""
        return ProfsParam1(self.nptprof, self.iouvext, self.iouvint, self.iba,
                             ProfilPrefs.nbpbf, ProfilPrefs.pourcentbf)

    def isSuspect(self):
        suspect = [
                    self.nptext  <= 0,
                    self.nptint  <= 0,
                    self.iouvint <= 0,
                    self.iouvext <= 0,
                    self.nptret  <= 0,#ouverture a l'extrados
                    self.nptouv  <  0,#recouvrement
                    self.iouvint <  self.iouvext,#recouvrement
                    self.nptouv  >  self.nptprof,
                    self.iba     >= self.nptprof,
                    self.iba     <= 0
                   ]
        s = False
        for k, v in enumerate(suspect) :
            s = s or v
        return s

    def isForbiden(self):
        forbiden = [
                    self.nptext  <= 0,
                    self.nptint  <= 0,
                    self.nptouv  <  0,
                    self.iouvint <= 0,
                    self.iouvext <= 0,
                    self.nptret  <= 0,
                    self.iba     <= 0,
                    self.iouvint <  self.iouvext,
                    self.nptouv  >  self.nptprof,
                    self.iba     >= self.nptprof,
                   ]
        f = False
        for v in forbiden :
            f = f or v
        return f
    @property
    def nptext(self):
        return self.iouvext + 1

    @property
    def nptint(self):
        return self.nptprof - self.iouvint

    @property
    def nptexttheorique(self):
        return self.iba + 1

    @property
    def nptinttheorique(self):
        return self.nptprof - self.iba

    @property
    def nptouv(self):#nb de points dans l'ouverture NON compris les bouts
        return self.nptprof - self.nptext - self.nptint

    @property
    def nptret(self):#nb points retour, Y COMPRIS le point BA et ouvext
        return self.iouvext - self.iba + 1

    @property
    def iouverture(self):
        return self.iouvext, self.iouvint
    @iouverture.setter
    def iouverture(self, (ie, ii)):
        u"""
        Afin d'éviter les erreurs, on ne peut modifier les valeurs de iouvext et iouvint que via
        >>> iouverture = (ke, ki)
        """
#         debug( iouvext=ie, iouvint=ii)
        if ie < self.iba or ii < self.iba :
            rdebug('Attention, ouverture (partiellement) sur l\'extrados, l\'echantillonnage est impossible (TODO)')
            self._iouvext, self._iouvint = ie, ii
        elif 0 <= ie <= ii <= self.nptprof :
            self._iouvext, self._iouvint = ie, ii
        else:
            rdebug('Valeurs incorrectes, on devrait avoir 0 <= iouvext(=%d) <= iouvint(=%d) <= nptprof(=%d)'%(ie, ii, self.nptprof))
            self._iouvext, self._iouvint = ie, ii

#             raise ValueError('Valeurs incorrectes, on devrait avoir 0 <= iouvext(=%d) <= iouvint(=%d) <= nptprof(=%d)'%(ie, ii, self.nptprof))
    @property
    def iouvext(self):
        return self._iouvext
    @iouvext.setter
    def iouvext(self, k):
        raise NotImplementedError("On ne peut pas modifier iouvext seul.\
         On doit modifier les deux points de l'ouverture en utilisant 'iouverture = (ke, ki)'.")
    @property
    def iouvint(self):
        return self._iouvint
    @iouvint.setter
    def iouvint(self, k):
        raise NotImplementedError("On ne peut pas modifier iouvint seul.\
         On doit modifier les deux points de l'ouverture en utilisant 'iouverture = (ke, ki)'.")
    @property
    def nptprof(self):
        return self._nptprof
    @nptprof.setter
    def nptprof(self, npt):
        raise ValueError("c'est ici !!! : npt=%d"%npt)
    @property
    def iba(self):
#         rdebug( "BUG : iba peut varier d'un profil a l'autre")
        return self._iba
    @iba.setter
    def iba(self, k):
        self._iba = k
#         raise NotImplementedError
############################################################
    @property
    def copy(self):
        return ProfsParamNew(self.nptprof, self.iouvext, self.iouvint, self.iba)
    @property
    def info(self):
        return ['<%s>'%self.__class__.__name__,
                '  nb points profil             = %d'%(self.nptprof),
                '  nb points extrados tissu     = %d'%(self.nptext),
                '  nb points extrados theorique = %d'%(self.nptexttheorique),
                '  indice BA                    = %d'%(self.iba),
                '  nb points retour             = %d'%(self.nptret),
                '  indices ouverture            = (%d,%d)'%self.iouverture,
                '  nb points ouverture          = %d'%(self.nptouv),
                '  nb points intrados tissu     = %d'%(self.nptint),
                '  nb points intrados theorique = %d'%(self.nptinttheorique),
                ]

    def shortinfo(self):
#         return ' <%s> [nptprof=%d, iouverture=%s, iba=%d]'%(self.__class__.__name__, self.nptprof, str(self.iouverture), self.iba, )
        return ' [nptprof=%d, iouverture=%s, iba=%d]'%(self.nptprof, str(self.iouverture), self.iba, )

    def __str__(self):
        return self.shortinfo()
        return '\n'.join(self.info)

    @property
    def dump(self):
        return self.toDump()

    def toDump(self):
        todump = {'nptprof': self.nptprof, 'iouvext':self.iouvext,
                  'iouvint':self.iouvint, 'iba':self.iba}
        return todump

    def load(self,dump):
        self.__init__(**dump)
class ProfsParam1(object):
    u'''
    Attributs :
    ---------
    * Modifiables :
        self.nptprof    => nb de points profil
        self.iouvext, self.iouvint => index ouverture
        self.iba        => index BA = le point le plus éloigné du BF
        self.nbpbf      => nb de pts réservés au BF pour pinces
        self.pourcentbf => % de corde au BF occupé par les nbpbf points réservés

    * NON modifiables : les autres
        self.nptint  => nb points intrados
        self.nptext  => nb points extrados
        self.nptret  => nb points retour Y COMPRIS le point BA (BA -> iouvext)
        self.copy (property) fournit un clone de self.

    Méthodes & exemples :
    -------------------
    >>> pp1 = ProfsParam1(10, 6, 7)
    >>> pp2 = ProfsParam1(11, 5, 7)
    >>> pp1 == pp2 #tester si deux profparam sont identiques
    False
    >>> pp1 != pp2 #tester si deux profparam sont différents
    True
    >>> print pp1
    nb points [ext+int] : 16=[10+6]; iouverture=(9,10); iba=7, nb points retour=3
    '''
#     def __init__(self, npe=0, npi=0, iba=0):
    def __init__(self, nptprof=0, iouvext=0, iouvint=0, iba=0,
                 nbpbf=0, pourcentbf=0.0):
        """Suivant la date de création du projet, et l'endroit d'ou il est invoqué,
        l'instanciation peut être faite de différentes manières :
        - ProfsParam1(pp) ou pp peut être
            * un dict {'nptprof':nptprof, 'iouvext':iouvext, 'iouvint':iouvint,
                       'iba':iba, 'nbpbf':nbpbf, 'pourcentbf':pourcentbf}
            * une tuple pp=(nptprof, iouvext, iouvint, iba, [nbpbf, [pourcentbf]])
            * une liste pp=[nptprof, iouvext, iouvint, iba, [nbpbf, [pourcentbf]]]
            * un str ou unicode ou QString
        - ProfsParam1(nptprof, [iouvext, [iouvint, [iba, [nbpbf, [pourcentbf]]]]])
        """
        if isinstance(nptprof, (tuple,list)):
            raise TypeError('Attention, il faut probablement ecrire pp=ProfsParam1(*liste), au lieu de pp=ProfsParam1(liste)')
        elif isinstance(nptprof, (ProfsParam, ProfsParam1)) :
            raise TypeError('Attention, il faut probablement ecrire pp=x, au lieu de pp=ProfsParam1(x)')
        elif isinstance(nptprof, dict) :
            raise TypeError('Attention, il faut probablement ecrire pp=ProfsParam1(**x), au lieu de pp=ProfsParam1(x)')
#         elif isinstance(nptprof, (str, QString, unicode)) :
#             raise TypeError('Attention, il faut probablement ecrire pp=ProfsParam1(**x), au lieu de pp=ProfsParam1(*x)')
        self._nptprof = nptprof
        self._iba = iba
        self.iouverture = ((iouvext,iouvint))
        try :
            self.nbpbf = nbpbf
        except AttributeError :
            self.nbpbf = ProfilPrefs.nbpbf
        try :
            self.pourcentbf = pourcentbf
        except AttributeError :
            self.pourcentbf = ProfilPrefs.pourcentbf
#         self._iouvext = iouvext
#         self._iouvint = iouvint

    def __eq__(self, other):
        return  self.nptprof    == other.nptprof and\
                self.iouvext    == other.iouvext and\
                self.iouvint    == other.iouvint and\
                self.iba        == other.iba     and\
                self.nbpbf      == other.nbpbf   and\
                self.pourcentbf == other.pourcentbf

    def __ne__(self, other):
        return not self.__eq__(other)

    def castToProfsParam1(self):
        """pour compatibilite ascendante"""
        return self

    def isSuspect(self):
        suspect = [
                    self.nptext  <= 0,
                    self.nptint  <= 0,
                    self.iouvint <= 0,
                    self.iouvext <= 0,
                    self.nptret  <= 0,#ouverture a l'extrados
                    self.nptouv  <  0,#recouvrement
                    self.iouvint <  self.iouvext,#recouvrement
                    self.nptouv  >  self.nptprof,
                    self.iba     >= self.nptprof,
                    self.iba     <= 0,
                    self.nbpbf   >= self.nptext,
                    self.nbpbf   >= self.nptint
                   ]
        msgs = [
                'nptext=%-2d  <= 0'%self.nptext,
                'nptint=%-2d  <= 0'%self.nptint,
                'iouvint=%-2d <= 0'%self.iouvint,
                'iouvext=%-2d <= 0'%self.iouvext,
                'nptret=%-2d  <= 0,#ouverture a l\'extrados'%self.nptret,
                'nptouv=%-2d  <  0,#recouvrement'%self.nptouv,
                'iouvint=%-2d <  self.iouvext=%d,#recouvrement'%(self.iouvint,self.iouvext),
                'nptouv=%-2d  >  self.nptprof=%d'%(self.nptouv,self.nptprof),
                'iba=%-2d     >= self.nptprof=%d'%(self.iba,self.nptprof),
                'iba=%-2d     <= 0'%self.iba,
                'nbpbf=%-2d   >= self.nptext=%d'%(self.nbpbf,self.nptext),
                'nbpbf=%-2d   >= self.nptint=%d'%(self.nbpbf,self.nptint)
               ]
#         print self
        for val, msg in zip(suspect, msgs):#suspect :
#             print "suspect=%s, msg='%s'"%(val,msg)
            if val :
#                 rdebug("Profparam suspect : %s"%msg)
                return msg
        return ''

    def isForbiden(self):
        forbiden = [
                    self.nptext  <= 0,
                    self.nptint  <= 0,
                    self.nptouv  <  0,
                    self.iouvint <= 0,
                    self.iouvext <= 0,
                    self.nptret  <= 0,
                    self.iba     <= 0,
                    self.iouvint <  self.iouvext,
                    self.nptouv  >  self.nptprof,
                    self.iba     >= self.nptprof,
                    self.nbpbf   >= self.nptext,
                    self.nbpbf   >= self.nptint
                   ]
        f = False
        for v in forbiden :
            f = f or v
        return f
    @property
    def nptext(self):
        return self.iouvext + 1

    @property
    def nptint(self):
        return self.nptprof - self.iouvint

    @property
    def nptexttheorique(self):
        return self.iba + 1

    @property
    def nptinttheorique(self):
        return self.nptprof - self.iba

    @property
    def nptouv(self):#nb de points dans l'ouverture NON compris les bouts
        return self.nptprof - self.nptext - self.nptint

    @property
    def nptret(self):#nb points retour, Y COMPRIS le point BA et ouvext
        return self.iouvext - self.iba + 1

    @property
    def iouverture(self):
        return self.iouvext, self.iouvint
    @iouverture.setter
    def iouverture(self, (ie, ii)):
        u"""
        Afin d'éviter les erreurs, on ne peut modifier les valeurs de iouvext et iouvint que via
        >>> iouverture = (ke, ki)
        """
#         debug( iouvext=ie, iouvint=ii)
        if ie < self.iba or ii < self.iba :
            rdebug('Attention, ouverture (partiellement) sur l\'extrados, l\'echantillonnage est impossible (TODO)')
#             self._iouvext, self._iouvint = ie, ii
        elif 0 <= ie <= ii <= self.nptprof :
            pass
#             self._iouvext, self._iouvint = ie, ii
        else:
            msg = 'Valeurs incorrectes, on devrait avoir 0 <= iouvext(=%d) <= iouvint(=%d) <= nptprof(=%d)'%(ie, ii, self.nptprof)
            msg = msg+'\nTODO:en cours de saisie, ce message ne devrait pas s\'afficher'
        self._iouvext, self._iouvint = ie, ii
#             raise ValueError('Valeurs incorrectes, on devrait avoir 0 <= iouvext(=%d) <= iouvint(=%d) <= nptprof(=%d)'%(ie, ii, self.nptprof))
    @property
    def iouvext(self):
        return self._iouvext
    @iouvext.setter
    def iouvext(self, k):
        raise NotImplementedError("On ne peut pas modifier iouvext seul.\
         On doit modifier les deux points de l'ouverture en utilisant 'iouverture = (ke, ki)'.")
    @property
    def iouvint(self):
        return self._iouvint
    @iouvint.setter
    def iouvint(self, k):
        raise NotImplementedError("On ne peut pas modifier iouvint seul.\
         On doit modifier les deux points de l'ouverture en utilisant 'iouverture = (ke, ki)'.")
    @property
    def nptprof(self):
        return self._nptprof
    @nptprof.setter
    def nptprof(self, npt):
        raise ValueError("c'est ici !!! : npt=%d"%npt)
    @property
    def iba(self):
#         rdebug( "BUG : iba peut varier d'un profil a l'autre")
        return self._iba
    @iba.setter
    def iba(self, k):
        self._iba = k
#         raise NotImplementedError
############################################################
    @property
    def copy(self):
        return ProfsParam1(self.nptprof, self.iouvext,
                           self.iouvint, self.iba,
                           self.nbpbf, self.pourcentbf)
    @property
    def info(self):
        return ['<%s>'%self.__class__.__name__,
                '  nb points profil             = %d'%(self.nptprof),
                '  nb points extrados tissu     = %d'%(self.nptext),
                '  nb points extrados theorique = %d'%(self.nptexttheorique),
                '  indice BA                    = %d'%(self.iba),
                '  nb points retour             = %d'%(self.nptret),
                '  indices ouverture            = (%d,%d)'%self.iouverture,
                '  nb points ouverture          = %d'%(self.nptouv),
                '  nb points intrados tissu     = %d'%(self.nptint),
                '  nb points intrados theorique = %d'%(self.nptinttheorique),
                '  nb points reserves BF        = %d'%(self.nbpbf),
                '  pourcentage reserve BF       = %d'%(self.pourcentbf),
                ]

    def shortinfo(self):
#         return ' <%s> [nptprof=%d, iouverture=%s, iba=%d]'%(self.__class__.__name__, self.nptprof, str(self.iouverture), self.iba, )
        return ' [nptprof=%d, iouverture=%s, iba=%d, nbpbf=%d]'%(self.nptprof, str(self.iouverture), self.iba, self.nbpbf)

    def __str__(self):
        return self.shortinfo()
        return '\n'.join(self.info)

    @property
    def dump(self):
        return self.toDump()

    def toDump(self):
        todump = {'nptprof': self.nptprof, 'iouvext':self.iouvext,
                  'iouvint':self.iouvint, 'iba':self.iba,
                  'nbpbf':self.nbpbf, 'pourcentbf':self.pourcentbf}
        return todump

    def load(self,dump):
        self.__init__(**dump)


class ProfsParam(object):
    '''
    Attributs :
    ---------
    * Modifiables
        self.nptext => nb points extrados (BF -> Ouverture)
        self.nptint => nb points intrados (Ouverture -> BF)
        self.iba    => index BA
        self.iouvext => index ouverture extrados
        self.iouvint => index ouverture intrados == 1 + self.iouvext

    * NON modifiables :
        self.nptprof => nb points profil
        self.nptret  => nb points retour (BA -> iouvext)
        self.copy (property) fournit un clone de self.

    Méthodes & exemples :
    -------------------
    >>> pp1 = ProfsParam(10, 6, 7)
    >>> pp2 = ProfsParam(11, 5, 7)
    >>> pp1 == pp2 #tester si deux profparam sont identiques
    False
    >>> pp1 != pp2 #tester si deux profparam sont différents
    True
    >>> print pp1
    nb points [ext+int] : 16=[10+6]; iouverture=(9,10); iba=7, nb points retour=3
    '''
#     def __init__(self, npe=0, npi=0, iba=0):
#     def __init__(self,*args,**dump):
    def __init__(self, nptext=0, nptint=0, iba=0):
        # Valeurs par défaut : appel du type ProfsParam()
        self.nptext=nptext
        self.nptint=nptint
        self.iba=iba
        #appel du type # ProfsParam(1,2,3)
#         if len(args)>0 :
#             try :
#                 self.nptext=args[0]
#                 self.nptint=args[1]
#                 self.iba=args[2]
#             except IndexError :
#                 pass
        #Appels du type ProfsParam(npe=1, npi=2, iba=3)
#         try : self.nptext=dump.pop('npe')
#         except KeyError :pass
#         try : self.nptint=dump.pop('npi')
#         except KeyError : pass
#         try : self.iba=dump.pop('iba')
#         except KeyError : pass
#         debug(None, nptext=self.nptext, nptint=self.nptint, nba=self.iba)
    def __eq__(self,other):
        return  self.nptext==other.nptext and\
                self.nptint==other.nptint and\
                self.iba==other.iba
    def __ne__(self,other):
        return not self.__eq__(other)

    def castToProfsParam1(self):
        """retourne self sous forme d'un ProfsParam1"""
        return ProfsParam1(self.nptprof, self.iouvext, self.iouvint, self.iba,
                             ProfilPrefs.nbpbf, ProfilPrefs.pourcentbf)

    def _getnptprof(self):#nb points profil
        return self.nptext+self.nptint
    nptprof=property(_getnptprof)

    @property
    def iouverture(self):
        return self.iouvext, self.iouvint

    def _setiouvext(self,i):
        delta=i-self.iouvext
        self.nptint-=delta
        self.nptext+=delta

    def _getiouvext(self):#index ouverture extrados
        return self.nptext-1
    iouvext=property(_getiouvext,_setiouvext)

    def _setouvint(self,i):
        self.ouvext=i-1

    def _getiouvint(self):#index ouverture intrados
        return self.nptext
    iouvint=property(_getiouvint,_setouvint)

    def _getnptret(self):#nb points retour, Y COMPRIS le point BA
        return self.iouvext-self.iba+1
    nptret=property(_getnptret)

    def _getCopy(self):
        return ProfsParam(self.nptext,self.nptint,self.iba)
    copy=property(_getCopy)

    def _getInfos(self):
        return [
                'nb points [ext+int] : %d=[%d+%d]'%(self.nptprof,self.nptext,self.nptint),
                'nb points retour=%d, iba=%d'%(self.nptret,self.iba),
                'iouverture=(%d,%d)'%(self.iouvext,self.iouvint),
                ]
    info=property(_getInfos)
    def __str__(self):
        return '[nptprof=%d, iouverture=%s, iba=%d] : <%s>'%(self.nptprof, self.iouverture, self.iba, self.__class__.__name__)
        return '; '.join(self.info)

    def toDump(self):
        todump={'npe':self.nptext,'npi':self.nptint,'iba':self.iba}
        return todump
    def load(self,dump):
        self.__init__(**dump)

if __name__=='__main__' :
    pp = ProfsParam1(nptprof=10, iouvext=4, iouvint=8, iba=5, nbpbf=5, pourcentbf=15)
    pp = ProfsParam1(nptprof=50, iouvext=30, iouvint=31, iba=27, nbpbf=20, pourcentbf=15)
    print pp
    dp = pp.toDump()
    pp.load(dp)
    ppp = ProfsParam1(**dp)
    print ppp==pp
#     pp._nptprof = 60
#     print pp
#     pp.iouverture = (30, 35)
#     print pp
    import cPickle
    ppkl = cPickle.dumps(pp, protocol=0)
#     print ppkl
    ppp = cPickle.loads(ppkl)
    print 'pp==ppp ?',pp==ppp
    print 20*"="
#     p = ProfsParam1(nptprof=86, iouvext=46, iouvint=53, iba=40)
#     print 'p : ', p
#     print 'p.copy : ',p.copy
#     print 'p==p.copy ? :', p==p.copy
#     p.iouverture = 46, 80
#     print p
#     app = QApplication(sys.argv)
#     gpp = GuiProfParam(parent=None, profparam=pp)
#     print gpp.exec_()
#     print gpp.old
#     print '     ===>'
#     print gpp.profparam

#     sys.exit(app.exec_())
    exit()

