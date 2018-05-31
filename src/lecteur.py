#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Created on 17 mai 2012

@author: puiseux
'''
import sys
from path import Path
import numpy as np
from utilitaires import (find,findAll,findAllLines,_findRel,whoami,debug,rdebug)
from config import VALIDATION_DIR
TAGS=[]

class Lecteur(file):
    """
    Un fichier avec des repérages de tags, dans bookmark
    l'argument filename est obligatoire
    """
    def __init__(self,filename,readlines=True):

#        print>>sys.stderr, "%%%", filename
        file.__init__(self, filename, mode='r')
        self.bookmark={}
        self.position=None
        self.filename=Path(filename)
        if readlines :
            self.lines=self.readlines()
            self.close()

    def __len__(self):
        return len(self.lines)

    def lire(self):
        sys.stdout.flush()
        rdebug(" methode 'virtuelle' non implementee")
        return None
    run=lire
    
    @property
    def classname(self):
        return self.__class__.__name__

    @property
    def name(self):
        return self.filename.name

    @property
    def ext(self):
        return self.filename.ext

    def __str__(self):
#        return str(self.bookmark)
        return '\n'.join(['keys : ',str(self.bookmark.keys()),'values :' ,str(self.bookmark.values())])

    def __getitem__(self,*args,**kargs):
        return self.bookmark.__getitem__(*args,**kargs)

    def keys(self):
        return self.bookmark.keys()
    def items(self):
        return self.bookmark.items()

    def markLines(self,tags,deb=0,fin=-1,order=False):
        '''Comme self.mark, mais seules sont prises en compte les lignes EGALES à tags,
        (au lieu des lignes CONTENANT tags)'''
        if order :
            self.orderMark(tags,deb,fin)
        else :
            for tag in tags :
                self.bookmark[tag]=self.findAllLines(tag,deb,fin)
        return self.bookmark


    def markFirstOccurenceOfTags(self,tags,deb=0,fin=-1,order=False):
        u"""Trouve le numero de ligne (dans lignes) de la PREMIERE occurence des chaines de caracteres donnees dans tags.
        Alimente ensuite le dictionnaire self.bookmark
        parametres :
        ----------
        -tags=une liste de tags (chaines de caracteres) à reperer dans le fichier
        -deb, fin=numeros de lignes entre lesquels est limitée la recherche.
        """
        if order :
            self.orderMark(tags,deb,fin)
        else :
            for tag in tags :
                self.bookmark[tag]=self.findAll(tag,deb,fin)
        return self.bookmark
    mark = markFirstOccurenceOfTags

    def markFirstOccurenceOfOrderedTags(self,tags,deb=0,fin=-1):
        """
        Comme markFirstOccurenceOfTags, mais dans l'ordre des tags (pour optimiser la recherche)
        parametres:
        ----------
            tags, une liste de chaines de caracteres.
        fonctionnement:
        --------------
        Parcours de  self.lignes[deb:fin] pour trouver le numero de ligne
        de la premiere ocurrence de chaque tag,
        Les tags sont recherches DANS L'ORDRE donné par tags :
            si la première occurrence de tags[i] est trouvée en ligne n,
            alors tags[i+1] est recherché dans self.lignes[n:fin]
        """

        for tag in tags :
            n=self._findRel(tag,deb,fin)
            if n>=0 :
                self.bookmark[tag]=n+deb
                deb=n
        return self.bookmark
    orderMark = markFirstOccurenceOfOrderedTags

    def findAll(self,tag,first=0,last=-1):
        """
        Retourne une liste triée des numeros de ligne de
        TOUTES les occurences de tag dans lines[first:last]
        """
        linesn=findAll(tag,self.lines,first,last)
        try : self.position=linesn[-1]
        except IndexError : self.position=None
        return linesn

    def findAllLines(self,tag,first=0,last=-1):
        linesn=findAllLines(tag,self.lines,first,last)
        try : self.position=linesn[-1]
        except IndexError : self.position=None
        return linesn

    def _findRel(self,tag,first=0,last=-1):
        """
        Trouve tag dans lines[first:]
        ATTENTION, retourne le numéro de ligne RELATIF à first
        Pour avoir le numéro de ligne absolu, il faut lui rajouter 'first'
        Par contre, self.position est mis a jour correctement (numéro absolu)
        """
        line_number=_findRel(tag,self.lines,first,last)
        self.position=first+line_number
        return line_number

    def find(self,tag,first=0,last=-1):
        """
        Trouve tag dans lines[first:]
        ATTENTION, retourne le numero de ligne relatif à 'first' ???
        """
        self.position=find(tag,self.lines,first,last)
        return self.position

    def goToNextTag(self,tag):
        """Place le curseur à la prochaine valeur de tag ou sur None"""
        try :
            return self.find(tag,1+self.position)
        except TypeError :
            self.position=None
            return None
    findNext = goToNextTag

    def goToTag(self,tag,n=0):
        """Se place à la n-ieme occurence de tag (si existe)"""
        try : self.position=self.bookmark[tag][n]
        except (KeyError, IndexError) : self.position=None
        return self.position
    goto = goToTag

#     def seek(self,n):
#         """Se place à la n-ieme ligne"""
#         self.position=n
#         return self.position

    def advance(self,n=1):
        """Avance de n"""
        if self.position is not None : self.position+=n
        return self.position

    def step(self):
        """Avance de 1"""
        return self.advance()

#    def bloc(self, first, last):
#        try : return self.lines[first:last]
#        except IndexError : return None
#
    @property
    def current(self):
        try : return self.lines[self._position]
        except (KeyError,TypeError) : return None
#     current=property(_current)
    current_line=current

    def context(self,n=3):
        """Renvoit le contexte (2n lignes) : n ligne avant et n-1 lignes apres la position actuelle"""
        try : return self.lines[self._position-n:self._position+n]
        except TypeError : return 'Pas de contexte (fin de fichier-EOF ?)'

    @property
    def position(self):
        return self._position
    @position.setter
    def position(self,n):
        if 0 <= n < len(self.lines) :
            self._position=n
        else :
            self._position=None
    current_line_number = position
#     position=property(_getPosition,_setPosition)

    def _setTags(self,tags):
        for tag in tags :
            self.bookmark[tag]=None

    def _getTags(self):
        return sorted(self.bookmark.keys())
    tags=property(_getTags,_setTags)

    def rewind(self):
        self.position=0

class LecteurGnuplot(Lecteur):
    """
    Parametres:
    ----------
    filename : str ou unicode ou QString, nom de fichier existant
    patch : booleen
        - True (défaut) => la méthode lire() retourne un PatchCmarc
        - False => la méthode lire() retourne un np.ndarray
    dim : int (3 par défaut) nombre de coordonnées par point

    Méthodes :
    --------
    lire() : retourne un PatchCmarc si le paramètre patch est True, un np.ndarray sinon.
        Le PatchCmarc est shapé a (nn,n,3) = nb nervures, nb points par nervure, dimension
    goodLine()
    """
    def __init__(self,filename,patch=True):
        Lecteur.__init__(self,filename,readlines=True)
        self.patch=patch
        self.dim=3
#        self.lire()

    def badLine(self,num):
        return not self.goodLine(num)

    def goodLine(self,num):
        line=self.lines[num].strip()
        try :
            point=[float(w) for w in line.split()]
            if len(point)==self.dim :
                return True
            else :
                return False
        except ValueError:
            return False

    def lire(self):
        """
        Fichier gnuplot :
        ===============
                # commentaire
                # encore commentaire

                1. 2. 3.
                1. 2. 3.
                1. 2. 3.
                1. 2. 3.
                #nouvelle nervure

                1. 2. 3.
                1. 2. 3.
                1. 2. 3.
                1. 2. 3.

                #autre nervure.
                # Les lignes vides ou commentaires contigues sont
                # considérées comme une seule ligne séparant deux nervures.
                1. 2. 3.
                1. 2. 3.
                1. 2. 3.
                1. 2. 3.

        Fonctionnement
        ==============
            1. petite gymnastique pour trouver le nombre de nervures.
            2. lecture brutale par numpy
            3. reshape et calcul connectivités.
            Retourne un PatchCMARC
        """
        for nl,line in enumerate(self.lines) :
            if self.goodLine(nl) : break
        firstline=nl

        NN=0
        for nl,line in enumerate(self.lines[firstline:]) :
            if self.badLine(nl) and self.goodLine(nl+1) :
                NN+=1

        #print NN
        self.points = np.loadtxt(self.filename)
#         print self.points
        if not self.patch : return self.points
        self.points.shape = NN, -1, 3
#         self.voile = PatchCMARC(self.points,name=self.filename.name)
        return self.points

class LecteurMatlab(Lecteur):pass

class LecteurSommaire(Lecteur):
    '''
    lit toute les lignes qui ne commencent pas par un #, retourne la liste des lignes sous forme liste de string.
    '''
    def __init__(self,filename,readlines=True):
        super(LecteurSommaire,self).__init__(filename,readlines=True)

    def lire(self):
        lines=[line.strip() for line in self.lines if line.strip() and line[0]!='#']
        self.lines=lines
        return self.lines


if __name__=="__main__":
#    exit()
#    filename = "aero.dat"
#    filename = "elastic.dat"
#    filename =  Path(datadir,'DeuxNervures.txt')
    # from vtk_nervures.grids import PatchCMARC

    print '*******Lecteur******'
    filename=Path(VALIDATION_DIR,'dxf','cot_d03.DXF')
    f=Lecteur(filename)
    f.mark(['LINE'])
    print f.bookmark
#     exit()
#     print '*******LecteurSommaire******'
#     filename=Path(VALIDATION_DIR,'lissage4.dat')
#     f=LecteurSommaire(filename)
#     print f.lire()
#     filename=Path(VALIDATION_DIR,'lissage3.dat')
#     f=LecteurSommaire(filename)
#     print f.lire()
    print '*******LecteurGnuplot******'
    datadir=Path(VALIDATION_DIR,'simple')
    filename=Path(VALIDATION_DIR,'simple','simple2d.gnu')
    filename=Path(VALIDATION_DIR,'1.gnu')
    f=LecteurGnuplot(filename,patch=False)
    print 'filename=%s'%filename
    print f.points#lire()
    exit()
    filename=Path(VALIDATION_DIR,'simple','simple.gnu')
    f=LecteurGnuplot(filename,patch=False)
    print f.lire()
    f=LecteurGnuplot(filename)
    print '***',f.lire()
    filename=Path(VALIDATION_DIR,'simple','simple.aerodynamique3d.cmo')
    print '*******Lecteur******'
    f=Lecteur(filename)

    f.mark(["toto","STREAMLINE","DOWNSTREAM"])
    for tag,value in sorted(f.bookmark.items()) :
        print tag,' : ',value
    print f.position
    next_=f.findNext("WAKE")
    if next_ is not None : print f.lines[next_]
    f.mark(["STREAMLINE","DOWNSTREAM","toto"])
    for tag,value in sorted(f.bookmark.items()) :
        print tag,' : ',value
#    print f.voile
#
    tag="WAKE NUMBER"
    next_=f.find(tag)
    while next_ is not None :
        print f.position,next_,f.lines[next_]
        next_=f.findNext(tag)
    print f.current
    print f.context()




