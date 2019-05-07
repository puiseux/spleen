#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe xxx
Description :
@module : lecteurdxf
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 11 mars 2013
'''
import numpy as np
from lecteur import LecteurGnuplot,LecteurMatlab
from lecteurdxf import LecteurDXFNervures,LecteurDXF0, LecteurDXF1
from utilitairesdivers import (Path,trace,my2dPlot,rdebug)
import cPickle


def pointsFromFile(filename):
    filename = unicode(filename)
    filename = Path(filename)
    ext = filename.ext.lower()
#     debug(filename=filename, ext=ext)
    if ext in ('.dxf',) :
        # raise NotImplementedError
        lecteur = LecteurDXFNervures(filename)
        polygon = lecteur.lire()
        return polygon
    elif ext in ('.gnu',) :
        # raise NotImplementedError
        lecteur = LecteurGnuplot(filename,patch=False)
        lecteur.dim = 2
        lecteur.patch = False
        polygon = lecteur.lire()
        return polygon
    elif ext in ('.pts',) :
        raise NotImplementedError
        # lecteur = LecteurNervures(filename)
        # polygon = lecteur.lire()
        # if lecteur.genre.lower()=='voile' :
        #     raise IOError(whoami()+u"%s n'est pas un polygon2d mais 3d (voile?, profil3d?)"%filename.name)
        # elif lecteur.dimension == 3 :
        #     polygon = lecteur.points[:,:-1]
        # return polygon
    elif ext in ('.spl',) :
        """Une spline"""
        with open(filename,'r') as f :
            dump = eval(f.read())
        if isinstance(dump, dict):
            if dump.has_key('cpoints') :
                return np.asarray(dump['cpoints'])
            elif dump.has_key('dmodel') :
                return np.asarray(dump['dmodel']['cpoints'])
            else :
                raise IOError(u"je ne sais pas extraire les points de ce fichier : %s"%filename.name)
        else :
            raise IOError(u"%s est un fichier '.spl' ; il devrait contenir un dictionnaire (dict Python) bien constitue."%filename.name)
    elif ext in ('.pkl',) :
        u"""peut-etre Une spline format pkl"""
        with open(filename,'r') as f :
            dump = cPickle.load(f)
        if isinstance(dump, dict):
            if dump.has_key('cpoints') :
                return np.asarray(dump['cpoints'])
            elif dump.has_key('dmodel') :
                return np.asarray(dump['dmodel']['cpoints'])
            else :
                raise IOError(u"je ne sais pas extraire les points de ce fichier : %s"%filename.name)
        else :
            raise IOError(u"%s est un fichier '.pkl' ; je ne parviens pas à en extraire des points 2d."%filename.name)
    else :
#         rdebug()
        raise IOError(u'Je ne sais pas lire ce fichier "%s"'%filename)
        # return np.zeros((0,2))

def pointsFrom(x):#,close_=False):
    u'''
    Paramètre x:
        - un np.ndarray((n,2))
        - un chaine de caractères evaluable par eval(x)
        - un nom de fichier
    Retourne :
    --------
    - points : np.ndarray, de shape (N,2) si close_==False
    '''
    # debug(x)
    if x is None :
        return np.zeros((0,2))
    elif isinstance(x,(Path,str)):
        return pointsFromFile(x)#Fichier
    elif isinstance(x, (str, unicode)) :
        return pointsFrom(eval(x))#chaine de caracteres p.ex.'[1,2,5,12]'
    elif isinstance(x, (list, tuple)) :
        npa = np.asarray(x)
#         debug(shape=npa.shape)
        sh = npa.shape
        if len(sh) != 2 or sh[1] != 2 :
            npa.shape=(len(npa)/2,2)
        return npa
    elif isinstance(x, np.ndarray) :
        return x
    else :
        raise TypeError (u"pointsFrom() : cas non pris en charge")


class LecteurUniversel(object):
    """
    Pour lire les .dxf, .pts, .gnu et .m (matlab ou octave)
    Le constructeur lit (self.lire()) le fichier et alimente self.points
    """
    def __init__(self,filename,readlines=True):
        super(LecteurUniversel,self).__init__()
        rdebug(u'LecteurUniversel obsolète. Utiliser le lecteur adapté \n(LecteurDXF, LecteurDXF0, LecteurDXF1,LecteurNervures,LecteurSVG...\n')
        filename=Path(filename)
        ext=filename.ext.lower()
        self.ext=ext
        self.filename=filename
        if 'dxf' in ext :
            try :
                self.lecteur=LecteurDXF0(filename)
                self.lire()
            except (KeyError, AssertionError) :
                try :
                    self.lecteur=LecteurDXF1(filename)
                    self.lire()
                except (KeyError,AssertionError) :
                    try :
                        self.lecteur=LecteurDXFNervures(filename)
                        self.lire()
                    except (KeyError,AssertionError) :
                        trace(self, 'impossible de lire %s'%filename.name)

        elif 'pts' in ext :
            raise NotImplementedError
            # self.lecteur=LecteurNervures(filename)
            # self.lire()
        elif 'gnu' in ext :
            try :
                self.lecteur=LecteurGnuplot(filename,patch=True)
                self.lire()
            except ValueError :
                self.lecteur=LecteurGnuplot(filename,patch=False)
                self.lire()

        elif 'mtl' in ext or 'm' in ext or 'oct' in ext:
            self.lecteur=LecteurMatlab(filename)
            self.lire()
        elif 'svg' in ext:
            from utilitaires.lecteursvg import LecteurSVG#On ne peut pas faire l'import dans l'en-tête
            self.lecteur=LecteurSVG(filename)
            self.lire()

        else :
#             self.lecteur = None
#             self.points = []
            trace(self,'extension = %s non prise en charge'%self.ext)

    @property
    def points(self):
        return self.lecteur.points

    def lire(self):
        try : self.lecteur.lire()
        except AttributeError : return

    def plot(self):
        try : self.lecteur.plot()
        except AttributeError :
            shape=self.points.shape
            if len(shape)==3 :# 3D : (nbnerv, nbp, 3)
                trace(self,'3d : shape=%s'%str(shape))
#                 trace(self, self.points[0,:,0:-1].shape)
                my2dPlot((self.points[0,:,0:-1],self.points[1,:,0:-1]))
            elif len(shape)==2 : #2D
                trace(self,'2d : shape=%s'%str(shape))
                my2dPlot((self.points,))



###### tests #####

if __name__=="__main__":
    from config import VALIDATION_DIR
    print '*******Lecteur DXF******'
    filenames =[
                Path(VALIDATION_DIR,'decorations','EuroSport.svg'),
                Path(VALIDATION_DIR,'dxf','cot_d03.DXF'),
                Path(VALIDATION_DIR,'dxf','DecoITVBi.DXF'),
                Path(VALIDATION_DIR,'dxf','teste01.dxf'),
                Path(VALIDATION_DIR,'faial2.dxf'),#3d, je sais pas lire...
                Path(VALIDATION_DIR,'dxf','Faial2v5A.DXF'),
                Path(VALIDATION_DIR,'dxf','blocjonc.dxf'),
                Path(VALIDATION_DIR,'formatpts','diamir.pts'),
                Path(VALIDATION_DIR,'formatpts','ARBIZON.PTS'),
                Path(VALIDATION_DIR,'1.gnu'),
                Path(VALIDATION_DIR,'simple','simple.gnu'),
                Path(VALIDATION_DIR,'simple','simple2d.gnu'),
                ]
    for filename in filenames :
        print '****************',filename.name,'****************'
        fl=LecteurUniversel(filename)
        print type(fl.lecteur).__name__
        print str(fl.lecteur)
#     fl.lire()
        try :
            print 'points.shape : ', fl.points.shape
            fl.plot()
        except Exception as msg:
            print type(msg).__name__, msg
            print '*** Impossible de tracer'
            pass

    exit()
