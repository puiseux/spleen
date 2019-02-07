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
import sys
import numpy as np
from lecteur import Lecteur
from utilitaires.utilitairesdivers import (find,findAll,_findRel,whoami,Path,trace,my2dPlot)
# from config import VALIDATION_DIR
from utilitaires.utilitairesdivers import debug
from dxfgrabber import dxfentities
# from geomdl import BSpline, utilities
from matplotlib import pyplot as plot
# from model.basicobjects.splinesimple import NSplineSimple
TAGS=['PTSPROFIL','CLOITROU','CLOISON','BLOCK','ENDBLK','SECTION','ENTITIES','POLYLINE','VERTEX','SEQEND','POINT','ENDSEC','EOF','CIRCLE','LINE']
POINT_TAGS=('10','20','30')
ENDLINE_TAGS=('11','21','31')
import dxfgrabber

class LecteurDXF(Lecteur):
    u"""
    Permet de lire un fichier dxf.
    Chaque fonction lireXXX renvoit un une liste de polylignes,
    ou un tableau de points ou autre, voir chaque méthodes
    """
    def __init__(self,filename,readlines=True):
        super(LecteurDXF,self).__init__(filename,readlines)
        self.markLines(TAGS)
        """
        parametres:
        ----------
        filename = (obligatoire) Un nom de fichier avec des repérages de tags, dans bookmark

        """
#        self.clearBookmark()
#        self.ignobleBricolage()

#    def ignobleBricolage(self):
#        '''
#        Je transforme les tuples de self.bookmark en liste
#        pour faire des pop() au fur et a mesure de leur lecture
#        (voir trouveNext())
#        '''
#
#        for key, value in self.bookmark.iteritems():
#            self.bookmark[key] = list(value)
#    def clearBookmark(self):
#        '''
#        On ne garde en bookmark que les numéros de ligne ou tag est present seul sur la ligne
#        par exemple une ligne comme 'TOTO_EST_CONTENT' figure dans self.bookmark si 'TOTO' est dans les tags recherchés.
#        On élimine ici une telle ligne. On ne garde que les lignes qui contiennent un tag et rien d'autre.
#        '''
#        for key, lignes in self.bookmark.iteritems():
#            lignes = list(lignes)
#            for l in lignes :
#                if self.lines[l].strip() == key : pass
#                else :
#                    trace(self, key, l, self.lines[l])
#                    lignes.remove(l)
#            self.bookmark[key] = lignes
#
    def lireElement(self,tags=POINT_TAGS):
        element=[]
        ltags=list(tags)
#        print type(self.current)
        while self.current.strip()!='0' :
            if self.current.strip() in ltags :
                ltags.remove(self.current.strip())
                self.step()
                element.append(float(self.current))
            self.step()
        if len(element)!=len(tags) :
#             trace(self, u'Lecture élément dxf sur %s : problème, ligne %d, tags=%s'%(self.filename.name, self.position, tags)
            return None
        return element

    def lireVertex(self):
        return self.lireElement()

    def lireCircle(self):
        return self.lireElement(tags=('10','20','30','40'))

    def lireDXFLine(self):
#        print POINT_TAGS+ENDLINE_TAGS
        return self.lireElement(tags=POINT_TAGS+ENDLINE_TAGS)

    def lireDXFLines(self):
        dxflines=[]
        while self.gotoNext('LINE') :
            line=self.lireDXFLine()
            if line is not None :
                dxflines.append(line)
        self.dxflines=np.asarray(dxflines)

        return self.dxflines

    def lirePolyline(self):
        '''
        self.current doit valoir 'POLYLIGNE'
        '''
        polyline=[]
        polyline.append(self.lireVertex())

        fin=self.trouveNext('SEQEND')
#         trace(self,'POLYLIGNE', self.position, 'SEQEND', fin)
        while self.gotoNext('VERTEX') is not None and self.position<fin :
            polyline.append(self.lireVertex())
        return np.asarray(polyline)

    def lirePolylines(self):
        self.rewind()
        self.polylines=[]
        for numlig in self.bookmark['POLYLINE'] :
            self.position=numlig
            self.polylines.append(self.lirePolyline())
        return self.polylines

    lirePoint=lireVertex

    def lirePoints(self):
        dim=len(self.bookmark['POINT'])
        points=np.zeros((dim,3))
        n=0
        while self.gotoNext('POINT') is not None :
            point=self.lirePoint()
            points[n]=point
            n+=1
        self.points=points
        return points

    def lire(self):
        self.lirePolylines()
        self.lirePoints()
        self.lireDXFLines()

    run=lire

    def trouve(self,tag,n=0):
        """ne change pas de position
        Trouve la n-ieme occurence de tag ou None, """
        try :
            return self.bookmark[tag][n]
        except (KeyError,IndexError) :
            return None

    def trouveNext(self,tag):
        """ne change pas de position
        Trouve la prochaine occurence de tag ou None"""
        pos=self.position
        position=None
        for p in self.bookmark[tag] :
            if p>pos :
                position=p
                break
        return position

    def goto(self,tag,n=0):
        """Se place à la n-ieme occurence de tag (si existe)"""
        self.position=self.trouve(tag,n)
        return self.position

    def gotoNext(self,tag):
        """Se place à la prochaine occurence de tag (si existe)"""
        self.position=self.trouveNext(tag)
        return self.position


class LecteurDXFNervures(LecteurDXF):
    """
    Un fichier avec des repérages de tags, dans bookmark
    l'argument filename est obligatoire
    Uniquement pour les fichiers comme 'Faial2v5A.DXF'
    """
    def __init__(self,filename,readlines=True):
        super(LecteurDXF,self).__init__(filename,readlines)
        self.markLines(['PTSPROFIL','CLOITROU','CLOISON','ENTITIES','POINT'])

    def lire(self):
        self.goto('ENTITIES')
        self.points=self.lirePoints()[:-1,:-1]
        return self.points

    def plot(self):
        my2dPlot((self.points,))

class LecteurDXF0(object):
    def __init__(self,filename):
        self.dwg=dxfgrabber.readfile(filename)
        self.entities=self.dwg.entities.get_entities()
        self.polylines=[entity for entity in self.entities if entity.dxftype=='POLYLINE']
        self.lines=[entity for entity in self.entities if entity.dxftype=='LINE']
        self.dxfpoints=[entity for entity in self.entities if entity.dxftype=='POINT']
        self.vertex=[entity for entity in self.entities if entity.dxftype=='VERTEX']

    def _getPoints(self):
        points=[]
        for polyline in self.polylines :
#             trace(self, len(polyline))
            points+=[point for point in polyline.points]
        return np.asarray(points)

    def lire(self):
        pass
    points=property(_getPoints)

class LecteurDXF1(object):
    def __init__(self,filename):
        self.dxf = dxfgrabber.readfile(filename)
#         debug(self.dxf)
        self.entities=self.dxf.entities.get_entities()
        self.polylines=[entity for entity in self.entities if entity.dxftype=='POLYLINE']
        self.lines=[entity for entity in self.entities if entity.dxftype=='LINE']
        self.dxfpoints=[entity for entity in self.entities if entity.dxftype=='POINT']
        self.vertex=[entity for entity in self.entities if entity.dxftype=='VERTEX']
        self.splines=[entity for entity in self.entities if entity.dxftype=='SPLINE']
        self.types = set([entity.dxftype for entity in self.dxf.entities])
    @property
    def points(self):
        points=[]
        for polyline in self.polylines :
#             trace(self, len(polyline))
            points+=[point for point in polyline.points]
        return np.asarray(points)

    def lire(self):
        pass
    

###### tests #####

if __name__=="__main__":
    print '*******Lecteur DXF1******'
    filename = Path(VALIDATION_DIR,'decorations','Bringhen_D_cmyk gross.dxf')
    fl=LecteurDXF1(filename)
    dxf = fl.dxf
    fl.lire()
    #print fl.points#, fl.vertex, fl.lines, fl.polylines
    print("DXF version: {}".format(dxf.dxfversion))
    print(dxf.header)
    print(dxf.layers)
    types = [entity.dxftype for entity in dxf.entities]

    print(set(types))
    polylines = [entity for entity in dxf.entities if entity.dxftype==u'POLYLINE']
    lwpolylines = [entity for entity in dxf.entities if entity.dxftype==u'LWPOLYLINE']
    splines = [entity for entity in dxf.entities if entity.dxftype==u'SPLINE']
#     pol = polylines[0]
#     for p in pol : print p
#     spl = splines[0]
#     for p in spl.control_points : print(p)
    points = []
    from matplotlib import pyplot as plot
    for k, pol in enumerate(lwpolylines) :
        pp = np.asarray([p[:-1] for p in pol])
        plot.plot(pp[:,0],pp[:,1],'r')
    
    for k, pol in enumerate(polylines) :
        pp = np.asarray([p[:-1] for p in pol])
        plot.plot(pp[:,0],pp[:,1],'g')
    
    for k, spl in enumerate(splines) :
#         if k>0 : continue
        print('\n***   k=%d'%k)
        print('degree',spl.degree)
        print('tangent',spl.start_tangent,spl.end_tangent)
        print('fit_points',spl.fit_points)
        print('knots',spl.knots)
        print('normal_vector',spl.normal_vector)
        print('control_points',spl.control_points)
        print('flags',spl.flags)
        pp = np.asarray([p[:-1] for p in spl.control_points])
#         ppspl = NSplineSimple(points=spl.control_points, methode=('ius',3))
#         print ('ppspline:',ppspl)
        plot.plot(pp[:,0],pp[:,1],'b')
    plot.show()

    exit()
    print '*******Lecteur DXF0******'
    filename = Path(VALIDATION_DIR,'dxf','blocjonc.dxf')
    filename = Path(VALIDATION_DIR,'decorations','Bringhen_D_cmyk gross.dxf')
    fl=LecteurDXF0(filename)
    fl.lire()
    print fl.points#, fl.vertex, fl.lines, fl.polylines
    exit()
    print '*******Lecteur DXF******'
    filename=Path(VALIDATION_DIR,'dxf','cot_d03.DXF')
    filename=Path(VALIDATION_DIR,'dxf','teste01.dxf')
    filename=Path(VALIDATION_DIR,'faial2.dxf')#3d, je sais pas lire...
    filename=Path(VALIDATION_DIR,'dxf','Faial2v5A.DXF')

    fl=LecteurDXFNervures(filename)
    fl.lire()
    print fl.points.shape
    fl.plot()#????? marche pas...
    exit()
    f=LecteurDXF(filename)
#    f.lire()
    for key,value in f.bookmark.iteritems() :
        print key,' : ',len(value)
#        if key is 'LINE' : print value
#    print f
    f.rewind()
    if 1 :
        print f.goto('ENTITIES')
        XY0=f.lirePoints()[:-1,:-1]

    if 0 :
        print f.gotoNext('LINE')
        print f.lireDXFLine()
        print f.position
        print f.lireElement(POINT_TAGS)
        print f.position
        print f.gotoNext('LINE')
        print f.gotoNext('LINE')
        f.rewind()
        print f.lireDXFLines()
#        exit()
    if 0 :
        f.rewind()
        f.gotoNext('CIRCLE')
        print f.lireCircle()

        f.gotoNext('CIRCLE')

        print f.lireCircle()
        f.gotoNext('CIRCLE')
        print f.lireCircle()
        f.lirePoints()

    if 0 :
        f.rewind()
        print f.trouveNext('VERTEX')
        print f.gotoNext('VERTEX')
        print f.position
        print f.gotoNext('VERTEX')
        print f.position
        print f.lireVertex()
        f.rewind()
        print f.gotoNext('POLYLINE')
        print f.lirePolyline()
    if 0 :pass
#    f.rewind()
#    print f.lirePolylines()
#    print f.lirePoints()
#        f.rewind()
#        f.lire()

#    XY0 = f.polylines[0][:,:-1]
#    XY1 = f.polylines[1][:,:-1]
#    XY2 = f.polylines[2][:,:-1]
#    XY3 = f.points[:,:-1]
#    XY0 = f.dxflines[:,0:2]
#    print f.dxflines.shape, XY0.shape
#    print XY0[:,1]
#    XY5 = f.dxflines[:,3:-1]
    if 1 :
        my2dPlot((XY0,))#), XY2, XY3,XY4,XY5))
