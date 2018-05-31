#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures
Classe xxx
Description :
@module : lecteursvg
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 11 mars 2013
'''
import sys
import numpy as np
from lecteur import Lecteur
import utilitaires
from utilitaires import (my2dPlot,Path)
from config import VALIDATION_DIR
from utilitaires.utilitaires import debug, rdebug, hardScale, p2a
# from dxfgrabber import dxfentities
# from geomdl import BSpline, utilities
from matplotlib import pyplot as plt
from model.basicobjects.splinesimple import NSplineSimple
from model.basicobjects.polyline import NPolyLine, NPolygone
from svgpathtools import (svg2paths, svg2paths2, wsvg, kinks, smoothed_path,
                          disvg, Line, CubicBezier, QuadraticBezier, Arc,
                     )
from svgpathtools import Path as SVGPath
from time import sleep
from PyQt4.QtCore import QRectF
from gui.graphicsbase.graphicscommon import p2t

def continuousPathToArray(path, nbp=20):
    u"""
    Transforme un path svg CONTINU en np.ndarray de points.
    Paramètres:
    ----------
        - path est un path svg, continu, c'est à dire une liste
            de segments svg (Line, CubicBezier, QuadraticBezier, Arc)
            dont les segments vérifient :
            path[i-1].end == path[i].start, i=1,...,len(path)-1
        - nbp est le nombre de points à echantillonner sur la chaque segment
            de type autre que Line,
            pour les Lines, le nb de points est 2 (start et end points)
            valeur par défaut nbp = 20

    Retourne :
    --------
        un tableau np.ndarray contenant un echantillon raisonable de points de path.
    """
    assert (path.iscontinuous()),"Path non continu"
    points = []
    for seg in path :
        if isinstance(seg, Line) :
            points.append(seg.start)
        else :
            delta = 1.0/nbp
            T = np.arange(0, 1, delta)
            points.extend([seg.point(t) for t in T])
    points.append(seg.end)
    points = [(p.real, p.imag) for p in points]
    return points

def _continuousPathToSpline(path, nbp=20, methode=('ius',1)):
    u"""
    Transforme un path svg CONTINU en NSplineSimple.
    Paramètres:
    ----------
        - path est un path svg, continu, c'est à dire une liste
            de segments svg (Line, CubicBezier, QuadraticBezier, Arc)
            dont les segments vérifient :
            path[i-1].end == path[i].start, i=1,...,len(path)-1
        - nbp est le nombre de points de contrôle sur la chaque segment de type autre que Line,
            pour les Lines, le nb de points est 2 (start et end points)
            valeur par défaut nbp = 20
        - methode est le type de spline retourné,
            valeur par défaut=('ius',1) c'est à dire un Polyligne

    Retourne : une spline du type precisé par methode, qui interpole les points de path.
    --------
    """
    points = continuousPathToArray(path, nbp)
    return NSplineSimple(points=points, methode=methode)

def continuousPathToPolyligne(path, nbp=20):
    points = continuousPathToArray(path, nbp)
    return NPolyLine(points=points)

def _pathToSplines(path, nbp=20, methode=('ius',1)):
    u"""
    Transforme un path svg quelconque en autant de NSplineSimple
        que path contient de sous-paths continus.
    Paramètres:
    ----------
        - path est un path svg, c'est à dire une liste de segments svg
            (Line, CubicBezier, QuadraticBezier, Arc)

        - nbp est le nombre de points de contrôle des splines résultantes
            sur les segments de type autre que Line,
            pour les segments de type Line, le nb de points est 2 (start et end points)
            valeur par défaut nbp = 20

        - methode est le type des splines retournées,
            valeur par défaut=('ius',1) c'est à dire des Polylignes

    Retourne :
    --------
    une liste de splines du type precisé par methode,
    chaque spline interpole les points du sous-path continu.
    """
    cont_paths = path.continuous_subpaths()
    return [_continuousPathToSpline(path, nbp, methode) for path in cont_paths]

def pathToArray(path, nbp=20):
    u"""
    Echantillonne un path svg quelconque en autant de np.ndarray(N,2)
        que path contient de sous-paths continus.
    Paramètres:
    ----------
        - path est un path svg, c'est à dire une liste de segments svg
            (Line, CubicBezier, QuadraticBezier, Arc)

        - nbp est le nombre de points de contrôle des splines résultantes
            sur les segments de type autre que Line,
            pour les segments de type Line, le nb de points est 2 (start et end points)
            valeur par défaut nbp = 20

    Retourne :
    --------
    une liste de tableaux de points (np.ndarray(N,2)),
    chaque tableau de points est un échantillon de points du sous-path continu.
    """
    cont_paths = path.continuous_subpaths()
    return [continuousPathToArray(path, nbp) for path in cont_paths]

def pathToPolylignes(path, nbp=20):
    return [continuousPathToPolyligne(path, nbp)
            for path in path.continuous_subpaths()]

class LecteurSVG(Lecteur):
    u"""
    Permet de lire un fichier svg et de formatter son contenu sous une forme
    utilisable par Axile.
    Essentiellement un fichier SVG contient des "paths" SVG, c'est à dire
    des balises XML :
        <path> id="..." style="..." d="..." </path>
    où
    - id = un identificateur
    - style = couleur, remplissage,...
    - d = description du path sous forme condensée
    Ces paths sont lus et stockés en objets de type SVGPath (= svgpathtools.Path)
        consistant essentiellement en une liste de Line, QuadraticBezier,
        CubicBezier, et Arc et les attributs (id, style, d)
    Les paths ne sont pas forcément continus.
    Un path (de type SVGPath) non continu peut être décomposé en une liste de
        subpaths (de type SVGPath) continus avec la
        méthode SVGPath.continuous_subpaths().
    A chaque path continu, on associe, pour Axile (décoration ou autre) :
    - un tableau de points np.ndarray (méthode points()) contenant
        - soit des points de contrôle pour une spline
        - soit un échantillonnage "raisonnable" de points du path
    - le style du path
    - l'id du path
    """
    def __init__(self,filename):
        super(LecteurSVG, self).__init__(filename)
        self.paths, self.attributes, self.svg_attributes = svg2paths2(filename)
#         print len(self.paths), type(self.paths)#, len(self.attributes), len(self.svg_attributes)
#         print self.paths[0]
        self._parseAttributes()
#         print self.attributes
#         print self.svg_attributes
        #
    def __getitem__(self, k):
        """Le k-ieme path"""
        return self.paths[k]

    def __len__(self):
        """Nombre de paths"""
        return len(self.paths)

    def id(self,k):
        u"""nom du path"""
        return self.attributes[k][u'id']
    #name=id
    def fill(self,k):
        u"""couleur remplissage"""
        return self.attributes[k][u'fill']
    couleur=fill
    def fills(self):
        u"""Liste des couleurs des self.paths (continus ou non)"""
#         return [self.fill(k) for k in range(len(self))]
        return [attr[u'fill'] for attr in self.attributes]
    couleurs=fills
#     def fill_opacity(self,k):pass
#     def fill_rule(self,k):pass
    def stroke(self,k):
        u"""couleur contour"""
        return self.attributes[k][u'stroke']
    contour=stroke
    def strokes(self):
        u"""Liste des couleurs contours, des self.paths (continus ou non)"""
#         return [self.stroke(k) for k in range(len(self))]
        return [attr[u'stroke'] for attr in self.attributes]
    contours=strokes

    def resize(self, scale, center, translate=None):
        """translation, puis homothétie de rapport scale=sx,sy, de centre center"""
        raise NotImplementedError, 'TODO'
        for path in self.paths :
            pass
    def _parseAttributes(self):
        u"""On ne garde que
        - 'id' => nom
        - 'fill' => couleur du fond
        - 'stroke' => couleur du tracé
        Pour le style d'un path donné, on a parfois une proprieté unique
        style="fill:#1c63b7;fill-opacity:1;fill-rule:nonzero;stroke:none"
        parfois plusieurs propriétés sous forme
        fill="green" fill-opacity="0.5" stroke="black" stroke-width="2"
        """
        for k, dattr in enumerate(self.attributes) :
            new_attr = {u'id':u'path',u'fill':u'none',u'stroke':u'none',u'style':u'none',u'd':u'none'}#defaut
            if dattr.has_key(u'style') :
                style = dattr[u'style']
                for f in style.split(';'):
                    e = f.split(':')
                    if e[0].strip() == 'fill':
                        new_attr[u'fill'] = e[1].strip()
                    elif e[0].strip() == 'stroke':
                        new_attr[u'stroke'] = e[1].strip()

            for key in (u'id', u'fill', u'stroke', u'd', u'style') :
                if dattr.has_key(key) :
                    new_attr[key] = dattr[key]
            self.attributes[k] = new_attr

    def boundingRect(self, k=None):
        if k is not None :
            return self[k].bbox()
        else :
            Ls = []
            for path in self.paths :
                if path :
                    Ls.extend(path._segments)
            P = SVGPath(*Ls)
            return P.bbox()
#             bbs = [path.bbox() for path in self.paths]
#             xmins, xmaxs, ymins, ymaxs = list(zip(*bbs))
#             xmin = min(xmins)
#             xmax = max(xmaxs)
#             ymin = min(ymins)
#             ymax = max(ymaxs)
#             #(0.0, 3953.78, 0.0, 670.039)
#             return xmin, xmax, ymin, ymax

    def points(self, k=None, nbp=5):
        u"""
        La liste self.paths contient N paths Pi=self[i], 0 <= i < N.
        Chaque path Pi peut être continu ou non. S'il n'est pas continu,
            on le décompose en subpaths continus PC[j], 0 <= j < 0,1,2...Nc.
            De sorte que l'on ne traite que de paths continus
        Chaque path, Pi ou PCj est constitué de "segments" au sens svg, de type
            Line, QuadraticBezier, CubicBezier, Arc.
        Pour un path continu, les segments S[l] 0 <= l <Ns sont contigus i.e.
            S[l-1].end == S[l].start pour 0<l<n
        À chaque (sub)path CONTINU, on associe un polyligne
            constitué de points échantillonnés sur le path sous forme d'un
            tableau de points np.ndarray de shape (xx,2).
        Pour une Line on prélève ses deux points start et stop
        Pour un segment de type QuadraticBezier, CubicBezier, Arc on prélève
            nbp points régulièrement espacés sur le segment.


        Retourne
        - si k>=0, la liste des np.ndarray des subpaths de self[k]
        - si k==-1, la liste à plat de tous les np.ndarray() de tous les subpaths.
        - si k==None, la même chose que pour k==-1, sauf que la liste n'est pas à plat,
            elle a une longueur de len(self), et son i-ème élément est
            la liste de tous les subpath de self[i]
        """
        if k is None :
            points = [pathToArray(path, nbp) for path in self.paths]
#         elif k == -1 :
#             points = []
#             for path in self.paths :
#                 points.extend(pathToArray(path, nbp))
        else :#k entier
            points = pathToArray(self[k], nbp)
        return points


    def polylignes(self, k=None, nbp=5):
        u"""
        Le k-ième path sous forme d'une liste de polylignes.
        le k-ieme path peut être continu ou non.
        Il est composé de "segments" de type
            Line, QuadraticBezier, CubicBezier, Arc.
        S'il a nc composantes continues, le nombre de polylignes retourné est nc
        nbp est le nombre de points de contrôle pour chaque "segment" de type
        QuadraticBezier, CubicBezier, et Arc (i.e. autre que Line)
        si k==-1, retourne liste à plat de tous les polylignes.
        si k==None, retourne la liste de tous les paths de self,
        où chaque path est sous la forme d'une liste de polylignes continus
        """
        if k is None :
            polylignes = [pathToPolylignes(path, nbp) for path in self.paths]
#         elif k == -1 :
#             polylignes = []
#             for path in self.paths :
#                 polylignes.extend(pathToPolylignes(path, nbp))
        else :#k entier
            polylignes = pathToPolylignes(self[k], nbp)
        return polylignes


#     def spline(self, k, nbp=100, methode=('ius',1)):
#         u"""Le k-ième path sous forme de spline (methode)
#         nbp est le nombre de points de controle de la spline
#         """
#         path = self[k]
#         delta = 1.0/nbp
#         T = np.arange(0, 1+delta, delta)
#         points = [path.point(t) for t in T]
#         points = [(p.real, p.imag) for p in points]
# #         print points
#         return NSplineSimple(points=points, methode=methode)

    def lire(self):
        pass

    @property
    def info(self):
        Nbsubp = sum(len(path.continuous_subpaths()) for path in self.paths)
        xmin, xmax, ymin, ymax = self.boundingRect(None)
        infos=[
                u"<%s>"%self.classname,
                u'%20s = '%'name'+'%s'%self.name,
                u'%20s = '%'nb paths'+'%d'%len(self),
                u'%20s = '%'nb subpath continus'+'%d'%Nbsubp,
                u'%20s = '%'bbox '+'%g, %g, %g, %g'%(xmin, xmax, ymin, ymax) + ' (xmin, xmax, ymin, ymax)'
#                 u'%20s = '%'role'+'%s'%self.role,
#                 u'%20s = '%'closed'+'%s'%self.isClosed(),
#                 u'%20s = '%'nb pts controle'+"%d"%len(self.cpoints),
#                 u'%20s = '%'methode'+"%s"%str(self.methode),
#                 u'%20s = '%'precision'+"%s"%str(self.precision),
#                 u'%20s = '%'mode echantillonage'+"%s"%self.mode,
#                 u'%20s = '%'nb pts echantillon'+"%s"%self.nbpe,
#                 u'%20s = '%'nb updates'+"%s"%self.nbupdate,
                ]
        D = []
        infos.append(u'%20s :  %5s, %10s, %10s, %11s, %23s'%('les subpaths','num', 'nom', 'couleur', 'nb segments', 'nb de subpaths continus'))#+'%s'%D)
        Nbsubp = 0
        for k, (path, attr) in enumerate(zip(self.paths, self.attributes)) :
            name = attr[u'id']
            color = self.couleur(k)
            contour = self.contour(k)
            nbseg = len(path)
            nbsubp = len(path.continuous_subpaths())
            Nbsubp += nbsubp
            infos.append(u'%20s    %5d, %10s, %10s, %11d, %23d'%(' ',k, name, color, nbseg, nbsubp))#+'%s'%D)

        return infos

    def resizedPoints(self, scale, center, translate=None):
        u"""Retourne les points de self redimensionnés
        - par translation puis
        - par homothétie de centre center, de rapport scale[0] ,scale[1]
        H(X+T) = S.(X+T) + (Id-S).C; S=matrice diagonale(scale0, scale1), C=center
        """
#         debug(scale=scale, center=center, translate=translate)
        Points = []
        for subpaths in self.points() : #subpaths = liste de np.ndarray des subpaths continus
            P = []
            for points in subpaths:#points=un path continu=un ndarray
#                 debug(points)
                points = np.asarray(points)
                if translate is not None :
                    points += translate
                hardScale(points, scale, center)
                P.append(points)
            Points.append(P)
        return Points

    def resizedPointsToFitInRect(self, rect):
        """Resize et translate pour que le svg soit entièrement contenu
        dans le rectangle rect (de cotés parallèles aux axes)
        rect est un tuple (xmin, xmax, ymin, ymax) ou bien un QRectF"""
        sw, sh, sc = self.width(), self.height(), np.asarray(self.center())
        if isinstance(rect, QRectF) :
            w, h, c  = rect.width(), rect.height(), np.asarray(p2t(rect.center()))
        else :
            (xmin, xmax, ymin, ymax) = rect
#             w, h, c  = xmax-xmin, ymax-ymin, 0.5*np.asarray([xmax+xmin, ymax-ymin])
            qrect = QRectF(xmin, ymin, xmax-xmin, ymax-ymin)
            return self.resizedPointsToFitInRect(qrect)
        s = min(w/sw, h/sh)

        return self.resizedPoints(scale=(s,s), center=c, translate=c-sc)

    def center(self, k=None):
        xm, xM, ym, yM = self.boundingRect(k)
        return (xm+self.width(k)/2, ym+self.height(k)/2)

    def width(self, k=None):
        xm, xM, ym, yM = self.boundingRect(k)
        return xM-xm

    def height(self, k=None):
        xm, xM, ym, yM = self.boundingRect(k)
        return yM-ym

    def __str__(self):
        return '\n ==> '.join(self.info)

###### tests #####

if __name__=="__main__":
    def testPathToPolyligne():
        continuous_path =  SVGPath(
            CubicBezier(start=(540.305+223.059j), control1=(541.074+230.176j), control2=(542.559+237.961j), end=(544.609+245.051j)),
            Line(start=(544.609+245.051j), end=(631.121+550.996j)),
            Line(start=(631.121+550.996j), end=(726.211+550.996j)),
            Line(start=(726.211+550.996j), end=(643.008+257.996j)),
            CubicBezier(start=(643.008+257.996j), control1=(640.977+250.93j), control2=(639.629+244.414j), end=(638.965+237.961j)),
            CubicBezier(start=(638.965+237.961j), control1=(636.172+212.078j), control2=(648.184+197.805j), end=(693.516+197.805j)),
            CubicBezier(start=(693.516+197.805j), control1=(745.254+197.805j), control2=(766.211+218.555j), end=(776.504+254.133j)),
            Line(start=(776.504+254.133j), end=(860.754+550.996j)),
            Line(start=(860.754+550.996j), end=(955.859+550.996j)),
            Line(start=(955.859+550.996j), end=(869.336+245.051j)),
            CubicBezier(start=(869.336+245.051j), control1=(848.594+172.637j), control2=(806.074+113.762j), end=(674.746+113.762j)),
            CubicBezier(start=(674.746+113.762j), control1=(570.637+113.762j), control2=(534.199+166.152j), end=(540.305+223.059j))
            )

        discontinuous_path =  SVGPath(CubicBezier(start=(1648.6+430.684j), control1=(1651.47+457.266j), control2=(1635.71+473.398j), end=(1589.79+473.398j)),
                                      CubicBezier(start=(1589.79+473.398j), control1=(1540.63+473.398j), control2=(1512.46+452.051j), end=(1502.33+417.781j)),
                                      Line(start=(1502.33+417.781j), end=(1455.62+253.477j)),
                                      CubicBezier(start=(1455.62+253.477j), control1=(1454.41+248.309j), control2=(1453.24+243.105j), end=(1452.73+238.586j)),
                                      CubicBezier(start=(1452.73+238.586j), control1=(1449.93+212.703j), control2=(1465.2+197.805j), end=(1511.74+197.805j)),
                                      CubicBezier(start=(1511.74+197.805j), control1=(1560.88+197.805j), control2=(1588.76+216.57j), end=(1599.2+253.477j)),
                                      Line(start=(1599.2+253.477j), end=(1645.92+417.781j)),
                                      CubicBezier(start=(1645.92+417.781j), control1=(1647.08+422.285j), control2=(1648.21+426.863j), end=(1648.6+430.684j)),
                                      CubicBezier(start=(1355.23+228.242j), control1=(1355.92+234.695j), control2=(1357.37+241.836j), end=(1359.35+248.309j)),
                                      Line(start=(1359.35+248.309j), end=(1409.47+426.152j)),
                                      CubicBezier(start=(1409.47+426.152j), control1=(1432.25+506.406j), control2=(1486.96+557.488j), end=(1608.54+557.488j)),
                                      CubicBezier(start=(1608.54+557.488j), control1=(1710.74+557.488j), control2=(1752.37+499.27j), end=(1746.04+440.43j)),
                                      CubicBezier(start=(1746.04+440.43j), control1=(1745.44+434.578j), control2=(1743.96+427.441j), end=(1742.02+420.996j)),
                                      Line(start=(1742.02+420.996j), end=(1692.08+245.051j)),
                                      CubicBezier(start=(1692.08+245.051j), control1=(1669.06+162.938j), control2=(1615.23+113.762j), end=(1487.85+113.762j)),
                                      CubicBezier(start=(1487.85+113.762j), control1=(1386.29+113.762j), control2=(1348.96+170.012j), end=(1355.23+228.242j))
                                      )
        polc = continuousPathToPolyligne(continuous_path)
        polc.plot(plt)
        pols = pathToPolylignes(discontinuous_path)
        for pol in pols :
            pts = pol.cpoints
            name = pol.name
            plt.plot(pts[:,0], pts[:,1], '.-',label=name)
        plt.legend()
        plt.show()

    def testLecteurSVG():
        points = np.asarray([[0.,0],[1,0],[1,1],[0,1]])
        print hardScale(points, [0.5,0.5], [2,1], [10,10])
        print '*******Lecteur SVG******'
        filename = Path(VALIDATION_DIR,'dxf','blocjonc.dxf')
        filename = Path(VALIDATION_DIR,'decorations','Bringhen_D_cmyk gross.dxf')
        filename = Path(VALIDATION_DIR,'dxf','cot_d03.DXF')
        filename = Path(VALIDATION_DIR,'dxf','teste01.dxf')
        filename = Path(VALIDATION_DIR,'faial2.dxf')#3d, je sais pas lire...
        filename = Path(VALIDATION_DIR,'dxf','Faial2v5A.DXF')

        filename = utilitaires.Path(VALIDATION_DIR,'decorations','Bringhen_D_cmyk gross.svg')
        filename = utilitaires.Path(VALIDATION_DIR,'decorations','EuroSport.svg')
        svg = LecteurSVG(filename)
        print svg

#         Polylignes = svg.polylignes()#liste de liste de np.array
        debug('bbox[1]', svg.boundingRect(1))
        debug('bbox[None]', svg.boundingRect())
        from matplotlib import pyplot as plt
        for k,polylignes in enumerate(svg.polylignes()) :
            color = '#777777'#svg.couleur(k)
            if color == '#ffffff' : color = '#000000'
            for j, pol in enumerate(polylignes) :
                pts = pol.cpoints
                plt.plot(pts[:,0], pts[:,1], color)#, label=str((k,j)))

#         xm, xM, ym, yM = svg.boundingRect()
#         bbcenter = (-0.5*(xm+xM),-0.5*(ym+yM))
#         rect = QRectF(xm,ym,xM-xm,yM-ym)#rectangle d'encombrement svg
        sw, sh, sc = svg.width(), svg.height(), np.asarray(svg.center())
        rect = -10000, -1000, 12000, 3000
        qrect = QRectF(*rect)#rectangle d'encombrement voile
        rect = xm, xM, ym, yM = qrect.left(), qrect.right(), qrect.top(), qrect.bottom()

        vw, vh, vc  = qrect.width(), qrect.height(), np.asarray(p2t(qrect.center()))

        print'(sw,sh)', (sw,sh)
        print'(vw,vh)', (vw,vh)
        s = min(vw/sw, vh/sh)
        plt.plot([xm, xM, xM, xm, xm],[ym, ym, yM, yM, ym], 'r.-')
        Points = svg.resizedPoints(scale=(s,s), center=vc, translate=vc-sc)
        Points = svg.resizedPointsToFitInRect(qrect)
        for k,P in enumerate(Points) :
            color = svg.couleur(k)
            if color == '#ffffff' : color = '#000000'
            for j, pts in enumerate(P) :
                plt.plot(pts[:,0], pts[:,1], color)#, label='transl')

        plt.legend()
        plt.show()

    ###TESTS####
#     testPathToPolyligne()
    testLecteurSVG()

