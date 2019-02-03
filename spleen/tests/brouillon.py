#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on Feb 1, 2019

@author: puiseux
'''
from numpy import (asarray, min, max, argmin, linspace,sqrt,ndarray,zeros,ones)
from numpy.linalg import norm
from splinesimple import NSplineSimple
from matplotlib import pyplot as plt
from scipy.optimize._minimize import minimize_scalar
from utilitaires.lecteurs import pointsFrom
plt.rcParams["figure.figsize"] = (20,10)
from utilitaires.utilitaires import debug, XY, rdebug, dist2
NUM = 62

def projection(S, p, discret=0, t0=0.0, t1=1.0, debog=None):
    u"""
    Calcule le projeté du point p sur la portion S01 de la spline S
        S01 = {S(t), t dans I=[t0,t1]},
        i.e. le min de la fonction phi : s-> dist(p,S(s)), t dans I.
    Pour cela,
    - si discret=0, on appelle scipy.minimize_scalar() sur l'intervalle I.
        Ce faisant, on risque fort de ne pas tomber sur le bon min,
        lorsqu'il y a plusieurs minima locaux.
        Typiquement, si S est un profil (fermé) et p=centre de gravité de
        S, la courbe t -> dist(p,S(t)) présente 2 minima locaux.
        On doit donc localiser le min absolu en discrétisant la spline.
        C'est ce qui est fait si discret>0
    - si discret>0, on localise le min absolu en discrétisant la spline
        => S.dpoints = S(T) ou T contient 1+discret pas de temps
        linéairement répartis, i.e. discret intervalles de largeur dt.
        # on cherche tm=le min (discret) des distances dist(p, S.dpoints)
        # on raffine ensuite sur l'intervalle [tm-dt, tm+dt]
    :param p : tuple, list ou ndarray((2,)), les coord. du point p.
    :param discret : int, le nombre de points de discrétisation
    :param t0,t1: float, float, l'intervalle de recherche du min.
        on doit avoir 0.0 <= t0 < t1 <= 1.0
    :return tt, dd, nev, msg : float, float, int, str
        - tt dans [t0,t1] est le temps sur S du projeté p' = S(tt)
        - dd est la distance euclidienne dist(p,p')
        - nev le nombre d'évaluations de la fonction phi(s)
        - msg un message fournit par minimize_scalar()
    """
    if discret>0 :
        S.nbpd = discret#_dpoints est effacé et recalculé à la demande
        eps = 1.0/(discret-1)
        #les distances de p aux points de dpoints
        D = norm(S.dpoints-p,axis=1)
        idx = argmin(D)
        d = D[idx]
#         debug(d=d,idx=idx)
        T = linspace(0,1,discret)
        #t est le temps du point dpoints[idx] sur le polygone S.dpoints
        t = T[idx]
        #raffinement
        t0, t1 = max([0,t-eps]), min([1,t+eps])
        tt, dd, nev, msg = projection(S, p, discret=0, t0=max([0,t-eps]), t1=min([1,t+eps]))
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
        res = minimize_scalar(lambda t: (a-S.sx(t))**2 + (b-S.sy(t))**2,
                              bounds=(t0, t1),
                              method='bounded',
                              options={'xatol':1.0e-9})
        return res.x, sqrt(res.fun), res.nfev, res.message

def projectionObj(S, obj, discret=0):
    u"""
    Calcule la projection de l'objet obj sur la spline S par appels répétés à
        projection(S,obj[k],...)
    :param obj: NSplineSimple ou ndarray((n,2),dtype=float) ou bien
                un seul point (float, float) a projeter sur S.
        # si obj et une NSplineSimple on projete les points de obj.cpoints sur S
        # si obj est un ndarray((n,2)) on projete les points de obj sur S
    :param discret: int, un nombre de points de discrétisation pour le calcul
        des projections.
    :return tprojs, dists, nevs, msgs: cf doc scipy.optimize.OptimizeResult
        - tprojs est la liste des parametres t des (c)points de obj sur S
        - dists est la liste des distances à S des projetés
        - nevs est la liste des nb d'évaluation de phi
        - msgs sont des messages de convergence donnés par scipy.minimize_scalar
    La fonction phi(t) (distance de p à S, lorsque p est un point (x,y))
    peut admettre plusieurs minima locaux.
    Il faut les trouver tous et les comparer entre eux pour obtenir le vrai
    minimum.
    C'est ce que fait l'appel à projete()
    """
    if isinstance(obj, (tuple, list, ndarray)) and len(obj)==2 :
        return projection(S, p=obj, discret=discret)
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
            debog = False if k != NUM else k
            tprojs[k], dists[k], nevs[k], msgs = projection(S, p, discret, debog=debog)
        return tprojs, dists, nevs, msgs

def elaguer(S, eps=0.5, replace=False):
    u"""
    On cherche une spline se avec un minimum de points de contrôle
    et qui soit à une distance de S.cpoints < eps (en ‰ de la longueur de S)

    La distance(S.cpoints, se) est le plus grand écart entre
    la spline calculée et les points de contrôle de la spline S.
    autrement dit le max des distances d(S.cpoints[k], se), k=0,1,...
    où d(P, se) est le min de la fonction t -> norme(P-se(t)).
    Voir la methode S.distanceSplineToObj().
    On discrétise finement S (une seule fois) -> tableau de points D0
    On initialise se partant des quelques points de contrôle dont
    S.cpoints[0] et S.cpoints[-1].
    On discrétise finement se (à chaque ajout de point) -> tableau de points D1
    puis on rajoute des points de contrôle à se jusqu'a obtention de la précision désirée (eps)
    Pour rajouter un point, on repère la distance point à point de D0 et D1,
    disons qu'elle est au i-eme point et on insère dans se un point de contrôle
    que l'on positionne exactement en D0[i].

    [ne marche pas =>] Quand la précision désirée est atteinte, on met
    un petit coup d'optimisation pour améliorer la position des points
    de contrôle de se.

    :param eps: float, la spline resultante est à une distance eps au maximum de S.
    :param replace: bool, si True, la spline élaguée remplace S.
    :return: se, pm, (n0,n1)

        - se est la spline NSplineSimple élaguée
        - pm = S.pourMille(d) est la précisison en valeur relative :
            la spline se ne s'éloigne pas à plus de pm ‰ de S.
        - n0, n1 = nb de points de contrôle avant et après élagage

    """
    debog=True
    if len(S) < 10 :
        rdebug(u"Élagage inutile, %d points de contrôle seulement"%len(S))
        return S, S.pourMille(0),(len(S),len(S))
    nd = 50
    c0 = S.cpoints.copy()
    t, m = S.methode
    n = len(S)
    if t == 'cubic' and m == 'periodic' :#il faut au moins 3 points et c[0] = c[-1]
        seen = [0, n/4, n/2, 3*n/4, n-1]
        debug('cubic, periodic, 5 points',seen=seen)
        c1 = pointsFrom(S[seen])
    elif dist2(c0[0], c0[-1]) < 1.0e-6 :
        seen = [0, n/3, 2*n/3, n-1]
        c1 = pointsFrom(S[seen])
        debug('ferme, 4 points',seen=seen)
    else:
        seen = [0,n-1]
        c1 = pointsFrom(S[seen])
        debug('ouvert, 2 points',seen=seen)
    se = NSplineSimple(cpoints=c1,
                       methode=S.methode,
                       nbpd=S.nbpd,
                       name='%s-elaguee'%S.name,
                       mode=S.mode)
    if debog :
        X0, Y0 = XY(c0)
        more = [(X0, Y0, 'k.-','c0'),]
        texts = [(X0[0],Y0[0]-0.05,u'%d'%0),(X0[-1],Y0[-1]-0.05,u'%d'%(len(X0)-1))]
        texts+= [(c1[0,0],c1[0,1],u'se:%d'%0),(c1[-1,0],c1[-1,1],u'se:%d'%(len(c1)-1))]
    T1, D01, _, _ = projectionObj(se, c0, discret=nd)#T1, D01 de longueur len(c0)
    while len(seen)<len(c0):
#     while len(seen)<18:
        #indices du-des point-s ou la distance entre les deux splines est max
        debug(seen=seen)
        D01[seen] = 0.0#On ne veut pas remettre 2 fois le même point (=>point double)
        d = max(D01)#le point de c0 le plus eloigné de se
        imaxd = (D01 == d).nonzero()
        dm = S.pourMille(sqrt(d))
        idx0 = imaxd[0][0]
#         if idx0 == 58 :
#             debug('D%d = %s'%(idx0,D01.tolist()))
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
#             print u'    idx0  = %d'%idx0
            print u'    insertion dans se en position i =',i#, ' ; dist=',D[idx],' ; maxD=',d
#             print '    dist=maxD ?',D01[idx0]==d
#             print '    *** les distances [d(c0[i],se) i=0...len(c0) ](len(se) = %d points) à c0'%len(se), D01.tolist()
            pj = se(t1)
            more += [((pos0[0],pj[0]),(pos0[1],pj[1]),'go-',u'à ajouter en position %d'%(i))]
            texts += [(pos0[0],pos0[1],u'%d'%(i))]
#                 debug(more=more)
            se.plot(plt,
                    titre=u'Spline élaguée : \nseen=%s, \ndist=%.1g ‰'%(str(seen),d),
                    more=more,texts=texts, show=True)
            more = more[:-1]
#                 texts=texts[:-2]
            c1 = se.cpoints
            texts = [(c1[0,0],c1[0,1],u'se:%d'%0),(c1[-1,0],c1[-1,1],u'se:%d'%(len(se)))]

#         try :
        se.insertPoint(pos0, i)
#         except ValueError as msg :#Impossible d'inserer le point (point double ?)
#             debug(u'Ca devrait pas arriver, je tente autre chose', msg, pos0=pos0)
#             raise msg
#             Td = linspace(0,1,nd)
#             d0 = asarray(S(Td))
#             d1 = asarray(se(Td))
#             ad = norm(d0-d1, 2, axis=1)#[1:-1]#ecart point à point des deux splines discrétisées
#             mad = (ad == max(ad)).nonzero()#indice des points ou la distance entre les deux splines est max
#             idx = mad[0]
#             t = Td[idx][0]
#             pos0 = d0[idx][0]
#             try :
#                 se.insertPoint(pos0)
#             except ValueError as msg :
#                 rdebug(msg, pos0=pos0)
#                 rdebug(u'Precision non atteinte, iteration %d : dist = %.2e '%(len(seen),d))
#                 break
#                     c1 = se.cpoints.copy()

#             d = [distancePointSpline(c0[i], se, t0=t).fun for i, t in enumerate(T0)]
#         D01 = distanceSplineToObj(se, c0, discret=nd)[1]
        T1, D01, _, _ = projectionObj(se, c0, discret=nd)
        dm = S.pourMille(d)
        debug(u'dist-%d = %.2g‰ '%(len(seen),dm))
        c1 = se.cpoints.copy()
        if dm<eps : break
    if len(se) == len(S) :#même nb de points : on garde la spline initiale.
        se.cpoints = S.cpoints.copy()
        n0 = len(S)
        n1 = len(se)
        return se, d,(n0,n1)
    #ici on peut mettre un petit coup d'optimisation pour replacer les points de contrôle de se
    #ca ne marche pas du tout !!!
    n0 = len(S)
    n1 = len(se)
    debug('Apres ELAGAGE : dist = %.2g mm/m ; %d => %d '%(dm,n0,n1))

    if replace :
        S.cpoints = se.cpoints
    return se, d, (n0,n1)

#celui-ci marche (ouvert
ouvert = asarray([[2.619476, 0.0], [2.58009, 0.0068919999999999815], [2.557583, 0.010827999999999977], [2.530067, 0.015635999999999983], [2.496937, 0.021418999999999994], [2.461583, 0.02757799999999999], [2.422409, 0.03436299999999998], [2.378257, 0.04185999999999998], [2.3387450000000003, 0.048644999999999994], [2.291149, 0.05678799999999998], [2.234816, 0.066379], [2.159633, 0.07910099999999998], [2.079646, 0.09254499999999999], [1.9955470000000002, 0.10657999999999998], [1.905313, 0.12152399999999997], [1.811423, 0.13694800000000001], [1.71244006, 0.15306999999999998], [1.6140466, 0.168954], [1.516619, 0.18454500000000001], [1.414023, 0.20081299999999996], [1.311207, 0.21696299999999996], [1.206928, 0.23318699999999998], [1.1125500000000001, 0.247734], [1.003597, 0.26343], [0.905398, 0.275535], [0.808524, 0.28580099999999997], [0.7151310000000001, 0.29162899999999997], [0.6254, 0.29386999999999996], [0.5399, 0.29217], [0.45913000000000004, 0.28627199999999997], [0.38362000000000007, 0.276041], [0.31381000000000014, 0.26148899999999997], [0.25014000000000003, 0.24278], [0.19300000000000006, 0.22022999999999998], [0.14275000000000015, 0.194295], [0.09970000000000012, 0.16553299999999999], [0.06411000000000011, 0.13456199999999996], [0.036190000000000166, 0.10200599999999999], [0.016129999999999978, 0.068438], [0.0040400000000000436, 0.034329], [0.0, 0.0], [0.0015500000000001624, -0.012182000000000012], [0.0048200000000000465, -0.023260000000000017], [0.010990000000000055, -0.034644200000000014], [0.02076000000000011, -0.04379510000000002], [0.03133000000000008, -0.05087460000000002], [0.041340000000000154, -0.05733490000000001], [0.13478000000000012, -0.11178900000000001], [0.14381, -0.11709790000000002], [0.15451000000000015, -0.12239003000000001], [0.16610999999999998, -0.12685575000000002], [0.1810400000000001, -0.13177934000000002], [0.21382000000000012, -0.13938238], [0.25636000000000014, -0.14548950000000002], [0.31609, -0.1508921], [0.3831200000000001, -0.1552542], [0.4642600000000001, -0.15996470000000002], [0.5438800000000001, -0.16315420000000003], [0.6317600000000001, -0.16668880000000003], [0.730441, -0.1705661], [0.8221980000000001, -0.1715721], [0.8962490000000001, -0.17146740000000002], [1.036013, -0.16904200000000003], [1.115243, -0.16583710000000002], [1.216719, -0.16066940000000002], [1.3206470000000001, -0.1539469]])#, [2.619476, 0.0]])
#celui-ci ne marche pas (c0[32] et c0[58]...)
ferme = asarray([[2.619476, 0.0], [2.58009, 0.0068919999999999815], [2.557583, 0.010827999999999977], [2.530067, 0.015635999999999983], [2.496937, 0.021418999999999994], [2.461583, 0.02757799999999999], [2.422409, 0.03436299999999998], [2.378257, 0.04185999999999998], [2.3387450000000003, 0.048644999999999994], [2.291149, 0.05678799999999998], [2.234816, 0.066379], [2.159633, 0.07910099999999998], [2.079646, 0.09254499999999999], [1.9955470000000002, 0.10657999999999998], [1.905313, 0.12152399999999997], [1.811423, 0.13694800000000001], [1.71244006, 0.15306999999999998], [1.6140466, 0.168954], [1.516619, 0.18454500000000001], [1.414023, 0.20081299999999996], [1.311207, 0.21696299999999996], [1.206928, 0.23318699999999998], [1.1125500000000001, 0.247734], [1.003597, 0.26343], [0.905398, 0.275535], [0.808524, 0.28580099999999997], [0.7151310000000001, 0.29162899999999997], [0.6254, 0.29386999999999996], [0.5399, 0.29217], [0.45913000000000004, 0.28627199999999997], [0.38362000000000007, 0.276041], [0.31381000000000014, 0.26148899999999997], [0.25014000000000003, 0.24278], [0.19300000000000006, 0.22022999999999998], [0.14275000000000015, 0.194295], [0.09970000000000012, 0.16553299999999999], [0.06411000000000011, 0.13456199999999996], [0.036190000000000166, 0.10200599999999999], [0.016129999999999978, 0.068438], [0.0040400000000000436, 0.034329], [0.0, 0.0], [0.0015500000000001624, -0.012182000000000012], [0.0048200000000000465, -0.023260000000000017], [0.010990000000000055, -0.034644200000000014], [0.02076000000000011, -0.04379510000000002], [0.03133000000000008, -0.05087460000000002], [0.041340000000000154, -0.05733490000000001], [0.13478000000000012, -0.11178900000000001], [0.14381, -0.11709790000000002], [0.15451000000000015, -0.12239003000000001], [0.16610999999999998, -0.12685575000000002], [0.1810400000000001, -0.13177934000000002], [0.21382000000000012, -0.13938238], [0.25636000000000014, -0.14548950000000002], [0.31609, -0.1508921], [0.3831200000000001, -0.1552542], [0.4642600000000001, -0.15996470000000002], [0.5438800000000001, -0.16315420000000003], [0.6317600000000001, -0.16668880000000003], [0.730441, -0.1705661], [0.8221980000000001, -0.1715721], [0.8962490000000001, -0.17146740000000002], [1.036013, -0.16904200000000003], [1.115243, -0.16583710000000002], [1.216719, -0.16066940000000002], [1.320647, -0.1539469], [ 1.75359   , -0.10263127], [ 2.186533  , -0.05131563], [2.619476, 0.0]])
c0 = ouvert
s0 = NSplineSimple(cpoints=c0, methode=('cubic', 'not-a-knot'), name=u's0', nbpd=100)
se = NSplineSimple(**{'name': u's0-elaguee', 
                      'classename': 'NSplineSimple', 
                      'methode': ('cubic', 'not-a-knot'), 
                      'role': 'NSplineSimple', 
                      'nbpd': 100, 
                      'cpoints': [[2.619476, 0.0], [0.0040400000000000436, 0.034329], [1.3206470000000001, -0.1539469]]})
D62 = asarray([1.5924604453904028, 0.6968182821515951, 0.24823339559579746, 0.6459610340321058, 0.9656523857649552, 1.1018692940222883, 1.0468273993282031, 0.7988838722961408, 0.35861342019467385, 0.28503399095548265])
T62 = asarray([0.0, 0.1111111111111111, 0.2222222222222222, 0.3333333333333333, 0.4444444444444444, 0.5555555555555556, 0.6666666666666666, 0.7777777777777777, 0.8888888888888888, 1.0])
D62 = asarray([1.5924604453904028, 1.501904484164042, 1.4130746709431556, 1.3259903626188043, 1.2406752819244435, 1.157158799515338, 1.0754776862036406, 0.9956785441628345, 0.9178212353917319, 0.8419838012316887, 0.7682696525772014, 0.696818282151595, 0.6278215343426847, 0.5615487642931251, 0.4983862890768003, 0.438899505717411, 0.38392892494254566, 0.33472835678200785, 0.2931236603633685, 0.26156429344323645, 0.24273648838859918, 0.23842624817197627, 0.2482333955957976, 0.26959311093990856, 0.2991415454982086, 0.33392765257059426, 0.37179857498491603, 0.4113008194488482, 0.4514741549107515, 0.49168079488994376, 0.5314897091519653, 0.5706031891991561, 0.6088109150058918, 0.6459610340321059, 0.6819416615590411, 0.716668814437721, 0.7500783782310689, 0.7821206505797568, 0.8127565625284021, 0.8419550140960456, 0.8696909636470616, 0.8959440362741039, 0.9206974954865575, 0.9439374731643753, 0.9656523857649558, 0.9858324866599741, 1.0044695192159525, 1.021556445302336, 1.03708723088877, 1.0510566752925883, 1.0634602741203745, 1.0742941084523314, 1.0835547546396893, 1.091239210423775, 1.0973448340785323, 1.1018692940222883, 1.104810526906993, 1.1061667026221433, 1.1059361949807587, 1.104117557111181, 1.1007095007797147, 1.0957108790287624, 1.0891206716436905, 1.0809379730669928, 1.0711619824668042, 1.0597919957435644, 1.046827399328203, 1.032267665691812, 1.016112350554801, 0.9983610918579823, 0.9790136106451138, 0.9580697141145108, 0.9355293012380238, 0.9113923715362328, 0.885659037864853, 0.8583295444491171, 0.8294042919635135, 0.7988838722961404, 0.766769116930027, 0.733061164907153, 0.6977615596205667, 0.6608723891190814, 0.6223964939132313, 0.5823377827528813, 0.5407017272275966, 0.49749616466688035, 0.45273265815510266, 0.4064289212337135, 0.35861342019467335, 0.3093348213571847, 0.25868344031369106, 0.20684696164850253, 0.1542853334838542, 0.1024588442133435, 0.05838243039687396, 0.05561321046811833, 0.10061340159161396, 0.15838864619813148, 0.22042644140857262, 0.28503399095548265])
T62 = asarray([0.0, 0.010101010101010102, 0.020202020202020204, 0.030303030303030304, 0.04040404040404041, 0.05050505050505051, 0.06060606060606061, 0.07070707070707072, 0.08080808080808081, 0.09090909090909091, 0.10101010101010102, 0.11111111111111112, 0.12121212121212122, 0.13131313131313133, 0.14141414141414144, 0.15151515151515152, 0.16161616161616163, 0.17171717171717174, 0.18181818181818182, 0.19191919191919193, 0.20202020202020204, 0.21212121212121213, 0.22222222222222224, 0.23232323232323235, 0.24242424242424243, 0.25252525252525254, 0.26262626262626265, 0.27272727272727276, 0.2828282828282829, 0.29292929292929293, 0.30303030303030304, 0.31313131313131315, 0.32323232323232326, 0.33333333333333337, 0.3434343434343435, 0.3535353535353536, 0.36363636363636365, 0.37373737373737376, 0.38383838383838387, 0.393939393939394, 0.4040404040404041, 0.4141414141414142, 0.42424242424242425, 0.43434343434343436, 0.4444444444444445, 0.4545454545454546, 0.4646464646464647, 0.4747474747474748, 0.48484848484848486, 0.494949494949495, 0.5050505050505051, 0.5151515151515152, 0.5252525252525253, 0.5353535353535354, 0.5454545454545455, 0.5555555555555556, 0.5656565656565657, 0.5757575757575758, 0.5858585858585859, 0.595959595959596, 0.6060606060606061, 0.6161616161616162, 0.6262626262626263, 0.6363636363636365, 0.6464646464646465, 0.6565656565656566, 0.6666666666666667, 0.6767676767676768, 0.686868686868687, 0.696969696969697, 0.7070707070707072, 0.7171717171717172, 0.7272727272727273, 0.7373737373737375, 0.7474747474747475, 0.7575757575757577, 0.7676767676767677, 0.7777777777777778, 0.787878787878788, 0.797979797979798, 0.8080808080808082, 0.8181818181818182, 0.8282828282828284, 0.8383838383838385, 0.8484848484848485, 0.8585858585858587, 0.8686868686868687, 0.8787878787878789, 0.888888888888889, 0.8989898989898991, 0.9090909090909092, 0.9191919191919192, 0.9292929292929294, 0.9393939393939394, 0.9494949494949496, 0.9595959595959597, 0.9696969696969697, 0.9797979797979799, 0.98989898989899, 1.0])
plt.figure()
plt.plot(T62,D62)
plt.xlabel('t')
plt.ylabel('d')
ce = se.cpoints
T, D, _, _ = projectionObj(se, s0, discret=100)
t, d, nev, msg = projection(se, s0[62], discret=100, t0=0, t1=1)
print t, d, nev, msg 
pjs = se(T)
X0,Y0 = XY(c0)
Xd,Yd = XY(se.dpoints)
X,Y = XY(pjs)
#
# # more = [(X,Y,'g*',u'projetés'),(X0,Y0,'go','uc0')]
more = [([X[k],X0[k]], [Y[k],Y0[k]], 'bo-', '%d'%k) for k in range(len(s0)-5,len(s0))]
more+= [(X0,Y0,'ro-','c0')]
more+= [(Xd,Yd,'bo-','se.dpoints')]
se.plot(plt,more=more)
exit()
se,d,(n0,n1) = elaguer(s0, 0.5)
se.nbpd=1000
se.plot(plt,more=[(X0,Y0,'ro-','c0')])
exit()
x,y= 0,0
idx = 58

p = c0[idx]
nd = se.nbpd

t, d, ev, msg = projection(se, p, discret=nd)
debug('==> c0[%d]'%(idx))
print "    discret = %-6d     ; t=%-6.3g, d=%-6.3g, ev=%2d, msg='%s'"%(nd, t, d, ev, msg)
T = linspace(0,1,6)
TT =zip(T[:-1],T[1:])
for t0, t1 in TT :
    t, d, ev, msg = projection(se, p, t0=t0, t1=t1)
    print '    (t0, t1)=(%4.3g,%4.3g) ;'%(t0,t1), "t=%-6.3g, d=%-6.3g, ev=%2d, msg='%s'"%(t, d, ev, msg)

T = linspace(0,1,nd)
dp = se(T)# identique à dp = se.dpoints
dists = asarray([norm(d-p) for d in se(T)])
# print norm(D01-dists)
plt.figure(1)
plt.plot(T,dists)
plt.xlabel('t')
plt.ylabel(u'y=dist(p,se(t))')
plt.title(u'distances de p=c0[%d] à se'%idx)
# plt.show()

# exit()

# plt.figure(2)
# plt.plot(Dd, 'r*-', label=u'D discret=100')
# plt.plot(D, 'b*-', label=u'D discret=0')
# plt.xlabel(u'num. point de c0')
# plt.ylabel(u'distance a se et le temps du projeté sur se')
# # plt.legend()
# plt.title(u'les distances des points de c0 a se')
# # plt.figure(3)
# plt.plot(Td, 'r.-', label=u'T discret=100')
# plt.plot(T, 'b.-', label=u'T discret=0')
# # plt.xlabel(u'num. point de c0')
# # plt.ylabel(u't')
# plt.legend()
# plt.title(u'les temps dans se des projetés des c0[i] sur se')
#
# plt.figure(4)
# plt.plot(D19,'r.-', label='D19')
# plt.plot(D20,'g.-', label='D20')
# plt.plot(D21,'b.-', label='D21')
# plt.title(u'Les distances c0[19,20,21] a se.dpoints (%d points)'%len(se.dpoints))
# plt.xlabel(u'num. point de se.dpoints')
# plt.legend()
more=[(X,Y, 'k.-','c0'),(x,y,'mo','ajout')],
se.plot(plt, show=False)
plt.show()
