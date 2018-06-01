#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Axile -- Outil de conception/simulation de parapentes Nervures

Classe : LecteurData
Description : lecteur générique pour fichier paramètres.
@module : inout.lecteurdata
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 18 janv. 2013
'''
import sys,os
from lecteur import Lecteur
from utilitaires import whoami,trace
from config import VALIDATION_DIR,DATA_DIR

class LecteurData(Lecteur):
    """
    Lecteur generique pour fichiers parametres
    - Pour chaque ligne, tout ce qui est après '#' est ignoré
    - Les lignes vides (ou avec que des espaces) sont ignorées
    La méthode Lire(filename) retourne un dictionnaire des variables lues sur filename.
    Exemple, pour un fichier 'machin.dat' ainsi constitué (les < xxx > ne font pas partie du fichier):
    -----------------------------------------------------
    # commentaire
    # encore commentaire
    <blank line>
    var1 1. 2. 3.                # liste de reels
    var2 2                       # un entier
    var3 'chaine de caracteres'  # une chaine
    var3 "chaine de caracteres"  # une chaine
    var5 <expression>            # toute expression Python évaluable par eval, dans le contexte
    var4 variable illisible      # ...
    <fin de fichier machin.dat>

    L'instruction de lecture est :
    >>> d = LecteurData("machin.dat").lire()
    d est un dictionnaire.

    Pour instancier ces variables, dans un objet 'toto' :
    >>> for key, value in d.iteritems() :
    ...     setattr(toto, key, value)

    C'est tout.
    """
    def __init__(self,filename):
        Lecteur.__init__(self,filename,readlines=True)

    def split(self,line):
        """Partage au premier espace trouvé"""
        line=line.strip()
        comm=line.find('#')
        if comm>=0 :
            return self.split(line[:comm])
        try : b=line.index(' ')
        except ValueError : return '' #pas de blanc, pas de valeur
        return line[:b].strip(),line[b:].strip()

#    def correctiveAction(self, word):
#        """
#        Voir self.value(), en cas de dictionnaires imbriqués, il faut instancier les valeurs du dictionnaire
#        comme des variables locales
#        exemple : on a lu un dictionnaires self.dict['d'] qui vaut {'a':1}
#        et on tombe sur une ligne
#        d1 {'aa':d}
#        word vaut donc "{'aa':d}" qui n'est pas évaluable car d est inconnu dans le contexte,
#        on ne connaît que self.dict['d’] qui vaut {'a':1}
#        on instancie self.dict['d’] comme une variable locale d
#        on réévalue word
#        """

    def value(self,word=""):
        if len(word)==0 :
            return None
        try :
            return eval(word)
        except (NameError) as msg :
            try : #peut-etre des dictionnaires imbriqués
                for key,value in self.dict.iteritems() :
                    if isinstance(value, (str, unicode)) :
                        exec('%s="%s"'%(key,str(value)))
                    else :
                        exec('%s=%s'%(key,str(value)))
                return eval(word)
            except (ValueError,NameError) as msg:
                print>>sys.stderr,whoami(self)," NameError::valeur illisible : %s"%word,msg
        except (ValueError) as msg :
            print>>sys.stderr,whoami(self)," ValueError::valeur illisible : %s"%word,msg
        return None


    def lire(self):
        """
        Fichier data :
        ===============
        # commentaire
        # encore commentaire

        var1 1. 2. 3.                # liste de reels
        var2 2                       # un entier
        var3 'chaine de caracteres'  # une chaine
        var3 "chaine de caracteres"  # une chaine
        var5 <expression>            # toute expression Python évaluable par eval, dans le contexte
        var4 variable illisible      # ...

        """
        self.dict={}
        for line in self.lines :
            line=line.strip()
            if not line : continue
            try :
                key,value=self.split(line)
            except ValueError :
                continue
            self.dict[key]=self.value(value)

        return self.dict

if __name__=="__main__":
    from path import Path
    class A():pass
    datadir=Path(DATA_DIR,'preferences')
    files = sorted(datadir.listdir('*.dat'))
    files = (Path(datadir, 'animation.dat'),)
    for filename in files :
        fn = filename.name
        print '\n'.join([len(fn)*'=',fn,len(fn)*'='])
        d = LecteurData(Path(datadir,filename)).lire()
        a = A()
        for key,value in sorted(d.iteritems()) :
            print "%20s = %-20s => %10s"%(key,str(value),type(value))
            setattr(a,key,value)

        print 21*"-"+'\n'+u"variables instanciées"+'\n'+21*"-"
        for key, val in sorted(a.__dict__.iteritems()) :
            print u"attribut %s=%s"%(key,val)
    # exit()
    print files
    filename = "profil.dat"
    lecteur=LecteurData(Path(DATA_DIR,"preferences",filename))
    p,n=lecteur.lire()['defaultprofil']
    print type(p),type(n)
    print p, n
    # exit()
    filename="aero.dat"
    f=LecteurData(Path(datadir,filename))
    d=f.lire()
#    Exemple d'instanciation dans une variable 'a'
    # class A():pass
    # a=A()
    for key,value in sorted(d.iteritems()) :
        print "%10s = %10s : %10s"%(key,str(value),type(value))
        setattr(a,key,type(value)(value))

    print u"variable instanciée"
    for key, val in sorted(a.__dict__.iteritems()) :
        print u"attribut %s=%r"%(key,val)
    exit()
