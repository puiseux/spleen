#!/usr/bin/python
#-*-coding: utf-8 -*-
u"""
Pierre PUISEUX
Octobre 96, Mai 2012, juillet 2016
"""
import sys
sys.path.append("/Users/puiseux/Documents/GitHub/spleen/spleen")
for p in sorted(sys.path) : print p
from path import Path
u"""
- Le module 'config.py' EST sur le dépot central,
- Chaque utilisateur a sa propre configuration, qui est lue sur 'configSpleen.txt'
    lequel 'configSpleen.txt' N'EST PAS sur le dépot car il est propre à chaque utilisateur.
- Un 'configSpleen.txt'  minimal est le suivant :

        ####################################################################
        #Configuration des répertoires utiles à Axile :
        #   ROOT_DIR : répertoire d'installation (path complet)
        #   RUNS_DIR : répertoire de stockage des projets (path complet)
        #Ce fichier doit se trouver sous <ROOT_DIR>/dist
        ####################################################################
        ROOT_DIR = /Volumes/developpement/axile-shape-edit
        RUNS_DIR = /Volumes/developpement/axile-shape-edit/workspace
        VTK_ENABLED = True ou False
        ####################################################################

- ROOT_DIR est le répertoire d'installation et
- RUNS_DIR est le répertoire de stockage des projets
- VTK_ENABLED = False si vtk5 n'est pas dispo, dans ce cas, la vue 3d est désactivée.
- Remarques :
  --------
    - Dans le fichier 'configAxile.txt' toutes les lignes commencant par '#' sont ignorées
    - hors d'une ligne commentaire le signe '=' ne doit être utilisé que pour indiquer
        la valeur d'une variable (comme ROOT_DIR)
    - Une variable simple, peut être évaluée :
    - une ligne qui ne commence pas par '#' est séparée en 2 morceaux par le signe '='.
        1. A gauche le nom de la variable,
        2. à droite sa valeur, et c'est tout.
            des blancs sont autorisés SAUF dans le nom de variable.
        NOM_VARIABLE =    1                => OK
        NOM_VARIABLE2= True                => OK
        NOM_VARIABLE3 = [1, 2]             => OK
        MON_PI = 3.14   la valeur de pi    => NON évalué car des caracteres autres que espaces à la suite de la valeur
    - le fichier confAxile.txt doit se trouver dans le répertoire <ROOT_DIR>/dist/guimain
    - les paths doivent être écrits sous la forme
        ROOT_DIR = Le/path/vers/le/repertoire
        ROOT_DIR = 'Le/path/vers/le/repertoire'
        ROOT_DIR = 'Le\path\vers\le\repertoire' etc...

"""
dirname = Path(__file__).dirname
path_ = Path(dirname, "confSpleen.txt")

with open(path_) as file_:
    for line in file_:
        line = line.strip()
        if not line or line[0] == '#' : continue
        if '='  in line :
            words = line.split('=')
            if 'ROOT_DIR' in words[0]:
                ROOT_DIR = Path(words[1].strip().strip("'").strip('"'))#lstrip removes whitespaces at the beginning
            if 'RUNS_DIR' in words[0]:
                RUNS_DIR = Path(words[1].strip().strip("'").strip('"'))
            if 'VTK_ENABLED' in words[0] :
                VTK_ENABLED = eval(words[1])
WORK_DIR = RUNS_DIR
SOURCES_DIR = Path(ROOT_DIR,'sources')
# sys.path.append(SOURCES_DIR)
DATA_DIR = Path(ROOT_DIR,'data')
VALIDATION_DIR = Path(ROOT_DIR,'validation')
TRASH_DIR = Path(RUNS_DIR,'trash')
TEST_MODE = False
print 'config.py = OK'
print "ROOT_DIR       = '%s'"%ROOT_DIR
print "RUNS_DIR       = '%s'"%RUNS_DIR
print "SOURCES_DIR    = '%s'"%SOURCES_DIR
print "DATA_DIR       = '%s'"%DATA_DIR
print "VALIDATION_DIR = '%s'"%VALIDATION_DIR
print "VTK_ENABLED    = '%s'"%VTK_ENABLED

