#!/usr/local/bin/python2.7
#coding: utf-8
from setuptools import setup
import spleen
setup(name='spleen',
      version=spleen.__version__,
      description=u'splines 2d paramétriques d\'interpolation et d\'ajustement',
      url='',
      author='Pierre Puiseux',
      author_email='pierre@puiseux.name',
      license='None',
      packages=find_packages(),
      zip_safe=False,
      classifiers=[
        "Programming Language :: Python",
        "Development Status :: 1 - Planning",
        "License :: None",
        "Natural Language :: French",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2.7",
        "Topic :: Mathematics, spline",
    ],)
