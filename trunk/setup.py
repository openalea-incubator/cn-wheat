# -*- coding: latin-1 -*-
"""
Notes:

- use setup.py develop when tracking in-development code
- when removing modules or data files from the project, run setup.py clean --all and delete any obsolete .pyc or .pyo.

"""
import ez_setup
ez_setup.use_setuptools()

import sys
from setuptools import setup, find_packages

import cnwheat

if sys.version_info < (2, 7):
    print('ERROR: CN-Wheat requires at least Python 2.7 to run.')
    sys.exit(1)
    
if sys.version_info >= (3, 0):
    print('WARNING: CN-Wheat has not been tested with Python 3.')

setup(
    name = "CN-Wheat",
    version=cnwheat.__version__,
    packages = find_packages(),
    
    install_requires = ['numpy>=1.7.2', 'pandas>=0.14.0', 'scipy>=0.12.1'], 
    include_package_data = True, 
    
    # metadata for upload to PyPI
    author = "C.Chambon, R.Barillot",
    author_email = "camille.chambon@grignon.inra.fr, romain.barillot@grignon.inra.fr",
    description = "Model of CN distribution for wheat",
    long_description = "Modèle de distribution spatiale de l'azote et du carbone chez le blé",
    license = "", # TODO
    keywords = "", # TODO
    url = "https://sourcesup.renater.fr/projects/cn-wheat/",
    download_url = "", # TODO    
)
