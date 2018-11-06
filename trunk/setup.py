# -*- coding: latin-1 -*-
"""

    setup
    ~~~~~
    
    Setup script for installation.
    
    See README.md for installing procedure.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
    .. seealso:: Barillot et al. 2016.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import ez_setup
import pkg_resources

ez_setup.use_setuptools()

import sys
from setuptools import setup, find_packages

import cnwheat

if sys.version_info < (2, 7):
    print('ERROR: CN-Wheat requires at least Python 2.7 to run.')
    sys.exit(1)

if sys.version_info >= (3, 0):
    print('WARNING: CN-Wheat has not been tested with Python 3.')

pkg_resources.require('numpy>=1.15.1', 'pandas>=0.23.4', 'scipy>=1.1.0',
                      'Respi-Wheat')

setup(
    name = "CN-Wheat",
    version=cnwheat.__version__,
    packages = find_packages(),
    include_package_data=True,
    author = "R.Barillot, C.Chambon and B.Andrieu",
    author_email = "romain.barillot@inra.fr, camille.chambon@inra.fr, bruno.andrieu@inra.fr",
    description = "CN-Wheat is a model of CN distribution for wheat",
    long_description = """CN-Wheat is a Functional-Structural Plant Model which 
simulates the distribution of carbon and nitrogen into wheat culms in relation 
to photosynthesis, N uptake, metabolite turnover, root exudation and tissue death.""",
    license = "CeCILL-C",
    keywords = "functional-structural plant model, wheat, ode, system integration, scipy", 
    url = "https://sourcesup.renater.fr/projects/cn-wheat/",
    download_url = "https://sourcesup.renater.fr/frs/download.php/latestzip/2088/CN-Wheat-Stable-latest.zip",
)
