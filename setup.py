# -*- coding: latin-1 -*-
'''
Notes:

- recommend the users to use setup.py develop when tracking in-development code, and tell them that this needs to be run after every update or commit
- when removing modules or data files from the project, do remind the users to run setup.py clean --all and delete any obsolete .pyc or .pyo.


'''
import ez_setup
ez_setup.use_setuptools() # TODO: check that it works with the current installed version of setuptools

from setuptools import setup, find_packages

setup(
    name = "CN-Wheat",
    version = "0.0.1",
    package_dir = {'': 'src'},
    packages = find_packages('src'),
    
    install_requires = ['numpy', 'pandas', 'scipy', 'matplotlib'], 
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
