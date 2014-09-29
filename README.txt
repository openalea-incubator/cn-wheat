====== CN-Wheat ======

**Authors** : C.Chambon and R.Barillot

**Institutes** : INRA

**Type** : Python

**Status** : Release candidate

**Version** : 0.0.1

**License** : TODO

**URL** : https://sourcesup.renater.fr/projects/cn-wheat/


===== About =====

CN-Wheat is a model of CN distribution for wheat.


=== Description ===

Modèle de distribution spatiale de l'azote et du carbone chez le blé


=== Content ===

The CN-Wheat package contains : 

* an implementation of the model,
* post processing which can be applied on CN-Wheat output. 


=== Requirements ===

* Python >= 2.7, http://www.python.org/
* NumPy >= 1.7.2, http://www.numpy.org/
* SciPy >= 0.12.1, http://docs.scipy.org/
* Pandas >= 0.14.0, http://pandas.pydata.org/
* Matplotlib >= 1.2.0, http://matplotlib.sourceforge.net/
* To run the tests: Nose >= 1.3.0, http://nose.readthedocs.org/
* To build the documentation: Sphinx >= 1.1.3, http://sphinx-doc.org/


=== Get the sources ===

svn checkout https://subversion.renater.fr/cn-wheat

This creates the directory ``cn-wheat``.


=== Install ===

Go to directory "cn-wheat" and run: python setup.py install (or python setup.py develop) 


=== Run the tests (after installation) ===

Go to directory "cn-wheat" and run: nosetests


=== Build the documentation ===

Go to directory "cn-wheat" and run: python setup.py build_sphinx

Then open the file :download:`doc/_build/html/index.html <_build/html/index.html>` 
with your web browser.


=== Mailing list ===

cn-wheat@groupes.renater.fr


=== Bug reports ===

https://sourcesup.renater.fr/tracker/?atid=5057&group_id=1515&func=browse

