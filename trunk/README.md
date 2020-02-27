# CN-Wheat

This is the Read Me file of the model *CN-Wheat*, a model of CN distribution for wheat.

*CN-Wheat* is a Functional-Structural Plant Model which simulates the distribution 
of carbon and nitrogen into wheat culms in relation to photosynthesis, 
N uptake, metabolite turnover, root exudation and tissue death. 

This model was first produced as part of the project BreedWheat over the last 
three years, through the Investment for the Future programme managed by the 
French Research National Agency (ANR-10-BTBR-03). The aim of the project BreedWheat was 
to improve the competiveness of the French wheat breeding sector, through the 
definition/ identification of ideotypes, parameters of interest maximizing grain yield 
and quality under sustainable agricultural systems and climate scenarios. 

These researches lead to the publication of project report, and two articles:

* Barillot, R., Chambon, C., & Andrieu, B. (2016). CN-Wheat, a functional–structural model 
  of carbon and nitrogen metabolism in wheat culms after anthesis. I. Model description. 
  Annals of Botany, 118(5), 997‑1013. https://doi.org/10.1093/aob/mcw143 
* and Barillot, R., Chambon, C., & Andrieu, B. (2016). CN-Wheat, a functional–structural 
  model of carbon and nitrogen metabolism in wheat culms after anthesis. II. Model evaluation. 
  Annals of Botany, 118(5), 1015‑1031. https://doi.org/10.1093/aob/mcw144

## 1. Getting Started

These instructions will get you a copy of *CN-Wheat* up and running on your local 
machine.

### 1.1 Prerequisites

To install and use *CN-Wheat*, you need first to install the dependencies.

*CN-Wheat* has been tested on Windows 10 64 bit and Linux Fedora 24 64 bit.
 
#### 1.1.1 Install the dependencies on Windows 10 64 bit

1. Install Python  

    * go to https://www.python.org/downloads/windows/download, 
    * click on "Latest Python 2 Release [...]", 
    * download "Windows x86-64 MSI installer" and install it selecting the following options:
        * install for all users,
        * default destination directory,
        * install all subfeatures, including subfeature "Add python.exe to Path".

2. Install NumPy:  

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy, 
    * download `NumPy+MKL` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `NumPy+MKL` for Python 2 64 bit,
        * install `NumPy+MKL` from the downloaded wheel file.  
          For example, if you downloaded file "numpy‑1.13.1+mkl‑cp27‑cp27m‑win_amd64.whl", 
          type: `pip install "numpy‑1.13.1+mkl‑cp27‑cp27m‑win_amd64.whl"`.

3. Install SciPy  

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy, 
    * download `SciPy` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `SciPy` for Python 2 64 bit,
        * install `SciPy` from the downloaded wheel file.  
          For example, if you downloaded file "scipy‑0.19.1‑cp27‑cp27m‑win_amd64.whl", 
          type: `pip install "scipy‑0.19.1‑cp27‑cp27m‑win_amd64.whl"`.

4. Install Pandas  

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#pandas, 
    * download `Pandas` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Pandas` for Python 2 64 bit,
        * install `Pandas` from the downloaded wheel file.  
          For example, if you downloaded file "pandas‑0.20.3‑cp27‑cp27m‑win_amd64.whl", 
          type: `pip install "pandas‑0.20.3‑cp27‑cp27m‑win_amd64.whl"`.
          
5. Install Matplotlib

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib, 
    * download `Matplotlib` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Matplotlib` for Python 2 64 bit,
        * install `Matplotlib` from the downloaded wheel file.  
          For example, if you downloaded file "matplotlib‑2.0.2‑cp27‑cp27m‑win_amd64.whl", 
          type: `pip install "matplotlib‑2.0.2‑cp27‑cp27m‑win_amd64.whl"`.
          
6. Install Sphinx

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#misc, 
    * download `Sphinx` for Python 2,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Sphinx` for Python 2,
        * install `Sphinx` from the downloaded wheel file.  
          For example, if you downloaded file "Sphinx‑1.6.3‑py2.py3‑none‑any.whl", 
          type: `pip install "Sphinx‑1.6.3‑py2.py3‑none‑any.whl"`.
          
7. Install Nose

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#misc, 
    * download `Nose` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Nose` for Python 2 64 bit,
        * install `Nose` from the downloaded wheel file.  
          For example, if you downloaded file "nose‑1.3.7‑py2‑none‑any.whl", 
          type: `pip install "nose‑1.3.7‑py2‑none‑any.whl"`.
          
8. Install Coverage

    * go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#coverage, 
    * download `Coverage` for Python 2 64 bit,
    * install it using `pip` installer: 
        * open a command line interpreter,
        * go to the directory where you saved `Coverage` for Python 2 64 bit,
        * install `Coverage` from the downloaded wheel file.  
          For example, if you downloaded file "coverage‑4.4.1‑cp27‑cp27m‑win_amd64.whl", 
          type: `pip install "coverage‑4.4.1‑cp27‑cp27m‑win_amd64.whl"`.
          
9. Install Respi-Wheat

    * download the lastest public release of model *Respi-Wheat* from https://sourcesup.renater.fr/frs/download.php/latestzip/2087/Respi-Wheat-Stable-latest.zip 
      and install it:
        * unzip it: you should obtain a zip file `respi-wheat_*.zip`,
        * unzip the zip file `respi-wheat_*.zip`: you should obtain a folder `respi-wheat`,
        * open a command line interpreter and go to folder `respi-wheat`,
        * run command: `python setup.py install --user`.

On Windows 10 64 bit, *CN-Wheat* has been tested with the following versions of the dependencies:  

* Python 2.7.13 64 bit,
* NumPy+MKL 1.13.1 64 bit,
* SciPy 0.19.1 64 bit,
* Pandas 0.20.3 64 bit,
* Matplotlib 2.0.2 64 bit,
* Sphinx 1.6.3,
* Nose 1.3.7,
* Coverage 4.4.1 64 bit.

#### 1.1.2 Install the dependencies on Linux Fedora 24 64 bit

To install the dependencies on Linux Fedora 24 64 bit:

* open a terminal,
* run this command with superuser privileges: `dnf -y install python2 python2-numpy python2-scipy python2-pandas python2-matplotlib python2-sphinx python2-nose python2-coverage`
* download the lastest public release of model *Respi-Wheat* from https://sourcesup.renater.fr/frs/download.php/latestzip/2087/Respi-Wheat-Stable-latest.zip and install it:
    * unzip file `Respi-Wheat-Stable-latest.zip`: you should obtain a zip file `respi-wheat_*.zip`,
    * unzip file `respi-wheat_*.zip`: you should obtain a folder `respi-wheat`,
    * go to folder `respi-wheat`,
    * run command: `python setup.py install --user`.

On Linux Fedora 24 64 bit, *CN-Wheat* has been tested with the following versions of the dependencies:  

* Python 2.7.13 64 bit,
* NumPy 1.11.0 64 bit,
* SciPy 0.16.1 64 bit,
* Pandas 0.18.0 64 bit,
* Matplotlib 1.5.2rc2 64 bit,
* Sphinx 1.4.8,
* Nose 1.3.7,
* Coverage 4.4.1 64 bit.


### 1.2 Installing

__Note__: We suppose you already installed the dependencies for your operating system. Otherwise follow these [instructions](prerequisites "Prerequisites").

You can install *CN-Wheat* either in "install" or "develop" mode.

#### 1.2.1 Install *CN-Wheat* in "install" mode

Install *CN-Wheat* in "install" mode if you're not going to develop, edit or debug 
it, i.e. you just want to used it as third party package.

To install *CN-Wheat* in "end-user" mode:

* open a command line interpreter,
* go to your local copy of project *CN-Wheat*,
* run command: `python setup.py install --user`.

#### 1.2.2 Install *CN-Wheat* in "develop" mode

Install *CN-Wheat* in "develop" mode if you want to get *CN-Wheat* installed and then 
be able to frequently edit the code and not have to re-install *CN-Wheat* to have the 
changes to take effect immediately.

To install *CN-Wheat* in "develop" mode:

* open a command line interpreter,
* go to your local copy of project *CN-Wheat*,
* run command: `python setup.py develop --user`.

### 1.3 Running

__Note__: We suppose you already installed the model. Otherwise follow these [instructions](installing "Installing").

To run a simulation example, compute post-processing and generate graphs for validation:

* open a command line interpreter,
* go to the directory `example/` of your local copy of project *CN-Wheat*,
* run command: `python main.py`.

See the user guide for a step by step explanation of how to set and run model *CN-Wheat*.

## 2. Reading the docs

To build the user and reference guides:

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the top directory of your local copy of the project,
* run this command: `python setup.py build_sphinx`,
* and direct your browser to file `doc/_build/html/index.html`.

## 3. Testing

The automated test permits to verify that the model implementation accurately 
represents the developer’s conceptual description of the model and its solution.

The automated test:

* initializes the model from input data in CSV files,
* runs the model on 2 steps, forcing the photosynthesis and senescence parameters 
  before each run of the model,
* concatenate the outputs of the model in dataframes, with one dataframe per topological scale,
* write the outputs dataframes to CSV files,
* compare actual to expected outputs,
* raise an error if actual and expected outputs are not equal up to a given tolerance.     

To run the automated test with coverage report:

* install the model (see [Installation of the model](installing "Installing")), 
* open a command line interpreter,
* go to the directory `test` of your local copy of the project,
* and run this command: `nosetests --with-coverage --cover-package=cnwheat test_cnwheat.py`.

The automated test does not verify the validity of the model, i.e. it doesn't permit 
to determine the degree to which the model is an accurate representation of the 
real world from the perspective of the intended uses of the model.  
To help verifying the validity of the model, use the plotting tools implemented 
in module `cnwheat.tools`.   

## Deployment

*CN-Wheat* can be coupled with other ecophysiological models, to simulate the interaction 
between CN distribution and (for example) leaves elongation, photosynthesis, growth, 
senescence, light interception and topology of wheat crops.  
Please contact <cn-wheat-request@groupes.renater.fr> for more information about the 
possibility of coupling and integrate *CN-Wheat* with other ecophysiological models. 

## Built With

* [Python](http://www.python.org/), [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/), 
  [Pandas](http://pandas.pydata.org/), [Respi-Wheat](https://sourcesup.renater.fr/projects/respi-wheat): 
  implementation and deployment of the model,
* [Matplotlib](http://matplotlib.org/): generation of graphs to help validating the model, 
* [Sphinx](http://sphinx-doc.org/): building of the documentation, 
* [Nose](http://nose.readthedocs.org/): run of the automated tests,
* [Coverage](http://nedbatchelder.com/code/coverage/): coverage of code testing.

## Contributing

First, send an email to <cn-wheat-request@groupes.renater.fr> to be added to the project.  

Then,
 
* check for open issues or open a fresh issue to start a discussion around a
  feature idea or a bug: https://sourcesup.renater.fr/tracker/?group_id=1515.
* If you feel uncomfortable or uncertain about an issue or your changes, feel
  free to email <cn-wheat@groupes.renater.fr>.

## Contact

For any question, send an email to <cn-wheat-request@groupes.renater.fr>.

## Versioning

We use an SVN repository on [SourceSup](https://sourcesup.renater.fr) for 
versioning: https://sourcesup.renater.fr/projects/cn-wheat/.  
If you need an access to the current in development version of the model, please send 
an email to <cn-wheat-request@groupes.renater.fr>.

## Authors

See file [AUTHORS](AUTHORS) for details

## License

This project is licensed under the CeCILL-C License - see file [LICENSE](LICENSE) for details

## Acknowledgments

The research leading these results has received funding through the 
Investment for the Future programme managed by the Research National Agency 
(BreedWheat project ANR-10-BTBR-03).
