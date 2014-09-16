# -*- coding: latin-1 -*-
'''
Created on 2 sept. 2014

@author: cchambon, adapted from http://www.landmap.ac.uk/index.php/Learning-Materials/Python-Scripting/6.4-Fitting-linear-equations#sthash.wDZ5zBrD.dpuf
'''

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def plot_linear_regression(ModelMaker_array, python_array, plot_filepath=None):
    '''Perform a linear regression of `ModelMaker_array` vs `python_array` 
    and create a plot showing the fit against the original data.
    If `plot_filepath` is not None, save the plot to a PNG file.
    
    This function permits to compare the output computed using ModelMaker and the 
    output computed using :func:`scipy.integrate.odeint`.
    
    :Parameters:
    
        - `ModelMaker_array` (:class:`numpy.ndarray`) - The output computed using ModelMaker.
        
        - `python_array` (:class:`numpy.ndarray`) - The output computed using :func:`scipy.integrate.odeint`.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. 
        
    '''
    # Perform fit
    (aCoeff, bCoeff, rVal, pVal, stdError) = stats.linregress(ModelMaker_array, python_array)
    
    # Use fits to predict Python output for a range of diameters
    ModelMaker_samples_array = np.linspace(min(ModelMaker_array), max(ModelMaker_array), 1000)
    python_predict_array = aCoeff * ModelMaker_samples_array + bCoeff

    # Create a string, showing the form of the equation (with fitted coefficients) and r squared value.
    # Coefficients are rounded to two decimal places.
    equation = str(round(aCoeff,2)) + 'x + ' + str(round(bCoeff,2)) + ' (r$^2$ = ' + str(round(rVal**2,2)) + ')'

    # Plot fit against original data
    plt.plot(ModelMaker_array, python_array,'.')
    plt.plot(ModelMaker_samples_array, python_predict_array)
    plt.xlabel('ModelMaker')
    plt.ylabel('Python')
    plt.legend(['ModelMaker vs Python', equation])

    # Save plot
    if plot_filepath is not None:
        plt.savefig(plot_filepath, dpi=200, format='PNG')

