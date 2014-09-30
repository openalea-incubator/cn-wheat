# -*- coding: latin-1 -*-
"""
    cnwheat.post_processing
    ~~~~~~~~~~~~~~~~~~~~~~~

    Post processings to apply on cnwheat output.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.txt.
    :license: TODO, see LICENSE.txt for details.
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def plot_linear_regression(x_array, y_array, x_label='x', y_label='y', plot_filepath=None):
    '''Perform a linear regression of `x_array` vs `y_array` 
    and create a plot showing the fit against the original data.
    If `plot_filepath` is not None, save the plot to a PNG file.
    
    This is derived from http://www.landmap.ac.uk/index.php/Learning-Materials/Python-Scripting/6.4-Fitting-linear-equations#sthash.wDZ5zBrD.dpuf, 
    which is: Copyright TODO
    
    :Parameters:
    
        - `x_array` (:class:`numpy.ndarray`) - The x.
        
        - `y_array` (:class:`numpy.ndarray`) - The y.
        
        - `x_label` (:class:`str`) - The label of the axis 'x'. Default is 'x'.
        
        - `y_label` (:class:`str`) - The label of the axis 'y'. Default is 'y'.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. 
            If `None`, do not save the plot.
        
    '''
    # Perform fit
    (aCoeff, bCoeff, rVal, pVal, stdError) = stats.linregress(x_array, y_array)
    
    # Use fits to predict y output for a range of diameters
    x_samples_array = np.linspace(min(x_array), max(x_array), 1000)
    y_predict_array = aCoeff * x_samples_array + bCoeff

    # Create a string, showing the form of the equation (with fitted coefficients) and r squared value.
    # Coefficients are rounded to two decimal places.
    equation = 'y = {} x + {} (R$^2$ = {})'.format(round(aCoeff,2), round(bCoeff,2), round(rVal**2,2))
    
    plt.figure()
    
    # Plot fit against original data
    plt.plot(x_array, y_array,'.')
    plt.plot(x_samples_array, y_predict_array)
    plt.title('{} vs {}'.format(x_label, y_label))
    
    x_label = 'x = {}'.format(x_label)
    plt.xlabel(x_label)
    y_label = 'y = {}'.format(y_label)
    plt.ylabel(y_label)
    
    plt.legend(['x vs y', equation])

    # Save plot
    if plot_filepath is not None:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        
        
def plot_column(x, y, x_label='', y_label='', title='', matplotlib_kwargs={}, plot_filepath=None):
    '''Plot y = f(x), using `matplotlib_kwargs` to set the `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `plot_filepath` is not None, save the plot to a PNG file.
    
    :Parameters:
    
        - `x` (:class:`numpy.ndarray`) - The x.
        
        - `y` (:class:`numpy.ndarray`) - The y.
        
        - `title` (:class:`str`) - the title of the plot.
        
        - `matplotlib_kwargs` (:class:`dict`) - The `kwargs` of :func:`matplotlib.pyplot.plot`. 
            If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. If `None`, do not save the plot.
        
    '''
    matplotlib_kwargs['label'] = y_label
    plt.figure()
    plt.plot(x, y, **matplotlib_kwargs)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.title(title)
    
    # Save plot
    if plot_filepath is not None:
        plt.savefig(plot_filepath, dpi=200, format='PNG')


def plot_columns(dataframe, title='', column_to_matplotlib_kwargs={}, plot_filepath=None):
    '''Plot the columns in `column_to_matplotlib_kwargs`, with `x_label='t'` and 
    `y_labels=column_to_matplotlib_kwargs.keys()`.
    For each `column` in `column_to_matplotlib_kwargs`, use `kwargs=column_to_matplotlib_kwargs[column]` 
    to set the `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `column_to_matplotlib_kwargs` is empty, then plot all the columns of `dataframe` and use
    the default `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `plot_filepath` is not None, save the plot to a PNG file.
    
    :Parameters:
    
        - `dataframe` (:class:`pandas.DataFrame`) - A dataframe with a column 't'.
        
        - `title` (:class:`str`) - the title of the plot.
        
        - `column_to_matplotlib_kwargs` (:class:`dict`) - A dictionary which keys are 
            the name of the columns to plot, and values are the `kwargs` of :func:`matplotlib.pyplot.plot`. 
            If `column_to_matplotlib_kwargs` is empty, then plot all the columns of `dataframe` and use
            the default `kwargs` of :func:`matplotlib.pyplot.plot`.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. If `None`, do not save the plot.
        
    '''
    if len(column_to_matplotlib_kwargs) == 0:
        column_to_matplotlib_kwargs = dict.fromkeys(dataframe.columns, {})
    
    x_label = 't'
    x = dataframe[x_label]

    plt.figure()
    
    for (column, matplotlib_kwargs) in column_to_matplotlib_kwargs.items():
        y = dataframe[column]
        if 'label' not in matplotlib_kwargs:
            matplotlib_kwargs['label'] = column
        plt.plot(x, y, **matplotlib_kwargs)
        
    plt.xlabel(x_label)
    plt.ylabel('See the legend')
    plt.legend()
    plt.title(title)
    
    # Save plot
    if plot_filepath is not None:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
         
         

         
         
         
         
    
    