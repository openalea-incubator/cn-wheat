# -*- coding: latin-1 -*-
"""
    Linear regression tool
    ~~~~~~~~~~~~~~~~~~~~~~

    Perform a linear regression of two arrays and create a plot showing
    the fit against the original data.
    Display the plot or save it to a PNG file.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def plot_linear_regression(x_array, y_array, x_label='x', y_label='y', plot_filepath=None):
    """Perform a linear regression of `x_array` vs `y_array`
    and create a plot showing the fit against the original data.
    If `plot_filepath` is not None, save the plot to a PNG file. Otherwise display the plot.

    This is derived from http://www.landmap.ac.uk/index.php/Learning-Materials/Python-Scripting/6.4-Fitting-linear-equations#sthash.wDZ5zBrD.dpuf,
    which is: Copyright TODO

    :Parameters:

        - `x_array` (:class:`numpy.ndarray`) - The x.

        - `y_array` (:class:`numpy.ndarray`) - The y.

        - `x_label` (:class:`str`) - The label of the axis 'x'. Default is 'x'.

        - `y_label` (:class:`str`) - The label of the axis 'y'. Default is 'y'.

        - `plot_filepath` (:class:`str`) - The file path to save the plot in.
            If `None`, do not save the plot.

    :Examples:

    >>> import pandas as pd
    >>> modelmaker_output_df = pd.read_csv('modelmaker_output.csv') # 'modelmaker_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # 'cnwheat_output.csv' must contain at least the column 'Sucrose_Phloem'
    >>> plot_linear_regression(modelmaker_output_df.Sucrose_Phloem,
                               cnwheat_output_df.Sucrose_Phloem,
                               x_label='modelmaker_{}'.format('Sucrose_Phloem'),
                               y_label='cnwheat_{}'.format('Sucrose_Phloem'),
                               plot_filepath='compare.png')

    """
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
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()

