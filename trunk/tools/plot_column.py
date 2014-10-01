# -*- coding: latin-1 -*-
"""
    Plot series tool
    ~~~~~~~~~~~~~~~~

    Plot a series of data and save the plot to a PNG file.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import matplotlib.pyplot as plt

def plot_series(x, y, x_label='', y_label='', title='', matplotlib_kwargs={}, plot_filepath=None):
    """Plot `y` = f(`x`), using `matplotlib_kwargs` to set the `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `plot_filepath` is not None, save the plot to a PNG file.
    
    :Parameters:
    
        - `x` (:class:`numpy.ndarray`) - The x.
        
        - `y` (:class:`numpy.ndarray`) - The y.
        
        - `title` (:class:`str`) - the title of the plot.
        
        - `matplotlib_kwargs` (:class:`dict`) - The `kwargs` of :func:`matplotlib.pyplot.plot`. 
            If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. If `None`, do not save the plot.
        
    """
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
        
        