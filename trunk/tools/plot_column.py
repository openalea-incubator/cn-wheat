# -*- coding: latin-1 -*-
"""
    Plot series tool
    ~~~~~~~~~~~~~~~~

    Plot a series of data. Display the plot or save it to a PNG file.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import matplotlib.pyplot as plt

def plot_series(x, y, x_label='', y_label='', title='', matplotlib_kwargs={}, plot_filepath=None):
    """Plot `y` = f(`x`), using `matplotlib_kwargs` to set the `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `plot_filepath` is not None, save the plot to a PNG file. Otherwise display the plot.
    
    :Parameters:
    
        - `x` (:class:`numpy.ndarray`) - The x.
        
        - `y` (:class:`numpy.ndarray`) - The y.
        
        - `x_label` (:class:`str`) - The *x* label of the plot.
        
        - `y_label` (:class:`str`) - The *y* label of the plot.
        
        - `title` (:class:`str`) - the title of the plot.
        
        - `matplotlib_kwargs` (:class:`dict`) - The `kwargs` of :func:`matplotlib.pyplot.plot`. 
            If `matplotlib_kwargs` is empty, then use the default `kwargs` of :func:`matplotlib.pyplot.plot`.
        
        - `plot_filepath` (:class:`str`) - The file path to save the plot in. If `None`, do not save the plot.
        
    :Examples:        

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose_Phloem'
    >>> plot_series(x=cnwheat_output_df.t,
                    y=cnwheat_output_df.Conc_Sucrose_Phloem,
                    x_label='t', 
                    y_label='Conc_Sucrose_Phloem',
                    title='{} = f({})'.format('Conc_Sucrose_Phloem', 't'), 
                    matplotlib_kwargs={'color': 'green', 'linestyle': 'solid', 'marker': 'o', 'markerfacecolor': 'blue', 'markersize': 12},
                    plot_filepath='Conc_Sucrose_Phloem.png')
        
    """
    if 'label' not in matplotlib_kwargs:
        matplotlib_kwargs['label'] = y_label
    plt.figure()
    plt.plot(x, y, **matplotlib_kwargs)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.title(title)
    
    # Save plot
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()
        
