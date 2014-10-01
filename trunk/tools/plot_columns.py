# -*- coding: latin-1 -*-
"""
    Plot dataframe tool
    ~~~~~~~~~~~~~~~~~~~

    Plot each column of a dataframe using the column 't' as x, and save the 
    plot to a PNG file.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import matplotlib.pyplot as plt

def plot_dataframe(dataframe, title='', column_to_matplotlib_kwargs={}, plot_filepath=None):
    """Plot the columns in `column_to_matplotlib_kwargs`, with `x_label='t'` and 
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
        
    """
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
        
        