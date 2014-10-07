# -*- coding: latin-1 -*-
"""
    Plot dataframe tool
    ~~~~~~~~~~~~~~~~~~~

    Plot columns of a dataframe. Display the plot or save it to a PNG file.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import matplotlib.pyplot as plt

def plot_dataframe(dataframe, x_label='', y_label='', title='', column_to_matplotlib_kwargs={}, plot_filepath=None):
    """Plot the columns in `column_to_matplotlib_kwargs`, with `x_label='t'` and
    `y_labels=column_to_matplotlib_kwargs.keys()`.
    For each `column` in `column_to_matplotlib_kwargs`, use `kwargs=column_to_matplotlib_kwargs[column]`
    to set the `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `column_to_matplotlib_kwargs` is empty, then plot all the columns of `dataframe` and use
    the default `kwargs` of :func:`matplotlib.pyplot.plot`.
    If `plot_filepath` is not None, save the plot to a PNG file. Otherwise display the plot.

    :Parameters:

        - `dataframe` (:class:`pandas.DataFrame`) - A dataframe with a column 't'.

        - `x_label` (:class:`str`) - The *x* label of the plot.

        - `y_label` (:class:`str`) - The *y* label of the plot.

        - `title` (:class:`str`) - the title of the plot.

        - `column_to_matplotlib_kwargs` (:class:`dict`) - A dictionary which keys are
            the name of the columns to plot, and values are the `kwargs` of :func:`matplotlib.pyplot.plot`.
            If `column_to_matplotlib_kwargs` is empty, then plot all the columns of `dataframe` and use
            the default `kwargs` of :func:`matplotlib.pyplot.plot`.

        - `plot_filepath` (:class:`str`) - The file path to save the plot in. If `None`, do not save the plot.

    :Examples:

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # 'cnwheat_output.csv' must contain at least the columns 't', 'SUCROSE_lamina' and 'SUCROSE_phloem'
    >>> plot_dataframe(cnwheat_output_df,
                       x_label='t',
                       y_label='SUCROSE',
                       title='{} = f({})'.format('SUCROSE', 't'),
                       column_to_matplotlib_kwargs={'SUCROSE_lamina': {'color': 'green', 'linestyle': 'solid', 'marker': 'o', 'markerfacecolor': 'blue', 'markersize': 12},
                                                    'SUCROSE_phloem': {'color': 'red', 'linestyle': 'dashed', 'marker': '*', 'markerfacecolor': 'yellow', 'markersize': 16}},
                       plot_filepath='SUCROSE.png')

    """
    if len(column_to_matplotlib_kwargs) == 0:
        column_to_matplotlib_kwargs = dict.fromkeys(dataframe.columns, {})

    x = dataframe[x_label]

    plt.figure()
    ax = plt.subplot(111)

    for (column, matplotlib_kwargs) in column_to_matplotlib_kwargs.items():
        y = dataframe[column]
        if 'label' not in matplotlib_kwargs:
            matplotlib_kwargs['label'] = column
        ax.plot(x, y, **matplotlib_kwargs)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(loc='upper center',  bbox_to_anchor=(0.5, 1.15))
    ax.set_title(title)

    # Save plot
    if plot_filepath is None:
        plt.show()
    else:
        plt.savefig(plot_filepath, dpi=200, format='PNG')
    plt.close()
