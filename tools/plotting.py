# -*- coding: latin-1 -*-
"""
    plotting
    ~~~~~~~~

    Tools to plot the outputs of CN-Wheat.

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

import types

import matplotlib.pyplot as plt

from cnwheat.simulation import CNWheat

def plot(outputs, x_name, y_name, x_label='', y_label='', title='', filters={}, plot_filepath=None):
    """Plot `outputs`, with `x=x_name` and `y=y_name`.
    
    The general algorithm is: 
    
        * find the scale of `outputs` and keep only the needed columns, 
        * apply `filters` to `outputs` and make groups according to the scale,
        * plot each group as a new line,
        * save or display the plot.
    
    :Parameters:

        - `outputs` (:class:`pandas.DataFrame`) - The outputs of CN-Wheat.

        - `x_name` (:class:`str`) - x axis of the plot.
        
        - `y_name` (:class:`str`) - y axis of the plot.

        - `x_label` (:class:`str`) - The x label of the plot. Default is ''.

        - `y_label` (:class:`str`) - The y label of the plot. Default is ''.

        - `title` (:class:`str`) - the title of the plot. Default is ''.
        
        - `filters` (:class:`dict`) - A dictionary whose keys are the columns of 
          `outputs` for which we want to apply a specific filter. 
          The value associated to each key is a criteria that the rows of `outputs` 
          must satisfy to be plotted. The values can be either one value or a list of values. 
          If no value is given for any column, then all rows are plotted (default).

        - `plot_filepath` (:class:`str`) - The file path to save the plot. 
          If `None`, do not save the plot but display it.

    :Examples:

    >>> import pandas as pd
    >>> cnwheat_output_df = pd.read_csv('cnwheat_output.csv') # in this example, 'cnwheat_output.csv' must contain at least the columns 't' and 'Conc_Sucrose'.
    >>> plot(cnwheat_output_df, x_name = 't', y_name = 'Conc_Sucrose', x_label='Time (Hour)', y_label=u'[Sucrose] (µmol g$^{-1}$ mstruct)', title='{} = f({})'.format('Conc_Sucrose', 't'), filters={'plant': 1, 'axis': 'MS', 'organ': 'Lamina', 'element': 1})

    """
    
    # finds the scale of `outputs`
    group_keys = [key for key in CNWheat.ELEMENTS_INDEXES if key in outputs and key != x_name and key != y_name]
    
    # keep only the needed columns (to make the grouping faster)
    outputs = outputs[group_keys + [x_name, y_name]]
    
    # apply filters to outputs
    for key, value in filters.iteritems():
        if key in outputs:
            # convert to list if needed
            try:
                _ = iter(value)
            except TypeError:
                values = [value]
            else:
                values = value
                # handle strings too
                if isinstance(values, types.StringTypes):
                    values = [values]
            # select data from outputs
            outputs = outputs[outputs[key].isin(values)]
            
    # makes groups according to the scale
    outputs_grouped = outputs.groupby(group_keys)
    
    # plots each group as a new line
    plt.figure()
    ax = plt.subplot(111)
    
    for outputs_group_name, outputs_group in outputs_grouped:
        line_label = '_'.join([str(key) for key in outputs_group_name])
        ax.plot(outputs_group[x_name], outputs_group[y_name], **{'label': line_label})
        
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(prop={'size':10}, framealpha=0.5)
    ax.set_title(title)

    if plot_filepath is None:
        # display the plot
        plt.show()
    else:
        # save the plot
        plt.savefig(plot_filepath, dpi=200, format='PNG')
        plt.close()
