# -*- coding: latin-1 -*-

'''
Created on 28 juil. 2014

@author: cchambon
'''

from __future__ import division # use "//" to do integer division

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d

# parameters
ALPHA_AXIS = 1
ALPHA_LAMINA = 1
DELTA_DSTORAGE = 0.0001
K_EAR = 2
K_STORAGE = 20
K_SUCROSE = 0.66
LAMINA_AREA = 0.00346
MSTRUCT_AXIS = 2.08
MSTRUCT_LAMINA = 0.14
R = 20000000
VMAX_EAR = 0.026
VMAX_STORAGE = 2
VMAX_SUCROSE = 1

# variables
def photosynthesis(An, lamina_area):
    return An * lamina_area


def D_storage(storage):
    '''Flow from STORAGE to SUCROSE_lamina
    '''
    return max(0, DELTA_DSTORAGE * (storage/(MSTRUCT_LAMINA*ALPHA_LAMINA)))  


def S_storage(triosesp):
    '''Flow from TRIOSESP to STORAGE
    '''
    return ((max(0, triosesp)/(MSTRUCT_LAMINA*ALPHA_LAMINA)) * VMAX_STORAGE) / ((max(0, triosesp)/(MSTRUCT_LAMINA*ALPHA_LAMINA)) +K_STORAGE)
    
    
def S_sucrose(triosesp):
    '''Flow from TRIOSESP to SUCROSE_lamina
    '''
    return ((max(0,triosesp)/(MSTRUCT_LAMINA*ALPHA_LAMINA)) * VMAX_SUCROSE) / ((max(0, triosesp)/(MSTRUCT_LAMINA*ALPHA_LAMINA)) +K_SUCROSE)

def func(y, t0, An_linear_interpolation):
    '''Compute the derivative of y at t0.
    '''
    storage_i = y[0]
    sucrose_lamina_i = y[1]
    triosesp_i = y[2]
    # the model equations (see Barillot et al. 2014)
    storage = (S_storage(triosesp_i) - D_storage(storage_i)) * (MSTRUCT_LAMINA*ALPHA_LAMINA) # µmol of C
    sucrose_lamina = (S_sucrose(triosesp_i) + D_storage(storage_i)) * (MSTRUCT_LAMINA*ALPHA_LAMINA) # µmol of C
    triosesp = photosynthesis(An_linear_interpolation(t0), LAMINA_AREA) - (S_sucrose(triosesp_i) + S_storage(triosesp_i)) * (MSTRUCT_LAMINA*ALPHA_LAMINA) # µmol of C
    return [storage, sucrose_lamina, triosesp]
    
    
def run(start_time, stop_time, number_of_output_steps, initial_conditions, An):
    '''
    Compute STORAGE, SUCROSE_lamina and TRIOSESP between `start_time` and `stop_time`, for y(t=0) = `initial_conditions`. 
    
    The photosynthesis is computed from the assimilation `An` and LAMINA_AREA.
    
    Parameters
    ----------
    start_time: int
        The starting of the time grid.
    stop_time: int
        The end of the time grid.
    number_of_output_steps: float, optional
        Number of time points for which to solve for y. 
    initial_conditions: list
        List of initial conditions: STORAGE(t=0), SUCROSE_lamina(t=0), TRIOSESP(t=0). 
    An: pandas.Series
        The assimilation of the leaf at each hour t (µmol m-2 s-1).
     
    Returns
    -------
    out : pandas.DataFrame
        Dataframe containing the value of Photosynthesis, STORAGE, SUCROSE_lamina and TRIOSESP 
        for each desired time.
    
    '''
    An_linear_interpolation = interp1d(An.index, An, bounds_error=False, 
                                       fill_value=An[An.last_valid_index()]) # odeint will evaluate An_linear_interpolation at time values past the last requested time: fill them by the last value of An
    t = np.linspace(start_time, stop_time, number_of_output_steps)
    soln = odeint(func, initial_conditions, t, (An_linear_interpolation,))
    return pd.DataFrame.from_items([('t', t), ('Photosynthesis', photosynthesis(An[t],LAMINA_AREA)),
                                    ('STORAGE', soln[:, 0]), ('SUCROSE_lamina', soln[:, 1]),
                                    ('TRIOSESP', soln[:, 2])])
    
    


