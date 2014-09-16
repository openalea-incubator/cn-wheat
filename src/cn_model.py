# -*- coding: latin-1 -*-

'''
Created on 28 juil. 2014

@author: cchambon
'''

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d


def func(y, t0, lamina, An_linear_interpolation):
    '''Compute the derivative of y at t0.
    '''
    STORAGE = y[0]
    SUCROSE_lamina = y[1]
    TRIOSESP = y[2]
    
    # variable
    Photosynthesis = lamina.calculate_Photosynthesis(An_linear_interpolation(t0))
    
    # flows
    D_storage = lamina.calculate_D_storage(STORAGE)
    S_storage = lamina.calculate_S_storage(TRIOSESP)
    S_sucrose = lamina.calculate_S_sucrose(TRIOSESP)
    
    # compartments
    STORAGE_derivative = lamina.calculate_STORAGE_derivative(S_storage, D_storage)
    SUCROSE_lamina_derivative = lamina.calculate_SUCROSE_lamina_derivative(S_sucrose, D_storage)
    TRIOSESP_derivative = lamina.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
    
    return [STORAGE_derivative, SUCROSE_lamina_derivative, TRIOSESP_derivative]
    
    
def run(start_time, stop_time, number_of_output_steps, lamina_1):
    '''
    Compute STORAGE, SUCROSE_lamina and TRIOSESP between `start_time` and `stop_time`, for y(t=0) = `initial_conditions`. 
    
    The Photosynthesis is computed from the Assimilation `An` and LAMINA_AREA.
    
    Parameters
    ----------
    start_time: int
        The starting of the time grid.
    stop_time: int
        The end of the time grid.
    number_of_output_steps: float
        Number of time points for which to solve for y. 
    initial_conditions: list
        List of initial conditions: STORAGE(t=0), SUCROSE_lamina(t=0), TRIOSESP(t=0). 
    An: pandas.Series
        The Assimilation of the lamina at each hour t (µmol m-2 s-1). Due to odeint 
        integration solver, `An` must provide Assimilation value at t=stop_time+1.   
     
    Returns
    -------
    out : pandas.DataFrame
        Dataframe containing the value of Photosynthesis, STORAGE, SUCROSE_lamina and TRIOSESP 
        for each desired time.
        
    References
    ----------
    .. [HNW93] Barillot et al. 2014.
    
    '''
    An_linear_interpolation = interp1d(lamina_1.Assimilation.index, lamina_1.Assimilation)
    t = np.linspace(start_time, stop_time, number_of_output_steps)
    soln = odeint(func, lamina_1.get_initial_conditions(), t, (lamina_1, An_linear_interpolation))
    STORAGE = soln[:, 0]
    SUCROSE_lamina = soln[:, 1]
    TRIOSESP = soln[:, 2]
    return pd.DataFrame.from_items([('t', t), 
                                    # variables
                                    ('Photosynthesis', lamina_1.calculate_Photosynthesis(lamina_1.Assimilation[t])),
                                    # flows
                                    ('D_storage', map(lamina_1.calculate_D_storage, STORAGE)),
                                    ('S_storage', map(lamina_1.calculate_S_storage, TRIOSESP)),
                                    ('S_sucrose', map(lamina_1.calculate_S_sucrose, TRIOSESP)),
                                    # compartments
                                    ('STORAGE', STORAGE), 
                                    ('SUCROSE_lamina', SUCROSE_lamina),
                                    ('TRIOSESP', TRIOSESP)])
    
    


