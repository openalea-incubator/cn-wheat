# -*- coding: latin-1 -*-

'''
Created on 28 juil. 2014

@author: cchambon
'''

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d

import organ


def calculate_all_derivatives(y, t, organs, phloem):
    '''Compute the derivative of y at t.
    '''
    y_iter = iter(y)
    y_derivatives = []
    
    if phloem is not None:
        SUCROSE_phloem = y_iter.next()
    else:
        SUCROSE_phloem = 0
        
    for organ_ in organs:
        
        if isinstance(organ_, organ.Lamina):
            STORAGE = y_iter.next()
            SUCROSE_lamina = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = organ_.calculate_Photosynthesis(t)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE_lamina, SUCROSE_phloem)
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_lamina_derivative = organ_.calculate_SUCROSE_lamina_derivative(S_sucrose, D_storage, organ_.Loading_sucrose)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([STORAGE_derivative, SUCROSE_lamina_derivative, TRIOSESP_derivative])
            
        elif isinstance(organ_, organ.Ear):
            Ear_value = y_iter.next()
            # flows
            organ_.Unloading = organ_.calculate_Unloading(SUCROSE_phloem)
            # compartments derivatives
            Ear_value_derivative = organ_.calculate_Ear_value_derivative(organ_.Unloading)
            y_derivatives.append(Ear_value_derivative)
            
    if phloem is not None:
        SUCROSE_phloem_derivative = phloem.calculate_SUCROSE_phloem_derivative(organs)
        y_derivatives.insert(0, SUCROSE_phloem_derivative)
            
    return y_derivatives
    
    
def run(start_time, stop_time, number_of_output_steps, organs, phloem=None):
    '''
    Compute CN exchanges between laminae, ear and phloem.
    
    The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps. 
    
    Parameters
    ----------
    start_time: int
        The starting of the time grid.
    stop_time: int
        The end of the time grid.
    number_of_output_steps: float
        Number of time points for which to compute the CN exchanges in the system. 
    organs: list
        List of organs. Each organ can be organ.Lamina or organ.Ear. 
    phloem: organ.Phloem or None
        The phloem of the system. 
        None if there is no phloem in the system (default).
     
    Returns
    -------
    out : pandas.DataFrame
        Dataframe containing the CN exchanges between laminae, ear and phloem 
        for each desired time.
        
    References
    ----------
    .. Barillot et al. 2014.
    
    '''
    
    # interpolate the Assimilation and construct the list of initial conditions
    initial_conditions = []
    
    if phloem is not None:
        initial_conditions.extend(phloem.get_initial_conditions())
        
    for organ_ in organs:
        if isinstance(organ_, organ.Lamina):
            organ_.An_linear_interpolation = interp1d(organ_.Assimilation.index, organ_.Assimilation)
        initial_conditions.extend(organ_.get_initial_conditions())
            
    t = np.linspace(start_time, stop_time, number_of_output_steps)
    
    soln = odeint(calculate_all_derivatives, initial_conditions, t, (organs, phloem))
    
    soln_iter = iter(soln.T)
    
    result_items = [('t', t)]
    
    if phloem is not None:
        SUCROSE_phloem = soln_iter.next()
        variables = [('Conc_Sucrose_phloem', map(phloem.calculate_Conc_Sucrose_phloem, SUCROSE_phloem))]
        compartments = [('SUCROSE_phloem', SUCROSE_phloem)]
        result_items.extend(variables + compartments)
    else:
        SUCROSE_phloem = np.zeros_like(t)
        
    for organ_ in organs:
        if isinstance(organ_, organ.Lamina):
            STORAGE = soln_iter.next()
            SUCROSE_lamina = soln_iter.next()
            TRIOSESP = soln_iter.next()
            variables = [('Photosynthesis%s' % organ_.name, organ_.calculate_Photosynthesis(t)),
                         ('Conc_TriosesP%s' % organ_.name, organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         ('Conc_Storage%s' % organ_.name, organ_.calculate_Conc_Storage(STORAGE)),
                         ('Conc_Sucrose_lamina%s' % organ_.name, organ_.calculate_Conc_Sucrose_lamina(SUCROSE_lamina))]
            flows = [('D_storage%s' % organ_.name, map(organ_.calculate_D_storage, STORAGE)),
                     ('S_storage%s' % organ_.name, map(organ_.calculate_S_storage, TRIOSESP)),
                     ('S_sucrose%s' % organ_.name, map(organ_.calculate_S_sucrose, TRIOSESP)),
                     ('Loading_sucrose%s' % organ_.name, map(organ_.calculate_Loading_sucrose, SUCROSE_lamina, SUCROSE_phloem))]
            compartments = [('STORAGE%s' % organ_.name, STORAGE), 
                            ('SUCROSE_lamina%s' % organ_.name, SUCROSE_lamina),
                            ('TRIOSESP%s' % organ_.name, TRIOSESP)]
        elif isinstance(organ_, organ.Ear):
            Ear_value = soln_iter.next()
            variables = [('Dry_mass_ear', organ_.calculate_Dry_mass_ear(Ear_value))]
            flows = [('Unloading', map(organ_.calculate_Unloading, SUCROSE_phloem))]
            compartments = [('Ear', Ear_value)]
        result_items.extend(variables + flows + compartments)
            
    return pd.DataFrame.from_items(result_items)


