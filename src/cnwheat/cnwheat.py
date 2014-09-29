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


def _calculate_all_derivatives(y, t, phloem, organs):
    '''Compute the derivative of `y` at `t`.
    '''
    y_iter = iter(y)
    y_derivatives = []

    SUCROSE_phloem = y_iter.next()

    for organ_ in organs:

        if isinstance(organ_, organ.Lamina):
            STORAGE = y_iter.next()
            SUCROSE = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = organ_.calculate_Photosynthesis(t)
            Rd = organ_.calculate_Rd(Photosynthesis)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE, SUCROSE_phloem)
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_derivative = organ_.calculate_SUCROSE_derivative(S_sucrose, D_storage, organ_.Loading_sucrose, Rd)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([STORAGE_derivative, SUCROSE_derivative, TRIOSESP_derivative])

        elif isinstance(organ_, organ.Chaff):
            STORAGE = y_iter.next()
            SUCROSE = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = organ_.calculate_Photosynthesis(t)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE, SUCROSE_phloem)
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_derivative = organ_.calculate_SUCROSE_derivative(S_sucrose, D_storage, organ_.Loading_sucrose)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([STORAGE_derivative, SUCROSE_derivative, TRIOSESP_derivative])

        elif isinstance(organ_, (organ.Internode, organ.Peduncle, organ.Sheath)):
            FRUCTAN = y_iter.next()
            STORAGE = y_iter.next()
            SUCROSE = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = organ_.calculate_Photosynthesis(t)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE, SUCROSE_phloem)
            Regul_Sfructan = organ_.calculate_Regul_Sfructan(organ_.Loading_sucrose)
            S_fructan = organ_.calculate_S_fructan(SUCROSE, Regul_Sfructan)
            D_fructan = organ_.calculate_D_fructan(SUCROSE, FRUCTAN)

            # compartments derivatives
            FRUCTAN_derivative = organ_.calculate_FRUCTAN_derivative(S_fructan, D_fructan)
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_derivative = organ_.calculate_SUCROSE_derivative(S_sucrose, D_storage, S_fructan, D_fructan, organ_.Loading_sucrose)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([FRUCTAN_derivative, STORAGE_derivative, SUCROSE_derivative, TRIOSESP_derivative])

        elif isinstance(organ_, organ.Grains):
            STORAGE = y_iter.next()
            STRUCTURE = y_iter.next()
            # needed variables
            RGR_structure = organ_.calculate_RGR_structure(SUCROSE_phloem)
            # flows
            organ_.Loading_sucrose_structure = organ_.calculate_Loading_sucrose_structure(t, STRUCTURE, RGR_structure)
            organ_.Loading_sucrose_storage = organ_.calculate_Loading_sucrose_storage(t, SUCROSE_phloem)
            organ_.STRUCTURE = STRUCTURE
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(organ_.Loading_sucrose_storage, organ_.STRUCTURE)
            STRUCTURE_derivative = organ_.calculate_STRUCTURE_derivative(organ_.Loading_sucrose_structure)
            y_derivatives.extend([STORAGE_derivative, STRUCTURE_derivative])

        elif isinstance(organ_, organ.Roots):
            Sucrose = y_iter.next()
            # flows
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE_phloem)
            # compartments derivatives
            Sucrose_derivative = organ_.calculate_Sucrose_derivative(organ_.Loading_sucrose)
            y_derivatives.extend([Sucrose_derivative])

    SUCROSE_phloem_derivative = phloem.calculate_SUCROSE_derivative(organs)
    y_derivatives.insert(0, SUCROSE_phloem_derivative)

    return y_derivatives


def run(start_time, stop_time, number_of_output_steps, phloem, organs):
    '''
    Compute CN exchanges between `organs` and `phloem`.

    The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

    :Parameters:

        - `start_time` (:class:`int`) - The starting of the time grid.

        - `stop_time` (:class:`int`) - The end of the time grid.

        - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.

        - `phloem` (:class:`cnwheat.organ.Phloem`) - The phloem of the system.

        - `organs` (:class:`list`) - List of :class:`cnwheat.organ.Organ`. Does not include the phloem.

    :Returns:
        Dataframe containing the CN exchanges between the *organs* and *phloem*
        for each desired time.

    :Returns Type:
        :class:`pandas.DataFrame`

    .. seealso:: Barillot et al. 2014.

    '''

    # interpolate the Assimilation and construct the list of initial conditions
    initial_conditions = []

    initial_conditions.extend(phloem.get_initial_conditions())

    for organ_ in organs:
        try:
            organ_.An_linear_interpolation = interp1d(organ_.Assimilation.index, organ_.Assimilation)
        except:
            pass
        initial_conditions.extend(organ_.get_initial_conditions())

    t = np.linspace(start_time, stop_time, number_of_output_steps)

    soln = odeint(_calculate_all_derivatives, initial_conditions, t, (phloem, organs))

    soln_iter = iter(soln.T)

    # construct the table of the results
    result_items = [('t', t)]

    SUCROSE_phloem = soln_iter.next()
    variables = [(('Conc_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_Conc_Sucrose(SUCROSE_phloem)),
                 (('Conc_C_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_Conc_C_Sucrose(SUCROSE_phloem))]
    compartments = [(('SUCROSE_%s' % phloem.name).rstrip('_'), SUCROSE_phloem)]
    result_items.extend(variables + compartments)

    for organ_ in organs:
        if isinstance(organ_, organ.Lamina):
            STORAGE = soln_iter.next()
            SUCROSE = soln_iter.next()
            TRIOSESP = soln_iter.next()
            Photosynthesis = map(organ_.calculate_Photosynthesis, t)
            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), Photosynthesis),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Storage(STORAGE)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Sucrose(SUCROSE)),
                         (('Rd_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Rd, Photosynthesis))]
            flows = [(('D_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_storage, STORAGE)),
                     (('S_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_storage, TRIOSESP)),
                     (('S_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_sucrose, TRIOSESP)),
                     (('Loading_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose, SUCROSE, SUCROSE_phloem))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('SUCROSE_%s' % organ_.name).rstrip('_'), SUCROSE),
                            (('TRIOSESP_%s' % organ_.name).rstrip('_'), TRIOSESP)]
        elif isinstance(organ_, organ.Chaff):
            STORAGE = soln_iter.next()
            SUCROSE = soln_iter.next()
            TRIOSESP = soln_iter.next()
            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), organ_.calculate_Photosynthesis(t)),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Storage(STORAGE)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Sucrose(SUCROSE))]
            flows = [(('D_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_storage, STORAGE)),
                     (('S_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_storage, TRIOSESP)),
                     (('S_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_sucrose, TRIOSESP)),
                     (('Loading_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose, SUCROSE, SUCROSE_phloem))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('SUCROSE_%s' % organ_.name).rstrip('_'), SUCROSE),
                            (('TRIOSESP_%s' % organ_.name).rstrip('_'), TRIOSESP)]
        elif isinstance(organ_, (organ.Internode, organ.Peduncle, organ.Sheath)):
            FRUCTAN = soln_iter.next()
            STORAGE = soln_iter.next()
            SUCROSE = soln_iter.next()
            TRIOSESP = soln_iter.next()
            Loading_sucrose = map(organ_.calculate_Loading_sucrose, SUCROSE, SUCROSE_phloem)
            Regul_Sfructan = map(organ_.calculate_Regul_Sfructan, Loading_sucrose)
            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), organ_.calculate_Photosynthesis(t)),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Storage(STORAGE)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Sucrose(SUCROSE)),
                         (('Conc_Fructan_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Fructan(FRUCTAN)), 
                         (('Regul_Sfructan_%s' % organ_.name).rstrip('_'), Regul_Sfructan)]
            flows = [(('D_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_storage, STORAGE)),
                     (('S_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_storage, TRIOSESP)),
                     (('S_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_sucrose, TRIOSESP)),
                     (('Loading_sucrose_%s' % organ_.name).rstrip('_'), Loading_sucrose),
                     (('D_fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_fructan, SUCROSE_phloem, FRUCTAN)),
                     (('S_fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_fructan, SUCROSE_phloem, Regul_Sfructan))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('SUCROSE_%s' % organ_.name).rstrip('_'), SUCROSE),
                            (('TRIOSESP_%s' % organ_.name).rstrip('_'), TRIOSESP),
                            (('FRUCTAN_%s' % organ_.name).rstrip('_'), FRUCTAN)]
        elif isinstance(organ_, organ.Grains):
            STORAGE = soln_iter.next()
            STRUCTURE = soln_iter.next()
            RGR_structure = map(organ_.calculate_RGR_structure, SUCROSE_phloem)
            variables = [(('Dry_mass_%s' % organ_.name).rstrip('_'), organ_.calculate_Dry_mass(STRUCTURE, STORAGE)),
                         (('RGR_structure_%s' % organ_.name).rstrip('_'), RGR_structure)]
            flows = [(('Loading_sucrose_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose_storage, t, SUCROSE_phloem)),
                     (('Loading_sucrose_structure_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose_structure, t, STRUCTURE, RGR_structure))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('STRUCTURE_%s' % organ_.name).rstrip('_'), STRUCTURE)]
        elif isinstance(organ_, organ.Roots):
            Sucrose = soln_iter.next()
            variables = [(('Dry_mass_%s' % organ_.name).rstrip('_'), organ_.calculate_Dry_mass(Sucrose))]
            flows = [(('Loading_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose, SUCROSE_phloem))]
            compartments = [(('Sucrose_%s' % organ_.name).rstrip('_'), Sucrose)]
        result_items.extend(variables + flows + compartments)

    return pd.DataFrame.from_items(result_items)

