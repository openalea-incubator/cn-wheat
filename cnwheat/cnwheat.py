# -*- coding: latin-1 -*-
"""
    cnwheat.cnwheat
    ~~~~~~~~~~~~~~~

    Front-end to run the model.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

from __future__ import division # use "//" to do integer division

import sys

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d

import organ
import photosynthesis

class CNWheatError(Exception): pass
class CNWheatInputError(CNWheatError): pass
class CNWheatRunError(CNWheatError): pass

class CNWheat(object):
    """
    Compute CN exchanges in wheat architecture defined by lamina, sheaths, internodes, peduncles, a chaff, a phloem, roots and grains.
    This class permits to initialize and run the model.
    
    :Parameters:
    
        - `organs` (:class:`list`) - List of :class:`cnwheat.organ.Organ`.
        
        - `meteo` (:class:`pandas.DataFrame`) - a :class:`pandas.DataFrame` which index 
          is time in hours and which columns are "Tac", "hs" and "Ca". Due to the 
          solver of :func:`scipy.integrate.odeint`, `meteo` must provide data for t = `stop_time` + 1.
      
    """
    
    def __init__(self, organs, meteo):
        
        self.organs = organs #: the organs used in the model
        self.meteo = meteo #: the meteo data used in the model
        
        # interpolate meteo data
        
        #: A dictionary which contains, for each data of meteo data, a function which 
        #: permits to find the value of new points using interpolation. Keys 
        #: are the name of meteo data (:class:`str`). Values are :class:`scipy.interpolate.interp1d`.
        self.meteo_interpolations = {} 

        for column in meteo.columns:
            self.meteo_interpolations[column] = interp1d(meteo.index, meteo[column])
            
        # interpolate the PAR and construct the list of initial conditions
        self.initial_conditions = [] #: the initial conditions of the compartments in the organs
        self.organs_without_phloem = [] #: the list of organs without the phloem
    
        for organ_ in organs:
            if isinstance(organ_, organ.PhotosyntheticOrgan):
                organ_.PAR_linear_interpolation = interp1d(organ_.PAR.index, organ_.PAR)
            if isinstance(organ_, organ.Phloem):
                self.phloem = organ_ # the phloem
                self.initial_conditions = organ_.get_initial_conditions() + self.initial_conditions
            else:
                self.organs_without_phloem.append(organ_)
                self.initial_conditions.extend(organ_.get_initial_conditions())
        
        self.progressbar = ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver 
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar 
                
                
    def run(self, start_time, stop_time, number_of_output_steps, photosynthesis_computation_interval=0, odeint_mxstep=5000, show_progressbar=False):
        """
        Compute CN exchanges between organs :attr:`organs`, using interpolated meteo data :attr:`meteo_interpolations`
        The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.
    
        :Parameters:
    
            - `start_time` (:class:`int`) - The starting of the time grid.
    
            - `stop_time` (:class:`int`) - The end of the time grid.
    
            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.
    
            - `photosynthesis_computation_interval` (:class:`int`) - The interval which defined the time steps
              at which photosynthesis is computed. For example, if `photosynthesis_computation_interval` = 4,
              then photosynthesis is computed at t=0, t=4, t=8, ...
              This permits to save computation time for large simulation.
              If `photosynthesis_computation_interval`=0 (the default), then photosynthesis 
              is computed for each time step demanded by the solver.
    
            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0, then `mxstep` is determined by the solver.
              The default value ( `5000` ) normally permits to solve the current model. User should increased this value if a more complex model is defined
              and if this model make the integration failed.
              
            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: DO NOT show the progress bar.
    
        :Returns:
            Dataframe containing the CN exchanges between organs
            for each desired time.
    
        :Returns Type:
            :class:`pandas.DataFrame`
    
        .. warning:: due to the solver of :func:`scipy.integrate.odeint`, `meteo` must provide data for t = `stop_time` + 1.
                     For the same reason, the attribute `PAR` of each organ of `organs` must also provide data for t = `stop_time` + 1.
                     This is automatically checked by the current function.
    
        .. seealso:: Barillot et al. 2014.
    
        """
        
        self._check_inputs_consistency(start_time, stop_time)
        
        t = np.linspace(start_time, stop_time, number_of_output_steps)
        
        self.show_progressbar = show_progressbar
        if self.show_progressbar:
            self.progressbar.set_t_max(stop_time)
        
        self._clear_organs()
        
        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, (photosynthesis_computation_interval,), full_output=True, mxstep=odeint_mxstep)
    
        if not set(infodict['mused']).issubset([1,2]): # I'm not sure if this test is robust or not... Beware especially when scipy is updated.
            raise CNWheatRunError("Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'.")
    
        cnwheat_output_df = self._format_solver_output(t, soln)
        
        return cnwheat_output_df
        
    
    def _check_inputs_consistency(self, start_time, stop_time):
        """
        Check the consistency of meteo and photosynthesis data. 
        Raise an exception if meteo or photosynthesis data is not valid. 
        """
        
        # check the consistency of meteo
        lowest_t = self.meteo.first_valid_index()
        if start_time < lowest_t:
            raise CNWheatInputError('the lowest t ({}) in meteo data is greater than start_time ({}).'.format(lowest_t, start_time))
    
        solver_upper_boundary = stop_time + 1
        highest_t = self.meteo.last_valid_index()
        if highest_t < solver_upper_boundary:
            raise CNWheatInputError("""the highest t ({}) in meteo data is lower than stop_time + 1 = {}.
                            scipy.integrate.odeint requires the highest t to be equal or
                            greater than stop_time + 1""".format(highest_t, solver_upper_boundary))
    
        # check the consistency of the PAR
        for organ_ in self.organs_without_phloem:
            if isinstance(organ_, organ.PhotosyntheticOrgan):
                lowest_t = organ_.PAR.first_valid_index()
                if start_time < lowest_t:
                    raise CNWheatInputError('the lowest t ({}) in the PAR of {} is greater than start_time ({}).'.format(lowest_t, organ_.name, start_time))
                highest_t = organ_.PAR.last_valid_index()
                if highest_t < solver_upper_boundary:
                    raise CNWheatInputError("""the highest t ({}) in the PAR of {} is lower than stop_time + 1 = {}.
                                    scipy.integrate.odeint requires the highest t to be equal or
                                    greater than stop_time + 1""".format(highest_t, organ_.name, solver_upper_boundary))


    def _clear_organs(self):
        for organ_ in self.organs_without_phloem:
            if isinstance(organ_, organ.PhotosyntheticOrgan):
                organ_.photosynthesis_mapping.clear()
    

    def _calculate_all_derivatives(self, y, t, photosynthesis_computation_interval):
        """Compute the derivative of `y` at `t`.
        
        :meth:`_calculate_all_derivatives` is passed as **func** argument to 
        :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
        :meth:`_calculate_all_derivatives` is called automatically by 
        :func:`scipy.integrate.odeint <scipy.integrate.odeint>`.
        
        First call to :meth:`_calculate_all_derivatives` uses `y` = **y0** and 
        `t` = **t** [0], where **y0** and **t** are arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
        
        Following calls to :meth:`_calculate_all_derivatives` use `t` in [min( **t** ), max( **t** ) + 1] where 
        **t** is an argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`. `y` is 
        computed automatically by the solver.
        
        :Parameters:
    
            - `y` (:class:`list`) - The current y values. The y values of :attr:`phloem` 
              must appear first. `y` is automatically set by 
              :func:`scipy.integrate.odeint`. User does not have control over `y`. 
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.odeint`, `y` = **y0** 
              where **y0** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `y` is 
              computed automatically by the solver.
    
            - `t` (:class:`float`) - The current t at which we want to compute the derivatives. 
              `t` is automatically set by :func:`scipy.integrate.odeint`. 
              User does not have control over `t`.
              At first call to :meth:`_calculate_all_derivatives` :func:`scipy.integrate.odeint`, 
              `t` = **t** [0], where **t** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs 
              to the interval [min( **t** ), max( **t** ) + 1], where **t** is an 
              argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
    
            - `photosynthesis_computation_interval` (:class:`int`) - The interval which defined the time steps
              at which photosynthesis is computed. For example, if `photosynthesis_computation_interval` = 4,
              then photosynthesis is computed at t=0, t=4, t=8, ...
              This permits to save computation time for large simulation.
              If `photosynthesis_computation_interval` = 0 (the default), then photosynthesis is computed 
              for each time step demanded by the solver.
    
        :Returns:
            The derivatives of `y` at `t`. 
    
        :Returns Type:
            :class:`list`
    
    
        """
        y_iter = iter(y)
        y_derivatives = []
    
        sucrose_phloem = y_iter.next()
    
        for organ_ in self.organs_without_phloem:
    
            if isinstance(organ_, organ.PhotosyntheticOrgan):
                storage = y_iter.next()
                sucrose = y_iter.next()
                triosesP = y_iter.next()
                fructan = y_iter.next()
                # calculate the time step at which we want to compute the photosynthesis
                if photosynthesis_computation_interval == 0: # i.e. we want to compute the photosynthesis for each time step demanded by the solver
                    t_inf = t
                else:
                    t_inf  = t // photosynthesis_computation_interval * photosynthesis_computation_interval
                # calculate the photosynthesis of organ_ only if it has not been already calculated at t_inf
                if t_inf not in organ_.photosynthesis_mapping:
                    PAR_t_inf = organ_.PAR_linear_interpolation(t_inf)
                    Tac_t_inf = self.meteo_interpolations['Tac'](t_inf)
                    hs_t_inf = self.meteo_interpolations['hs'](t_inf)
                    Ca_t_inf = self.meteo_interpolations['Ca'](t_inf)
                    An = photosynthesis.PhotosynthesisModel.calculate_An(t_inf, PAR_t_inf, Tac_t_inf, Ca_t_inf, hs_t_inf)
                    organ_.photosynthesis_mapping[t_inf] = organ_.calculate_photosynthesis(t_inf, An)
                photosynthesis_ = organ_.photosynthesis_mapping[t_inf]
                # flows
                d_storage = organ_.calculate_d_storage(storage)
                s_storage = organ_.calculate_s_storage(triosesP)
                s_sucrose = organ_.calculate_s_sucrose(triosesP)
                organ_.loading_sucrose = organ_.calculate_loading_sucrose(sucrose, sucrose_phloem)
                regul_s_fructan = organ_.calculate_regul_s_fructan(organ_.loading_sucrose)
                d_fructan = organ_.calculate_d_fructan(sucrose, fructan)
                s_fructan = organ_.calculate_s_fructan(sucrose, regul_s_fructan)
                # compartments derivatives
                storage_derivative = organ_.calculate_storage_derivative(s_storage, d_storage)
                sucrose_derivative = organ_.calculate_sucrose_derivative(s_sucrose, d_storage, organ_.loading_sucrose, s_fructan, d_fructan)
                triosesP_derivative = organ_.calculate_triosesP_derivative(photosynthesis_, s_sucrose, s_storage)
                fructan_derivative = organ_.calculate_fructan_derivative(s_fructan, d_fructan)
                y_derivatives.extend([storage_derivative, sucrose_derivative, triosesP_derivative, fructan_derivative])
    
            elif isinstance(organ_, organ.Grains):
                storage = y_iter.next()
                structure = y_iter.next()
                # needed variables
                RGR_structure = organ_.calculate_RGR_structure(sucrose_phloem)
                # flows
                organ_.unloading_sucrose_structure = organ_.calculate_unloading_sucrose_structure(t, structure, RGR_structure)
                organ_.unloading_sucrose_storage = organ_.calculate_unloading_sucrose_storage(t, sucrose_phloem)
                organ_.structure = structure
                # compartments derivatives
                storage_derivative = organ_.calculate_storage_derivative(organ_.unloading_sucrose_storage, organ_.structure)
                structure_derivative = organ_.calculate_structure_derivative(organ_.unloading_sucrose_structure)
                y_derivatives.extend([storage_derivative, structure_derivative])
    
            elif isinstance(organ_, organ.Roots):
                sucrose = y_iter.next()
                # flows
                organ_.unloading_sucrose = organ_.calculate_unloading_sucrose(sucrose_phloem)
                # compartments derivatives
                sucrose_derivative = organ_.calculate_sucrose_derivative(organ_.unloading_sucrose)
                y_derivatives.extend([sucrose_derivative])
    
        sucrose_phloem_derivative = self.phloem.calculate_sucrose_derivative(self.organs_without_phloem)
        y_derivatives.insert(0, sucrose_phloem_derivative)
        
        if self.show_progressbar:
            self.progressbar.update(t)
    
        return y_derivatives
    
    
    def _format_solver_output(self, t, solver_output):
        """
        Create a :class:`pandas.DataFrame` wich columns are: 
            
            * the time grid `t`,
            * the output of the solver `solver_output`,
            * and intermediate and post-processed variables usefull for debug and validation 
        """
        
        solver_output_iter = iter(solver_output.T)
    
        result_items = [('t', t)]
    
        sucrose_phloem = solver_output_iter.next()
        variables = [('Conc_Sucrose_{}'.format(self.phloem.name).rstrip('_'), self.phloem.calculate_conc_sucrose(sucrose_phloem)),
                     ('Conc_C_Sucrose_{}'.format(self.phloem.name).rstrip('_'), self.phloem.calculate_conc_c_sucrose(sucrose_phloem))]
        compartments = [('Sucrose_{}'.format(self.phloem.name).rstrip('_'), sucrose_phloem)]
        result_items.extend(variables + compartments)
    
        for organ_ in self.organs_without_phloem:
            if isinstance(organ_, organ.PhotosyntheticOrgan):
                storage = solver_output_iter.next()
                sucrose = solver_output_iter.next()
                triosesP = solver_output_iter.next()
                fructan = solver_output_iter.next()
                loading_sucrose = map(organ_.calculate_loading_sucrose, sucrose, sucrose_phloem)
                regul_s_fructan = map(organ_.calculate_regul_s_fructan, loading_sucrose)
    
                An = np.array(map(photosynthesis.PhotosynthesisModel.calculate_An,
                                  t,
                                  organ_.PAR_linear_interpolation(t),
                                  self.meteo_interpolations['Tac'](t),
                                  self.meteo_interpolations['Ca'](t),
                                  self.meteo_interpolations['hs'](t)))
    
                variables = [('Photosynthesis_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_photosynthesis, t, An)),
                             ('Photosynthesis_Surfacic_Rate_{}'.format(organ_.name).rstrip('_'), An),
                             ('Conc_TriosesP_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_triosesP(triosesP)),
                             ('Conc_Storage_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_storage(storage)),
                             ('Conc_Sucrose_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_sucrose(sucrose)),
                             ('Conc_Fructan_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_fructan(fructan)),
                             ('Regul_S_Fructan_{}'.format(organ_.name).rstrip('_'), regul_s_fructan)]
                flows = [('D_Storage_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_d_storage, storage)),
                         ('S_Storage_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_storage, triosesP)),
                         ('S_Sucrose_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_sucrose, triosesP)),
                         ('Loading_Sucrose_{}'.format(organ_.name).rstrip('_'), loading_sucrose),
                         ('D_Fructan_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_d_fructan, sucrose, fructan)),
                         ('S_Fructan_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_fructan, sucrose, regul_s_fructan))]
                compartments = [('Storage_{}'.format(organ_.name).rstrip('_'), storage),
                                ('Sucrose_{}'.format(organ_.name).rstrip('_'), sucrose),
                                ('TriosesP_{}'.format(organ_.name).rstrip('_'), triosesP),
                                ('Fructan_{}'.format(organ_.name).rstrip('_'), fructan)]
            elif isinstance(organ_, organ.Grains):
                storage = solver_output_iter.next()
                structure = solver_output_iter.next()
                RGR_structure = map(organ_.calculate_RGR_structure, sucrose_phloem)
                variables = [('Dry_Mass_{}'.format(organ_.name).rstrip('_'), organ_.calculate_dry_mass(structure, storage)),
                             ('RGR_Structure_{}'.format(organ_.name).rstrip('_'), RGR_structure)]
                flows = [('Unloading_Sucrose_Storage_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose_storage, t, sucrose_phloem)),
                         ('Unloading_Sucrose_Structure_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose_structure, t, structure, RGR_structure))]
                compartments = [('Storage_{}'.format(organ_.name).rstrip('_'), storage),
                                ('Structure_{}'.format(organ_.name).rstrip('_'), structure)]
            elif isinstance(organ_, organ.Roots):
                sucrose = solver_output_iter.next()
                variables = [('Dry_Mass_{}'.format(organ_.name).rstrip('_'), organ_.calculate_dry_mass(sucrose))]
                flows = [('Unloading_Sucrose_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose, sucrose_phloem))]
                compartments = [('Sucrose_{}'.format(organ_.name).rstrip('_'), sucrose)]
            result_items.extend(variables + flows + compartments)
    
        return pd.DataFrame.from_items(result_items)
    
    
class ProgressBarError(Exception): pass
    
class ProgressBar(object):
    """
    Display a console progress bar.
    """
    
    def __init__(self, bar_length=20, title='', block_character='#', uncomplete_character='-'):
        if bar_length <= 0:
            raise ProgressBarError('bar_length <= 0')
        self.bar_length = bar_length #: the number of blocks in the progress bar. MUST BE GREATER THAN ZERO !
        self.t_max = 1 #: the maximum t that the progress bar can display. MUST BE GREATER THAN ZERO !
        self.block_interval = 1 #: the time interval of each block. MUST BE GREATER THAN ZERO !
        self.last_upper_t = 0 #: the last upper t displayed by the progress bar 
        self.progress_mapping = {} #: a mapping to optimize the refresh rate
        self.title = title #: the title to write on the left side of the progress bar
        self.block_character = block_character #: the character to represent a block
        self.uncomplete_character = uncomplete_character #: the character to represent the uncompleted part of the progress bar 
        
    def set_t_max(self, t_max):
        """"Set :attr:`t_max` and update other attributes accordingly.
        """
        if t_max <= 0:
            raise ProgressBarError('t_max <= 0')
        self.t_max = t_max
        self.block_interval = self.t_max / self.bar_length
        self.last_upper_t = 0
        self.progress_mapping.clear()
        
    def update(self, t):
        """Update the progress bar if needed. 
        """
        t = min(t, self.t_max)
        if t < self.last_upper_t:
            return
        else:
            self.last_upper_t = t
        t_inf = t // self.block_interval * self.block_interval
        if t_inf not in self.progress_mapping:
            progress = t / self.t_max
            block = int(round(self.bar_length * progress))
            text = "\r{0}: [{1}] {2:>5d}%".format(self.title, 
                                                  self.block_character * block + self.uncomplete_character * (self.bar_length - block), 
                                                  int(progress*100))
            self.progress_mapping[t_inf] = text
            sys.stdout.write(self.progress_mapping[t_inf])
            sys.stdout.flush()
            