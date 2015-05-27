# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to couple CN-Wheat and Farquhar-Wheat.

    You must first install :mod:`cnwheat` and :mod:`farquharwheat` (and add them to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''

import os
import time, datetime
import profile, pstats

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cnwheat import simulation
from cnwheat import model as cnwheat_model
from farquharwheat import model as photosynthesis_model
from cnwheat import tools
from senescwheat.model import SenescenceModel

t0 = time.time()

INPUTS_DIRPATH = 'inputs'

PAR_FILENAME = 'PAR.csv'

GREEN_AREA_FILENAME = 'green_area.csv'

METEO_FILENAME = 'meteo.csv'

OUTPUTS_DIRPATH = 'outputs'

GRAPHS_DIRPATH = 'graphs' # GRAPHS_DIRPATH must be an existing directory

PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
PHYTOMERS_OUTPUTS_FILENAME = 'phytomers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

OUTPUTS_PRECISION = 6

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # set to logging.DEBUG for debugging or logging.INFO for INFO level

tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL)


def compute_CN_distrib(run_simu=True, make_graphs=True):

    plants_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PLANTS_OUTPUTS_FILENAME)
    axes_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME)
    phytomers_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PHYTOMERS_OUTPUTS_FILENAME)
    organs_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ORGANS_OUTPUTS_FILENAME)
    elements_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)

    if run_simu:

        population = cnwheat_model.Population()

        plant = cnwheat_model.Plant(index=1)
        population.plants.append(plant)

        axis = cnwheat_model.Axis(axis_type=cnwheat_model.Axis.Types.MAIN_STEM, index=0)
        plant.axes.append(axis)

        axis.grains = cnwheat_model.Grains(starch=0, structure=4000, proteins=170)

        axis.roots = cnwheat_model.Roots(mstruct=0.504, Nstruct=0.01, sucrose=900, nitrates=250, amino_acids=60)
        axis.soil = cnwheat_model.Soil(volume=0.00875, nitrates=1875) # 1875 �mol of nitrates in this volume is almost equivalent to 3 g m-2 of N fertilizer

        axis.phloem = cnwheat_model.Phloem(sucrose=0, amino_acids=0)

        # Phytomer 1
        phytomer1 = cnwheat_model.Phytomer(index=1)

        phytomer1.lamina = cnwheat_model.Lamina()
        lamina_element = cnwheat_model.LaminaElement(area=0.00346, green_area=0.00346, mstruct=0.14, Nstruct=0.00102, width= 0.018, height=0.6,
                                        starch=0, sucrose=252, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=16, proteins=380)
        phytomer1.lamina.exposed_element = lamina_element

        phytomer1.sheath = cnwheat_model.Sheath()
        sheath_element = cnwheat_model.SheathElement(area=0.0006, green_area=0.0006, mstruct=0.103, Nstruct=0.00068, width=0.0011, height=0.5,
                                        starch=0, sucrose=185, triosesP=0, fructan=0,
                                        nitrates=0 , amino_acids=12, proteins=130)
        phytomer1.sheath.exposed_element = sheath_element

        # Internode enclosed
        phytomer1.internode = cnwheat_model.Internode()
        internode_enclosed_element = cnwheat_model.InternodeElement(area=0.001129, green_area=0.001129, mstruct=0.1415, Nstruct=0.00064, width=0.00257, height=0.3,
                                              starch=0, sucrose=255, triosesP=0, fructan=0,
                                              nitrates=0, amino_acids=17, proteins=66)
        phytomer1.internode.enclosed_element = internode_enclosed_element

        # Internode exposed
        internode_exposed_element = cnwheat_model.InternodeElement(area=0.000371, green_area=0.000371, mstruct=0.0465, Nstruct=0.00021, width=0.00257, height=0.4,
                                              starch=0, sucrose=84, triosesP=0, fructan=0,
                                              nitrates=0, amino_acids=5, proteins=22)
        phytomer1.internode.exposed_element = internode_exposed_element

        axis.phytomers.append(phytomer1)

        # Phytomer 2
        phytomer2 = cnwheat_model.Phytomer(index=2)

        phytomer2.lamina = cnwheat_model.Lamina()
        lamina_element = cnwheat_model.LaminaElement(area=0.0034, green_area=0.0034, mstruct=0.09, Nstruct=0.00083, width= 0.014, height=0.38,
                                        starch=0, sucrose=162, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=10, proteins=210)
        phytomer2.lamina.exposed_element = lamina_element

        phytomer2.sheath = cnwheat_model.Sheath()
        sheath_element = cnwheat_model.SheathElement(area=0.0005, green_area=0.0005, mstruct=0.069, Nstruct=0.00021, width=0.00091, height=0.3,
                                        starch=0, sucrose=124, triosesP=0, fructan=0,
                                        nitrates=0 , amino_acids=8, proteins=47)
        phytomer2.sheath.exposed_element = sheath_element

        phytomer2.internode = cnwheat_model.Internode()
        internode_element = cnwheat_model.InternodeElement(area=0.0004, green_area=0.0004, mstruct=0.18, Nstruct=0.00033, width=0.00099, height=0.18,
                                              starch=0, sucrose=324, triosesP=0, fructan=0,
                                              nitrates=0, amino_acids=21, proteins=20)
        phytomer2.internode.enclosed_element = internode_element

        axis.phytomers.append(phytomer2)

        # Phytomer 3
        phytomer3 = cnwheat_model.Phytomer(index=3)

        phytomer3.lamina = cnwheat_model.Lamina()
        lamina_element = cnwheat_model.LaminaElement(area=0.00228, green_area=0.00228, mstruct=0.05, Nstruct=0.00053, width= 0.0125, height=0.24,
                                        starch=0, sucrose=90, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=6, proteins=85)
        phytomer3.lamina.exposed_element = lamina_element

        phytomer3.sheath = cnwheat_model.Sheath()
        sheath_element = cnwheat_model.SheathElement(area=0.0004, green_area=0.0004, mstruct=0.043, Nstruct=0.00011, width=0.00051, height=0.18,
                                        starch=0, sucrose=77, triosesP=0, fructan=0,
                                        nitrates=0 , amino_acids=5, proteins=13)
        phytomer3.sheath.exposed_element = sheath_element

        phytomer3.internode = cnwheat_model.Internode()
        internode_element = cnwheat_model.InternodeElement(area=0.00025, green_area=0.00025, mstruct=0.154, Nstruct=0.00014, width=0.00093, height=0.08,
                                              starch=0, sucrose=277, triosesP=0, fructan=0,
                                              nitrates=0, amino_acids=18, proteins=20)
        phytomer3.internode.enclosed_element = internode_element

        axis.phytomers.append(phytomer3)

        # Phytomer 4 (reproductive)
        phytomer4 = cnwheat_model.Phytomer(index=4)

        # Enclosed peduncle
        phytomer4.peduncle = cnwheat_model.Peduncle()
        peduncle_enclosed_element = cnwheat_model.PeduncleElement(area=0.00159, green_area=0.00159, mstruct=0.170, Nstruct=0.00086, width= 0.00349, height=0.65,
                                            starch=0, sucrose=306, triosesP=0, fructan=0, nitrates=0,
                                            amino_acids=20, proteins=120)
        phytomer4.peduncle.enclosed_element = peduncle_enclosed_element

        # Exposed peduncle
        peduncle_exposed_element = cnwheat_model.PeduncleElement(area=0.00081, green_area=0.00081, mstruct=0.087, Nstruct=0.00044, width= 0.00349, height=0.5,
                                            starch=0, sucrose=156, triosesP=0, fructan=0, nitrates=0,
                                            amino_acids=10, proteins=61)
        phytomer4.peduncle.exposed_element = peduncle_exposed_element
        axis.phytomers.append(phytomer4)

        # Phytomer 5 (reproductive)
        phytomer5 = cnwheat_model.Phytomer(index=5)
        phytomer5.chaff = cnwheat_model.Chaff()
        chaff_element = cnwheat_model.ChaffElement(area=0.00075, green_area=0.00075, mstruct=0.21, Nstruct=0.00107, width=0.00265, height= 0.7, starch=0,
                                      sucrose=378, triosesP=0, fructan=0, nitrates=0, amino_acids=25,
                                      proteins=260)
        phytomer5.chaff.exposed_element = chaff_element
        axis.phytomers.append(phytomer5)

        # Get PAR data
        PAR_filepath = os.path.join(INPUTS_DIRPATH, PAR_FILENAME)
        PAR_df = pd.read_csv(PAR_filepath)
        PAR_grouped = PAR_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

        # Get green area data
        green_area_filepath = os.path.join(INPUTS_DIRPATH, GREEN_AREA_FILENAME)
        green_area_df = pd.read_csv(green_area_filepath)
        green_area_grouped = green_area_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

        # get meteo data
        meteo_df = tools.read_t_data(INPUTS_DIRPATH, METEO_FILENAME)

        # create the simulation
        simulation_ = simulation.Simulation()
        
        # initialize the simulation from population
        simulation_.initialize(population)

        # run the models
        start_time = 0
        stop_time = 960 # 960
        photosynthesis_model_ts = 2
        cn_model_ts = 1 #241

        all_plants_df_list = []
        all_axes_df_list = []
        all_phytomers_df_list = []
        all_organs_df_list = []
        all_elements_df_list = []

        for t_photosynthesis_model in xrange(start_time, stop_time, photosynthesis_model_ts):
            # run the model of photosynthesis and update the population
            for plant in population.plants:
                plant_index = plant.index
                for axis in plant.axes:
                    axis_id = axis.id

                    # Root growth and senescence
                    mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth = SenescenceModel.calculate_roots_mstruct_growth(axis.roots.sucrose, axis.roots.mstruct, 3600) #3600 is temporary, will be moved
                    axis.roots.mstruct_C_growth = mstruct_C_growth
                    axis.roots.Nstruct_N_growth = Nstruct_N_growth
                    mstruct_senescence, Nstruct_senescence = SenescenceModel.calculate_roots_senescence(axis.roots.mstruct, axis.roots.Nstruct, 3600)
                    axis.roots.mstruct_senescence = mstruct_senescence
                    delta_mstruct, delta_Nstruct = SenescenceModel.calculate_delta_mstruct_roots(mstruct_growth, Nstruct_growth, mstruct_senescence, Nstruct_senescence)
                    axis.roots.mstruct += delta_mstruct
                    axis.roots.Nstruct += delta_Nstruct

                    for phytomer in axis.phytomers:
                        phytomer_index = phytomer.index
                        for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                            if organ is None:
                                continue
                            organ_type = organ.__class__.__name__
                            for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                                if element is None:
                                    continue

                                # Senescence
                                group_id = (t_photosynthesis_model, plant_index, axis_id, phytomer_index, organ_type, element_type)
                                new_green_area, relative_delta_green_area = SenescenceModel.calculate_relative_delta_green_area(green_area_grouped, group_id, element.green_area)
                                element.green_area = new_green_area
                                element.relative_delta_green_area = relative_delta_green_area
                                new_mstruct, new_Nstruct = SenescenceModel.calculate_delta_mstruct_shoot(relative_delta_green_area, element.mstruct, element.Nstruct)
                                element.mstruct = new_mstruct
                                element.Nstruct = new_Nstruct
                                new_SLN = SenescenceModel.calculate_surfacic_nitrogen(element.nitrates, element.amino_acids, element.proteins, element.Nstruct, new_green_area)
                                element.surfacic_nitrogen = new_SLN

                                # PAR and photosynthesis
                                PAR = PAR_grouped.get_group((t_photosynthesis_model, plant_index, axis_id, phytomer_index, organ_type, element_type)).PAR.values[0]
                                Ag, An, Rd, Tr, Ts, gs = photosynthesis_model.Model.calculate_An(element.surfacic_nitrogen, element.width, element.height, PAR, meteo_df['air_temperature'][t_photosynthesis_model],
                                    meteo_df['ambient_CO2'][t_photosynthesis_model], meteo_df['humidity'][t_photosynthesis_model],
                                    meteo_df['Wind'][t_photosynthesis_model], organ_type)
                                element.Ag = Ag
                                element.An = An
                                element.Rd = Rd
                                element.Tr = Tr
                                element.Ts = Ts
                                element.gs = gs

            for t_cn_model in xrange(t_photosynthesis_model, t_photosynthesis_model + photosynthesis_model_ts, cn_model_ts):
                # run the model of CN exchanges ; the population is internally updated by the model of CN exchanges
                simulation_.run(start_time=t_cn_model, stop_time=t_cn_model+cn_model_ts, number_of_output_steps=cn_model_ts+1)
                
                # format the outputs to Pandas dataframes
                all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = simulation_.format_outputs()
            
                all_plants_df_list.append(all_plants_df)
                all_axes_df_list.append(all_axes_df)
                all_phytomers_df_list.append(all_phytomers_df)
                all_organs_df_list.append(all_organs_df)
                all_elements_df_list.append(all_elements_df)

        global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
        global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
        global_plants_df.to_csv(plants_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
        global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
        global_axes_df.to_csv(axes_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
        global_phytomers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
        global_phytomers_df.to_csv(phytomers_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
        global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
        global_organs_df.to_csv(organs_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
        global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
        global_elements_df.to_csv(elements_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

        execution_time = int(time.time()-t0)
        print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

    ######POST-PROCESSING##
    if make_graphs:
        x_name = 't'
        x_label='Time (Hour)'

        # Photosynthetic organs
        ph_elements_output_df = pd.read_csv(elements_outputs_filepath)

        graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (�mol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (�mol m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration (mm H$_{2}$0 h$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (�mol C h$^{-1}$)', 'Ts': u'Temperature surface (�C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
                           'Conc_TriosesP': u'[TriosesP] (�mol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (�mol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (�mol g$^{-1}$ mstruct)',
                           'Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (�mol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'SLN': u'Surfacic nitrogen content (g m$^{-2}$)',
                           'Nitrates_import': u'Total nitrates imported (�mol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (�mol N h$^{-1}$)',
                           'S_Amino_Acids': u'[Rate of amino acids synthesis] (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (�mol N g$^{-1}$ mstruct h$^{-1}$)',
                           'Loading_Sucrose': u'Loading Sucrose (�mol C sucrose g$^{-1}$ mstruct h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (�mol N amino acids g$^{-1}$ mstruct h$^{-1}$)',
                           'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (�mol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (�mol C h$^{-1}$)', 'R_residual': u'Respiration residual (�mol C h$^{-1}$)',
                           'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)', 'remob_starch_senescence': u'Remobilized C from starch (�mol h$^{-1}$)', 'remob_fructan_senescence': u'Remobilized C from fructan (�mol h$^{-1}$)', 'remob_proteins_senescence': u'Remobilized N from prot (�mol h$^{-1}$)'}


        for org_ph in (['Lamina'], ['Sheath'], ['Internode'], ['Peduncle', 'Chaff']):
            for variable_name, variable_label in graph_variables_ph_elements.iteritems():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                tools.plot_cnwheat_ouputs(ph_elements_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org_ph},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # Roots, grains and phloem
        organs_output_df = pd.read_csv(organs_outputs_filepath)

        graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                            'Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (�mol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                            'Uptake_Nitrates':u'Nitrates uptake (�mol h$^{-1}$)', 'Unloading_Sucrose': u'Unloaded_Sucrose (�mol C sucrose g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (�mol N AA g$^{-1}$ mstruct h$^{-1}$)',
                            'S_Amino_Acids': u'Rate of amino acids synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (�mol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (�mol N AA h$^{-1}$)',
                            'R_Nnit_upt': u'Respiration nitrates uptake (�mol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (�mol C h$^{-1}$)', 'R_residual': u'Respiration residual (�mol C h$^{-1}$)',
                            'R_grain_growth_struct': u'Respiration grain structural growth (�mol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (�mol C h$^{-1}$)',
                            'R_growth': u'Growth respiration of roots (�mol C h$^{-1}$)', 'Nstruct_N_growth': u'Growth of Nstruct (�mol N h$^{-1}$)', 'mstruct_growth': u'growth of structural dry mass (�mol C g$^{-1}$ mstruct h$^{-1}$)', 'mstruct_senescence': u'Structural dry mass lost by root senescence (g h$^{-1}$)','mstruct': u'Structural mass (g)',
                            'C_exudation': u'Carbon lost by root exudation (�mol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (�mol N g$^{-1}$ mstruct h$^{-1}$'}

        for org in (['Roots'], ['Grains'], ['Phloem']):
            for variable_name, variable_label in graph_variables_organs.iteritems():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # Soil
        graph_name = 'Conc_Nitrates_Soil.PNG'
        tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = 'Conc_Nitrates',
                              x_label=x_label,
                              y_label=u'[Nitrates] (�mol m$^{-3}$)',
                              filters={'organ': 'Soil'},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)


        # Integrated graphs
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)

        N_compartments = ['nitrates', 'amino_acids', 'proteins']
        t = ph_elements_output_df['t'].unique()

        # Photosynthetic organs
        ph_elements_output_df_group = ph_elements_output_df.groupby('t').sum()
        ax1.plot(t, ph_elements_output_df_group[N_compartments[0]]*14*1E-3, label='Nitrates ph', color='g', linestyle ="-.", linewidth=2)
        ax1.plot(t, ph_elements_output_df_group[N_compartments[1]]*14*1E-3, label='Amino acids ph', color='g', linestyle =":", linewidth=2)
        ax1.plot(t, ph_elements_output_df_group[N_compartments[2]]*14*1E-3, label='Proteins ph', color='g', linestyle ="--", linewidth=2)
        Nstruct_ph = ph_elements_output_df_group['Nstruct']*1E3
        ax1.plot(t, Nstruct_ph, label='Nstruct ph', color='g', linestyle ="-", linewidth=2)

        # Roots
        roots_output_df_group = organs_output_df[organs_output_df['organ']=='Roots'].groupby('t').sum()
        ax1.plot(t, roots_output_df_group[N_compartments[0]]*14*1E-3, label='Nitrates roots', color='k', linestyle ="-.", linewidth=2)
        ax1.plot(t, roots_output_df_group[N_compartments[1]]*14*1E-3, label='Amino acids roots', color='k', linestyle =":", linewidth=2)
        Nstruct_roots = roots_output_df_group['Nstruct']*1E3
        ax1.plot(t, Nstruct_roots, label='Nstruct roots', color='k', linestyle ="-", linewidth=2)

        # Phloem
        AA_phlo = organs_output_df[organs_output_df['organ']=='Phloem'].groupby('t').sum()[N_compartments[1]]*14*1E-3
        ax1.plot(t, AA_phlo, label='Amino acid phloem', color='b', linestyle =":", linewidth=2)

        # Grains
        prot_grains = organs_output_df[organs_output_df['organ']=='Grains'].groupby('t').sum()[N_compartments[2]]*14*1E-3
        ax1.plot(t, prot_grains, label='Proteins grains', color='y', linestyle ="--", linewidth=2)

        # Total
        N_tot = ph_elements_output_df_group[N_compartments].sum(axis=1)*14*1E-3 + Nstruct_ph + roots_output_df_group[N_compartments].sum(axis=1)*14*1E-3 + Nstruct_roots + AA_phlo + prot_grains
        ax1.plot(t, N_tot, label='Total', color='r', linewidth=2)
        ax2.plot(t, ph_elements_output_df_group[N_compartments].sum(axis=1)*14*1E-3 + Nstruct_ph, label='N ph', color='g', linewidth=2)
        ax2.plot(t, roots_output_df_group[N_compartments].sum(axis=1)*14*1E-3 + Nstruct_roots, label='N roots', color='k', linewidth=2)
        ax2.plot(t, AA_phlo, label='N phloem', color='b', linewidth=2)
        ax2.plot(t, prot_grains, label='N grains', color='y', linewidth=2)
        ax2.plot(t, N_tot, label='Total', color='r', linewidth=2)

        # Formatting
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.set_ylabel('Total tiller N mass (mg)')
        ax2.set_ylabel('Total tiller N mass (mg)')
        ax2.set_xlabel('Day')
        plt.tight_layout(rect=[0, 0, 0.67, .95])
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Cumulative_N.PNG'), dpi=200, format='PNG')
        plt.close()

        # Respiration plots
        ph_elements_output_df['day'] = ph_elements_output_df['t']//24+1
        organs_output_df['day'] = organs_output_df['t']//24+1
        days = ph_elements_output_df['day'].unique()

        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].sum()*12E-9/140E-4
        R_Nnit_red = ph_elements_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/140E-4 + organs_output_df.groupby('day')['R_Nnit_red'].sum()*12E-9/140E-4
        R_residual = ph_elements_output_df.groupby('day')['R_residual'].sum()*12E-9/140E-4 + organs_output_df.groupby('day')['R_residual'].sum()*12E-9/140E-4
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].sum()*12E-9/140E-4
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].sum()*12E-9/140E-4
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].sum()*12E-9/140E-4
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(days, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(days, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(days, R_residual, label='R_residual')
        ax1.plot(days, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(days, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(days, R_grain_growth_starch, label='R_grain_growth_starch')

        # Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Total tiller respiration (kg C m$^{-2}$ d$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Respiration_total.PNG'), dpi=200, format='PNG')
        plt.close()

        # 2nd plot
        R_phloem_loading = ph_elements_output_df.groupby('day')['R_phloem_loading'].mean()
        R_Nnit_red = ph_elements_output_df.groupby('day')['R_Nnit_red'].mean() + organs_output_df.groupby('day')['R_Nnit_red'].mean()
        R_residual = ph_elements_output_df.groupby('day')['R_residual'].mean() + organs_output_df.groupby('day')['R_residual'].mean()
        R_Nnit_upt = organs_output_df.groupby('day')['R_Nnit_upt'].mean()
        R_grain_growth_struct =  organs_output_df.groupby('day')['R_grain_growth_struct'].mean()
        R_grain_growth_starch =  organs_output_df.groupby('day')['R_grain_growth_starch'].mean()
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(days, R_phloem_loading, label='R_phloem_loading')
        ax1.plot(days, R_Nnit_red, label='R_Nnit_red')
        ax1.plot(days, R_residual, label='R_residual')
        ax1.plot(days, R_Nnit_upt, label='R_Nnit_upt')
        ax1.plot(days, R_grain_growth_struct, label='R_grain_growth_struct')
        ax1.plot(days, R_grain_growth_starch, label='R_grain_growth_starch')

        # Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Mean hourly tiller respiration (�mol C h$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Respiration_total2.PNG'), dpi=200, format='PNG')
        plt.close()

        # Daily N uptake plot
        organs_output_df['day'] = organs_output_df['t']//24+1
        fig, ax1 = plt.subplots(1, sharex=True)
        days = organs_output_df['day'].unique()
        daily_N_uptk = organs_output_df.groupby('day')['Uptake_Nitrates'].aggregate(np.sum)
        ax1.plot(days, daily_N_uptk)

        # Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Daily nitrates uptake (�mol N day$^{-1}$)')
        ax1.set_xlabel('Day')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Uptake_Nitrates_Roots_daily.PNG'), dpi=200, format='PNG')
        plt.close()


        # Total dry mass accumulation
        C_MOLAR_MASS = 12
        N_MOLAR_MASS = 14

        ## TriosesP
        TRIOSESP_MOLAR_MASS_C_RATIO = 0.21
        C_triosesP = ph_elements_output_df.groupby('t')['triosesP'].aggregate(np.sum)
        triosesP = (C_triosesP*1E-6*C_MOLAR_MASS)/TRIOSESP_MOLAR_MASS_C_RATIO
        ## Sucrose
        SUCROSE_MOLAR_MASS_C_RATIO = 0.42
        C_sucrose = ph_elements_output_df.groupby('t')['sucrose'].aggregate(np.sum) + organs_output_df.groupby('t')['sucrose'].aggregate(np.sum)
        sucrose = (C_sucrose*1E-6*C_MOLAR_MASS)/SUCROSE_MOLAR_MASS_C_RATIO
        ## Starch
        HEXOSE_MOLAR_MASS_C_RATIO = 0.4
        C_starch = ph_elements_output_df.groupby('t')['starch'].aggregate(np.sum) + organs_output_df.groupby('t')['starch'].aggregate(np.sum)
        starch = (C_starch*1E-6*C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO
        ## Fructans
        C_fructan = ph_elements_output_df.groupby('t')['fructan'].aggregate(np.sum)
        fructan = (C_fructan*1E-6*C_MOLAR_MASS)/HEXOSE_MOLAR_MASS_C_RATIO
        ## Structural dry mass
        RATIO_C_MSTRUCT = 0.384
        mstruct = ph_elements_output_df.groupby('t')['mstruct'].aggregate(np.sum)/RATIO_C_MSTRUCT + organs_output_df.groupby('t')['mstruct'].aggregate(np.sum)# supposed to include Nstruct also?
        ## Nitrates
        NITRATES_MOLAR_MASS_N_RATIO = 0.23
        N_nitrates = ph_elements_output_df.groupby('t')['nitrates'].aggregate(np.sum) + organs_output_df.groupby('t')['nitrates'].aggregate(np.sum)
        nitrates = (N_nitrates*1E-6*N_MOLAR_MASS)/NITRATES_MOLAR_MASS_N_RATIO
        ## Amino acids
        AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145
        N_amino_acids = ph_elements_output_df.groupby('t')['amino_acids'].aggregate(np.sum) + organs_output_df.groupby('t')['amino_acids'].aggregate(np.sum)
        amino_acids = (N_amino_acids*1E-6*N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO
        ## Proteins
        N_proteins = ph_elements_output_df.groupby('t')['proteins'].aggregate(np.sum) + organs_output_df.groupby('t')['proteins'].aggregate(np.sum)
        proteins = (N_proteins*1E-6*N_MOLAR_MASS)/AMINO_ACIDS_MOLAR_MASS_N_RATIO

        total_dry_mass = triosesP + sucrose + starch + fructan + mstruct + nitrates + amino_acids + proteins
        fig, ax1 = plt.subplots(1, sharex=True)
        ax1.plot(t, total_dry_mass)
        # Formatting
        ax1.legend(prop={'size':10}, framealpha=0.5)
        ax1.set_ylabel(u'Total dry mass (g)')
        ax1.set_xlabel('Hour')
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Total_dry_mass.PNG'), dpi=200, format='PNG')
        plt.close()

if __name__ == '__main__':
    compute_CN_distrib(run_simu=True, make_graphs=True)
##    # Profiling
##    filename = 'profile.pstats'
##    profile.run('compute_CN_distrib(make_graphs=True)', filename)
##    stats = pstats.Stats(filename)
##    stats.sort_stats('time')