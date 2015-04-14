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

import logging
import logging.config
import json

import pandas as pd
import matplotlib.pyplot as plt

from cnwheat import simulation
from cnwheat import model as cnwheat_model
from farquharwheat import model as photosynthesis_model
from cnwheat import tools

t0 = time.time()

INPUTS_DIRPATH = 'inputs'

PAR_FILENAME = 'PAR.csv'

METEO_FILENAME = 'meteo.csv'

OUTPUTS_DIRPATH = 'outputs'

PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
PHYTOMERS_OUTPUTS_FILENAME = 'phytomers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

OUTPUTS_PRECISION = 6

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # set to logging.DEBUG for debugging

def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_compartments=False, log_derivatives=False, log_model=False):
    """Setup logging configuration.

    :Parameters:

        - `config_filepath` (:class:`str`) - the file path of the logging
          configuration.

        - `level` (:class:`int`) - the global level of the logging. Use either
          `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or
          `logging.CRITICAL`.

        - `log_compartments` (:class:`bool`) - if `True`, log the values of the compartments.
          `False` otherwise.

        - `log_derivatives` (:class:`bool`) - if `True`, log the values of the derivatives.
          `False` otherwise.

        - `log_model` (:class:`bool`) - if `True`, log the messages from :mod:`cnwheat.model`.
          `False` otherwise.

    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    logging.getLogger('cnwheat.compartments').disabled = not log_compartments # set to False to log the compartments
    logging.getLogger('cnwheat.derivatives').disabled = not log_derivatives # set to False to log the derivatives
    logging.getLogger('cnwheat.model').disabled = not log_model # set to False to log the derivatives

setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL)


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


population = cnwheat_model.Population()
 
plant = cnwheat_model.Plant(index=1)
population.plants.append(plant)
 
axis = cnwheat_model.Axis(axis_type=cnwheat_model.Axis.Types.MAIN_STEM, index=0)
plant.axes.append(axis)
 
axis.grains = cnwheat_model.Grains(starch=0, structure=4000, proteins=170)
 
axis.roots = cnwheat_model.Roots(mstruct=0.504, Nstruct=0.01, sucrose=0, nitrates=0, amino_acids=0)
 
axis.phloem = cnwheat_model.Phloem(sucrose=0, amino_acids=0)
 
# Phytomer 1
phytomer1 = cnwheat_model.Phytomer(index=1)
 
phytomer1.lamina = cnwheat_model.Lamina()
lamina_element = cnwheat_model.LaminaElement(area=0.00346, mstruct=0.14, Nstruct=0.00102, width= 0.018, height=0.6,
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                amino_acids=0, proteins=380)
phytomer1.lamina.exposed_element = lamina_element
 
phytomer1.sheath = cnwheat_model.Sheath()
sheath_element = cnwheat_model.SheathElement(area=0.0006, mstruct=0.103, Nstruct=0.00068, width=0.0011, height=0.5,
                                starch=0, sucrose=0, triosesP=0, fructan=0,
                                nitrates=0 , amino_acids=0, proteins=130)
phytomer1.sheath.exposed_element = sheath_element
 
# Internode enclosed
phytomer1.internode = cnwheat_model.Internode()
internode_enclosed_element = cnwheat_model.InternodeElement(area=0.0012, mstruct=0.148, Nstruct=0.00067, width=0.00257, height=0.3,
                                      starch=0, sucrose=0, triosesP=0, fructan=0,
                                      nitrates=0, amino_acids=0, proteins=20)
phytomer1.internode.enclosed_element = internode_enclosed_element
 
# Internode exposed
internode_exposed_element = cnwheat_model.InternodeElement(area=0.0003, mstruct=0.04, Nstruct=0.00018, width=0.00257, height=0.4,
                                      starch=0, sucrose=0, triosesP=0, fructan=0,
                                      nitrates=0, amino_acids=0, proteins=90)
phytomer1.internode.exposed_element = internode_exposed_element
 
axis.phytomers.append(phytomer1)
 
# Phytomer 2
phytomer2 = cnwheat_model.Phytomer(index=2)
 
phytomer2.lamina = cnwheat_model.Lamina()
lamina_element = cnwheat_model.LaminaElement(area=0.0034, mstruct=0.09, Nstruct=0.00083, width= 0.014, height=0.38,
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                amino_acids=0, proteins=210)
phytomer2.lamina.exposed_element = lamina_element
 
phytomer2.sheath = cnwheat_model.Sheath()
sheath_element = cnwheat_model.SheathElement(area=0.0005, mstruct=0.069, Nstruct=0.00021, width=0.00091, height=0.3,
                                starch=0, sucrose=0, triosesP=0, fructan=0,
                                nitrates=0 , amino_acids=0, proteins=47)
phytomer2.sheath.exposed_element = sheath_element
 
phytomer2.internode = cnwheat_model.Internode()
internode_element = cnwheat_model.InternodeElement(area=0.0004, mstruct=0.18, Nstruct=0.00033, width=0.00099, height=0.18,
                                      starch=0, sucrose=0, triosesP=0, fructan=0,
                                      nitrates=0, amino_acids=0, proteins=20)
phytomer2.internode.exposed_element = internode_element
 
axis.phytomers.append(phytomer2)
 
# Phytomer 3
phytomer3 = cnwheat_model.Phytomer(index=3)
 
phytomer3.lamina = cnwheat_model.Lamina()
lamina_element = cnwheat_model.LaminaElement(area=0.00228, mstruct=0.05, Nstruct=0.00053, width= 0.0125, height=0.24,
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                amino_acids=0, proteins=85)
phytomer3.lamina.exposed_element = lamina_element
 
phytomer3.sheath = cnwheat_model.Sheath()
sheath_element = cnwheat_model.SheathElement(area=0.0004, mstruct=0.043, Nstruct=0.00011, width=0.00051, height=0.18,
                                starch=0, sucrose=0, triosesP=0, fructan=0,
                                nitrates=0 , amino_acids=0, proteins=13)
phytomer3.sheath.exposed_element = sheath_element
 
phytomer3.internode = cnwheat_model.Internode()
internode_element = cnwheat_model.InternodeElement(area=0.00025, mstruct=0.154, Nstruct=0.00014, width=0.00093, height=0.08,
                                      starch=0, sucrose=0, triosesP=0, fructan=0,
                                      nitrates=0, amino_acids=0, proteins=20)
phytomer3.internode.exposed_element = internode_element
 
axis.phytomers.append(phytomer3)
 
# Phytomer 4 (reproductive)
phytomer4 = cnwheat_model.Phytomer(index=4)
 
# Enclosed peduncle
phytomer4.peduncle = cnwheat_model.Peduncle()
peduncle_enclosed_element = cnwheat_model.PeduncleElement(area=0.00155, mstruct=0.168, Nstruct=0.00085, width= 0.00349, height=0.65,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=30)
phytomer4.peduncle.enclosed_element = peduncle_enclosed_element
 
# Exposed peduncle
peduncle_exposed_element = cnwheat_model.PeduncleElement(area=0.00085, mstruct=0.089, Nstruct=0.00045, width= 0.00349, height=0.5,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=180)
phytomer4.peduncle.exposed_element = peduncle_exposed_element
axis.phytomers.append(phytomer4)
 
# Phytomer 5 (reproductive)
phytomer5 = cnwheat_model.Phytomer(index=5)
phytomer5.chaff = cnwheat_model.Chaff()
chaff_element = cnwheat_model.ChaffElement(area=0.00075, mstruct=0.21, Nstruct=0.00107, width=0.00265, height= 0.7, starch=0,
                              sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0,
                              proteins=260)
phytomer5.chaff.exposed_element = chaff_element
axis.phytomers.append(phytomer5)
 
# Get PAR data
PAR_filepath = os.path.join(INPUTS_DIRPATH, PAR_FILENAME)
PAR_df = pd.read_csv(PAR_filepath)
PAR_grouped = PAR_df.groupby(simulation.CNWheat.ELEMENTS_INDEXES)
 
# get meteo data
meteo_df = read_t_data(INPUTS_DIRPATH, METEO_FILENAME)
 
# initialize the model of CN exchanges
cnwheat_ = simulation.CNWheat(population=population)
 
# run the models
start_time = 0
stop_time = 10 # 960
photosynthesis_model_ts = 2
cn_model_ts = 1 #241
 
all_plants_df_list = []
all_axes_df_list = []
all_phytomers_df_list = []
all_organs_df_list = []
all_elements_df_list = []
 
for t_photosynthesis_model in xrange(start_time, stop_time, photosynthesis_model_ts):
    # update the population
    population.t = t_photosynthesis_model
    # run the model of photosynthesis and update the population
    for plant in population.plants:
        plant_index = plant.index
        for axis in plant.axes:
            axis_id = axis.id
            for phytomer in axis.phytomers:
                phytomer_index = phytomer.index
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    organ_type = organ.__class__.__name__
                    for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                        if element is None:
                            continue
                        PAR = PAR_grouped.get_group((t_photosynthesis_model, plant_index, axis_id, phytomer_index, organ_type, element_type)).PAR.values[0]
                        An, Rd, Tr, Ts, gs = photosynthesis_model.PhotosynthesisModel.calculate_An(element.surfacic_nitrogen, element.width, element.height, PAR, meteo_df['air_temperature'][t_photosynthesis_model],
                            meteo_df['ambient_CO2'][t_photosynthesis_model], meteo_df['humidity'][t_photosynthesis_model],
                            meteo_df['Wind'][t_photosynthesis_model], organ_type)
                        element.An = An
                        element.Rd = Rd
                        element.Tr = Tr
                        element.Ts = Ts
                        element.gs = gs
 
    for t_cn_model in xrange(t_photosynthesis_model, t_photosynthesis_model + photosynthesis_model_ts, cn_model_ts):
        # update the population
        population.t = t_cn_model
        # run the model of CN exchanges ; the population is internally updated by the model of CN exchanges
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = cnwheat_.run(start_time=t_cn_model, stop_time=t_cn_model+cn_model_ts, number_of_output_steps=cn_model_ts+1)
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
        all_elements_df_list.append(all_elements_df)
 
global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
global_plants_df.drop_duplicates(subset=simulation.CNWheat.PLANTS_INDEXES, inplace=True)
plants_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PLANTS_OUTPUTS_FILENAME)
global_plants_df.to_csv(plants_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
 
global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
global_axes_df.drop_duplicates(subset=simulation.CNWheat.AXES_INDEXES, inplace=True)
axes_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME)
global_axes_df.to_csv(axes_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
 
global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
global_phytomers_df.drop_duplicates(subset=simulation.CNWheat.PHYTOMERS_INDEXES, inplace=True)
phytomers_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PHYTOMERS_OUTPUTS_FILENAME)
global_phytomers_df.to_csv(phytomers_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
 
global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
global_organs_df.drop_duplicates(subset=simulation.CNWheat.ORGANS_INDEXES, inplace=True)
organs_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ORGANS_OUTPUTS_FILENAME)
global_organs_df.to_csv(organs_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
 
global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
global_elements_df.drop_duplicates(subset=simulation.CNWheat.ELEMENTS_INDEXES, inplace=True)
elements_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)
global_elements_df.to_csv(elements_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
 
execution_time = int(time.time()-t0)
print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))


####POST-PROCESSING##

graphs_dirpath = 'graphs' # graphs_dirpath must be an existing directory
x_name = 't'
x_label='Time (Hour)'

# Photosynthetic organs
ph_elements_output_df = pd.read_csv(elements_outputs_filepath)

graph_variables_ph_elements = {'An': u'Net photosynthesis (�mol m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration (mm H$_{2}$0 h$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (�mol m$^{-2}$ s$^{-1}$)', 'Ts': u'Temperature surface (�C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
                   'Conc_TriosesP': u'[TriosesP] (�mol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (�mol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (�mol g$^{-1}$ mstruct)',
                   'Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (�mol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'SLN': u'Surfacic nitrogen content (g m$^{-2}$)',
                   'Nitrates_import': u'Total nitrates imported (�mol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (�mol N h$^{-1}$)',
                   'S_Amino_Acids': u'[Rate of amino acids synthesis] (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (�mol N g$^{-1}$ mstruct h$^{-1}$)',
                   'Loading_Sucrose': u'Loading Sucrose (�mol C sucrose g$^{-1}$ mstruct h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (�mol N amino acids g$^{-1}$ mstruct h$^{-1}$)',
                   'green_area': u'Green area (m$^{2}$)', 'R_phloem': u'Respiration phloem loading (�mol h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (�mol h$^{-1}$)', 'R_residual': u'Respiration residual (�mol h$^{-1}$)'}


for org_ph in (['Lamina'], ['Sheath'], ['Internode'], ['Peduncle', 'Chaff']):
    for variable_name, variable_label in graph_variables_ph_elements.iteritems():
        graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
        tools.plot_cnwheat_ouputs(ph_elements_output_df,
                      x_name = x_name,
                      y_name = variable_name,
                      x_label=x_label,
                      y_label=variable_label,
                      filters={'organ': org_ph},
                      plot_filepath=os.path.join(graphs_dirpath, graph_name),
                      explicit_label=False)

# Roots, grains and phloem
organs_output_df = pd.read_csv(organs_outputs_filepath)

graph_variables_ph_elements = {'Conc_Sucrose':u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                    'Conc_Nitrates_Soil':u'[Nitrates] (�mol m$^{-3}$)','Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (�mol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                    'Uptake_Nitrates':u'Nitrates uptake (�mol h$^{-1}$)', 'Unloading_Sucrose': u'Unloaded_Sucrose (�mol C sucrose g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (�mol N AA g$^{-1}$ mstruct h$^{-1}$)',
                    'S_Amino_Acids': u'Rate of amino acids synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (�mol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (�mol N AA h$^{-1}$)',
                    'R_Nnit_upt': u'Respiration nitrates uptake (�mol h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (�mol h$^{-1}$)', 'R_residual': u'Respiration residual (�mol h$^{-1}$)'}

for org in (['Roots'], ['Grains'], ['Phloem']):
    for variable_name, variable_label in graph_variables_ph_elements.iteritems():
        graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
        tools.plot_cnwheat_ouputs(organs_output_df,
                      x_name = x_name,
                      y_name = variable_name,
                      x_label=x_label,
                      y_label=variable_label,
                      filters={'organ': org},
                      plot_filepath=os.path.join(graphs_dirpath, graph_name),
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
Nstruct_ph = ph_elements_output_df_group['Nstruct']
ax1.plot(t, Nstruct_ph, label='Nstruct ph', color='g', linestyle ="-", linewidth=2)

# Roots
roots_output_df_group = organs_output_df[organs_output_df['organ']=='Roots'].groupby('t').sum()
ax1.plot(t, roots_output_df_group[N_compartments[0]]*14*1E-3, label='Nitrates roots', color='k', linestyle ="-.", linewidth=2)
ax1.plot(t, roots_output_df_group[N_compartments[1]]*14*1E-3, label='Amino acids roots', color='k', linestyle =":", linewidth=2)
Nstruct_roots = roots_output_df_group['Nstruct']
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
ax2.set_xlabel(x_label)
plt.tight_layout(rect=[0, 0, 0.67, .95])
plt.savefig(os.path.join(graphs_dirpath, 'Cumulative_N.PNG'), dpi=200, format='PNG')
plt.close()