# -*- coding: latin-1 -*-
"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test the CN-Wheat model.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    CSV files must contain only ASCII characters and ',' as separator.

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

import os

import numpy as np
import pandas as pd

from cnwheat import simulation
from cnwheat import model as cnwheat_model
from cnwheat import tools

import time, datetime

t0 = time.time()

INPUTS_DIRPATH = 'inputs'

PHOTOSYNTHESIS_DATA_FILENAME = 'photosynthesis_data.csv'
SENESCENCE_DATA_FILENAME = 'senescence_data.csv'

OUTPUTS_DIRPATH = 'outputs'

DESIRED_PLANTS_OUTPUTS_FILENAME = 'desired_plants_outputs.csv'
DESIRED_AXES_OUTPUTS_FILENAME = 'desired_axes_outputs.csv'
DESIRED_PHYTOMERS_OUTPUTS_FILENAME = 'desired_phytomers_outputs.csv'
DESIRED_ORGANS_OUTPUTS_FILENAME = 'desired_organs_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'

ACTUAL_PLANTS_OUTPUTS_FILENAME = 'actual_plants_outputs.csv'
ACTUAL_AXES_OUTPUTS_FILENAME = 'actual_axes_outputs.csv'
ACTUAL_PHYTOMERS_OUTPUTS_FILENAME = 'actual_phytomers_outputs.csv'
ACTUAL_ORGANS_OUTPUTS_FILENAME = 'actual_organs_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'


def test_run():

    population = cnwheat_model.Population()

    plant = cnwheat_model.Plant(index=1)
    population.plants.append(plant)

    axis = cnwheat_model.Axis()
    plant.axes.append(axis)

    axis.grains = cnwheat_model.Grains(starch=0, structure=4000, proteins=170)

    axis.roots = cnwheat_model.Roots(mstruct=0.504, Nstruct=0.01, sucrose=900, nitrates=250, amino_acids=60)

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

    # Get photosynthesis data
    photosynthesis_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_DATA_FILENAME)
    photosynthesis_data_df = pd.read_csv(photosynthesis_data_filepath)
    photosynthesis_data_grouped = photosynthesis_data_df.groupby(simulation.Simulation.ELEMENTS_INDEXES)
 
    # Get senescence and growth data
    senescence_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_DATA_FILENAME)
    senescence_data_df = pd.read_csv(senescence_data_filepath)
    senescence_data_grouped = senescence_data_df.groupby(simulation.Simulation.ELEMENTS_INDEXES)

    # initialize the simulator
    simulation_ = simulation.Simulation(population=population)

    start_time = 0
    stop_time = 48
    time_step = 4

    all_plants_df_list = []
    all_axes_df_list = []
    all_phytomers_df_list = []
    all_organs_df_list = []
    all_elements_df_list = []

    for t in xrange(start_time, stop_time, time_step):
        # update the population
        for plant in simulation_.population.plants:
            plant_index = plant.index
            for axis in plant.axes:
                axe_id = axis.id

                # Root growth and senescence
                group = senescence_data_grouped.get_group((t, plant_index, axe_id, 0, 'Roots', 'enclosed'))
                row_index = group.first_valid_index()
                axis.roots.mstruct_C_growth = group.mstruct_growth[row_index]
                axis.roots.Nstruct_N_growth = group.Nstruct_N_growth[row_index]
                axis.roots.mstruct_senescence = group.mstruct_senescence[row_index]
                axis.roots.mstruct = group.mstruct[row_index]
                axis.roots.Nstruct = group.Nstruct[row_index]

                for phytomer in axis.phytomers:
                    phytomer_index = phytomer.index
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        organ_type = organ.__class__.__name__
                        for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                            if element is None:
                                continue
                            group_photo = photosynthesis_data_grouped.get_group((t, plant_index, axe_id, phytomer_index, organ_type, element_type))
                            group_senesc = senescence_data_grouped.get_group((t, plant_index, axe_id, phytomer_index, organ_type, element_type))
                            row_index_photo = group_photo.first_valid_index()
                            row_index_sensc = group_senesc.first_valid_index()

                            # Senescence
                            element.green_area = group_senesc.green_area[row_index_sensc]
                            element.relative_delta_green_area = group_senesc.relative_delta_green_area[row_index_sensc]
                            element.mstruct = group_senesc.mstruct[row_index_sensc]
                            element.Nstruct = group_senesc.Nstruct[row_index_sensc]
                            element.surfacic_nitrogen = group_senesc.SLN[row_index_sensc]

                            element.Ag = group_photo.Ag[row_index_photo]
                            element.An = group_photo.An[row_index_photo]
                            element.Rd = group_photo.Rd[row_index_photo]
                            element.Tr = group_photo.Tr[row_index_photo]
                            element.Ts = group_photo.Ts[row_index_photo]
                            element.gs = group_photo.gs[row_index_photo]

        # run the model of CN exchanges ; the population is internally updated by the model of CN exchanges
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = simulation_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
        all_elements_df_list.append(all_elements_df)

    execution_time = int(time.time()-t0)
    print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_plants_df, DESIRED_PLANTS_OUTPUTS_FILENAME, ACTUAL_PLANTS_OUTPUTS_FILENAME)

    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_axes_df, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME)

    global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
    global_phytomers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_phytomers_df, DESIRED_PHYTOMERS_OUTPUTS_FILENAME, ACTUAL_PHYTOMERS_OUTPUTS_FILENAME)

    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_organs_df, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME)

    global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
    global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_elements_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)

if __name__ == '__main__':
    test_run()
