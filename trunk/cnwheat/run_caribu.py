from alinea.adel.astk_interface import AdelWheat


from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight
from alinea.caribu.label import encode_label
import pandas as pd

#: the columns which define the topology in the input/output dataframe
DATAFRAME_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']
CARIBU_OUTPUTS = ['Eabsm2']


def to_dataframe(g, res_by_id):
    """
    Convert outputs of run_caribu to Pandas dataframe.

    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`).
        - `res_by_id` (:class:`dict`) - The outputs of run_caribu sorted by id.

    :Returns:
        The outputs in a dataframe.

    :Returns Type:
        :class:`pandas.DataFrame`
    """

    ids = []
    for vid in res_by_id['label'].keys():
        ind = (int(g.index(g.complex_at_scale(vid, 1))), g.label(g.complex_at_scale(vid, 2)), int(g.index(g.complex_at_scale(vid, 3))), g.label(g.complex_at_scale(vid, 4)), g.label(vid))
        ids.append(ind)

    ids_df = pd.DataFrame(sorted(ids), columns=DATAFRAME_TOPOLOGY_COLUMNS)
    Eabsm2_list = [res_by_id['Eabsm2'][Eabsm2] for Eabsm2 in sorted(res_by_id['Eabsm2'].keys())] # Used in order to sort values by keys. TODO: must be a faster way
    data_df = pd.DataFrame({'Eabsm2': Eabsm2_list})
    df = pd.concat([ids_df, data_df], axis=1)
    df.reset_index(drop=True, inplace=True)
    return df

def run_caribu(g, adel_wheat):
    """
    Run Caribu from the MTG. Write 'Eabsm2' variable in the MTG and a dataframe.

    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`).
        - `adel_wheat` - The model ADEL-Wheat.

    :Returns:
        The outputs in a dataframe.

    :Returns Type:
        :class:`pandas.DataFrame`
    """

    geom = g.property('geometry')

    # Instance of CaribuScene defined with the pattern
    c_scene = CaribuScene(pattern = adel_wheat.domain)

    # Diffuse light sources
    sky = GenSky.GenSky()(1, 'soc', 4, 5) # (Energy, soc/uoc, azimuts, zenits)
    sky = GetLight.GetLight(sky)
    c_scene.addSources(sky)

    # Labels
    labels = []
    for vid in geom.keys():
        if g.class_name(vid) == 'HiddenElement':
            continue
        elif g.class_name(vid) in ('LeafElement1', 'LeafElement'):
            plant_id = g.index(g.complex_at_scale(vid, scale=1))
            label = encode_label(opt_id=1, opak=1, plant_id=plant_id)[0]
        elif g.class_name(vid) == 'StemElement':
            plant_id = g.index(g.complex_at_scale(vid, scale=1))
            label = encode_label(opt_id=1, opak=0, plant_id=plant_id)[0]
        else:
            label = encode_label()[0]
            print 'Warning: unknown element type {}, vid={}'.format(g.class_name(vid), vid)
        labels.append(label)

    # Add scene to CaribuScene
    idmap = c_scene.add_Shapes(geom, canlabels = labels)

    # Run Caribu
    output = c_scene.runCaribu(direct = True, infinity=True)
##    # Visualisation
##    c_scene.plot(output=output)

    # Aggregation of results
    res_by_id = c_scene.output_by_id(output, idmap)

    # Update MTG
    if not 'Eabsm2' in g.properties():
        g.add_property('Eabsm2')
    g.property('Eabsm2').update(res_by_id['Eabsm2'])

    df = to_dataframe(g, res_by_id)

    return df