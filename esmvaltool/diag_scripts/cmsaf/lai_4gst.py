"""
ESMValTool diagnostic for GLOBMAP LAI 4GST.

xxxxxxxxxxxxx add text here
"""

import logging

import iris
import iris.coord_categorisation as icc

import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

ARGMAX = iris.analysis.Aggregator("argmax", np.argmax, units_func=lambda units: 1)

def _get_input_cubes(metadata):
    """Load the data files into cubes.

    Based on the hydrology diagnostic.

    Inputs:
    metadata = List of dictionaries made from the preprocessor config

    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
    inputs = {}
    ancestors = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        inputs[short_name] = cube
        ancestors[short_name] = [filename]

    return inputs, ancestors


def _make_plots(lst_diff_data, lst_diff_data_low, lst_diff_data_high, config):
    """Create and save the output figure.

    
    """
     


def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.
    xxxxxxxxxxxxxxxxxx update this to GLOBMAP LAI stuff
    Inputs:
    attributes = dictionary of ensembles/models used, the region bounds
                 and years of data used.
    ancestor_files = list of data files used by the diagnostic.

    Outputs:
    record = dictionary of provenance records.
    """
    caption = "Timeseries of ESA CCI LST difference to mean of "\
        + "model ensembles calculated over region bounded by latitude "\
        + "{lat_south} to {lat_north}, longitude {lon_west} to {lon_east} "\
        + "and for model/ensembles {ensembles}. "\
        + "Shown for years {start_year} to {end_year}.".format(**attributes)

    record = {
        'caption': caption,
        'statistics': ['mean', 'stddev'],
        'domains': ['reg'],
        'plot_types': ['times'],
        'authors': ['king_robert'],
        # 'references': [],
        'ancestors': ancestor_files
    }

    return record

def calc_4gst(obs_mean):
    """
    data_dict : the loaded_data dictionary
    """

    for this_year in [1992]: # make this cover all years

        this_year_cube = obs_mean.extract(iris.Constraint(year=this_year))

        print(this_year_cube)
        max_loc = this_year_cube.collapsed('time', ARGMAX)
        print('###############')
        print(max_loc.data)
        print(len(this_year_cube.coord('time').points))
        print(0/0)
    return None

def _diagnostic(config):
    """Perform the .......

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------
    xxxxxxxxxxxxxxxxxxxxxxxx
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()

    loaded_data = {}
    ancestor_list = [] # what is this actually used for?
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes
        # ancestor_list.append(ancestors['ts'][0])

    print(loaded_data)
    # OBS loaded_data['GLOBMAP_LAI']['lairpk']


    # note LST data needed some time coords changing
    # CMIP data had 360 day calendar, CCI data has 365 day calendar

    
    # obs lai, becareful with hard coded names here!!!11
    obs_lai = loaded_data['GLOBMAP_LAI']['lairpk']
    icc.add_day_of_year(obs_lai, 'time')
    icc.add_year(obs_lai, 'time')
    obs_lai_mean = obs_lai.collapsed(['latitude','longitude'], iris.analysis.MEAN)
    print(obs_lai_mean)
    print(obs_lai_mean.data)
    print(obs_lai_mean.coord('day_of_year'))
    

    calc_4gst(obs_lai_mean)


    # Plotting
    # _make_plots(lst_diff_cube, lst_diff_cube_low, lst_diff_cube_high, config)

    # Provenance
    # Get this information form the data cubes
    # add this back in using LST diagnostic as template

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
