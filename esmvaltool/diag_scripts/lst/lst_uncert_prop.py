"""
ESMValTool diagnostic for ESA CCI LST V3 data - Uncertainity Propagation
"""

# Expected keys for OBS LST data in loaded_data:
# ts_day
# ts_night
# lst_unc_loc_sfc_day 
# lst_unc_loc_sfc_night
# lst_unc_loc_atm_day 
# lst_unc_loc_atm_night
# lst_unc_sys_day 
# lst_unc_sys_night
# lst_unc_ran_day 
# lst_unc_ran_night

import logging
import time
import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)


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


def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.

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


def _diagnostic(config):
    """Perform the control for the ESA CCI LST diagnostic.

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------
    
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()
    loaded_data = {}
    ancestor_list = []

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        print(dataset, metadata)
        cubes, ancestors = _get_input_cubes(metadata)
        print(cubes, ancestors)
        loaded_data[dataset] = cubes

    # for now leave this so it appears in the log file
    print('###########################################')
    print(loaded_data)
    print('##########################################')


    # add calls to propagation equations here

    ts_day_mean = eq_arithmetic_mean(loaded_data['ESACCI-LST']['ts_day'])


    print('*********************')
    print(f'ts_day mean {ts_day_mean.data}')

# make function for each propagation equation
def eq_arithmetic_mean(cube):
    """Arithmetic mean of cube, across latitude and longitude
    ATBD eq 1
    """

    out_cube = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)

    return out_cube


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
