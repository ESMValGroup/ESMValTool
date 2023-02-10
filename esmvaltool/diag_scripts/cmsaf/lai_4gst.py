"""
ESMValTool diagnostic for GLOBMAP LAI 4GST.

xxxxxxxxxxxxx add text here
"""

import logging

import iris
import iris.coord_categorisation as icc
import iris.quickplot as qplt

import matplotlib.pyplot as plt
import numpy as np
import datetime
from scipy.stats import linregress

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

ARGMAX = iris.analysis.Aggregator("argmax",
                                  np.argmax,
                                  units_func=lambda units: 1)

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


# def _make_plots(lst_diff_data, lst_diff_data_low, lst_diff_data_high, config):
#    """Create and save the output figure.
#    """


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

    list_of_years = obs_mean.coord('year').points
    list_of_times = obs_mean.coord('time').points
    # use this to be able to index the whole timeseries

    times_of_year_max = obs_mean.aggregated_by('year', iris.analysis.MAX)

    qplt.plot(obs_mean)
    plt.savefig('1_obs_mean.png')
    plt.close()

    for i, this_year in enumerate(np.unique(list_of_years)):
        print(i, this_year)
        if i == 0: continue
        # this skips the first year as we dont have lai data from before to shift
        #times_wanted = iris.Constraint(time

        central_time = times_of_year_max[i].coord('time').points[0]

        early = central_time - 182
        late = central_time + 182

        point_list_1 = np.where(list_of_times < late)[0]
        point_list_2 = np.where(list_of_times > early)[0]
        wanted = np.intersect1d(point_list_1, point_list_2)
        

        X = obs_mean[wanted]

        qplt.plot(X, linewidth=4, color='k')

        cols = ['red', 'green', 'blue', 'magenta']

        this_year_points = X.coord('time').points
        print(this_year_points)
        this_regress = []

        # check for EVG first
        evgt = (np.max(X.data) - np.min(X.data)) / np.mean(X.data)

        for j in range(4):

            this_quarter_1 = np.where(this_year_points >= this_year_points[0] + j * 92)
            this_quarter_2 = np.where(this_year_points <= this_year_points[0] + (j + 1)* 92)
            print(this_quarter_2)
            quarter_points = np.intersect1d(this_quarter_1, this_quarter_2)
            print(quarter_points)
            data = X[quarter_points].data
            days = X[quarter_points].coord('day_of_year').points

            qplt.plot(X[quarter_points], '--', linewidth=2, color=cols[j])

            regress = linregress(days, data)

            if regress.slope <=0:
                col = 'red'
            else:
                col ='blue'

            qplt.plot(X[quarter_points] * 0, linewidth=10, color=col)

            if evgt < 0.25:
                col = 'green'
            else:
                col = 'black'

            qplt.plot(X * 0 + 10, linewidth=10, color=col)

            this_regress.append(np.sign(regress.slope))

        plt.title(f'{this_year} {this_regress}')

        plt.savefig(f'2_{this_year}.png')
        plt.close()

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
    obs_lai_mean = obs_lai.collapsed(['latitude','longitude'],
                                     iris.analysis.MEAN)

    calc_4gst(obs_lai_mean)

    # Provenance
    # Get this information form the data cubes
    # add this back in using LST diagnostic as template

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
