"""
ESMValTool diagnostic for ESA CCI LST data.

The code uses the all time average monthly data.
The ouptput is a timeseries plot of the mean differnce of
CCI LST to model ensemble average, with the ensemble spread
represented by a standard deviation either side of the mean.
"""

import logging

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


def _make_plots(lst_diff_data, lst_diff_data_low, lst_diff_data_high, config):
    """Create and save the output figure.

    The plot is a mean differnce with +/- one standard deviation
    of the model spread,

    Inputs:
    lst_diff_data = cube of the mean difference
    lst_diff_data_low = cube of the mean difference
                        with model minus standard deviation
    lst_diff_data_high = cube of the mean difference
                        with model plus standard deviation
    config = The config dictionary from the preprocessor

    Outputs:
    Saved figure
    """
    fig, ax = plt.subplots(figsize=(20, 15))

    ax.plot(lst_diff_data.data, color='black', linewidth=4)
    ax.plot(lst_diff_data_low.data, '--', color='blue', linewidth=3)
    ax.plot(lst_diff_data_high.data, '--', color='blue', linewidth=3)
    ax.fill_between(range(len(lst_diff_data.data)),
                    lst_diff_data_low.data,
                    lst_diff_data_high.data,
                    color='blue',
                    alpha=0.25)

    # make X ticks
    x_tick_list = []
    time_list = lst_diff_data.coord('time').units.num2date(
        lst_diff_data.coord('time').points)
    for item in time_list:
        if item.month == 1:
            x_tick_list.append(item.strftime('%Y %b'))
        elif item.month == 7:
            x_tick_list.append(item.strftime('%b'))
        else:
            x_tick_list.append('')

    ax.set_xticks(range(len(lst_diff_data.data)))
    ax.set_xticklabels(x_tick_list, fontsize=18, rotation=45)

    # make Y ticks
    y_lower = np.floor(lst_diff_data_low.data.min())
    y_upper = np.ceil(lst_diff_data_high.data.max())
    ax.set_yticks(np.arange(y_lower, y_upper + 0.1, 2))
    ax.set_yticklabels(np.arange(y_lower, y_upper + 0.1, 2), fontsize=18)
    ax.set_ylim((y_lower - 0.1, y_upper + 0.1))

    ax.set_xlabel('Date', fontsize=20)
    ax.set_ylabel('Difference / K', fontsize=20)

    ax.grid()

    lons = lst_diff_data.coord('longitude').bounds
    lats = lst_diff_data.coord('latitude').bounds

    ax.set_title('Area: lon %s lat %s' % (lons[0], lats[0]), fontsize=22)

    fig.suptitle('ESACCI LST - CMIP6 Historical Ensemble Mean', fontsize=24)

    plot_path = get_plot_filename('timeseries', config)
    plt.savefig(plot_path)
    plt.close('all')  # Is this needed?


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
    figures made by make_plots.
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config['input_data'].values()

    loaded_data = {}
    ancestor_list = []
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes
        ancestor_list.append(ancestors['ts'][0])

    # loaded data is a nested dictionary
    # KEY1 model ESACCI-LST or something else
    # KEY2 is ts, the surface temperature
    # ie loaded_data['ESACCI-LST']['ts'] is the CCI cube
    #    loaded_data['MultiModelMean']['ts'] is CMIP6 data, emsemble means
    #    similarly dor Std, see preprocessor

    # The Diagnostic uses CCI - MODEL

    # CMIP data had 360 day calendar, CCI data has 365 day calendar
    # Assume the loaded data is all the same shape
    loaded_data['MultiModelMean']['ts'].remove_coord('time')
    loaded_data['MultiModelMean']['ts'].add_dim_coord(
        loaded_data['ESACCI-LST']['ts'].coord('time'), 0)
    loaded_data['MultiModelStd_Dev']['ts'].remove_coord('time')
    loaded_data['MultiModelStd_Dev']['ts'].add_dim_coord(
        loaded_data['ESACCI-LST']['ts'].coord('time'), 0)

    # Make a cube of the LST difference, and with +/- std of model variation
    lst_diff_cube = loaded_data['ESACCI-LST']['ts'] - loaded_data[
        'MultiModelMean']['ts']
    lst_diff_cube_low = loaded_data['ESACCI-LST']['ts'] - (
        loaded_data['MultiModelMean']['ts'] +
        loaded_data['MultiModelStd_Dev']['ts'])
    lst_diff_cube_high = loaded_data['ESACCI-LST']['ts'] - (
        loaded_data['MultiModelMean']['ts'] -
        loaded_data['MultiModelStd_Dev']['ts'])

    # Plotting
    _make_plots(lst_diff_cube, lst_diff_cube_low, lst_diff_cube_high, config)

    # Provenance
    # Get this information form the data cubes
    data_attributes = {}
    data_attributes['start_year'] = lst_diff_cube.coord('time').units.num2date(
        lst_diff_cube.coord('time').points)[0].year
    data_attributes['end_year'] = lst_diff_cube.coord('time').units.num2date(
        lst_diff_cube.coord('time').points)[-1].year
    data_attributes['lat_south'] = lst_diff_cube.coord('latitude').bounds[0][0]
    data_attributes['lat_north'] = lst_diff_cube.coord('latitude').bounds[0][1]
    data_attributes['lon_west'] = lst_diff_cube.coord('longitude').bounds[0][0]
    data_attributes['lon_east'] = lst_diff_cube.coord('longitude').bounds[0][1]
    data_attributes['ensembles'] = ''

    for item in input_metadata:
        if 'ESACCI' in item['alias'] or 'MultiModel' in item[
                'alias'] or 'OBS' in item['alias']:
            continue
        data_attributes['ensembles'] += "%s " % item['alias']

    record = _get_provenance_record(data_attributes, ancestor_list)
    plot_file = get_plot_filename('timeseries', config)
    with ProvenanceLogger(config) as provenance_logger:
        provenance_logger.log(plot_file, record)


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
