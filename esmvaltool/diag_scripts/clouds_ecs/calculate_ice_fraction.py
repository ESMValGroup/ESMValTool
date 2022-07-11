"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import cf_units
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
    io,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)

LINE_LEGEND = {
    'ECS_high_hist': 'ECS_high',
    'ECS_med_hist': 'ECS_med',
    'ECS_low_hist': 'ECS_low',
}

LINE_COLOR = {
    'ECS_high_hist': 'royalblue',
    'ECS_high_scen': 'royalblue',
    'ECS_med_hist': 'green',
    'ECS_med_scen': 'green',
    'ECS_low_hist': 'orange',
    'ECS_low_scen': 'orange',
    'CMIP6': 'firebrick',
    'CMIP5': 'royalblue',
    'CMIP3': 'darkcyan',
    'OBS': 'black'
}

LINE_DASH = {
    'ECS_high_hist': 'solid',
    'ECS_high_scen': 'dashed',
    'ECS_med_hist': 'solid',
    'ECS_med_scen': 'dashed',
    'ECS_low_hist': 'solid',
    'ECS_low_scen': 'dashed',
    'CMIP6': 'solid',
    'CMIP5': 'solid',
    'CMIP3': 'solid',
    'OBS': 'solid'
}

def _get_provenance_record(caption):
    """Create a provenance record."""

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'bock_lisa',
        ],
        'references': [
            'acknow_project',
        ],
    }
    return record


def write_data(data, cfg):
    """Write netcdf file."""
    var_attr =  {'short_name': 'ice_fraction',
                  'long_name': 'ice fraction',
                  'units': cf_units.Unit('1')
                 }

    global_attrs = {'project': list(cfg['input_data'].values())[0]['project']}
    caption = 'Ice fraction'
    filename = var_attr['short_name']
    attributes = {}
    caption += '.'
    attributes.update(global_attrs)
    path = get_diagnostic_filename(filename, cfg)

    data_dict = {}
    data_dict["Temperature"] = data[0]
    data_dict["Ice fraction mean"] = data[1]
    data_dict["Ice fraction p5"] = data[2]
    data_dict["Ice fraction p95"] = data[3]

    io.save_scalar_data(data_dict,
                    path,
                    var_attr,
                    attributes=attributes)



def plot_icefrac(ice_frac, legend, cfg):
    """Create diagnostic data and plot it."""

    cube_label = legend
    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.plot(ice_frac[0], ice_frac[1], label=cube_label, color=line_color,
              linestyle=line_dash)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Ice fraction')

    logger.info("Plotting %s", legend)


def plot_errorband(x_axis, y_axis_p5, y_axis_p95, legend, cfg):
    """Create diagnostic data and plot it."""

    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.fill_between(x_axis, y_axis_p5, y_axis_p95, color=line_color,
                     linestyle=line_dash, alpha=.1)

    logger.info("Plotting %s", legend)


def create_bins(lower_bound, width, quantity):
    """Create_bins returns an equal-width (distance) partitioning."""

    bins = []
    for low in range(lower_bound, 
                     lower_bound + quantity*width + 1, width):
        bins.append((low, low+width))
    return bins


def find_bin(value, bins):
    
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1


def calculate_icefrac(input_data, cfg, description=None):
    """Calculate ice fraction."""
    logger.info("Calculating ice fraction.")
    msg = '' if description is None else f' for {description}'
    ancestors = []
    ice_frac = []
    ice_frac_all = []

    # ta bins
    ta_bins = create_bins(lower_bound=230,
                          width=5,
                          quantity=13)
    ta_axis = [np.mean(x) for x in ta_bins]

    print(ta_bins)
    print(ta_axis)

    # Iterate over all datasets and save ice fraction
    for dataset in select_metadata(input_data, short_name='ta'):
        dataset_name = dataset['dataset']
        logger.debug("Calculating ice fraction of dataset '%s'", dataset_name)
        cli_data = select_metadata(input_data,
                                    short_name='cli',
                                    dataset=dataset_name)
        clw_data = select_metadata(input_data,
                                    short_name='clw',
                                    dataset=dataset_name)
        if not cli_data:
            logger.debug(
                "No 'cli' data for '%s' available, skipping ice fracion "
                "calculation for it", dataset_name)
            continue
        elif not clw_data:
            logger.debug(
                "No 'clw' data for '%s' available, skipping ice fracion "
                "calculation for it", dataset_name)
            continue
        ta_cube = dataset['cube']
        cli_cube = cli_data[0]['cube']
        clw_cube = clw_data[0]['cube']
        ancestors.extend(dataset['ancestors'] + cli_data[0]['ancestors']
                         + clw_data[0]['ancestors'])

        cli_cube.data[cli_cube.data < 0.001] = np.nan
        clw_cube.data[clw_cube.data < 0.001] = np.nan

        # Calculate ice fraction
        ice_fraction = (cli_cube / (clw_cube + cli_cube))

        binned_weights = [[] for _ in range(len(ta_bins))]

        longitude = ta_cube.coord('longitude').points
        latitude = ta_cube.coord('latitude').points
        air_pressure = ta_cube.coord('air_pressure').points

        for ipres, pressure in enumerate(air_pressure):
            for ilat, lat in enumerate(latitude):
                for ilon, lon in enumerate(longitude):
                    bin_index = find_bin(ta_cube.data[ipres, ilat, ilon], ta_bins)
                    if (bin_index >= 0 and ice_fraction.data[ipres, ilat, ilon] > 0.):
                        #print(bin_index, ta_cube.data[ipres, ilat, ilon], 
                        #      ice_fraction.data[ipres, ilat, ilon], ipres, ilat, ilon)
                        binned_weights[bin_index].append(ice_fraction.data[ipres, ilat, ilon])

        ice_frac_all.append([np.mean(ibin) for ibin in binned_weights])

        print(ice_frac_all)

    ice_frac.append(ta_axis)
    ice_frac.append(np.average(ice_frac_all, axis=0))
    ice_frac.append(np.quantile(ice_frac_all, 0.05, axis=0))
    ice_frac.append(np.quantile(ice_frac_all, 0.95, axis=0))

    print(ice_frac)

    return ice_frac


def preprocess_data(cfg, year_idx=None):
    """Preprocess data."""
    input_data = deepcopy(list(cfg['input_data'].values()))

    project = input_data[0]['project']
    new_input_data = []
    for (var, var_data) in group_metadata(input_data, 'short_name').items():
        grouped_data = group_metadata(var_data, 'dataset')
        for (dataset_name, datasets) in grouped_data.items():
            logger.debug("Preprocessing '%s' data for dataset '%s'", var,
                         dataset_name)

            cube = iris.load_cube(datasets[0]['filename'])

            logger.debug("Running example computation")
            cube = iris.util.squeeze(cube)

            new_input_data.append({
                **datasets[0],
                'ancestors': [datasets[0]['filename']],
                'cube': cube,
            })

    return new_input_data


def set_default_cfg(cfg):
    """Set default values for cfg."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    return cfg


def main(cfg):
    """Run diagnostic."""
    cfg = set_default_cfg(cfg)
    input_data = preprocess_data(cfg) 

    ice_frac = calculate_icefrac(input_data, cfg)

    plot_icefrac(ice_frac, "ECS_high_hist", cfg)

    plot_errorband(ice_frac[0], ice_frac[2], ice_frac[3], "ECS_high_hist", cfg)

#    all_vars = list(group_metadata(input_data, 'short_name'))
#    vars_groups = group_metadata(input_data, 'short_name')

#    groups = group_metadata(input_data, 'variable_group', sort='dataset')
#    print(groups)

#    cubes = iris.cube.CubeList()

#    plt.figure(figsize=(8, 12))

#    for group_name in groups:
#        logger.info("Processing group %s", group_name)
#
#        for attributes in groups[group_name]:
#            logger.info("Loop dataset %s", attributes['dataset'])
#            if attributes['dataset'] == 'MultiModelMean':
#              logger.info("Processing dataset %s", attributes['dataset'])
#              input_file = attributes['filename']
#
#              cube = compute_diagnostic(input_file)
#              if cube.var_name == 'ta':
#                cube_ta = cube
#
#              cubes.append(cube)
#
#              plot_diagnostic(cube, group_name, plot_type, cfg)
#
#            elif attributes['dataset'] == 'MultiModelP5':
#              logger.info("Processing dataset %s", attributes['dataset'])
#              input_file = attributes['filename']
#              cube_p5 = compute_diagnostic(input_file)
#              cubes.append(cube_p5)
#
#            elif attributes['dataset'] == 'MultiModelP95':
#              logger.info("Processing dataset %s", attributes['dataset'])
#              input_file = attributes['filename']
#              cube_p95 = compute_diagnostic(input_file)
#              cubes.append(cube_p95)
#
#        #if cube_p5 and cube_p95:
#        if group_name != 'OBS':
#          plot_errorband(cube_p5, cube_p95, group_name, plot_type, cfg)
#
#    if plot_type == 'height':
#      plt.ylim(1000.,100.)
#      plt.yscale('log')
#      title = 'Vertical mean of ' + attributes['long_name']
#    elif plot_type == 'zonal':
#      title = 'Zonal mean of ' + attributes['long_name']
#    else:
#      title = attributes['long_name']
#
#    plt.title(title)
#    plt.legend(ncol=1)
#    plt.grid(True)
#
#    for group_name in cfg['group_by']:
#
#        logger.info("Processing group %s", group_name)
#
#        for attributes_1 in groups[group_name[0]]:
#            logger.info("Loop dataset %s", attributes_1['dataset'])
#            if attributes_1['dataset'] == 'MultiModelMean':
#              logger.info("Processing dataset %s", attributes_1['dataset'])
#              input_file_1 = attributes_1['filename']
#
#        for attributes_2 in groups[group_name[1]]:
#            logger.info("Loop dataset %s", attributes_2['dataset'])
#            if attributes_2['dataset'] == 'MultiModelMean':
#              logger.info("Processing dataset %s", attributes_2['dataset'])
#              input_file_2 = attributes_2['filename']
#
#        cube = compute_diff(input_file_1, input_file_2)
#
#        cubes.append(cube)
#
#        plot_diagnostic_diff(cube, group_name[0], plot_type, cfg)
#
#
    caption = ("Ice fraction")

    path = get_diagnostic_filename('ice_fraction', cfg)

    #provenance_record = get_provenance_record(
    #    attributes, ancestor_files=cfg['input_files'])

    # Provenance
    provenance_record = _get_provenance_record(caption)
    provenance_record['ancestors'] = cfg['input_files']
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(path, provenance_record)

    basename = 'ice_fraction'

    # Save the data used for the plot
    #save_data(basename, provenance_record, cfg, ice_frac)
    write_data(ice_frac, cfg)

    title = 'Ice fraction'

    plt.title(title)

    # And save the plot
    save_figure(basename, provenance_record, cfg)



if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
