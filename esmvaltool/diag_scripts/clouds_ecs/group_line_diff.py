"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
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

def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    #print(attributes)
    caption = ("Average {long_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def _get_cube_list(input_files):
    """Get :class:`iris.cube.CubeList` of input files."""
    cubes = iris.cube.CubeList()

    # Input files
    for filename in input_files:
        logger.info("Loading '%s'", filename)
        cube = _load_cube_with_dataset_coord(filename)
        cube.attributes['filename'] = filename
        cubes.append(cube)

    # Check metadata of cubes
    for cube in cubes:
        check_metadata(cube.attributes)

    return cubes



def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def compute_diff(filename1, filename2):
    """Compute difference between two cubes."""
    logger.debug("Loading %s", filename1)
    cube1 = iris.load_cube(filename1)
    cube2 = iris.load_cube(filename2)

    cube = cube2 - cube1
    cube.metadata = cube1.metadata
    return cube


def plot_diagnostic(cube, legend, plot_type, cfg):
    """Create diagnostic data and plot it."""

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
    else:
        cube_label = legend
        line_color = LINE_COLOR.get(legend, legend)
        line_dash = LINE_DASH.get(legend, legend)

        plt.subplot(211)

        if cube.var_name == 'pr':
          cube.units = cube.units / 'kg m-3'
          cube.data = cube.core_data() / 1000.
          cube.convert_units('mm day-1')
        elif cube.var_name == 'cli':
            cube.convert_units('g/kg')
        elif cube.var_name == 'clw':
            cube.convert_units('g/kg')

        if plot_type == 'height':
          cube.coord('air_pressure').convert_units('hPa')
          y_axis = cube.coord('air_pressure')
          qplt.plot(cube, y_axis, label=cube_label, color=line_color, 
                    linestyle=line_dash)
        else:
          lat = cube.coord('latitude')
          qplt.plot(lat, cube, label=cube_label, color=line_color,
                    linestyle=line_dash)

        logger.info("Plotting %s", legend)


def plot_diagnostic_diff(cube, legend, plot_type, cfg):
    """Create diagnostic data and plot it."""

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
    else:
        cube_label = LINE_LEGEND.get(legend, legend)
        line_color = LINE_COLOR.get(legend, legend)
        line_dash = LINE_DASH.get(legend, legend)

        plt.subplot(212)

        if cube.var_name == 'pr':
          cube.units = cube.units / 'kg m-3'
          cube.data = cube.core_data() / 1000.
          cube.convert_units('mm day-1')
        elif cube.var_name == 'cli':
            cube.convert_units('g/kg')
        elif cube.var_name == 'clw':
            cube.convert_units('g/kg')

        if plot_type == 'height':
          cube.coord('air_pressure').convert_units('hPa')
          y_axis = cube.coord('air_pressure')
          qplt.plot(cube, y_axis, label=cube_label, color=line_color, 
                    linestyle=line_dash)
        else:
          lat = cube.coord('latitude')
          qplt.plot(lat, cube, label=cube_label, color=line_color,
                    linestyle=line_dash)

        logger.info("Plotting %s", legend)


def plot_errorband(cube1, cube2, legend, plot_type, cfg):
    """Create diagnostic data and plot it."""

    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.subplot(211)

    if cube1.var_name == 'pr':
      cube1.units = cube1.units / 'kg m-3'
      cube1.data = cube1.core_data() / 1000.
      cube1.convert_units('mm day-1')
      cube2.units = cube2.units / 'kg m-3'
      cube2.data = cube2.core_data() / 1000.
      cube2.convert_units('mm day-1')
    elif cube1.var_name == 'cli':
        cube1.convert_units('g/kg')
        cube2.convert_units('g/kg')
    elif cube1.var_name == 'clw':
        cube1.convert_units('g/kg')
        cube2.convert_units('g/kg')

    if plot_type == 'height':
      cube1.coord('air_pressure').convert_units('hPa')
      cube2.coord('air_pressure').convert_units('hPa')
      y_axis = cube1.coord('air_pressure').points
      plt.fill_betweenx(y_axis, cube1.data, cube2.data, color=line_color,
                        linestyle=line_dash, alpha=.1)
    else:
      lat = cube1.coord('latitude').points
      plt.fill_between(lat, cube1.data, cube2.data, color=line_color,
                       linestyle=line_dash, alpha=.1)
    logger.info("Plotting %s", legend)


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    cfg.setdefault('filename_attach', 'base')
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    plot_type = cfg['plot_type']

    input_data = list(cfg['input_data'].values())
    all_vars = list(group_metadata(input_data, 'short_name'))

    groups = group_metadata(input_data, 'variable_group', sort='dataset')

    cubes = iris.cube.CubeList()

    plt.figure(figsize=(8, 12))

    for group_name in groups:
        logger.info("Processing variable %s", group_name)

        for attributes in groups[group_name]:
            logger.info("Loop dataset %s", attributes['dataset'])
            if attributes['dataset'] == 'MultiModelMean' or group_name == 'OBS':
              logger.info("Processing dataset %s", attributes['dataset'])
              input_file = attributes['filename']
              cube = compute_diagnostic(input_file)
              cubes.append(cube)

              plot_diagnostic(cube, group_name, plot_type, cfg)

            elif attributes['dataset'] == 'MultiModelP5':
              logger.info("Processing dataset %s", attributes['dataset'])
              input_file = attributes['filename']
              cube_p5 = compute_diagnostic(input_file)
              cubes.append(cube_p5)

            elif attributes['dataset'] == 'MultiModelP95':
              logger.info("Processing dataset %s", attributes['dataset'])
              input_file = attributes['filename']
              cube_p95 = compute_diagnostic(input_file)
              cubes.append(cube_p95)

        #if cube_p5 and cube_p95:
        if group_name != 'OBS':
          plot_errorband(cube_p5, cube_p95, group_name, plot_type, cfg)

    if plot_type == 'height':
      plt.ylim(1000.,100.)
      plt.yscale('log')
      plt.yticks([1000., 800., 600., 400., 300., 200., 100.], [1000, 800, 600, 400, 300, 200, 100])
      title = 'Vertical mean of ' + attributes['long_name']
    elif plot_type == 'zonal':
      title = 'Zonal mean of ' + attributes['long_name']
    else:
      title = attributes['long_name']

    plt.title(title)
    plt.legend(ncol=1)
    plt.grid(True)

    for group_name in cfg['group_by']:

        logger.info("Processing variable %s", group_name)

        for attributes_1 in groups[group_name[0]]:
            logger.info("Loop dataset %s", attributes_1['dataset'])
            if attributes_1['dataset'] == 'MultiModelMean':
              logger.info("Processing dataset %s", attributes_1['dataset'])
              input_file_1 = attributes_1['filename']

        for attributes_2 in groups[group_name[1]]:
            logger.info("Loop dataset %s", attributes_2['dataset'])
            if attributes_2['dataset'] == 'MultiModelMean':
              logger.info("Processing dataset %s", attributes_2['dataset'])
              input_file_2 = attributes_2['filename']

        cube = compute_diff(input_file_1, input_file_2)

        cubes.append(cube)

        plot_diagnostic_diff(cube, group_name[0], plot_type, cfg)


    if plot_type == 'height':
      plt.ylim(1000.,100.)
      plt.yscale('log')
      plt.yticks([1000., 800., 600., 400., 300., 200., 100.], [1000, 800, 600, 400, 300, 200, 100])
      plt.axvline(x=0, ymin=0., ymax=1., color='black', linewidth=3)
      title = 'Difference of vertical mean of ' + attributes['long_name']
    elif plot_type == 'zonal':
      plt.axhline(y=0, xmin=-90., xmax=90., color='black', linewidth=3)
      title = 'Difference of zonal mean of ' + attributes['long_name']
    else:
      title = attributes['long_name']

    plt.title(title)
    plt.legend(ncol=1)
    plt.grid(True)

    provenance_record = get_provenance_record(
        attributes, ancestor_files=cfg['input_files'])

    if plot_type == 'height':
      basename = 'level_diff_' + attributes['short_name'] + '_' + cfg['filename_attach']
    else:
      basename = 'zonal_diff_' + attributes['short_name'] + '_' + cfg['filename_attach']

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cubes)

    # And save the plot
    save_figure(basename, provenance_record, cfg)



if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
