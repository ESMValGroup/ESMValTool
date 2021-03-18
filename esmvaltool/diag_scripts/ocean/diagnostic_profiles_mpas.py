"""
Profile diagnostics.
====================

Diagnostic to produce figure of the profile over time from a cube.
These plost show cube value (ie temperature) on the x-axis, and depth/height
on the y axis. The colour scale is the time series.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, and depth component, but no
latitude or longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_profile:
      extract_volume:
        long1: 0.
        long2:  20.
        lat1:  -30.
        lat2:  30.
        z_min: 0.
        z_max: 3000.
      area_statistics:
        operator: mean


In order to add an observational dataset to the profile plot, the following
arguments are needed in the diagnostic script::

  diagnostics:
    diagnostic_name:
      variables:
        ...
      additional_datasets:
      - {observational dataset description}
      scripts:
        script_name:
          script: ocean/diagnostic_profiles.py
          observational_dataset: {observational dataset description}

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys

import numpy as np
import iris
import iris.coord_categorisation
import iris.exceptions
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvalcore.preprocessor._time import extract_time
from esmvalcore.preprocessor._regrid import extract_levels

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def determine_profiles_str(cube):
    """
    Determine a string from the cube, to describe the profile.

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    str
        Returns a string which describes the profile.
    """
    options = ['latitude', 'longitude']
    for option in options:
        coord = cube.coord(option)
        if len(coord.points) > 1:
            continue
        value = coord.points.mean()
        if option == 'latitude':
            return str(value) + ' N'
        if option == 'longitude':
            if value > 180.:
                return str(value - 360.) + ' W'
            return str(value) + ' E'
    return ''


def make_profiles_plots(
        cfg,
        metadata,
        filename,
        obs_metadata={},
        obs_filename='',
):
    """
    Make a profile plot for an individual model.

    The optional observational dataset can also be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    try:
        raw_times = diagtools.cube_time_to_float(cube)
    except iris.exceptions.CoordinateNotFoundError:
        return

    # Make annual or Decadal means from:
    if np.max(raw_times) - np.min(raw_times) < 20:
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
    else:
        cube = diagtools.decadal_average(cube)

    times_float = diagtools.cube_time_to_float(cube)
    time_0 = times_float[0]

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    cmap = plt.cm.get_cmap('jet')

    plot_details = {}
    for time_index, time in enumerate(times_float):
        if times_float[-1] == time_0:
            color = 'black'
        else:
            color = cmap((time - time_0) / (times_float[-1] - time_0))

        qplt.plot(cube[time_index, :], cube[time_index, :].coord('depth'),
                  c=color)

        plot_details[str(time_index)] = {'c': color, 'ls': '-', 'lw': 1,
                                         'label': str(int(time))}

    # Add observational data.
    if obs_filename:
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, metadata['short_name'])
        obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)

        obs_key = obs_metadata['dataset']
        qplt.plot(obs_cube, obs_cube.coord('depth'), c='black')

        plot_details[obs_key] = {'c': 'black', 'ls': '-', 'lw': 1,
                                 'label': obs_key}

    # Add title to plot
    title = ' '.join([
        metadata['dataset'],
        metadata['long_name'],
    ])
    plt.title(title)

    # Add Legend outside right.
    diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename).replace(
                '.nc', '_profile' + image_extention)
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='profile' + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def make_multi_model_profiles_plots(
        cfg,
        metadatas,
        short_name,
        obs_metadata={},
        obs_filename='',
        time_range = [2040., 2050.],
        figure_style = 'difference',
        fig = None,
        ax = None,

    ):
    """
    Make a profile plot for an individual model.

    The optional observational dataset can also be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
    # Load cube and set up units
    if fig is None:
        single_pane = True
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        single_pane = False

    cubes = {}
    for filename, metadata in metadatas.items():
        if short_name != metadata['short_name']:
            continue
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata['short_name'])
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')

        raw_times = diagtools.cube_time_to_float(cube)

        times_float = diagtools.cube_time_to_float(cube)
        dataset = metadata['dataset']
        scenario = metadata['exp']
        ensemble = metadata['ensemble']

        if scenario == 'historical':
            cube = extract_time(cube, 1950, 1, 1, 2000, 1, 1)
        else:
            cube = extract_time(cube, time_range[0], 1, 1, time_range[1], 1, 1)

        cube_mean = cube.copy().collapsed('time', iris.analysis.MEAN)
        cube_min = cube.copy().collapsed('time', iris.analysis.MIN)
        cube_max = cube.copy().collapsed('time', iris.analysis.MAX)

        cubes[(dataset, scenario, ensemble, 'mean')]  = cube_mean
        cubes[(dataset, scenario, ensemble, 'min')]  = cube_min
        cubes[(dataset, scenario, ensemble, 'max')]  = cube_max

    for (dataset, scenario, ensemble, metric), cube_mean in cubes.items():
        if metric !=  'mean': continue

        cube_min = cubes[(dataset, scenario, ensemble, 'min')]
        cube_max = cubes[(dataset, scenario, ensemble, 'max')]
        cube_hist = cubes[(dataset, 'historical', ensemble, 'mean')]

        color = diagtools.ipcc_colours[scenario]
        plot_details = {}


        if figure_style == 'difference':
            cube_mean = extract_levels(cube_mean, cube_hist.coord('depth').points, "nearest_horizontal_extrapolate_vertical")
            cube_min = extract_levels(cube_min, cube_hist.coord('depth').points, "nearest_horizontal_extrapolate_vertical")
            cube_max = extract_levels(cube_max, cube_hist.coord('depth').points, "nearest_horizontal_extrapolate_vertical")

            plt.plot(cube_mean.data - cube_hist.data, -1.*cube_hist.coord('depth').points,
                lw=2,
                c=color,
                label= scenario)

            plt.fill_betweenx(-1.*cube_hist.coord('depth').points,
                cube_min.data - cube_hist.data,
                x2=cube_max.data - cube_hist.data,
                alpha = 0.2,
                color=color)

        if figure_style == 'compare':
            cube_min = extract_levels(cube_min, cube_mean.coord('depth').points, "nearest_horizontal_extrapolate_vertical")
            cube_max = extract_levels(cube_max, cube_mean.coord('depth').points, "nearest_horizontal_extrapolate_vertical")
            plt.plot(cube_mean.data, -1.*cube_mean.coord('depth').points,
                lw=2,
                c=color, 
                label=scenario)

            plt.fill_betweenx(-1.*cube_mean.coord('depth').points,
                cube_min.data,
                x2=cube_max.data,
                alpha = 0.2,
                color=color)

        plot_details[scenario] = {'c': color, 'ls': '-', 'lw': 1,
                                             'label': scenario}

    # Add observational data.
    if obs_filename:
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, metadata['short_name'])
        obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)

        obs_key = obs_metadata['dataset']
        qplt.plot(obs_cube, obs_cube.coord('depth'), c='black')

        plot_details[obs_key] = {'c': 'black', 'ls': '-', 'lw': 1,
                                 'label': obs_key}

    time_str = '-'.join([str(t) for t in time_range])
    # Add title to plot
    title = ' '.join([
        short_name,
        figure_style,
        time_str
        ])

    plt.title(title)

    # Add Legend outside right.
    plt.legend()
    #diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    path = diagtools.get_image_path(
            cfg,
            metadata,
            prefix='_'.join(['multi_model', short_name, figure_style, str(time_str)]),
            suffix='profile' + image_extention,
        )

    if not single_pane:
        return fig, ax
    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def main(cfg):
    """
    Run the diagnostics profile tool.

    Load the config file, find an observational dataset filename,
    pass loaded into the plot making tool.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    metadatas = diagtools.get_input_files(cfg)
    obs_filename = ''
    obs_key = 'observational_dataset'
    obs_metadata = {}

    short_names = {}
    for fn, metadata in metadatas.items():
        short_names[metadata['short_name']] = True

    if obs_key in cfg:
        obs_filename = diagtools.match_model_to_key(obs_key,
                                                    cfg[obs_key],
                                                    metadatas)
        if obs_filename:
            obs_metadata = metadatas[obs_filename]
        else:
            obs_metadata = ''

    time_ranges = [[2040, 2050], [2090, 2100], [2050, 2100]]
    for time_range in time_ranges:
        for short_name in short_names.keys():
            for figure_style in ['compare', 'difference']:
                make_multi_model_profiles_plots(
                    cfg,
                    metadatas,
                    short_name,
                    figure_style=figure_style,
                    obs_metadata=obs_metadata,
                    obs_filename=obs_filename,
                    time_range=time_range,
                )
    return

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        for filename in sorted(metadatas.keys()):

            if filename == obs_filename:
                continue

            if metadatas[filename]['frequency'] == 'fx':
                continue

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            make_profiles_plots(cfg, metadatas[filename], filename,
                                obs_metadata=obs_metadata,
                                obs_filename=obs_filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
