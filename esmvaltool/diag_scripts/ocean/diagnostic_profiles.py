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
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        obs_key = 'observational_dataset'
        obs_filename = ''
        obs_metadata = {}

        if obs_key in cfg:
            obs_filename = diagtools.match_model_to_key(obs_key,
                                                        cfg[obs_key],
                                                        metadatas)
            if obs_filename:
                obs_metadata = metadatas[obs_filename]
            else:
                obs_metadata = ''
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
