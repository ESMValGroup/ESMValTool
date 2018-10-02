"""
Diagnostic Maps quad.

Diagnostic to produce an image showing four maps.
These plost show latitude vs longitude and the cube value is used as the colour
scale.
        model1              model 1 minus model2
        model2 minus obs    model1 minus obs

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be:
preprocessors:
  prep_map:
    extract_levels:
      levels:  [100., ]
      scheme: linear_extrap
    time_average:

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt

import numpy as np

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def match_moddel_to_key(model_type, cfg_dict, input_files_dict, ):
    """
    Match up the three models and observations dataset from the configs.

    This function checks that the control_model, exper_model and
    observational_dataset dictionairies from the recipe are matched with the
    input file dictionairy in the cfg metadata.
    """
    for input_file, intput_dict in input_files_dict.items():
        intersection = dict(intput_dict.items() & cfg_dict.items())
        if intersection == cfg_dict:
            return input_file
    logger.warning("Unable to match model: %s", model_type)
    return ''


def get_cube_range(cubes):
    """Determinue the minimum and maximum values of an array of cubes."""
    mins = []
    maxs = []
    for cube in cubes:
        mins.append(cube.data.min())
        maxs.append(cube.data.max())
    return [np.min(mins), np.max(maxs), ]


def get_cube_range_diff(cubes):
    """Determinue the largest deviation from zero in an array of cubes."""
    ranges = []
    for cube in cubes:
        ranges.append(np.abs(cube.data.min()))
        ranges.append(np.abs(cube.data.max()))
    return [-1. * np.max(ranges), np.max(ranges)]


def add_map_subplot(subplot, cube, n_points=15, title='', cmap=''):
    """Create a map subplot."""
    plt.subplot(subplot)
    qplt.contourf(cube, n_points, linewidth=0, cmap=plt.cm.get_cmap(cmap))
    plt.gca().coastlines()
    plt.title(title)


def multi_model_maps(
        cfg,
        input_files,
):
    """
    Make a simple map plot for an individual model.

    The cfg is the opened global config,
    input_files is the input files dictionairy
    filename is the preprocessing model file.
    """
    filenames = {}
    ctl_key = 'control_model'
    exp_key = 'exper_model'
    obs_key = 'observational_dataset'
    model_types = [ctl_key, exp_key, obs_key]
    for model_type in model_types:
        logger.debug(model_type, cfg[model_type])
        filenames[model_type] = match_moddel_to_key(model_type,
                                                    cfg[model_type],
                                                    input_files)

    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        cube = iris.load_cube(input_file)
        cube = diagtools.bgc_units(cube, input_files[input_file]['short_name'])

        cubes[model_type] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_type]:
            layers[layer] = True

    logger.debug('layers: %s', ', '.join(layers))
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    # ####
    # load names:
    exper = input_files[filenames[exp_key]]['dataset']
    control = input_files[filenames[ctl_key]]['dataset']
    obs = input_files[filenames[obs_key]]['dataset']
    long_name = cubes[exp_key][list(layers.keys())[0]].long_name

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:
        fig = plt.figure()
        fig.set_size_inches(9, 6)

        # Create the cubes
        cube221 = cubes[exp_key][layer]
        cube222 = cubes[exp_key][layer] - cubes[ctl_key][layer]
        cube223 = cubes[ctl_key][layer] - cubes[obs_key][layer]
        cube224 = cubes[exp_key][layer] - cubes[obs_key][layer]

        # create the z axis for plots 2, 3, 4.
        zrange = get_cube_range_diff([cube222, cube223, cube224])
        n_points = 15
        linspace = np.linspace(zrange[0], zrange[1], n_points, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(221, cube221, n_points=n_points, cmap='viridis',
                        title=exper)
        add_map_subplot(222, cube222, n_points=linspace, cmap='bwr',
                        title=' '.join([exper, 'minus', control]))
        add_map_subplot(223, cube223, n_points=linspace, cmap='bwr',
                        title=' '.join([control, 'minus', obs]))
        add_map_subplot(224, cube224, n_points=linspace, cmap='bwr',
                        title=' '.join([exper, 'minus', obs]))

        # Add overall title
        fig.suptitle(long_name, fontsize=14)

        # Determine image filename:
        fn_list = [long_name, exper, control, obs, str(layer)]
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def main(cfg):
    """
    Load the config file, and send it to the plot maker.

    The cfg is the opened global config.
    """
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )
        input_files = diagtools.get_input_files(cfg, index=index)
        # #####
        # Multi model time series
        multi_model_maps(
            cfg,
            input_files,
        )

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
