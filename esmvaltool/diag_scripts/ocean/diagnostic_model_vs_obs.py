"""
Diagnostic Model vs observations maps.

Diagnostic to produce an image showing four maps.
These plost show latitude vs longitude and the cube value is used as the colour
scale.
        model              obs
        model minus obs    model1 over obs

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
        print('\n------------\nmatch_moddel_to_key: intput_dict:', intput_dict)
        print('match_moddel_to_key: cfg_dict:', cfg_dict)
        try:
            intersection = (intput_dict.items() & cfg_dict.items())

        except: intersection = {}
        intersection = {i: j for i,j in intersection}
        print('match_moddel_to_key: intersection:', intersection)
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

def make_model_vs_obs_plots(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
    Make a model vs obs map plot

    The cfg is the opened global config,
    input_files is the input files dictionairy
    filename is the preprocessing model file.
    """
    filenames = {'model' : model_filename, 'obs': obs_filename}
    print('make_model_vs_obs_plots:', filenames)
    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        print('loading:', model_type, input_file)
        cube = iris.load_cube(input_file)
        #cube = diagtools.bgc_units(cube, metadatas[input_file]['short_name'])
        cubes[model_type] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_type]:
            layers[layer] = True

    logger.debug('layers: %s', ', '.join(layers))
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    # ####
    # load names:
    model = metadatas[filenames['model']]['dataset']
    obs = metadatas[filenames['obs']]['dataset']

    long_name = cubes['model'][list(layers.keys())[0]].long_name

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:
        fig = plt.figure()
        fig.set_size_inches(9, 6)

        # Create the cubes
        cube221 = cubes['model'][layer]
        cube222 = cubes['obs'][layer]
        cube223 = cubes['model'][layer] - cubes['obs'][layer]
        cube224 = cubes['model'][layer] / cubes['obs'][layer]

        # create the z axis for plots 2, 3, 4.
        zrange12 = get_cube_range([cube221, cube222])
        zrange3 = get_cube_range_diff([cube223])
        zrange4 = [0.1, 10.]

        n_points = 15
        linspace12 = np.linspace(zrange12[0], zrange12[1], n_points, endpoint=True)
        linspace3 = np.linspace(zrange3[0], zrange3[1], n_points, endpoint=True)
        logspace4 = np.logspace(zrange4[0], zrange4[1], 21, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(221, cube221, n_points=linspace12, cmap='viridis',
                        title=model)
        add_map_subplot(222, cube222, n_points=linspace12, cmap='viridis',
                        title=' '.join([obs, ]))
        add_map_subplot(223, cube223, n_points=linspace3, cmap='bwr',
                        title=' '.join([model, 'minus', obs]))
        add_map_subplot(224, cube224, n_points=logspace4, cmap='bwr',
                        title=' '.join([model, 'over', obs]))

        # Add overall title
        fig.suptitle(long_name, fontsize=14)

        # Determine image filename:
        fn_list = [long_name, model, obs, str(layer)]
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
            'metadata filename:\t%s, %s',
            index,
            metadata_filename,
        )
        metadatas = diagtools.get_input_files(cfg, index=index)

        model_type = 'observational_dataset'
        logger.debug('model_type: %s, %s', index, model_type,)
        logger.debug('cfg[model_type]: %s, %s', index, cfg[model_type])
        logger.debug('metadatas:  %s, %s', index, metadatas,)
        obs_filename = match_moddel_to_key('observational_dataset',
                                           cfg[model_type],
                                           metadatas)
        for filename in sorted(metadatas.keys()):

            if filename == obs_filename:
                continue
            if not os.path.exists(obs_filename):
                continue
            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # model vs obs plots
            make_model_vs_obs_plots(
                cfg,
                metadatas,
                filename,
                obs_filename)


    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
