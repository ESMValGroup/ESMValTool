"""
Model 1 vs Model 2 vs Observations diagnostics.
===============================================

Diagnostic to produce an image showing four maps, based on a comparison of two
differnt models results against an observational dataset. This process is
often used to compare a new iteration of a model under development against
a previous version of the same model. The four map plots are:

* Top left: model 1
* Top right: model 1 minus model 2
* Bottom left: model 2 minus obs
* Bottom right: model 1 minus obs

All four plots show latitude vs longitude and the cube value is used as the
colour scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_map:
      extract_levels:
        levels:  [100., ]
        scheme: linear_extrap
      climate_statistics:
        operator: mean

This diagnostic also requires the ``exper_model``, ``exper_model`` and
``observational_dataset`` keys in the recipe::

  diagnostics:
     diag_name:
       ...
       scripts:
         Global_Ocean_map:
           script: ocean/diagnostic_maps_quad.py
           exper_model:  {Model 1 dataset details}
           control_model: {Model 2 dataset details}
           observational_dataset: {Observational dataset details}

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk

"""
import logging
import os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def add_map_subplot(subplot, cube, nspace, title='', cmap=''):
    """
    Add a map subplot to the current pyplot figure.

    Parameters
    ----------
    subplot: int
        The matplotlib.pyplot subplot number. (ie 221)
    cube: iris.cube.Cube
        the iris cube to be plotted.
    nspace: numpy.array
        An array of the ticks of the colour part.
    title: str
        A string to set as the subplot title.
    cmap: str
        A string to describe the matplotlib colour map.

    """
    plt.subplot(subplot)
    qplot = qplt.contourf(cube, nspace, linewidth=0,
                          cmap=plt.cm.get_cmap(cmap))
    qplot.colorbar.set_ticks([nspace.min(),
                              (nspace.max() + nspace.min()) / 2.,
                              nspace.max()])

    plt.gca().coastlines()
    plt.title(title)


def multi_model_maps(
        cfg,
        input_files,
):
    """
    Make the four pane model vs model vs obs comparison plot.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    input_files: dict
        the metadata dictionairy

    """
    filenames = {}
    ctl_key = 'control_model'
    exp_key = 'exper_model'
    obs_key = 'observational_dataset'
    model_types = [ctl_key, exp_key, obs_key]
    for model_type in model_types:
        logger.debug(model_type, cfg[model_type])
        filenames[model_type] = diagtools.match_model_to_key(model_type,
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
        zrange1 = diagtools.get_cube_range([cube221, ])
        zrange2 = diagtools.get_cube_range_diff([cube222, cube223, cube224])

        linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
        linspace2 = np.linspace(zrange2[0], zrange2[1], 12, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(221, cube221, linspace1, cmap='viridis', title=exper)
        add_map_subplot(222, cube222, linspace2, cmap='bwr',
                        title=' '.join([exper, 'minus', control]))
        add_map_subplot(223, cube223, linspace2, cmap='bwr',
                        title=' '.join([control, 'minus', obs]))
        add_map_subplot(224, cube224, linspace2, cmap='bwr',
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

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

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
