"""
Diagnostic profile:

Diagnostic to produce images of the time development of a metric from
cubes. These plost show time on the x-axis and cube value (ie temperature) on
the y-axis.

Two types of plots are produced: individual model timeseries plots and
multi model time series plots. The inidivual plots show the results from a
single cube, even if this is a mutli-model mean made by the _multimodel.py
preproccessor. The multi model time series plots show several models
on the same axes, where each model is represented by a different line colour.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, no depth component, and no
latitude or longitude coordinates.

Some approproate preprocessors for a 3D+time field would be:

preprocessors:
  prep_timeseries_1:# For Global Volume Averaged
    average_volume:
      coord1: longitude
      coord2: latitude
      coordz: depth
  prep_timeseries_2: # For Global surface Averaged
    extract_levels:
      levels:  [0., ]
      scheme: linear_extrap
    average_area:
      coord1: longitude
      coord2: latitude


This tool is part of the ocean diagnostic tools package in the ESMValTool.


Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""

import logging
import os
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def timeplot(cube, **kwargs):
    """
    Make a time series plot

    Needed because iris version 1.13 fails due to the time axis.
    """
    if iris.__version__ > '2.0':
        qplt.plot(cube, kwargs)
    else:
        times = diagtools.timecoord_to_float(cube.coord('time'))
        plt.plot(times, cube.data, **kwargs)


def make_time_series_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple plot for an indivudual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)

        if multi_model:
            timeplot(cube_layer, label=metadata['dataset'], ls=':')
        else:
            timeplot(cube_layer, label=metadata['dataset'])

        # Add title, legend to plots
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Determine image filename:
        if multi_model:

            path = diagtools.get_image_path(
                cfg,
                metadata,
                prefix='MultiModel',
                suffix='_'.join(['timeseries', str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='timeseries_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_time_series(
        cfg,
        metadata,
):
    """
    Make a time series plot showing several models.

    This function makes a simple plot for an indivudual model.
    The cfg is the opened global config,
    metadata is the metadata dictionairy.
    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[filename] = cubes
        for layer in cubes:
            layers[layer] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:

        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('viridis')

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'

            if 'MultiModel' in metadata[filename]['dataset']:
                timeplot(
                    model_cubes[filename][layer],
                    c=color,
                    # label=metadata[filename]['dataset'],
                    ls=':',
                    lw=2.,
                )
                plot_details[filename] = {
                    'c': color,
                    'ls': ':',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }
            else:
                timeplot(
                    model_cubes[filename][layer],
                    c=color,
                    # label=metadata[filename]['dataset'])
                    ls='-',
                    lw=2.,
                )
                plot_details[filename] = {
                    'c': color,
                    'ls': '-',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }

            title = metadata[filename]['long_name']
            if layer != '':
                z_units = model_cubes[filename][layer].coords('depth')[0].units

        # Add title, legend to plots
        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels_',
                suffix='_'.join(['timeseries', str(layer) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        # Resize and add legend outside thew axes.
        plt.gcf().set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

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
            metadata_filename
        )

        metadatas = diagtools.get_input_files(cfg, index=index)

        #######
        # Multi model time series
        multi_model_time_series(
            cfg,
            metadatas,
        )

        for filename in sorted(metadatas.keys()):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            make_time_series_plots(cfg, metadatas[filename], filename)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
