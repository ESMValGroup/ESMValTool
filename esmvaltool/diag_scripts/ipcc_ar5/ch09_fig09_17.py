"""Time series diagnostics for ocean heat content.

Time series diagnostics produce figures of the time development of a
field from cubes. These plost show time on the x-axis and cube value
(ie temperature) on the y-axis.

Two types of plots are produced: individual model timeseries plots and
multi model time series plots. The inidivual plots show the results from a
single cube, even if this is a mutli-model mean made by the _multimodel.py
preproccessor. The multi model time series plots show several models
on the same axes, where each model is represented by a different line colour.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has a time component, no depth component, and no
latitude or longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_timeseries_1:# For Global Volume Averaged
      average_volume:
        coord1: longitude
        coord2: latitude
        coordz: depth

An approproate preprocessor for a 3D+time field at the surface would be::

    prep_timeseries_2: # For Global surface Averaged
      extract_levels:
        levels:  [0., ]
        scheme: linear_extrap
      average_area:
        coord1: longitude
        coord2: latitude

An approproate preprocessor for a 2D+time field would be::

    prep_timeseries_2: # For Global surface Averaged
      average_area:
        coord1: longitude
        coord2: latitude

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube

    """
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return

    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)


def make_time_series_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple time series plot for an indivudual model 1D cube.

    This tool loads the cube from the file, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.

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
        if layer != '':
            if cube_layer.coords('depth'):
                z_units = cube_layer.coord('depth').units
            else:
                z_units = ''
            title = ' '.join([title, '(', layer, str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')
        plt.ylabel(str(cube_layer.units))

        # Determine image filename:
        if multi_model:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                prefix='MultiModel',
                suffix='_'.join(['timeseries',
                                 str(layer) + image_extention]),
                metadata_id_list=[
                    'short_name', 'preprocessor', 'diagnostic', 'start_year',
                    'end_year'
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
        cubedic,
):
    """
    Make a time series plot showing several preprocesssed datasets.

    This tool loads several cubes from the files, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make the plot
    plot_details = {}
    cmap = plt.cm.get_cmap('viridis')

    # Plot each file in the group
    for index, filename in enumerate(sorted(cubedic)):
        if len(metadata) > 1:
            color = cmap(index / (len(metadata) - 1.))
        else:
            color = 'blue'

        if index == 0:
            filename0 = filename

        cube = cubedic[filename]

        if 'MultiModel' in metadata[filename]['dataset']:
            timeplot(
                cube,
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
                cube,
                c=color,
                # label=metadata[filename]['dataset'])
                ls='-',
                lw=2.,
            )
            plot_details[filename] = {
                'c': color,
                'ls': '-',
                'lw': 2.,
                'label': metadata[filename]['dataset'],
            }

    # Saving files:
    if cfg['write_plots']:
        path = diagtools.get_image_path(
            cfg,
            metadata[filename0],
            prefix='MultipleModels_',
            suffix='_'.join(['timeseries_0', image_extention]),
            metadata_id_list=[
                'short_name', 'preprocessor', 'diagnostic', 'start_year',
                'end_year'
            ],
        )

    # Resize and add legend outside the axes.
    plt.gcf().set_size_inches(9., 6.)
    diagtools.add_legend_outside_right(
        plot_details, plt.gca(), column_width=0.15)
    plt.ylabel(metadata[filename0]['short_name'])
    plt.title(metadata[filename0]['standard_name'])
    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()

    # Save each netcdf file in the group
    if cfg['write_netcdf']:
        metadata_id_list = [
            'project',
            'dataset',
            'mip',
            'exp',
            'ensemble',
            'short_name',
            'diagnostic',
            'start_year',
            'end_year',
        ]
        workdir = diagtools.folder(cfg['work_dir'])
        for index, filename in enumerate(sorted(cubedic)):
            cube = cubedic[filename]
            path_nc = workdir
            path_nc += '_'.join(
                [str(metadata[filename][b]) for b in metadata_id_list])
            path_nc = path_nc + '.nc'
            logger.info('path_nc = %s', path_nc)
            iris.save(cube, path_nc)


def global_sum(metadata, ):
    """
    Sum up the global heat content file to a 1-D time series.

    Parameters
    ----------
    metadata: dict
        The metadata dictioniary for a specific model.

    """
    cubedic = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])
        cube = cube.collapsed(['longitude', 'latitude', 'depth'],
                              iris.analysis.SUM)
        cube.remove_coord('day_of_month')
        cube.remove_coord('day_of_year')
        cube.remove_coord('month_number')
        #######
        # subtract reference year value 1971
        refvalue = cube.extract(iris.Constraint(year=1971))
        cube -= refvalue
        cubedic[filename] = cube
    return cubedic


def main(cfg):
    """
    Main routine of this diagnostic.

    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata-filename:')
        print(metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)
        #######
        # Sum up the cube over global domain
        cubedic = global_sum(metadatas)
    #######
    # Multi model time series
    multi_model_time_series(
        cfg,
        metadatas,
        cubedic,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
