"""
Diagnostic Maps:

Diagnostic to produce images of a map with coastlines from a cube.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

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

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
import matplotlib
matplotlib.use('Agg')  # noqa
import iris

from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import cartopy

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

# #Blues_r
# 0: (0.03137254901960784, 0.18823529411764706, 0.4196078431372549, 1.0),
# 0.025(0.03137254901960784, 0.21259515570934256, 0.45577854671280277, 1.0),
# 0.05 (0.03137254901960784, 0.23695501730103805, 0.4919492502883507, 1.0),

ice_cmap_dict = {'red': ((0., 0.0313, 0.0313),
                   (0.15, 0.0313, 1.),
                   (1., 1., 1.)),
         'green': ((0., 0.237, 0.237),
                   (0.15, 0.237,1.),
                   (1., 1., 1.)),
         'blue':  ((0., 0.456, 0.456),
                   (0.15, 0.456, 1.),
                   (1., 1., 1.))
        }
ice_cmap = matplotlib.colors.LinearSegmentedColormap('ice_cmap', ice_cmap_dict)#.reversed()


def calculate_area_time_series(cube,):
    """
    Calculate the area of unmasked cube cells.
    Requires a cube with two spacial dimensions. (no depth coordinate).

    Parameters
    ----------
        cube: iris.cube.Cube
            Original data

    Returns
    -------
        iris.cube.Cube
            collapsed cube, in units of m^2
    """
    data = []
    times = np.array(cube.coord('time').points.astype(float))
    for t, time in enumerate(times):
        icedata = cube[t].data
        icedata = np.ma.masked_where(icedata < 0.15,icedata)

        area = iris.analysis.cartography.area_weights(cube[t])
        totalArea = np.ma.masked_where(icedata.mask, area.data).sum()
        print(t,totalArea)
        data.append(totalArea)

    ######
    # Create a small dummy output array
    data = np.array(data)
    return times, data



def make_ts_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple time series plot for an individual model.

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

    # # Load threshold
    # threshold = float(cfg['threshold'])

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        #cube_layer.data = np.ma.masked_where(cube_layer.data <  threshold, cube_layer.data)
        for m,i in metadata.items():
                print(m,i)
        times, data = calculate_area_time_series(cube_layer)
        layer = str(layer)

        plt.plot(times, data,)

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['preprocessor']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', metadata['preprocessor'] + str(layer_index) + image_extention)
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='map_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()

def make_polar_map(
        cube,
        pole = 'North',
        cmap = 'Blues_r',
        zlim = [-0.001,100.001],
):
    """
    Make a polar map plot.

    The cube is the opened cube (two dimensional),
    pole is the polar region (North/South)
    cmap is the colourmap,
    zlim is the z limits (usually 0, 100)
    """
    fig = plt.figure()
    fig.set_size_inches(7, 7)

    if pole not in ['North', 'South']:
        logger.fatal('make_polar_map: hemisphere not provided.')


    if pole == 'North':# North Hemisphere
        ax1 = plt.subplot(111,projection=cartopy.crs.NorthPolarStereo())
        ax1.set_extent([-180, 180, 50, 90], cartopy.crs.PlateCarree())

    if pole == 'South':# South Hemisphere
        ax1 = plt.subplot(111,projection=cartopy.crs.SouthPolarStereo())
        ax1.set_extent([-180, 180, -90, -50], cartopy.crs.PlateCarree())

    qplt.contourf(cube, 20, vmim=zlim[0], vmax=zlim[1], cmap=cmap,linewidth=0, rasterized=True,)
    ax1.add_feature(cartopy.feature.LAND, zorder=10, facecolor = [0.8,0.8,0.8], )
    ax1.gridlines(linewidth=0.5, color='black', zorder=20, alpha=0.5, linestyle='--')#':',c='k',zorder=20,)
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('make_polar_map: Not able to add coastlines')
    return fig, ax1


def get_pole(cube):
    """ Return a hemisphere name as a string (Either North or South)."""
    margin = 5.
    print(np.min(cube.coord('latitude').points),np.max(cube.coord('latitude').points))
    if np.max(cube.coord('latitude').points) < 0. + margin: return 'South'
    if np.min(cube.coord('latitude').points) > 0. - margin: return 'North'
    logger.fatal('get_pole: Not able to determine hemisphere.')


def get_time_string(cube):
    """ Return a climatological season time string in the format "year season"."""
    season = cube.coord('clim_season').points
    year = cube.coord('year').points
    return str(year[0]) + ' ' + season[0].upper()


def make_map_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple map plot for an individual model.

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

    # Load threshold
    threshold = float(cfg['threshold'])

    # Making plots for each layer
    plot_types = ['Fractional cover', 'Ice Extent']
    plot_times = ['first', 'last'] #'decades']
    for plot_type, plot_time in product(plot_types, plot_times):
        for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
            layer = str(layer)

            if plot_type == 'Fractional cover':
                    cmap = 'Blues_r'
            if plot_type == 'Ice Extent':
                cmap = ice_cmap

            if plot_time == 'first':
                cube = cube_layer[0]
            if plot_time == 'last':
                cube = cube_layer[-1]

            # use cube to determine which hemisphere, season and year.
            pole = get_pole(cube)
            time_str = get_time_string(cube)

            # Make the polar map.
            fig, ax1 = make_polar_map(
                    cube,
                    pole = pole,
                    cmap = cmap)

            # Add title to plot
            title = ' '.join([metadata['dataset'], plot_type, time_str])
            if layer:
                title = ' '.join(
                    [title, '(', layer,
                     str(cube_layer.coords('depth')[0].units), ')'])
            plt.title(title)

            # Determine image filename:
            suffix  = '_'.join(['ortho_map', plot_type, time_str, str(layer_index) + image_extention])
            suffix = suffix.replace(' ', '')
            if multi_model:
                path = diagtools.folder(cfg['plot_dir']) + os.path.basename(filename)
                path = path.replace('.nc', suffix)
            else:
                path = diagtools.get_image_path(
                    cfg,
                    metadata,
                    suffix=suffix,
                )

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


        metadatas = diagtools.get_input_files(cfg, index=index)
        for filename in sorted(metadatas.keys()):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )
            ######
            # maps plots of individual models
            make_map_plots(cfg, metadatas[filename], filename)

            ######
            # time series plots o
            #make_ts_plots(cfg, metadatas[filename], filename)


    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
