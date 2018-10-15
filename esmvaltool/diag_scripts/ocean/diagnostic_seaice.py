"""
Diagnostic Maps.

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
import numpy as np
from itertools import product

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt

import cartopy

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def create_ice_cmap(threshold):
    """Create colour map with ocean blue below 15% and white above 15%."""
    threshold = threshold / 100.
    ice_cmap_dict = {'red': ((0., 0.0313, 0.0313),
                             (threshold, 0.0313, 1.),
                             (1., 1., 1.)),
                     'green': ((0., 0.237, 0.237),
                               (threshold, 0.237, 1.),
                               (1., 1., 1.)),
                     'blue':  ((0., 0.456, 0.456),
                               (threshold, 0.456, 1.),
                               (1., 1., 1.))}

    return matplotlib.colors.LinearSegmentedColormap('ice_cmap', ice_cmap_dict)


def calculate_area_time_series(cube, plot_type, threshold):
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
    times = diagtools.timecoord_to_float(cube.coord('time'))
    for time_itr, time in enumerate(times):
        icedata = cube[time_itr].data

        area = iris.analysis.cartography.area_weights(cube[time_itr])
        if plot_type.lower() == 'ice extent':
            # Ice extend is the area with more than 15% ice cover.
            icedata = np.ma.masked_where(icedata < threshold, icedata)
            total_area = np.ma.masked_where(icedata.mask, area.data).sum()
        if plot_type.lower() == 'ice area':
            # Ice area is cover * cell area
            total_area = np.sum(icedata * area)

        logger.debug('Calculating time series area: %s, %s, %s,',
                     time_itr, time, total_area)
        data.append(total_area)

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
    Make a ice extent time series plot for an individual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])
    cube = agregate_by_season(cube)

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # # Load threshold, pole, season.
    threshold = float(cfg['threshold'])
    pole = get_pole(cube)
    season = get_season(cube)

    # Making plots for each layer
    for plot_type in ['Ice Extent', 'Ice Area']:
        for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
            layer = str(layer)

            times, data = calculate_area_time_series(cube_layer,
                                                     plot_type,
                                                     threshold)

            plt.plot(times, data)

            # Add title to plot
            title = ' '.join([metadata['dataset'], pole, 'hemisphere',
                              season, plot_type])
            if layer:
                title = ' '.join(
                    [title, '(', layer,
                     str(cube_layer.coords('depth')[0].units), ')'])
            plt.title(title)

            # y axis label:
            plt.ylabel(' '.join([plot_type, 'm^2']))

            # Determine image filename:
            suffix = '_'.join(['ts', metadata['preprocessor'], season, pole,
                               plot_type, str(layer_index)])\
                     + image_extention
            suffix = suffix.replace(' ', '')
            if multi_model:
                path = diagtools.folder(
                    cfg['plot_dir']) + os.path.basename(filename)
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


def make_polar_map(
        cube,
        pole='North',
        cmap='Blues_r',
):
    """
    Make a polar map plot.

    The cube is the opened cube (two dimensional),
    pole is the polar region (North/South)
    cmap is the colourmap,
    """
    fig = plt.figure()
    fig.set_size_inches(7, 7)

    # ####
    # Set  limits, based on https://nedbatchelder.com/blog/200806/pylint.html

    if pole not in ['North', 'South']:
        logger.fatal('make_polar_map: hemisphere not provided.')

    if pole == 'North':  # North Hemisphere
        ax1 = plt.subplot(111, projection=cartopy.crs.NorthPolarStereo())
        ax1.set_extent([-180, 180, 50, 90], cartopy.crs.PlateCarree())

    if pole == 'South':  # South Hemisphere
        ax1 = plt.subplot(111, projection=cartopy.crs.SouthPolarStereo())
        ax1.set_extent([-180, 180, -90, -50], cartopy.crs.PlateCarree())

    linrange = np.linspace(0., 100., 21.)
    qplt.contourf(cube,
                  linrange,
                  cmap=cmap,
                  linewidth=0,
                  rasterized=True)
    plt.tight_layout()

    ax1.add_feature(cartopy.feature.LAND,
                    zorder=10,
                    facecolor=[0.8, 0.8, 0.8], )

    ax1.gridlines(linewidth=0.5,
                  color='black',
                  zorder=20,
                  alpha=0.5,
                  linestyle='--')
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('make_polar_map: Not able to add coastlines')
    return fig, ax1


def get_pole(cube):
    """Return a hemisphere name as a string (Either North or South)."""
    margin = 5.
    if np.max(cube.coord('latitude').points) < 0. + margin:
        return 'South'
    if np.min(cube.coord('latitude').points) > 0. - margin:
        return 'North'
    logger.fatal('get_pole: Not able to determine hemisphere.')
    return False


def get_time_string(cube):
    """Return a climatological season string in the format: "year season"."""
    season = cube.coord('clim_season').points
    year = cube.coord('year').points
    return str(int(year[0])) + ' ' + season[0].upper()


def get_year(cube):
    """Return the cube year as a string."""
    year = cube.coord('year').points
    return str(int(year))


def get_season(cube):
    """Return a climatological season time string."""
    season = cube.coord('clim_season').points
    return season[0].upper()


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
    cube = agregate_by_season(cube)

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention and threshold.
    image_extention = diagtools.get_image_format(cfg)
    threshold = float(cfg['threshold'])

    # Making plots for each layer
    plot_types = ['Fractional cover', 'Ice Extent']
    plot_times = [0, -1]
    for plot_type, plot_time in product(plot_types, plot_times):
        for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
            layer = str(layer)

            if plot_type == 'Fractional cover':
                cmap = 'Blues_r'
            if plot_type == 'Ice Extent':
                cmap = create_ice_cmap(threshold)

            cube = cube_layer[plot_time]

            # use cube to determine which hemisphere, season and year.
            pole = get_pole(cube)
            time_str = get_time_string(cube)

            # Make the polar map.
            fig, ax1 = make_polar_map(cube,
                                      pole=pole,
                                      cmap=cmap)

            # Add title to plot
            title = ' '.join([metadata['dataset'], plot_type, time_str])
            if layer:
                title = ' '.join([title, '(', layer,
                                  str(cube_layer.coords('depth')[0].units),
                                  ')'])
            plt.title(title)

            # Determine image filename:
            suffix = '_'.join(['ortho_map', plot_type, time_str,
                               str(layer_index)])
            suffix = suffix.replace(' ', '') + image_extention
            if multi_model:
                path = diagtools.folder(cfg['plot_dir'])
                path = path + os.path.basename(filename)
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


def agregate_by_season(cube):
    """
    Aggregate the cube into seasonal means.

    Note that it is not currently possible to do this in the preprocessor,
    as the seasonal mean changes the cube units.
    """
    if not cube.coords('clim_season'):
        iris.coord_categorisation.add_season(cube,
                                             'time',
                                             name='clim_season')
    if not cube.coords('season_year'):
        iris.coord_categorisation.add_season_year(cube,
                                                  'time',
                                                  name='season_year')
    return cube.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN)


def make_map_extent_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make an extent map plot showing several times for an individual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])
    cube = agregate_by_season(cube)

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold, pole and season
    threshold = float(cfg['threshold'])
    pole = get_pole(cube)
    season = get_season(cube)

    # Start making figure
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):

        fig = plt.figure()
        fig.set_size_inches(7, 7)

        if pole == 'North':  # North Hemisphere
            projection = cartopy.crs.NorthPolarStereo()
            ax1 = plt.subplot(111, projection=projection)
            ax1.set_extent([-180, 180, 50, 90], cartopy.crs.PlateCarree())

        if pole == 'South':  # South Hemisphere
            projection = cartopy.crs.SouthPolarStereo()
            ax1 = plt.subplot(111, projection=projection)
            ax1.set_extent([-180, 180, -90, -50], cartopy.crs.PlateCarree())

        ax1.add_feature(cartopy.feature.LAND,
                        zorder=10,
                        facecolor=[0.8, 0.8, 0.8])

        ax1.gridlines(linewidth=0.5,
                      color='black',
                      zorder=20,
                      alpha=0.5,
                      linestyle='--')

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('make_polar_map: Not able to add coastlines')

        times = np.array(cube.coord('time').points.astype(float))
        plot_desc = {}
        for time_itr, time in enumerate(times):
            cube = cube_layer[time_itr]
            line_width = 1
            color = plt.cm.jet(float(time_itr) / float(len(times)))
            label = get_year(cube)
            plot_desc[time] = {'label': label,
                               'c': [color, ],
                               'lw': [line_width, ],
                               'ls': ['-', ]}

            layer = str(layer)
            qplt.contour(cube,
                         [threshold, ],
                         colors=plot_desc[time]['c'],
                         linewidths=plot_desc[time]['lw'],
                         linestyles=plot_desc[time]['ls'],
                         rasterized=True)

        # Add legend
        legend_size = len(plot_desc.keys()) + 1
        ncols = int(legend_size / 25) + 1
        ax1.set_position([ax1.get_position().x0,
                          ax1.get_position().y0,
                          ax1.get_position().width * (1. - 0.1 * ncols),
                          ax1.get_position().height])

        fig.set_size_inches(7 + ncols * 1.2, 7)

        # Construct dummy plots.
        for i in sorted(plot_desc.keys()):
            plt.plot([], [],
                     c=plot_desc[i]['c'][0],
                     lw=plot_desc[i]['lw'][0],
                     ls=plot_desc[i]['ls'][0],
                     label=plot_desc[i]['label'],)

        legd = ax1.legend(loc='center left',
                          ncol=ncols,
                          prop={'size': 10},
                          bbox_to_anchor=(1., 0.5))
        legd.draw_frame(False)
        legd.get_frame().set_alpha(0.)

        # Add title to plot
        title = ' '.join([metadata['dataset'], ])
        if layer:
            title = ' '.join([title, '(', layer,
                              str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)

        # Determine image filename:
        suffix = '_'.join(['ortho_map', pole, season, str(layer_index)])
        suffix = suffix.replace(' ', '') + image_extention
        if multi_model:
            path = diagtools.folder(cfg['plot_dir'])
            path = path + os.path.basename(filename)
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
            # extent maps plots of individual models
            make_map_extent_plots(cfg, metadatas[filename], filename)

            ######
            # maps plots of individual models
            make_map_plots(cfg, metadatas[filename], filename)

            ######
            # time series plots o
            make_ts_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
