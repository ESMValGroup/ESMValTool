"""
Sea Ice Diagnostics.
====================

Diagnostic to produce a series of images which are useful for evaluating
the behaviour of the a sea ice model.

There are three kinds of plots shown here.
1. Sea ice Extent maps plots with a stereoscoic projection.
2. Maps plots of individual models ice fracrtion.
3. Time series plots for the total ice extent.

All three kinds of plots are made for both Summer and Winter in both the
North and Southern hemisphere.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

This diagnostic takes data from either North or South hemisphere, and
from either December-January-February or June-July-August. This diagnostic
requires the data to be 2D+time, and typically expects the data field to be
the sea ice cover.
An approproate preprocessor would be::

  preprocessors:
    timeseries_NHW_ice_extent: # North Hemisphere Winter ice_extent
      custom_order: true
      extract_time:
          start_year: 1960
          start_month: 12
          start_day: 1
          end_year: 2005
          end_month: 9
          end_day: 31
      extract_season:
        season: DJF
      extract_region:
        start_longitude: -180.
        end_longitude: 180.
        start_latitude: 0.
        end_latitude: 90.


Note that this recipe may not function on machines with no access to the
internet, as cartopy may try to download the shapefiles. The solution to
this issue is the put the relevant cartopy shapefiles on a disk visible to your
machine, then link that path to ESMValTool via the `auxiliary_data_dir`
variable. The cartopy masking files can be downloaded from::

    https://www.naturalearthdata.com/downloads/

Here, cartopy uses the 1:10, physical coastlines and land files::

        110m_coastline.dbf  110m_coastline.shp  110m_coastline.shx
        110m_land.dbf  110m_land.shp  110m_land.shx

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import itertools
import logging
import os
import sys

import cartopy
import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


# Note that this recipe may not function on machines with no access to
# the internet, as cartopy may try to download geographic files.


def create_ice_cmap(threshold=0.15):
    """
    Create colour map with ocean blue below a threshold and white above.

    Parameters
    ----------
    threshold: float
        The threshold for the line between blue and white.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap:
        The resulting colour map.

    """
    threshold = threshold / 100.
    ice_cmap_dict = {
        'red': ((0., 0.0313, 0.0313), (threshold, 0.0313, 1.), (1., 1., 1.)),
        'green': ((0., 0.237, 0.237), (threshold, 0.237, 1.), (1., 1., 1.)),
        'blue': ((0., 0.456, 0.456), (threshold, 0.456, 1.), (1., 1., 1.))
    }

    return matplotlib.colors.LinearSegmentedColormap('ice_cmap', ice_cmap_dict)


def calculate_area_time_series(cube, plot_type, threshold):
    """
    Calculate the area of unmasked cube cells.

    Requires a cube with two spacial dimensions. (no depth coordinate).

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube
    plot_type: str
        The type of plot: ice extent or ice area
    threshold: float
        The threshold for ice fraction (typically 15%)

    Returns
    -------
    numpy array:
        An numpy array containing the time points.
    numpy.array:
        An numpy array containing the total ice extent or total ice area.

    """
    data = []
    times = diagtools.cube_time_to_float(cube)
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

        logger.debug('Calculating time series area: %s, %s, %s,', time_itr,
                     time, total_area)
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
    Make a ice extent and ice area time series plot for an individual model.

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
    iris.coord_categorisation.add_year(cube, 'time')
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

            times, data = calculate_area_time_series(cube_layer, plot_type,
                                                     threshold)

            plt.plot(times, data)

            # Add title to plot
            title = ' '.join(
                [metadata['dataset'], pole, 'hemisphere', season, plot_type])
            if layer:
                title = ' '.join([
                    title, '(', layer,
                    str(cube_layer.coords('depth')[0].units), ')'
                ])
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
    Make a polar stereoscopic map plot.

    The cube is the opened cube (two dimensional),
    pole is the polar region (North/South)
    cmap is the colourmap,

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube
    pole: str
        The hemisphere
    cmap: str
        The string describing the matplotlib colourmap.

    Returns
    ----------
    matplotlib.pyplot.figure:
        The matplotlib figure where the map was drawn.
    matplotlib.pyplot.axes:
        The matplotlib axes where the map was drawn.

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

    linrange = np.linspace(0, 100, 21)
    qplt.contourf(cube, linrange, cmap=cmap, linewidth=0, rasterized=True)
    plt.tight_layout()

    try:
        ax1.add_feature(
            cartopy.feature.LAND,
            zorder=10,
            facecolor=[0.8, 0.8, 0.8],
        )
    except ConnectionRefusedError:
        logger.error('Cartopy was unable add coastlines due to  a '
                     'connection error.')
    ax1.gridlines(
        linewidth=0.5, color='black', zorder=20, alpha=0.5, linestyle='--')
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('make_polar_map: Not able to add coastlines')
    return fig


def get_pole(cube):
    """
    Figure out the hemisphere and returns it as a string (North or South).

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube

    Returns
    ----------
    str:
        The hemisphere (North or South)

    """
    margin = 5.
    if np.max(cube.coord('latitude').points) < 0. + margin:
        return 'South'
    if np.min(cube.coord('latitude').points) > 0. - margin:
        return 'North'
    logger.fatal('get_pole: Not able to determine hemisphere.')
    return False


def get_time_string(cube):
    """
    Return a climatological season string in the format: "year season".

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube

    Returns
    ----------
    str:
        The climatological season as a string

    """
    season = cube.coord('clim_season').points
    year = cube.coord('year').points
    return str(int(year[0])) + ' ' + season[0].upper()


def get_year(cube):
    """
    Return the cube year as a string.

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube

    Returns
    ----------
    str:
        The year as a string

    """
    year = cube.coord('year').points
    return str(int(year))


def get_season(cube):
    """
    Return a climatological season time string.

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube

    Returns
    ----------
    str:
        The climatological season as a string

    """
    season = cube.coord('clim_season').points
    return season[0].upper()


def make_map_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple map plot for an individual model.

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
    iris.coord_categorisation.add_year(cube, 'time')
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
    for plot_type, plot_time in itertools.product(plot_types, plot_times):
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
            make_polar_map(cube, pole=pole, cmap=cmap)

            # Add title to plot
            title = ' '.join([metadata['dataset'], plot_type, time_str])
            if layer:
                title = ' '.join([
                    title, '(', layer,
                    str(cube_layer.coords('depth')[0].units), ')'
                ])
            plt.title(title)

            # Determine image filename:
            suffix = '_'.join(
                ['ortho_map', plot_type, time_str,
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

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube

    Returns
    ----------
    iris.cube.Cube:
        Data Cube with the seasonal means

    """
    if not cube.coords('clim_season'):
        iris.coord_categorisation.add_season(cube, 'time', name='clim_season')
    if not cube.coords('season_year'):
        iris.coord_categorisation.add_season_year(
            cube, 'time', name='season_year')
    return cube.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN)


def make_map_extent_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make an extent map plot showing several times for an individual model.

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
    iris.coord_categorisation.add_year(cube, 'time')
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
        try:
            ax1.add_feature(
                cartopy.feature.LAND, zorder=10, facecolor=[0.8, 0.8, 0.8])
        except ConnectionRefusedError:
            logger.error('Cartopy was unable add coastlines due to  a '
                         'connection error.')

        ax1.gridlines(
            linewidth=0.5, color='black', zorder=20, alpha=0.5, linestyle='--')

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
        legend_size = len(plot_desc) + 1
        ncols = int(legend_size / 25) + 1
        ax1.set_position([
            ax1.get_position().x0,
            ax1.get_position().y0,
            ax1.get_position().width * (1. - 0.1 * ncols),
            ax1.get_position().height
        ])

        fig.set_size_inches(7 + ncols * 1.2, 7)

        # Construct dummy plots.
        for i in sorted(plot_desc):
            plt.plot(
                [],
                [],
                c=plot_desc[i]['c'][0],
                lw=plot_desc[i]['lw'][0],
                ls=plot_desc[i]['ls'][0],
                label=plot_desc[i]['label'],
            )

        legd = ax1.legend(
            loc='center left',
            ncol=ncols,
            prop={'size': 10},
            bbox_to_anchor=(1., 0.5))
        legd.draw_frame(False)
        legd.get_frame().set_alpha(0.)

        # Add title to plot
        title = ' '.join([
            metadata['dataset'],
        ])
        if layer:
            title = ' '.join([
                title, '(', layer,
                str(cube_layer.coords('depth')[0].units), ')'
            ])
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
    Load the config file and metadata, then pass them the plot making tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)
        for filename in sorted(metadatas):

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
