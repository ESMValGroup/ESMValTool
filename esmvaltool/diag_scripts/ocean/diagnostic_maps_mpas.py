"""
Maps diagnostics
================

Diagnostic to produce images of a map with coastlines from a cube.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

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
import cartopy.crs as ccrs

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvalcore.preprocessor._time import extract_time


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


long_name_dict = {
    'thetao': 'Temperature',
    'tos': 'Surface Temperature',
    'tob': 'Seafloor Temperature',
    'sos': 'Surface Salinity',
    'uo': 'Zonal Velocity',
    'vo': 'Meridional Velocity',
    'ph': 'Surface pH',
    'chl': 'Surface chlorophyll',
    'zos': 'Sea Surface Height',
    'no3': 'Dissolved Nitrate',
    'o2': 'Dissolved Oxygen',
    'intpp': 'Integrated Primary production'}

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
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

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

        qplt.contourf(cube_layer, 25, linewidth=0, rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join([
                title, '(', layer,
                str(cube_layer.coords('depth')[0].units), ')'
            ])
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', '_map_' + str(layer_index) + image_extention)
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


def make_map_contour(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple contour map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention and threshold.thresholds.
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    plot_details = {}
    colours = []
    thresholds = diagtools.load_thresholds(cfg, metadata)

    for itr, thres in enumerate(thresholds):
        if len(thresholds) > 1:
            colour = plt.cm.jet(float(itr) / float(len(thresholds) - 1.))
        else:
            colour = plt.cm.jet(0)
        label = str(thres) + ' ' + str(cube.units)
        colours.append(colour)
        plot_details[thres] = {'c': colour,
                               'lw': 1,
                               'ls': '-',
                               'label': label}

    linewidths = [1 for thres in thresholds]
    linestyles = ['-' for thres in thresholds]
    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)
        qplt.contour(cube_layer,
                     thresholds,
                     colors=colours,
                     linewidths=linewidths,
                     linestyles=linestyles,
                     rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')
        try:
            plt.gca().add_feature(cartopy.feature.LAND,
                                  zorder=10,
                                  facecolor=[0.8, 0.8, 0.8])
        except AttributeError:
            logger.warning('Not able to add coastlines')
        # Add legend
        diagtools.add_legend_outside_right(plot_details,
                                           plt.gca(),
                                           column_width=0.02,
                                           loc='below')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        depth_units = str(cube_layer.coords('depth')[0].units)
        if layer:
            title = '{} ({} {})'.format(title, layer, depth_units)
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = os.path.join(diagtools.folder(cfg['plot_dir']),
                                os.path.basename(filename))
            path = path.replace('.nc', '_contour_map_' + str(layer_index))
            path = path + image_extention
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='_contour_map_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_contours(
        cfg,
        metadata,
):
    """
    Make a contour map showing several models.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary.

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

    # Load threshold/thresholds.
    thresholds = diagtools.load_thresholds(cfg, metadata)

    # Make a plot for each layer and each threshold
    for layer, threshold in itertools.product(layers, thresholds):

        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('jet')
        land_drawn = False

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):

            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'
            linewidth = 1.
            linestyle = '-'

            # Determine line style for Observations
            if metadata[filename]['project'] in diagtools.get_obs_projects():
                color = 'black'
                linewidth = 1.7
                linestyle = '-'

            # Determine line style for MultiModel statistics:
            if 'MultiModel' in metadata[filename]['dataset']:
                color = 'black'
                linestyle = ':'
                linewidth = 1.4

            cube = model_cubes[filename][layer]
            qplt.contour(cube,
                         [threshold, ],
                         colors=[color, ],
                         linewidths=linewidth,
                         linestyles=linestyle,
                         rasterized=True)
            plot_details[filename] = {
                'c': color,
                'ls': linestyle,
                'lw': linewidth,
                'label': metadata[filename]['dataset']
            }

            if not land_drawn:
                try:
                    plt.gca().coastlines()
                except AttributeError:
                    logger.warning('Not able to add coastlines')
                plt.gca().add_feature(cartopy.feature.LAND,
                                      zorder=10,
                                      facecolor=[0.8, 0.8, 0.8])
                land_drawn = True

            title = metadata[filename]['long_name']
            if layer != '':
                z_units = model_cubes[filename][layer].coords('depth')[0].units
            units = str(model_cubes[filename][layer].units)

        # Add title, threshold, legend to plots
        title = ' '.join([title, str(threshold), units])
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
                suffix='_'.join(['_contour_map_',
                                 str(threshold),
                                 str(layer) + image_extention]),
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


def multi_model_map_figure(
        cfg,
        metadatas,
        short_name='',
        dataset='',
        ensemble='',
        figure_style = 'four_ssp',
        hist_time_range = [1990., 2000.],
        ssp_time_range = [2040., 2050.],
    ):
    """
    produce a monthly climatology figure.
    """
    central_longitude=-160.+3.5
    proj = ccrs.Robinson(central_longitude=central_longitude)

    fig = plt.figure()
    seq_cmap = 'viridis'
    div_cmap =' BrBG'
    if figure_style=='four_ssp':
        subplots = {221: 'ssp126', 222:'ssp245', 223:'ssp370', 224: 'ssp585'}
        subplot_style = {221: 'diff', 222: 'diff', 223: 'diff', 224: 'diff'}
        cmaps =  {221: div_cmap, 222:div_cmap, 223: div_cmap, 224: div_cmap}
    elif figure_style=='hist_and_ssp':
        subplots = {321: 'historical', 323: 'ssp126', 324: 'ssp245', 325: 'ssp370', 326: 'ssp585'}
        subplot_style = {321: 'mean', 323:'diff', 324: 'diff', 325: 'diff', 326: 'diff'}
        cmaps = {321: seq_cmap, 323:div_cmap, 324: div_cmap, 325: div_cmap, 326: div_cmap}
    elif figure_style in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
        subplots = {221: 'historical', 222:figure_style, 223:figure_style, 224: figure_style}
        subplot_style = {221:'mean',222: 'diff', 223: 'min_diff', 224: 'max_diff'}
        cmaps =  {221: seq_cmap, 222:div_cmap, 223: div_cmap, 224: div_cmap}
    else:
        assert 0
    cubes = {}
    for filename, metadata in metadatas.items():
        if short_name != metadata['short_name']:
            continue
        if dataset != metadata['dataset']:
            continue
        if ensemble != metadata['ensemble']:
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
            cube = extract_time(cube, hist_time_range[0], 1, 1, hist_time_range[1], 12, 31)
        else:
            cube = extract_time(cube, ssp_time_range[0], 1, 1, ssp_time_range[1], 12, 31)

        cube_mean = cube.copy().collapsed('time', iris.analysis.MEAN)
        cube_min = cube.copy().collapsed('time', iris.analysis.MIN)
        cube_max = cube.copy().collapsed('time', iris.analysis.MAX)

        cube_mean = cube_mean.intersection(longitude=(central_longitude-180., central_longitude+180.))
        cube_min = cube_min.intersection(longitude=(central_longitude-180., central_longitude+180.))
        cube_max = cube_max.intersection(longitude=(central_longitude-180., central_longitude+180.))

        cubes[(dataset, scenario, ensemble, 'mean')] = cube_mean
        cubes[(dataset, scenario, ensemble, 'min')] = cube_min
        cubes[(dataset, scenario, ensemble, 'max')] = cube_max

    if len(cubes.keys()) == 0:
        return

    for sbp, exp in subplots.items():
        ax = fig.add_subplot(subplot, projection=proj)
        sbp_style = subplot_style[sbp]

        if sbp_style in ['mean', 'min', 'max']:
            cube = cubes[(dataset, exp, ensemble, sbp_style)]
        if sbp_style == 'diff':
            cube = cubes[(dataset, exp, ensemble, 'mean')] - cubes[(dataset, 'historical', ensemble, 'mean')]
        if sbp_style == 'min_diff':
            cube = cubes[(dataset, exp, ensemble, 'min')] - cubes[(dataset, 'historical', ensemble, 'min')]
        if sbp_style == 'max_diff':
            cube = cubes[(dataset, exp, ensemble, 'max')] - cubes[(dataset, 'historical', ensemble, 'max')]

        qplot = iris.plot.contourf(
            cube,
            #nspace,
            linewidth=0,
            cmap=cmaps[sbp],
            extend='neither',)
            #zmin=plot_range[0],
            #zmax=plot_range[1])

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([exp, sbp_style])
        plt.title(title)

    suptitle = ' '.join([dataset, ensemble, long_name_dict[short_name],
                         '\n Historical', '-'.join([str(t) for t in hist_time_range]),
                         'vs SSP', '-'.join([str(t) for t in ssp_time_range]) ])

    pyplot.suptitle(suptitle)

    # save and close.
    path = diagtools.folder(cfg['plot_dir']+'/'+'_'.join([short_name, dataset, ensemble]))
    path += '_'.join([dataset, ensemble, short_name, figure_style, time_str])
    path += diagtools.get_image_format(cfg)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()


def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']
    short_names = {}
    datasets = {}
    ensembles = {}
    metadatas = diagtools.get_input_files(cfg, )
    for fn, metadata in metadatas.items():
        short_names[metadata['short_name']] = True
        datasets[metadata['dataset']] = True
        ensembles[metadata['ensemble']] = True

    for short_name, dataset, ensemble in itertools.product(short_names.keys(), datasets.keys(), ensembles.keys()):
        hist_time_ranges = [[1850, 2015], [1850, 1900], [1950, 2000], [1990, 2000], [1990, 2000]]
        ssp_time_ranges  = [[2015, 2100], [2050, 2100], [2050, 2100], [2040, 2050], [2090, 2100]]
        for hist_time_range, ssp_time_range in zip(hist_time_ranges, ssp_time_ranges):

            for figure_style in ['ssp126', 'ssp245', 'ssp370', 'ssp585', 'hist_and_ssp', 'four_ssp']:
                multi_model_map_figure(
                    cfg,
                    metadatas,
                    short_name=short_name,
                    ensemble=ensemble,
                    dataset=dataset,
                    figure_style=figure_style,
                    hist_time_range=hist_time_range,
                    ssp_time_range=ssp_time_range,
                )
                assert 0

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)
        thresholds = diagtools.load_thresholds(cfg, metadatas)

        if thresholds:
            #######
            # Multi model contour plots
            multi_model_contours(
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
            # Contour maps of individual model
            if thresholds:
                make_map_contour(cfg, metadatas[filename], filename)

            ######
            # Maps of individual model
            make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
