"""Model vs Observations maps Diagnostic.

Diagnostic to produce comparison maps of model(s) and data (if provided).
If observations are not provided, data maps for each model are drawn.

The image shows on top row observational data and the following subplot(s)
the comparison for each model by following the order of the recipe.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 2D + time field would be::

  preprocessors:
    prep_map:
      time_average:
      regrid:
        target_grid: 1x1
        scheme: linear

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: lovato_tomas
"""
import logging
import os
from pprint import pformat

import cartopy.crs as ccrs
import iris
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(plot_file, attributes, obsname, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    if obsname != '':
        caption = (
            "{long_name} bias for average between {start_year} and {end_year}".
            format(**attributes) + " against " + obsname + " observations.")
    else:
        caption = (
            "Average {long_name} between {start_year} and {end_year} ".format(
                **attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'map',
        'authors': [
            'lovato_tomas',
        ],
        'references': [
            'acknow_project',
        ],
        'plot_file': plot_file,
        'ancestors': ancestor_files,
    }
    return record


def add_map_plot(axs, plot_cube, cols):
    """
    Add a map in the current pyplot suplot.

    Parameters
    ----------
    axs: object
        The matplotlib.pyplot Axes object
    plot_cube: dictionary
        dictionary with data for plot defined in select_cubes
    cols: integer
        Number of columns in the multipanel plot
    """
    contour_lev = 13
    nspace = np.linspace(plot_cube['range'][0],
                         plot_cube['range'][1],
                         contour_lev,
                         endpoint=True)
    iris.plot.contourf(plot_cube['cube'],
                       nspace,
                       linewidth=0,
                       cmap=plt.cm.get_cmap(plot_cube['cmap']),
                       extend=plot_cube['extend'])

    axs.coastlines()
    gls = axs.gridlines(draw_labels=False, color='black', alpha=0.4)
    gls.ylocator = mticker.FixedLocator(np.linspace(-90., 90., 7))
    axs.set_title(plot_cube['title'], fontweight="bold", fontsize='large')

    if plot_cube['hascbar']:
        bba = (0., -0.1, 1, 1)
        if cols > 0:
            bba = ((0.5 - cols) * 1.1, -0.15, cols, 1)
        axins = inset_axes(
            axs,
            width="95%",
            height="6%",
            loc='lower center',
            bbox_to_anchor=bba,
            bbox_transform=axs.transAxes,
            borderpad=0,
        )
        cbar = plt.colorbar(orientation='horizontal', cax=axins)
        cbar.set_ticks(nspace[::2])


def make_subplots(cubes, layout, obsname, fig):
    """
    Realize subplots using cubes input data.

    Parameters
    ----------
    cubes: dict
        dictionary with data for plot defined in select_cubes
    layout : list
        subplot rows x cols organization
    obsname: string
        Observation data name
    fig: object
         The matplotlib.pyplot Figure object
    """
    proj = ccrs.Robinson(central_longitude=0)
    gsc = gridspec.GridSpec(layout[0], layout[1])
    row = 0
    col = 0
    for thename in cubes:
        axs = plt.subplot(gsc[row, col], projection=proj)
        add_map_plot(axs, cubes[thename], col)
        # next row & column indexes
        row = row + 1
        if row == layout[0]:
            row = 1 if obsname != '' else 0
            col = col + 1

    # Adjust subplots size & position
    plt.subplots_adjust(top=0.92,
                        bottom=0.08,
                        left=0.05,
                        right=0.95,
                        hspace=0.15,
                        wspace=0.15)

    # Vertically detach OBS plot and center
    if obsname != '':
        axs = fig.axes
        box = axs[0].get_position()
        shift = box.y0 * 0.05
        box.y0 = box.y0 + shift
        box.y1 = box.y1 + shift
        shift = box.x1 - box.x0
        box.x0 = 0.5 - shift * 0.5
        box.x1 = box.x0 + shift
        axs[0].set_position(box)


def load_cubes(filenames, obs_filename, metadata):
    """
    Organize data provided by recipe.

    Parameters
    ----------
    filenames: dict
        input files listed in the recipe
    obs_filename: str
        the preprocessed observations file.
    metadata: dict
        the input files dictionary
    """
    # check if observations are provided
    if obs_filename:
        obsname = metadata[obs_filename]['dataset']
        filenames.remove(obs_filename)
        filenames.insert(0, obs_filename)
    else:
        obsname = ''
        logger.info('Observations not provided. Plot each model data.')

    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for thename in filenames:
        logger.debug('loading: \t%s', thename)
        cube = iris.load_cube(thename)
        cube = diagtools.bgc_units(cube, metadata[thename]['short_name'])
        model_name = metadata[thename]['dataset']
        cubes[model_name] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_name]:
            layers[layer] = True

    logger.debug('layers: %s', layers)
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    return cubes, layers, obsname


def select_cubes(cubes, layer, obsname, metadata):
    """
    Create a dictionary of input layer data & metadata to plot.

    Parameters
    ----------
    cubes: list
        Input data iris cubes
    layer: list
        Data level to be plotted
    obsname: string
         Observation data name
    metadata: dict
        the first input file dictionary
    """
    plot_cubes = {}
    list_cubes = []

    for thename in cubes:
        plot_cubes[thename] = {
            'cube': cubes[thename][layer],
            'title': thename,
            'cmap': 'viridis',
            'range': None,
            'extend': 'neither',
            'hascbar': False
        }
        if (obsname != '') & (thename != obsname):
            plot_cubes[thename] = {
                'cube': cubes[thename][layer] - cubes[obsname][layer],
                'title': thename,
                'cmap': 'RdBu_r',
                'range': None,
                'extend': 'neither',
                'hascbar': False
            }
        if thename in [obsname, list(cubes.keys())[-1]]:
            plot_cubes[thename]['hascbar'] = True

        if thename == list(cubes.keys())[0]:
            cube = plot_cubes[thename]['cube']
            plot_cubes[thename][
                'title'] = thename + ' (' + cube.var_name + ') [' + str(
                    cube.units) + ']'
        list_cubes.append(plot_cubes[thename]['cube'])

    # get cubes data ranges
    mrange = diagtools.get_cube_range(list_cubes)
    if obsname != '':
        mrange = diagtools.get_cube_range([list_cubes[0]])
        drange = diagtools.get_cube_range(list_cubes[1:])

    # get user defined plot ranges
    user_range = {'maps': None, 'diff': None}
    if 'maps_range' in metadata:
        user_range['maps'] = metadata['maps_range']
    if 'diff_range' in metadata:
        user_range['diff'] = metadata['diff_range']

    # define contour levels using ranges
    for thename in cubes:
        mrange = mrange
        if user_range['maps']:
            mrange = user_range['maps']
            plot_cubes[thename]['extend'] = 'both'
        if (obsname != '') & (thename != obsname):
            mrange = drange
            if user_range['diff']:
                mrange = user_range['diff']
                plot_cubes[thename]['extend'] = 'both'
        plot_cubes[thename]['range'] = mrange

    return plot_cubes


def make_multiple_plots(cfg, metadata, obsname):
    """
    Produce multiple panel comparison maps of model(s) and data (if provided).

    If observations are not provided, plots of each model data are drawn.
    Put on top row observational data (if available) and in following subplots
    model difference (or data) organized in rows/cols using row/col layout.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the input files dictionary
    obsname: str
        the preprocessed observations file.
    """
    logger.debug('make_multiple_plots')
    # ####
    filenames = list(metadata.keys())
    varname = metadata[filenames[0]]['short_name']

    # plot setting
    layout = metadata[filenames[0]]['layout_rowcol']

    # load input data
    [cubes, layers, obsname] = load_cubes(filenames, obsname, metadata)

    if obsname != '':
        layout[0] = layout[0] + 1
    else:
        logger.info('Observations not provided. Plot each model data.')

    if len(filenames) > (layout[0] * layout[1]):
        raise ValueError(
            'Number of inputfiles is larger than layout scheme (rows x cols). '
            'Revise layout_rowcol size in recipe.')

    # Make a plot for each layer
    for layer in layers:

        fig = plt.figure()
        fig.set_size_inches(layout[1] * 4., layout[0] * 2. + 2.)

        # select cubes to plot
        plot_cubes = select_cubes(cubes, layer, obsname,
                                  metadata[filenames[0]])

        # create individual subplot
        make_subplots(plot_cubes, layout, obsname, fig)

        # Determine image filename:
        if obsname != '':
            plot_file = ['multimodel_vs', obsname, varname, str(layer), 'maps']
        else:
            plot_file = ['multimodel', varname, str(layer), 'maps']
        plot_file = '_'.join(plot_file)
        path = diagtools.folder(cfg['plot_dir']) + plot_file

        # Saving file:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path, dpi=200)

        # Provenance
        provenance_record = get_provenance_record(plot_file,
                                                  metadata[filenames[-1]],
                                                  obsname, filenames)
        logger.info("Recording provenance of %s:\n%s", plot_file,
                    pformat(provenance_record))
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(plot_file, provenance_record)

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
            'metadata filename:\t%s, %s',
            index,
            metadata_filename,
        )
        metadatas = diagtools.get_input_files(cfg, index=index)

        model_type = 'observational_dataset'
        logger.debug(
            'model_type: %s, %s',
            index,
            model_type,
        )
        logger.debug(
            'metadatas:  %s, %s',
            index,
            metadatas,
        )
        obs_filename = None
        if model_type in cfg.keys():
            obs_filename = diagtools.match_model_to_key(
                'observational_dataset', cfg[model_type], metadatas)
            if not os.path.isfile(obs_filename):
                logger.info('OBS file not found %s', obs_filename)

        make_multiple_plots(cfg, metadatas, obs_filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
