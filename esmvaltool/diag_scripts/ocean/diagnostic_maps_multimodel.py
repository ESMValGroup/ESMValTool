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

Author: lovato_tomas
"""
import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def add_map_plot(ax,
                 cube,
                 nspace,
                 cols,
                 title='',
                 cmap='',
                 extend='neither',
                 hascbar=False):
    """
    Add a map in the current pyplot suplot.

    Parameters
    ----------
    ax: object
        The matplotlib.pyplot Axes object
    cube: iris.cube.Cube
        the iris cube to be plotted.
    nspace: numpy.array
        An array of the ticks of the colour part.
    cols: integer
        Number of columns in the multipanel plot
    title: str
        A string to set as the subplot title.
    cmap: str
        A string to describe the matplotlib colour map.
    extend: str
        Contourf-coloring of values outside the levels range
    hascbar: logical
        Add colorbar to the subplot
    """
    iris.plot.contourf(
        cube,
        nspace,
        linewidth=0,
        cmap=plt.cm.get_cmap(cmap),
        extend=extend)

    ax.coastlines()
    gls = ax.gridlines(draw_labels=False, color='black', alpha=0.4)
    gls.ylocator = mticker.FixedLocator(np.linspace(-90., 90., 7))
    ax.set_title(title, fontweight="bold", fontsize='large')

    if hascbar:
        bba = (0., -0.1, 1, 1)
        if cols > 0:
            bba = ((0.5 - cols) * 1.1, -0.15, cols, 1)
        axins = inset_axes(
            ax,
            width="95%",
            height="6%",
            loc='lower center',
            bbox_to_anchor=bba,
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        cbar = plt.colorbar(orientation='horizontal', cax=axins)
        cbar.set_ticks(nspace[::2])


def adjust_subplot_spacing(fig, hasobs):
    """
    Adjust spacing to avoid subplots overlays and improve readability.

    Parameters
    ----------
    fig: object
         The matplotlib.pyplot Figure object
    hasobs : logical
         Check if obs are plotted
    """
    # Adjust subplots size & position
    plt.subplots_adjust(
        top=0.92, bottom=0.08, left=0.05, right=0.95, hspace=0.15, wspace=0.15)

    # Vertically detach OBS plot and center
    if hasobs:
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
        cube = diagtools.bgc_units(cube, cube.var_name)
        model_name = metadata[thename]['dataset']
        cubes[model_name] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_name]:
            layers[layer] = True

    logger.debug('layers: %s', layers)
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    return cubes, layers, obsname


def select_cubes(cubes, layer, obsname, user_range, clevels):
    """
    Create a dictionary of input layer data & metadata for plot.

    Parameters
    ----------
    cubes: list
        Input data iris cubes
    layer: list
        Data level to be plotted
    obsname: string
         Observation data name
    user_range: dict
         Plot ranges read from recipe
    clevels: integer
         Number of contour levels
    """
    plot_cubes = {}
    list_cubes = []

    for thename in cubes:
        plot_cubes[thename] = {
            'cube': cubes[thename][layer],
            'cmap': 'viridis',
            'nspace': None,
            'extend': 'neither'
        }
        if (obsname != '') & (thename != obsname):
            plot_cubes[thename] = {
                'cube': cubes[thename][layer] - cubes[obsname][layer],
                'cmap': 'RdBu_r',
                'nspace': None,
                'extend': 'neither'
            }
        list_cubes.append(plot_cubes[thename]['cube'])

    # get cubes data ranges
    mrange = diagtools.get_cube_range(list_cubes)
    if obsname != '':
        mrange = diagtools.get_cube_range([list_cubes[0]])
        drange = diagtools.get_cube_range(list_cubes[1:])

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
        plot_cubes[thename]['nspace'] = np.linspace(
            mrange[0], mrange[1], clevels, endpoint=True)

    return plot_cubes


def make_multiple_plots(cfg, metadata, obs_filename):
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
    obs_filename: str
        the preprocessed observations file.
    """
    logger.debug('make_multiple_plots')
    # ####
    filenames = list(metadata.keys())
    varname = metadata[filenames[0]]['short_name']
    user_range = {'maps': None, 'diff': None}
    if 'maps_range' in metadata[filenames[0]]:
        user_range['maps'] = metadata[filenames[0]]['maps_range']
    if 'diff_range' in metadata[filenames[0]]:
        user_range['diff'] = metadata[filenames[0]]['diff_range']

    # plot setting
    layout = metadata[filenames[0]]['layout_rowcol']
    proj = ccrs.Robinson(central_longitude=0)
    contour_lev = 13

    # load input data
    [cubes, layers, obsname] = load_cubes(filenames, obs_filename, metadata)

    if obsname != '':
        hasobs = True
        layout[0] = layout[0] + 1
    else:
        hasobs = False
        logger.info('Observations not provided. Plot each model data.')

    # Make a plot for each layer
    for layer in layers:

        fig = plt.figure()
        fig.set_size_inches(layout[1] * 4., layout[0] * 2. + 2.)

        # select cubes to plots
        plot_cubes = select_cubes(cubes, layer, obsname, user_range,
                                  contour_lev)

        # create subplots
        gsc = gridspec.GridSpec(layout[0], layout[1])
        row = 0
        col = 0
        for thename in plot_cubes:
            hascbar = False
            cube = plot_cubes[thename]['cube']
            model_name = thename

            if thename in [obsname, list(plot_cubes.keys())[-1]]:
                hascbar = True

            if thename == list(plot_cubes.keys())[0]:
                model_name = thename + ' (' + varname + ') [' + str(
                    cube.units) + ']'

            axs = plt.subplot(gsc[row, col], projection=proj)
            add_map_plot(
                axs,
                cube,
                plot_cubes[thename]['nspace'],
                col,
                cmap=plot_cubes[thename]['cmap'],
                title=model_name,
                extend=plot_cubes[thename]['extend'],
                hascbar=hascbar)

            # next row & column indexes
            row = row + 1
            if row == layout[0]:
                row = 1 if hasobs else 0
                col = col + 1

        adjust_subplot_spacing(fig, hasobs)

        # Determine image filename:
        fn_list = ['multimodel_vs', obsname, varname, str(layer), 'maps']
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path, dpi=200)

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
