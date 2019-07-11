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

Author: lova_to
"""
import logging
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import iris
import cartopy.crs as ccrs
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

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
        extend=extend,
        zmin=nspace.min(),
        zmax=nspace.max())

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
        the input files dictionairy
    obs_filename: str
        the preprocessed observations file.
    """
    logger.debug('make_multiple_plots')
    # ####
    filenames = list(metadata.keys())
    proj = ccrs.Robinson(central_longitude=0)
    # plot layout
    layout = metadata[filenames[0]]['layout_rowcol']

    # check if observations are provided
    hasobs = False
    obsname = ''
    if obs_filename != '':
        hasobs = True
        obsname = metadata[obs_filename]['dataset']
        filenames.remove(obs_filename)
        filenames.insert(0, obs_filename)
        layout[0] = layout[0] + 1
    else:
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

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    extend = 'neither'
    for layer in layers:

        fig = plt.figure()
        fig.set_size_inches(layout[1] * 4., layout[0] * 2. + 2.)

        # aggregate cubes
        name_cubes = []
        maps_cubes = []
        diff_cubes = []
        for thename in cubes:
            maps_cubes.append(cubes[thename][layer])
            name_cubes.append(thename)
            if (hasobs) & (thename != obsname):
                diff_cubes.append(cubes[thename][layer] -
                                  cubes[obsname][layer])

        # data ranges
        maps_range = diagtools.get_cube_range(maps_cubes)
        if hasobs:
            diff_range = diagtools.get_cube_range(diff_cubes)
            maps_range = diagtools.get_cube_range([cubes[obsname][layer]])

        if 'maps_range' in metadata[filenames[0]]:
            maps_range = metadata[filenames[0]]['maps_range']
            extend = 'both'
        if 'diff_range' in metadata[filenames[0]]:
            diff_range = metadata[filenames[0]]['diff_range']
            extend = 'both'

        # create subplots

        varunit = str(maps_cubes[0].units)
        varname = maps_cubes[0].var_name
        gsc = gridspec.GridSpec(layout[0], layout[1])
        yy = 0
        xx = 0
        clevels = 13
        for ii in range(len(maps_cubes)):
            hascbar = False
            cube = maps_cubes[ii]
            thename = name_cubes[ii]
            nspace = np.linspace(
                maps_range[0], maps_range[1], clevels, endpoint=True)
            cmap = 'viridis'

            if thename in [obsname, name_cubes[-1]]:
                hascbar = True

            if ii == 0:
                thename = thename + ' (' + varname + ') [' + varunit + ']'

            if (hasobs) & (ii > 0):
                cube = diff_cubes[ii - 1]
                nspace = np.linspace(
                    diff_range[0], diff_range[1], clevels, endpoint=True)
                cmap = 'RdBu_r'
                thename = thename

            axs = plt.subplot(gsc[xx, yy], projection=proj)
            add_map_plot(
                axs,
                cube,
                nspace,
                yy,
                cmap=cmap,
                title=thename,
                extend=extend,
                hascbar=hascbar)

            # next row & column indexes
            xx = xx + 1
            if xx == layout[0]:
                xx = 1 if hasobs else 0
                yy = yy + 1

        # Adjust subplots size & position
        plt.subplots_adjust(
            top=0.92,
            bottom=0.08,
            left=0.05,
            right=0.95,
            hspace=0.15,
            wspace=0.15)

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

        # Determine image filename:
        fn_list = ['multimodel_vs', obsname, varname, str(layer), 'maps']
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

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
        obs_filename = diagtools.match_model_to_key('observational_dataset',
                                                    cfg[model_type], metadatas)

        if not os.path.exists(obs_filename):
            logger.info('OBS file not found %s', obs_filename)
            obs_filename = ''

        make_multiple_plots(cfg, metadatas, obs_filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
