"""
Model vs Observations maps Diagnostic.
======================================

Diagnostic to produce comparison of model and data.
The first kind of image shows four maps and the other shows a scatter plot.

The four pane image is a latitude vs longitude figures showing:

* Top left: model
* Top right: observations
* Bottom left: model minus observations
* Bottom right: model over observations


The scatter plots plot the matched model coordinate on the x axis, and the
observational dataset on the y coordinate, then performs a linear
regression of those data and plots the line of best fit on the plot.
The parameters of the fit are also shown on the figure.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D + time field would be::

  preprocessors:
    prep_map:
      extract_levels:
        levels:  [100., ]
        scheme: linear_extrap
      climate_statistics:
        operator: mean
      regrid:
        target_grid: 1x1
        scheme: linear

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
import math

from matplotlib import pyplot
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt
import numpy as np
from scipy.stats import linregress

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def add_map_subplot(subplot, cube, nspace, title='',
                    cmap='', extend='neither', log=False):
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
    extend: str
        Contourf-coloring of values outside the levels range
    log: bool
        Flag to plot the colour scale linearly (False) or
        logarithmically (True)
    """
    plt.subplot(subplot)
    logger.info('add_map_subplot: %s', subplot)
    if log:
        qplot = qplt.contourf(
            cube,
            nspace,
            linewidth=0,
            cmap=plt.cm.get_cmap(cmap),
            norm=LogNorm(),
            zmin=nspace.min(),
            zmax=nspace.max())
        qplot.colorbar.set_ticks([0.1, 1., 10.])
    else:
        qplot = iris.plot.contourf(
            cube,
            nspace,
            linewidth=0,
            cmap=plt.cm.get_cmap(cmap),
            extend=extend,
            zmin=nspace.min(),
            zmax=nspace.max())
        cbar = pyplot.colorbar(orientation='horizontal')
        cbar.set_ticks(
            [nspace.min(), (nspace.max() + nspace.min()) / 2.,
             nspace.max()])

    plt.gca().coastlines()
    plt.title(title)


def make_model_vs_obs_plots(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
    Make a figure showing four maps and the other shows a scatter plot.

    The four pane image is a latitude vs longitude figures showing:

    * Top left: model
    * Top right: observations
    * Bottom left: model minus observations
    * Bottom right: model over observations

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        the input files dictionairy
    model_filename: str
        the preprocessed model file.
    obs_filename: str
        the preprocessed observations file.

    """
    filenames = {'model': model_filename, 'obs': obs_filename}
    logger.debug('make_model_vs_obs_plots filenames: %s', filenames)
    # ####
    # Load the data for each layer as a separate cube
    input_file = None
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        logger.debug('loading: \t%s, \t%s', model_type, input_file)
        cube = iris.load_cube(input_file)
        cube = diagtools.bgc_units(cube, metadata[input_file]['short_name'])
        cubes[model_type] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_type]:
            layers[layer] = True

    logger.debug('layers: %s', layers)
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    # ####
    # load names:
    model = metadata[filenames['model']]['dataset']
    obs = metadata[filenames['obs']]['dataset']

    long_name = cubes['model'][list(layers.keys())[0]].long_name
    units = str(cubes['model'][list(layers.keys())[0]].units)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:

        fig = plt.figure()
        fig.set_size_inches(9, 6)

        # Create the cubes
        cube221 = cubes['model'][layer]
        cube222 = cubes['obs'][layer]
        cube223 = cubes['model'][layer] - cubes['obs'][layer]
        cube224 = cubes['model'][layer] / cubes['obs'][layer]

        # create the z axis for plots 2, 3, 4.
        extend = 'neither'
        zrange12 = diagtools.get_cube_range([cube221, cube222])
        if 'maps_range' in metadata[input_file]:
            zrange12 = metadata[input_file]['maps_range']
            extend = 'both'
        zrange3 = diagtools.get_cube_range_diff([cube223])
        if 'diff_range' in metadata[input_file]:
            zrange3 = metadata[input_file]['diff_range']
            extend = 'both'

        cube224.data = np.ma.clip(cube224.data, 0.1, 10.)

        n_points = 12
        linspace12 = np.linspace(
            zrange12[0], zrange12[1], n_points, endpoint=True)
        linspace3 = np.linspace(
            zrange3[0], zrange3[1], n_points, endpoint=True)
        logspace4 = np.logspace(-1., 1., 12, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(
            221, cube221, linspace12, cmap='viridis', title=model,
            extend=extend)
        add_map_subplot(
            222, cube222, linspace12, cmap='viridis',
            title=' '.join([obs]),
            extend=extend)
        add_map_subplot(
            223,
            cube223,
            linspace3,
            cmap='bwr',
            title=' '.join([model, 'minus', obs]),
            extend=extend)
        if np.min(zrange12) > 0.:
            add_map_subplot(
                224,
                cube224,
                logspace4,
                cmap='bwr',
                title=' '.join([model, 'over', obs]),
                log=True)

        # Add overall title
        fig.suptitle(long_name + ' [' + units + ']', fontsize=14)

        # Determine image filename:
        fn_list = ['model_vs_obs', long_name, model, obs, str(layer), 'maps']
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path, dpi=200)

        plt.close()


def rounds_sig(value, sig=3):
    """
    Round a float to a specific number of sig. figs. & return it as a string.

    Parameters
    ----------
    value: float
       The float that is to be rounded.
    sig: int
        The number of significant figures.

    Returns
    ----------
    str:
        The rounded output string.

    """
    if value == 0.:
        return str(0.)
    if value < 0.:
        value = abs(value)
        return str(
            -1. * round(value, sig - int(math.floor(math.log10(value))) - 1))
    return str(round(value, sig - int(math.floor(math.log10(value))) - 1))


def add_linear_regression(plot_axes,
                          arr_x,
                          arr_y,
                          showtext=True,
                          add_diagonal=False,
                          extent=None):
    """
    Add a straight line fit to an axis.

    Parameters
    ----------
    plot_axes: matplotlib.pyplot.axes
        The matplotlib axes on which to plot the linear regression.
    arr_x: numpy.array
        The data for the x coordinate.
    arr_y: numpy array
        The data for the y coordinate.
    showtext: bool
        A flag to turn on or off the result of the fit on the plot.
    add_diagonal: bool
        A flag to also add the 1:1 diagonal line to the figure
    extent: list of floats
        The extent of the plot axes.
    """
    beta1, beta0, r_value, p_value, std_err = linregress(arr_x, arr_y)
    texts = [
        r'$\^\beta_0$ = ' + rounds_sig(beta0),
        r'$\^\beta_1$ = ' + rounds_sig(beta1),
        r'R = ' + rounds_sig(r_value),
        r'P = ' + rounds_sig(p_value),
        r'N = ' + str(int(len(arr_x)))
    ]
    thetext = '\n'.join(texts)

    if showtext:
        pyplot.text(
            0.04,
            0.96,
            thetext,
            horizontalalignment='left',
            verticalalignment='top',
            transform=plot_axes.transAxes)

    if extent is None:
        x_values = np.arange(arr_x.min(), arr_x.max(),
                             (arr_x.max() - arr_x.min()) / 20.)
        y_values = [beta0 + beta1 * a for a in x_values]
    else:
        minv = min(extent)
        maxv = max(extent)
        x_values = np.arange(minv, maxv, (maxv - minv) / 1000.)
        y_values = np.array([beta0 + beta1 * a for a in x_values])

        mask = (x_values < minv) + (y_values < minv) \
            + (x_values > maxv) + (y_values > maxv)
        x_values = np.ma.masked_where(mask, x_values)
        y_values = np.ma.masked_where(mask, y_values)

    pyplot.plot(x_values, y_values, 'k')

    if add_diagonal:
        axis = pyplot.gca().axis()
        step = (max(axis) - min(axis)) / 100.
        one_to_one = np.arange(min(axis), max(axis) + step, step)
        pyplot.plot(one_to_one, one_to_one, 'k--')


def make_scatter(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
    Makes Scatter plots of model vs observational data.

    Make scatter plot showing the matched model and observational data with the
    model data as the x-axis coordinate and the observational data as the
    y-axis coordinate. A linear regression is also applied to the matched
    data and the result of the fit is shown on the figure.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        the input files dictionairy
    model_filename: str
        the preprocessed model file.
    obs_filename: str
        the preprocessed observations file.
    """

    filenames = {'model': model_filename, 'obs': obs_filename}
    logger.debug('make_model_vs_obs_plots: \t%s', filenames)
    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        logger.debug('loading: \t%s, \t%s', model_type, input_file)
        cube = iris.load_cube(input_file)
        cube = diagtools.bgc_units(cube, metadata[input_file]['short_name'])
        cubes[model_type] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[model_type]:
            layers[layer] = True

    logger.debug('layers: %s', layers)
    logger.debug('cubes: %s', ', '.join(cubes))

    # ####
    # load names:
    model = metadata[filenames['model']]['dataset']
    obs = metadata[filenames['obs']]['dataset']

    long_name = cubes['model'][list(layers.keys())[0]].long_name

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:

        fig = plt.figure()
        fig.set_size_inches(7, 6)

        # Create the cubes
        model_data = np.ma.array(cubes['model'][layer].data)
        obs_data = np.ma.array(cubes['obs'][layer].data)

        mask = model_data.mask + obs_data.mask
        model_data = np.ma.masked_where(mask, model_data).compressed()
        obs_data = np.ma.masked_where(mask, obs_data).compressed()

        colours = 'gist_yarg'
        zrange = diagtools.get_array_range([model_data, obs_data])
        plotrange = [zrange[0], zrange[1], zrange[0], zrange[1]]

        x_scale = 'log'
        if np.min(zrange) * np.max(zrange) < -1:
            x_scale = 'linear'
        if np.min(zrange) < 0.:
            logger.info('Skip scatter for %s. Min is < 0', long_name)
            return

        pyplot.hexbin(
            model_data,
            obs_data,
            xscale=x_scale,
            yscale=x_scale,
            bins='log',
            # extent=np.log10(plotrange),
            gridsize=50,
            cmap=pyplot.get_cmap(colours),
            mincnt=0)
        cbar = pyplot.colorbar()
        cbar.set_label('log10(N)')

        pyplot.gca().set_aspect("equal")
        pyplot.axis(plotrange)

        add_linear_regression(
            pyplot.gca(),
            model_data,
            obs_data,
            showtext=True,
            add_diagonal=True,
            extent=plotrange)

        pyplot.title(long_name)
        pyplot.xlabel(model)
        pyplot.ylabel(obs)

        # Determine image filename:
        fn_list = [
            'model_vs_obs', long_name, model, obs,
            str(layer), 'scatter'
        ]
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
        for filename in sorted(metadatas.keys()):

            if filename == obs_filename:
                continue
            if not os.path.exists(obs_filename):
                continue
            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            # #####
            # model vs obs scatter plots
            make_scatter(cfg, metadatas, filename, obs_filename)

            # #####
            # model vs obs map plots
            make_model_vs_obs_plots(cfg, metadatas, filename, obs_filename)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
