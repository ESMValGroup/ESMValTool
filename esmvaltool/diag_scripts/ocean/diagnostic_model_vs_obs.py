"""
Model vs Observations maps Diagnostic
=====================================

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
      time_average:
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
import matplotlib
matplotlib.use('Agg')  # noqa
from matplotlib import ticker, pyplot
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt

import math
import numpy as np
from scipy.stats import linregress

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def add_map_subplot(subplot, cube, nspace, title='', cmap='', log=False):
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
    log: bool
        Flag to plot the colour scale linearly (False) or
        logarithmically (True)
    """

    plt.subplot(subplot)
    logger.info('add_map_subplot: %s', subplot)
    if log:
        qplot = qplt.contourf(cube, nspace, linewidth=0,
                              cmap=plt.cm.get_cmap(cmap),
                              norm=LogNorm(),
                              zmin=nspace.min(),
                              zmax=nspace.max())
        qplot.colorbar.set_ticks([0.1, 1., 10.])
    else:
        qplot = iris.plot.contourf(cube, nspace, linewidth=0,
                                   cmap=plt.cm.get_cmap(cmap),
                                   zmin=nspace.min(),
                                   zmax=nspace.max())
        cb = pyplot.colorbar(orientation='horizontal')
        cb.set_ticks([nspace.min(), (nspace.max() + nspace.min())/2.,
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

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers.keys():

        fig = plt.figure()
        fig.set_size_inches(9, 6)

        # Create the cubes
        cube221 = cubes['model'][layer]
        cube222 = cubes['obs'][layer]
        cube223 = cubes['model'][layer] - cubes['obs'][layer]
        cube224 = cubes['model'][layer] / cubes['obs'][layer]

        # create the z axis for plots 2, 3, 4.
        zrange12 = diagtools.get_cube_range([cube221, cube222])
        zrange3 = diagtools.get_cube_range_diff([cube223])

        cube224.data = np.ma.clip(cube224.data, 0.1, 10.)
        zrange4 = [0.1, 10.]

        n_points = 12
        linspace12 = np.linspace(zrange12[0], zrange12[1], n_points,
                                 endpoint=True)
        linspace3 = np.linspace(zrange3[0], zrange3[1], n_points,
                                endpoint=True)
        logspace4 = np.logspace(-1., 1., 12, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(221, cube221, linspace12, cmap='viridis',
                        title=model)
        add_map_subplot(222, cube222, linspace12, cmap='viridis',
                        title=' '.join([obs, ]))
        add_map_subplot(223, cube223, linspace3, cmap='bwr',
                        title=' '.join([model, 'minus', obs]))
        add_map_subplot(224, cube224, logspace4, cmap='bwr',
                        title=' '.join([model, 'over', obs]), log=True)

        # Add overall title
        fig.suptitle(long_name, fontsize=14)

        # Determine image filename:
        fn_list = ['model_vs_obs', long_name, model, obs, str(layer), 'maps']
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def round_sig(x, sig=3):
    """
    Function to round a float to a specific number of significant figures
    and return it as a string.

    Parameters
    ----------
    x: float
       The float that is to be rounded.
    sig: int
        The number of significant figures.

    Returns
    ----------
    str:
        The rounded output string.

    """
    if x == 0.:
        return str(0.)
    if x < 0.:
        return str(-1. * round(abs(x),
                               sig - int(math.floor(math.log10(abs(x)))) - 1))
    else:
        return str(round(x, sig - int(math.floor(math.log10(x))) - 1))


def add_linear_regression(ax, arr_x, arr_y, showtext=True, addOneToOne=False,
                          extent=None):
    """
    Add a straight line fit to an axis.

    Parameters
    ----------
    ax: matplotlib.pyplot.axes
        The matplotlib axes on which to plot the linear regression.
    arr_x: numpy.array
        The data for the x coordinate.
    arr_y: numpy array
        The data for the y coordinate.
    showtext: bool
        A flag to turn on or off the result of the fit on the plot.
    addOneToOne: bool
        A flag to also add a 1:1 line to the figure
    extent: list of floats
        The extent of the plot axes.
    """

    beta1, beta0, rValue, pValue, stdErr = linregress(arr_x, arr_y)
    texts = [r'$\^\beta_0$ = ' + round_sig(beta0),
             r'$\^\beta_1$ = ' + round_sig(beta1),
             r'R = ' + round_sig(rValue),
             r'P = ' + round_sig(pValue),
             r'N = '+str(int(len(arr_x)))]
    thetext = '\n'.join(texts)

    if showtext:
        pyplot.text(0.04, 0.96, thetext, horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes)

    if extent is None:
        fx = arange(arr_x.min(), arr_x.max(),
                    (arr_x.max() - arr_x.min()) / 20.)
        fy = [beta0 + beta1 * a for a in fx]
    else:
        minv = min(extent)
        maxv = max(extent)
        fx = np.arange(minv, maxv, (maxv - minv)/1000.)
        fy = np.array([beta0 + beta1 * a for a in fx])

        mask = (fx < minv) + (fy < minv) + (fx > maxv) + (fy > maxv)
        fx = np.ma.masked_where(mask, fx)
        fy = np.ma.masked_where(mask, fy)

    pyplot.plot(fx, fy, 'k')
    if addOneToOne:
        pyplot.plot(fx, fx, 'k--')


def make_scatter(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
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
    logger.debug('cubes: %s', ', '.join(cubes.keys()))

    # ####
    # load names:
    model = metadata[filenames['model']]['dataset']
    obs = metadata[filenames['obs']]['dataset']

    long_name = cubes['model'][list(layers.keys())[0]].long_name

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers.keys():

        fig = plt.figure()
        fig.set_size_inches(7, 6)

        # Create the cubes
        model_data = np.ma.array(cubes['model'][layer].data)
        obs_data = np.ma.array(cubes['obs'][layer].data)

        mask = model_data.mask + obs_data.mask
        model_data = np.ma.masked_where(mask, model_data).compressed()
        obs_data = np.ma.masked_where(mask, obs_data).compressed()

        colours = 'gist_yarg'
        zrange = get_array_range([model_data, obs_data])
        plotrange = [zrange[0], zrange[1], zrange[0], zrange[1]]

        hexbin = pyplot.hexbin(model_data,
                               obs_data,
                               # xscale='log',
                               # yscale='log',
                               bins='log',
                               # extent=np.log10(plotrange),
                               gridsize=50,
                               cmap=pyplot.get_cmap(colours),
                               mincnt=0)
        cb = pyplot.colorbar()
        cb.set_label('log10(N)')

        pyplot.gca().set_aspect("equal")
        pyplot.axis(plotrange)

        add_linear_regression(pyplot.gca(),
                              model_data, obs_data,
                              showtext=True,
                              addOneToOne=True,
                              extent=plotrange)

        pyplot.title(long_name)
        pyplot.xlabel(model)
        pyplot.ylabel(obs)

        # Determine image filename:
        fn_list = ['model_vs_obs', long_name, model, obs, str(layer),
                   'scatter']
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
        logger.debug('model_type: %s, %s', index, model_type,)
        logger.debug('metadatas:  %s, %s', index, metadatas,)
        obs_filename = diagtools.match_model_to_key('observational_dataset',
                                                    cfg[model_type],
                                                    metadatas)
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
            make_scatter(
                cfg,
                metadatas,
                filename,
                obs_filename)

            # #####
            # model vs obs map plots
            make_model_vs_obs_plots(
                cfg,
                metadatas,
                filename,
                obs_filename)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
