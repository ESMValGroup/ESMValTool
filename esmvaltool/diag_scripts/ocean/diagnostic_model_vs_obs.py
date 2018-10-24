"""
Diagnostic Model vs observations maps.

Diagnostic to produce an image showing four maps.
These plost show latitude vs longitude and the cube value is used as the colour
scale.
        model              obs
        model minus obs    model1 over obs

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

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def match_moddel_to_key(model_type, cfg_dict, input_files_dict, ):
    """
    Match up the three models and observations dataset from the configs.

    This function checks that the control_model, exper_model and
    observational_dataset dictionairies from the recipe are matched with the
    input file dictionairy in the cfg metadata.
    """
    for input_file, intput_dict in input_files_dict.items():
        print('\n------------\nmatch_moddel_to_key: intput_dict:', intput_dict)
        print('match_moddel_to_key: cfg_dict:', cfg_dict)
        #try:

        print (intput_dict, cfg_dict)
        intersect_keys = intput_dict.keys() & cfg_dict.keys()
        match = True
        for key in intersect_keys:
            if intput_dict[key] == cfg_dict[key]: continue
            match = False
        if match:
            return input_file
    logger.warning("Unable to match model: %s", model_type)
    return ''


def get_cube_range(cubes):
    """Determinue the minimum and maximum values of an array of cubes."""
    mins = []
    maxs = []
    for cube in cubes:
        mins.append(cube.data.min())
        maxs.append(cube.data.max())
    return [np.min(mins), np.max(maxs), ]


def get_array_range(arrays):
    """Determinue the minimum and maximum values of a list of arrays."""
    mins = []
    maxs = []
    for arr in arrays:
        mins.append(arr.min())
        maxs.append(arr.max())
    return [np.min(mins), np.max(maxs), ]



def get_cube_range_diff(cubes):
    """Determinue the largest deviation from zero in an array of cubes."""
    ranges = []
    for cube in cubes:
        ranges.append(np.abs(cube.data.min()))
        ranges.append(np.abs(cube.data.max()))
    print('get_cube_range_diff:', ranges, -1. * np.max(ranges), np.max(ranges))
    return [-1. * np.max(ranges), np.max(ranges)]


def add_map_subplot(subplot, cube, nspace, title='', cmap='', log=False):
    """Create a map subplot."""
    plt.subplot(subplot)
    print(subplot, cube.data.min(), cube.data.max(), nspace)

    if log:
        #if n_points[0] == 0.1 and n_points[-1] == 10.:
        #locator = ticker.LogLocator(numticks=12,) # subs=(10.,), numdecs=10, numticks=None)
           #print(locator.tick_values(), locator.view_limits())
        #qplt.contourf(cube, 12, linewidth=0, cmap=plt.cm.get_cmap(cmap), locator = locator)
        qplot = qplt.contourf(cube, nspace, linewidth=0, cmap=plt.cm.get_cmap(cmap), norm = LogNorm(), zmin = nspace.min(), zmax=nspace.max())
        #qplot.zmin(nspace.min())
        qplot.colorbar.set_ticks([0.1, 1., 10.])

    else:
        qplot = iris.plot.contourf(cube, nspace, linewidth=0, cmap=plt.cm.get_cmap(cmap), zmin = nspace.min(), zmax=nspace.max())
        cb = pyplot.colorbar(orientation='horizontal')
        ticks = [nspace.min(), (nspace.max() + nspace.min())/2., nspace.max()]
        print (subplot, ticks, nspace)
        cb.set_ticks(ticks)

    #qplot.zmin = nspace.min()
    #qplot.zmax = nspace.max()

    plt.gca().coastlines()
    plt.title(title)


def make_model_vs_obs_plots(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
    Make a model vs obs map plot

    The cfg is the opened global config,
    input_files is the input files dictionairy
    filename is the preprocessing model file.
    """
    filenames = {'model' : model_filename, 'obs': obs_filename}
    print('make_model_vs_obs_plots:', filenames)
    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        print('loading:', model_type, input_file)
        cube = iris.load_cube(input_file)
        #cube = diagtools.bgc_units(cube, metadatas[input_file]['short_name'])
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
        print (cube221)
        print (cube222)
        print ('\nmodel lat:',cube221.coord('latitude').coord_system, '\nobs lat:',cube222.coord('latitude').coord_system)
        #print (cube221.coord('longitude'), cube222.coord('longitude'))
        cube223 = cubes['model'][layer] - cubes['obs'][layer]
        cube224 = cubes['model'][layer] / cubes['obs'][layer]
        # cube223 = cubes['model'][layer].data - cubes['obs'][layer].data
        # cube224 = cubes['model'][layer].data / cubes['obs'][layer].data

        # create the z axis for plots 2, 3, 4.
        zrange12 = get_cube_range([cube221, cube222])
        zrange3 = get_cube_range_diff([cube223])

        cube224.data = np.ma.clip(cube224.data, 0.1, 10.)
        zrange4 = [0.1, 10.]

        n_points = 12
        linspace12 = np.linspace(zrange12[0], zrange12[1], n_points, endpoint=True)
        linspace3 = np.linspace(zrange3[0], zrange3[1], n_points, endpoint=True)
        logspace4 = np.logspace(-1., 1., 12, endpoint=True)
        #logspace4 = np.linspace(zrange4[0], zrange4[1], 21, endpoint=True)

        print('linspace3:', linspace3)
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
        fn_list = [long_name, model, obs, str(layer), 'maps']
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def round_sig(x, sig=3):
	"""
	:param x: a float
	:param sig: number of significant figures

	rounds a value to a specific number of significant figures.
	"""
	if x ==0.:	return str(0.)
	if x<0.:	return str(-1.* round(abs(x), sig-int(math.floor(math.log10(abs(x))))-1)	)
	else: 		return str(    round(x, sig-int(math.floor(math.log10(x)))-1))


def add_linear_regression(ax, x, y,showtext=True, addOneToOne=False, extent = [0,0,0,0]):
    """
    Adds a straight line fit to an axis.
    """
    def getLinRegText(ax, x, y, showtext=True):
        x = [a for a in x if (a is np.ma.masked)==False]
        y = [a for a in y if (a is np.ma.masked)==False]
        beta1, beta0, rValue, pValue, stdErr = linregress(x, y)
        thetext = r'$\^\beta_0$ = '+round_sig(beta0)        \
            + '\n'+r'$\^\beta_1$ = '+round_sig(beta1)    \
            + '\nR = '+ round_sig(rValue)        \
            + '\nP = '+round_sig(pValue)        \
            + '\nN = '+str(int(len(x)))
            #+ '\n'+r'$\epsilon$ = ' + round_sig(stdErr)    \
        if showtext: pyplot.text(0.04, 0.96,thetext ,
                     horizontalalignment='left',
                     verticalalignment='top',
                     transform = ax.transAxes)
        return beta1, beta0, rValue, pValue, stdErr

    b1, b0, rValue, pValue, stdErr = getLinRegText(ax, x, y, showtext =showtext)
    if extent == [0,0,0,0]:
        fx = arange(x.min(), x.max(), (x.max()-x.min())/20.)
        fy =[b0 + b1*a for a in fx]
    else:
        minv = min(extent)
        maxv = max(extent)
        fx = np.arange(minv, maxv, (maxv-minv)/1000.)
        fy = np.array([b0 + b1*a for a in fx])

        fx = np.ma.masked_where((fx<minv) + (fy < minv) + (fx>maxv) + (fy > maxv), fx)
        fy = np.ma.masked_where((fx<minv) + (fy < minv) + (fx>maxv) + (fy > maxv), fy)

    pyplot.plot(fx,fy, 'k')
    if addOneToOne: pyplot.plot(fx,fx, 'k--')


def make_scatter(
        cfg,
        metadata,
        model_filename,
        obs_filename):
    """
    Make a model vs obs scatter plot

    The cfg is the opened global config,
    input_files is the input files dictionairy
    filename is the preprocessing model file.
    """
    filenames = {'model' : model_filename, 'obs': obs_filename}
    print('make_model_vs_obs_plots:', filenames)
    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for model_type, input_file in filenames.items():
        print('loading:', model_type, input_file)
        cube = iris.load_cube(input_file)
        #cube = diagtools.bgc_units(cube, metadatas[input_file]['short_name'])
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

        colours = 'gist_yarg' # 'Greens'
        zrange = get_array_range([model_data, obs_data])
        print(model_data.shape, obs_data.shape)
        plotrange = [zrange[0], zrange[1], zrange[0], zrange[1]]


        hexbin = pyplot.hexbin(model_data,
                               obs_data,
                               #xscale='log',
                               #yscale='log',
                               bins='log',
                               #extent=np.log10(plotrange),
                               gridsize = 50,
                               cmap=pyplot.get_cmap(colours),
                               mincnt=0)
        cb = pyplot.colorbar()
        cb.set_label('log10(N)')
        #cb = pyplot.colorbar(ticks=[0, 1, 2, 3, 4, 5, 6, ],)
        #cb.set_ticklabels([r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$',])

        pyplot.gca().set_aspect("equal")
        pyplot.axis(plotrange)

        add_linear_regression(pyplot.gca(), model_data, obs_data, showtext =True, addOneToOne=True, extent=plotrange)

        pyplot.title(long_name)
        pyplot.xlabel(model)
        pyplot.ylabel(obs)

        # Determine image filename:
        fn_list = [long_name, model, obs, str(layer), 'scatter']
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

    The cfg is the opened global config.
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
        #print(cfg.keys(), '\n', metadatas[].keys())
        #logger.debug('cfg[model_type]: %s, %s', index, cfg[model_type])
        logger.debug('metadatas:  %s, %s', index, metadatas,)
        obs_filename = match_moddel_to_key('observational_dataset',
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

            ######
            # model vs obs scatter plots
            make_scatter(
                cfg,
                metadatas,
                filename,
                obs_filename)

            ######
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
