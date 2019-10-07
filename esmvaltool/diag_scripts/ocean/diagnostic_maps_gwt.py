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
import logging
import os
import sys
from itertools import product
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt
import cartopy
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def detrended_contourf(cube):
    """
    make the detrended contouf plot
    """
    cmap = 'RdBu_r'
    drange = diagtools.get_cube_range_diff([cube, ])
    dlinspace = np.linspace(drange[0], drange[1], 22, endpoint=True)
    ticks = [t for t in np.linspace(dlinspace.min(),dlinspace.max(), 7)]
    try:
        qplot = qplt.contourf(cube, dlinspace, cmap=cmap) # linewidth=0, rasterized=True,
        qplot.colorbar.set_ticks(ticks)
    except:
        print('Unable to plot cube:', cube)
        return False
    return True


def trended_contourf(cube, cmap='YlOrRd'):
    """
    make the detrended contouf plot
    """
    try:
        qplt.contourf(cube, 12, cmap=cmap) # linewidth=0, rasterized=True,
    except:
        print('Unable to plot cube:', cube)
        return False
    return True


def split_variable_groups(variable_group):
    """
    Split variable group into variable and experiment.
    """
    variable, exp, threshold = variable_group.split('_')
    if variable == 'tas':
        variable = 'Surface Temperature'
    exp = exp.upper()
    exp = ' '.join([exp[:3], exp[3], exp[4:]])
    if threshold == '15':
        threshold = '1.5'
    threshold += u'\N{DEGREE SIGN}'
    return variable, exp, threshold


def calc_ensemble_mean(cube_list):
    """
    Calculate the ensemble mean of a list of cubes.

    """
    if len(cube_list) ==1:
        return cube_list[0]

    cube_data = cube_list[0].data
    for c in cube_list[1:]:
        cube_data += c.data

    cube_data = cube_data/float(len(cube_list))
    ensemble_mean = cube_list[0]
    ensemble_mean.data = cube_data
    return ensemble_mean


def make_map_plots(
        cfg,
        metadata,
        cube,
        key,
        detrend,
        cmap='YlOrRd'
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
    detrend: bool
        key to show wherether the mean was subtracted.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join([metadata['dataset'], key, 'detrended'])+image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'single_plots'])+ suffix
    else:
        suffix = '_'.join([metadata['dataset'], key])+image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'single_plots'])+suffix

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube, cmap=cmap)
    if not success:
        print('Failed to make figure:', path)
        print('key:', key)
        plt.close()
        return

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # Add title to plot
    if detrend:
        title = ' '.join([metadata['dataset'], key, '- detrended'])
    else:
        title = ' '.join([metadata['dataset'], key, '- trend intact'])

    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def add_map_subplot(subplot, cube, nspace, title='', cmap=''):
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

    """
    plt.subplot(subplot)
    qplot = qplt.contourf(cube, nspace, linewidth=0,
                          cmap=plt.cm.get_cmap(cmap))
    qplot.colorbar.set_ticks([nspace.min(),
                              (nspace.max() + nspace.min()) / 2.,
                              nspace.max()])

    try: plt.gca().coastlines()
    except: pass
    plt.title(title)


def weighted_mean(cube):
    """
    Calculate the weighted mean.
    """
    return cube.data.mean()


def make_four_pane_map_plot(
        cfg,
        cube_ssp,
        cube_hist,
        cube_anomaly,
        cube_detrended,
        key,
        plot_type,
):
    """
    Make a four pane map plot for an individual.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.
    detrend: bool
        key to show wherether the mean was subtracted.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    suffix = '_'.join([plot_type, key]) + image_extention
    path = diagtools.folder([cfg['plot_dir'], 'four_pane_plots', plot_type]) + suffix

    fig = plt.figure()
    fig.set_size_inches(9, 6)

    # Create the cubes
    cube221 = cube_ssp
    cube222 = cube_hist
    cube223 = cube_anomaly
    cube224 = cube_detrended

    # create the z axis for plots 2, 3, 4.
    zrange1 = diagtools.get_cube_range([cube221, cube222 ])
    zrange3 = diagtools.get_cube_range_diff([cube223, ])
    zrange4 = diagtools.get_cube_range_diff([cube224, ])

    linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
    linspace3 = np.linspace(zrange3[0], zrange3[1], 12, endpoint=True)
    linspace4 = np.linspace(zrange4[0], zrange4[1], 12, endpoint=True)

    # Add the sub plots to the figure.
    add_map_subplot(221, cube221, linspace1, cmap='viridis',
                    title='SSP')
    add_map_subplot(222, cube222, linspace1, cmap='viridis',
                    title='Historical')
    add_map_subplot(223, cube223, linspace3, cmap='RdBu_r',
                    title='SSP minus historical')
    add_map_subplot(224, cube224, linspace4, cmap='RdBu_r',
                    title='Detrended SSP minus historical')

    # Add overall title
    fig.suptitle(key, fontsize=14)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_ensemble_map_plots(
        cfg,
        cube,
        variable_group,
        detrend,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for each variable_group.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    # Determine image filename:
    if detrend:
        suffix = '_'.join(['variable_group_ensembles', variable_group, 'Detrended']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'variable_group_ensembles']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', variable_group, 'Trend_intact']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'variable_group_ensembles']) + suffix

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube, cmap=cmap)
    if not success:
        print('Unable to make figure:', path)
        print('variable_group:', variable_group)
        plt.close()
        return

    #plt.gca().coastlines()

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # tas_ssp585_6
    variable, exp, threshold = split_variable_groups(variable_group)

    # Add title to plot
    if detrend:
        title = ' '.join([variable, '- ensemble mean of', exp, 'after', threshold, 'warming - detrended'])
    else:
        title = ' '.join([variable, '- ensemble mean of', exp, 'after', threshold, 'warming - trend intact'])
    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_threshold_ensemble_map_plots(
        cfg,
        cube,
        variable,
        threshold,
        detrend,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if detrend:
        suffix = '_'.join(['variable_group_ensembles', 'Detrended', threshold, variable]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'threshold_ensemble']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', 'Trend_intact', threshold, variable]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'threshold_ensemble']) + suffix

    # Making plots for each layer
    #try: qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)
    #except:
    #    print('Trying to make figure:', path)
    #    print('variable_group:', variable, threshold)
    #    print('Unable to plot cube:', cube)
    #    plt.close()
    #    return
    #    #qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

    # Making plots
    if detrend:
        success = detrended_contourf(cube)
    else:
        success = trended_contourf(cube,cmap=cmap)
    if not success:
        print('Unable to make figure:', path)
        print('variable_group:',threshold, variable)
        plt.close()
        return

    # plt.gca().coastlines()
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # Add title to plot
    if detrend:
        title = ' '.join(['Ensemble mean after', threshold, 'warming - detrended'])
    else:
        title = ' '.join(['Ensemble mean after', threshold, 'warming - trend intact'])

    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)
    plt.close()


def make_gwt_map_four_plots(cfg, ):
    """
    Make plots

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    do_single_plots: bool
        Make individual plots for each dataset.

    """
    metadatas = diagtools.get_input_files(cfg)
    do_single_plots=False
    do_variable_group_plots=True
    do_threshold_plots=True

    #print('\n', cfg.keys())

    files_dict = {}
    short_names = set()
    ensembles = set()
    exps = set()
    variable_groups = set()
    variables = set()
    thresholds = {}

    for fn, details in sorted(metadatas.items()):
        #print(fn, details.keys())
        short_names.add(details['short_name'])
        ensembles.add(details['ensemble'])
        exps.add(details['exp'])
        variable_groups.add(details['variable_group'])

        unique_key = (details['variable_group'], details['ensemble'])
        try:
            files_dict[unique_key].append(fn)
        except:
            files_dict[unique_key] = [fn, ]

    if len(short_names) ==1:
        short_name = list(short_names)[0]
    print(files_dict.keys())


    # lets look at  minus the historical
    hist_cubes = {variable_group:{} for variable_group in variable_groups}
    ssp_cubes = {variable_group:{} for variable_group in variable_groups}
    anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    detrended_anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    viable_keys = {}

    hist_cubes_fns = {}
    # Calculate the anomaly for each ensemble/threshold combination
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):
            # guess historical group name:
            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                continue

            print('Plotting:', ensemble, variable_group)
            variable, exp, threshold = split_variable_groups(variable_group)
            variables.add(variable)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[(historical_group, ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            cube = diagtools.bgc_units(cube, details['short_name'])

            # Check is historial cube is loaded already
            if fn_hist in hist_cubes_fns:
                cube_hist = hist_cubes_fns[fn_hist]
            else:
                cube_hist =  iris.load_cube( fn_hist)
                cube_hist = diagtools.bgc_units(cube_hist, details['short_name'])

            cube_anomaly = cube - cube_hist
            detrended_cube_anomaly = cube_anomaly - weighted_mean(cube_anomaly)

            ssp_cubes[variable_group][ensemble] = cube
            hist_cubes[variable_group][ensemble] = cube_hist
            anomaly_cubes[variable_group][ensemble] = cube_anomaly
            detrended_anomaly_cubes[variable_group][ensemble] = detrended_cube_anomaly
            viable_keys[(variable_group, ensemble)] = True
            try:
                thresholds[threshold].append([variable_group, ensemble])
            except:
                thresholds[threshold] = [[variable_group, ensemble], ]

            key = variable_group.replace('_',' ') + ' '+ensemble


    # single plots:
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):
            if not do_single_plots:
                continue
            if (variable_group, ensemble) not in viable_keys:
                print("Not in viable keys:", ensemble, variable_group)
                continue

            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                continue

            cube_ssp = ssp_cubes[variable_group][ensemble]
            cube_hist = hist_cubes[variable_group][ensemble]
            cube_anomaly = anomaly_cubes[variable_group][ensemble]
            cube_detrend = detrended_anomaly_cubes[variable_group][ensemble]
            key = ' '.join([short_name, variable_group, ensemble])

            make_four_pane_map_plot(
                cfg,
                cube_ssp,
                cube_hist,
                cube_anomaly,
                cube_detrend,
                key,
                'single_plots',
            )

    # Ensemble mean for each variable_group:
    for variable_group in sorted(variable_groups):
        # check to make plots.
        if not do_variable_group_plots:
            continue

        # guess historical group name:
        historical_group = variable_group[:variable_group.find('_')] + '_historical'
        if variable_group == historical_group:
            continue

        # Calculate ensemble means
        output_cubes = []
        for cubes in [ssp_cubes, hist_cubes, anomaly_cubes, detrended_anomaly_cubes]:
            cube_list = []
            for vari, cube in sorted(cubes[variable_group].items()):
                print(variable_group, vari)
                cube_list.append(cube)
            if cube_list == []:
                continue
            ensemble_mean = calc_ensemble_mean(cube_list)
            output_cubes.append(ensemble_mean)
        key = ' '.join([short_name, variable_group])

        if len(output_cubes) != 4:
            print('Problem with ', variable_group, 'four plots.', historical_group)
            assert 0

        make_four_pane_map_plot(
            cfg,
            output_cubes[0],
            output_cubes[1],
            output_cubes[2],
            output_cubes[3],
            key,
            'variable_group_ensembles',
            )

    # Ensemble mean for each threshold:
    for threshold, paths in sorted(thresholds.items()):
        if not do_threshold_plots: continue
        output_cubes = []
        for cubes in [ssp_cubes, hist_cubes, anomaly_cubes, detrended_anomaly_cubes]:
            cube_list = []
            for (variable_group, ensemble) in viable_keys.keys():
                var_g, exp, thres = split_variable_groups(variable_group)
                if thres != threshold:
                    continue
                cube_list.append(cubes[variable_group][ensemble])
            ensemble_mean = calc_ensemble_mean(cube_list)
            output_cubes.append(ensemble_mean)

        key = ' '.join([short_name, threshold])
        make_four_pane_map_plot(
            cfg,
            output_cubes[0],
            output_cubes[1],
            output_cubes[2],
            output_cubes[3],
            key,
            'threshold_ensemble',
            )




def make_gwt_map_plots(cfg, detrend = True,):
    """
    Make plots

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    detrend: bool
        Subtract the mean from the figure.
    do_single_plots: bool
        Make individual plots for each dataset.

    """
    metadatas = diagtools.get_input_files(cfg)
    do_single_plots=True
    do_variable_group_plots=True
    do_threshold_plots=True

    #print('\n', cfg.keys())

    files_dict = {}
    short_names = set()
    ensembles = set()
    exps = set()
    variable_groups = set()
    variables = set()
    thresholds = {}

    for fn, details in sorted(metadatas.items()):
        #print(fn, details.keys())
        short_names.add(details['short_name'])
        ensembles.add(details['ensemble'])
        exps.add(details['exp'])
        variable_groups.add(details['variable_group'])

        unique_key = (details['variable_group'], details['ensemble'])
        try:
            files_dict[unique_key].append(fn)
        except:
            files_dict[unique_key] = [fn, ]

    print(files_dict.keys())

    # lets look at  minus the historical
    anomaly_cubes = {variable_group:{} for variable_group in variable_groups}
    hist_cubes = {}
    # Calculate the anomaly for each ensemble/threshold combination
    for ensemble in sorted(ensembles):
        for variable_group in sorted(variable_groups):
            # guess historical group name:
            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                continue

            print('Plotting:', ensemble, variable_group)
            variable, exp, threshold = split_variable_groups(variable_group)
            variables.add(variable)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[(historical_group, ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            cube = diagtools.bgc_units(cube, details['short_name'])

            # Check is historial cube is loaded already
            if fn_hist in hist_cubes:
                cube_hist = hist_cubes[fn_hist]
            else:
                cube_hist =  iris.load_cube( fn_hist)
                hist_cubes[fn_hist] = cube_hist
                cube_hist = diagtools.bgc_units(cube_hist, details['short_name'])

            cube_anomaly = cube - cube_hist

            if detrend:
                cube_anomaly = cube_anomaly - weighted_mean(cube_anomaly)

            anomaly_cubes[variable_group][ensemble] = cube_anomaly
            try:
                thresholds[threshold].append([variable_group, ensemble])
            except:
                thresholds[threshold] = [[variable_group, ensemble], ]
            key = variable_group.replace('_',' ') + ' '+ensemble

            # Produce a plot of the anomaly.
            if do_single_plots:
                make_map_plots(cfg, details, cube_anomaly, key, detrend)

    historical_ensemble_mean={}
    # Ensemble mean for each variable_group:
    for variable_group in sorted(variable_groups):
        # check to make plots.
        if not do_variable_group_plots:
            continue

        # guess historical group name:
        historical_group = variable_group[:variable_group.find('_')] + '_historical'
        if variable_group == historical_group:
            continue

        cube_list = []
        for vari, cube in sorted(anomaly_cubes[variable_group].items()):
            print(variable_group, vari)
            cube_list.append(cube)

        if cube_list == []:
            continue
        ensemble_mean = calc_ensemble_mean(cube_list)
        make_ensemble_map_plots(cfg, ensemble_mean, variable_group, detrend)

    # Ensemble mean for each threshold:
    for variable in sorted(variables):
        if not do_threshold_plots: continue
        for threshold, paths in sorted(thresholds.items()):
            cube_list = []
            for [variable_group, ensemble] in sorted(paths):
                var_g, exp, threshold = split_variable_groups(variable_group)
                if var_g != variable:
                    continue
                cube_list.append(anomaly_cubes[variable_group][ensemble])

            if cube_list == []:
                continue
            ensemble_mean = calc_ensemble_mean(cube_list)
            make_threshold_ensemble_map_plots(cfg, ensemble_mean, variable, threshold, detrend, )


def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    make_gwt_map_four_plots(cfg)

    for detrend in [True, False, ]:
        make_gwt_map_plots(cfg, detrend)


    # for index, metadata_filename in enumerate(cfg['input_files']):
    #     logger.info(
    #         'metadata filename:\t%s',
    #         metadata_filename,
    #     )
    #
    #     metadatas = diagtools.get_input_files(cfg, index=index)
    #     #thresholds = diagtools.load_thresholds(cfg, metadatas)
    #
    #
    #     # for filename in sorted(metadatas.keys()):
    #     #
    #     #     logger.info('-----------------')
    #     #     logger.info(
    #     #         'model filenames:\t%s',
    #     #         filename,
    #     #     )
    #     #
    #     #     ######
    #     #     # Contour maps of individual model
    #     #     if thresholds:
    #     #         make_map_contour(cfg, metadatas[filename], filename)
    #     #
    #     #     ######
        #     # Maps of individual model
        #     make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
