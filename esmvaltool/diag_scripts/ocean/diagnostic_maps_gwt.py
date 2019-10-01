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

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


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

    # Load cube and set up units
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Making plots for each layer
    qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

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
        suffix = '_'.join(['variable_group_ensembles', metadata['dataset'], key, 'Detrended']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'variable_group_ensembles']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', metadata['dataset'], key, 'Trend_intact']) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'variable_group_ensembles']) + suffix

    # Making plots for each layer
    qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

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
        threshold,
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
    path = diagtools.folder([cfg['plot_dir'], 'threshold_ensemble_plots'])+'Ensemble_mean_'+threshold+image_extention
    if detrend:
        suffix = '_'.join(['variable_group_ensembles', metadata['dataset'], 'Detrended', threshold, ]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Detrended', 'threshold_ensemble']) + suffix
    else:
        suffix = '_'.join(['variable_group_ensembles', metadata['dataset'], 'Trend_intact', threshold, ]) + image_extention
        path = diagtools.folder([cfg['plot_dir'], 'Trend_intact', 'threshold_ensemble']) + suffix

    # Making plots for each layer
    qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

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


def calc_ensemble_mean(cube_list):
    """"
    """
    if len(cube_list)<=1:
        return cube_list[0]
    cube_data = cube_list[0].data
    for c in cube_list[1:]:
        cube_data += c.data
    cube_data = cube_data/float(len(cube_list))
    ensemble_mean = cube_list[0]
    ensemble_mean.data = cube_data
    return ensemble_mean


def make_gwt_map_plots(cfg, detrend = True, do_single_plots=False):
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
    #print('\n', cfg.keys())

    files_dict = {}
    short_names = set()
    ensembles = set()
    exps = set()
    variable_groups = set()
    thresholds = {}

    for fn, details in metadatas.items():
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

    # Calculate the anomaly for each ensemble/threshold combination
    for ensemble in ensembles:
        for variable_group in variable_groups:
            # guess historical group name:
            historical_group = variable_group[:variable_group.find('_')] +'_historical'
            if variable_group == historical_group:
                continue

            print('Plotting:', ensemble, variable_group)
            variable, exp, threshold = split_variable_groups(variable_group)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[(historical_group, ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            cube = diagtools.bgc_units(cube, details['short_name'])

            cube_hist =  iris.load_cube( fn_hist)
            cube_hist = diagtools.bgc_units(cube_hist, details['short_name'])

            cube.data = cube.data - cube_hist.data

            if detrend:
                cube.data = cube.data - cube.data.mean()

            anomaly_cubes[variable_group][ensemble] = cube
            try:
                thresholds[threshold].append([variable_group, ensemble])
            except:
                thresholds[threshold] = [[variable_group, ensemble], ]
            key = variable_group.replace('_',' ') + ' '+ensemble

            # Produce a plot of the anomaly.
            if do_single_plots:
                make_map_plots(cfg, details, cube, key, detrend)

    # Ensemble mean for each variable_group:
    for variable_group in variable_groups:
        # guess historical group name:
        historical_group = variable_group[:variable_group.find('_')] +'historical'
        if variable_group == historical_group:
            continue
        cube_list = []
        for vari, cube in anomaly_cubes[variable_group].items():
            print(variable_group, vari )
            cube_list.append(cube)

        if cube_list == []: continue
        ensemble_mean = calc_ensemble_mean(cube_list)
        make_ensemble_map_plots(cfg, ensemble_mean, variable_group)

    # Ensemble mean for each threshold:
    for threshold, paths in sorted(thresholds.items()):
        cube_list = []
        for [variable_group, ensemble] in paths:
            cube_list.append(anomaly_cubes[variable_group][ensemble])

        if cube_list == []: continue
        ensemble_mean = calc_ensemble_mean(cube_list)

        make_threshold_ensemble_map_plots(cfg, ensemble_mean, threshold, detrend)


def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    for detrend in [True, False]:
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
