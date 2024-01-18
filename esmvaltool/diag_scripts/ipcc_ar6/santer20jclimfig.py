#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Plot based on Santer et al. (2020).

###############################################################################
santer20jclim/santer20jclimfig.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP ans 4C project
###############################################################################

Description
-----------
    Water Vapor Path trends following Santer et al. (2020).

Configuration options
---------------------
sample_obs: optional, filter all data sets (netCDF file with 0 and 1 for
            used grid). The data sets must be interpolated to the same lat/lon
            grid as the filter.
            The filter must cover at least the time period used for the data.

###############################################################################

"""

import logging
import os
from collections import OrderedDict
from pprint import pformat
from cf_units import Unit
import iris
import iris.coord_categorisation as cat
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
# from scipy.stats import norm
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            plot, get_diagnostic_filename,
                                            get_plot_filename,
                                            group_metadata,
                                            select_metadata, run_diagnostic,
                                            variables_available,
                                            extract_variables)

logger = logging.getLogger(os.path.basename(__file__))


def _apply_filter(cfg, cube):
    """Apply filter from RSS Anomalies to all data and calculates mean."""
    if 'sample_obs' in cfg:
        filt = iris.load(os.path.join(cfg['auxiliary_data_dir'],
                                      cfg['sample_obs']))[0]

        for iii, dim_str in enumerate(['time', 'latitude', 'longitude']):
            filt = _fit_dim(filt, cube, dim_str, iii)

        cube.data = cube.data * filt.data

    coords = ('longitude', 'latitude')
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    cube_grid_areas = iris.analysis.cartography.area_weights(cube)
    cube = (cube.collapsed(coords, iris.analysis.MEAN,
                           weights=cube_grid_areas))

    return cube


def _fit_dim(filt, cube, dim_str, dim_nr):
    """Adjust array for given dimension."""
    dimfil = filt.coord(dim_str)
    dimcu = cube.coord(dim_str)
    if dim_str == 'time':
        start = dimcu.units.num2date(dimcu.points[0])
        end = dimcu.units.num2date(dimcu.points[-1])
        sdim = dimfil.nearest_neighbour_index(dimfil.units.date2num(start))
        edim = dimfil.nearest_neighbour_index(dimfil.units.date2num(end))
    else:
        start = dimcu.points[0]
        end = dimcu.points[-1]
        sdim = dimfil.nearest_neighbour_index(start)
        edim = dimfil.nearest_neighbour_index(end)

    if dim_nr == 0:
        filt = filt[sdim:edim + 1, :, :]
    elif dim_nr == 1:
        filt = filt[:, sdim:edim + 1, :]
    elif dim_nr == 2:
        filt = filt[:, :, sdim:edim + 1]

    return filt


def _calc_trend(cube_anom):
    """Calculate trends."""
    reg_var = stats.linregress(np.linspace(0, len(cube_anom.data) - 1,
                                           len(cube_anom.data)),
                               cube_anom.data)

    return reg_var.slope


def _calculate_anomalies(cube):
    """Remove annual cycle from time series."""
    c_n = cube.aggregated_by('month_number', iris.analysis.MEAN)
    month_number = cube.coord('month_number').points
    startind = np.where(month_number == 1)
    endind = np.where(month_number == 12)
    data_in = cube.data
    data_new = np.full(data_in.shape, 0.5)
    for iii in range((startind[0])[0], (endind[0])[-1], 12):
        data_new[iii:iii + 12] = ((cube.data[iii:iii + 12] - c_n.data) /
                                  c_n.data) * 100.0 * 120.0

    cube.data = data_new
    cube.units = Unit('percent')

    return cube


def _check_full_data(dataset_path, cube):
    """Check if cube covers time series from start year to end year."""
    # cstart_month = cube.coord('month_number').points[0]
    # cend_month = cube.coord('month_number').points[1]
    cstart_year = cube.coord('year').points[0]
    cend_year = cube.coord('year').points[-1]
    start_year = int(dataset_path['start_year'])
    end_year = int(dataset_path['end_year'])

    check = 0
    if start_year == cstart_year and end_year == cend_year:
        check = 1
        print("Full time series:")
        print(start_year, cstart_year, end_year, cend_year)
    else:
        print("FAILED:")
        print(start_year, cstart_year, end_year, cend_year)

    return check


def _get_hem_letter(lat):
    """Get S or N from Latitude."""
    shem = ''
    if lat < 0:
        shem = 'S'
    if lat > 0:
        shem = 'N'

    return shem


def _get_valid_datasets(input_data):
    """Get valid datasets list and number for each model."""
    number_of_subdata = OrderedDict()
    available_dataset = list(group_metadata(input_data, 'dataset'))
    valid_datasets = []
    period = {}
    period['start_year'] = []
    period['end_year'] = []
    period['span'] = []
    period['slat'] = []
    period['elat'] = []
    for dataset in available_dataset:
        meta = select_metadata(input_data, dataset=dataset)
        number_of_subdata[dataset] = float(len(meta))
        for dataset_path in meta:
            cube = iris.load(dataset_path['filename'])[0]
            cat.add_year(cube, 'time', name='year')
            if not _check_full_data(dataset_path, cube):
                number_of_subdata[dataset] = number_of_subdata[dataset] - 1
            else:
                valid_datasets.append(dataset_path['filename'])
                period['start_year'].append(int(dataset_path['start_year']))
                period['end_year'].append(int(dataset_path['end_year']))
                period['span'].append(int(dataset_path['end_year']) -
                                      int(dataset_path['start_year']) + 1)
                period['slat'].append(round(cube.coord('latitude').points[0]))
                period['elat'].append(round(cube.coord('latitude').points[-1]))
    if min(period['span']) == max(period['span']):
        period['common_span'] = str(min(period['span']))
    if min(period['start_year']) == max(period['start_year']) and\
       min(period['end_year']) == max(period['end_year']):
        period['common_period'] = str(min(period['start_year'])) + ' - ' +\
            str(min(period['end_year']))

    if min(period['slat']) == max(period['slat']) and\
       min(period['elat']) == max(period['elat']):
        shem = _get_hem_letter(min(period['slat']))
        ehem = _get_hem_letter(min(period['elat']))

        period['common_lats'] = str(abs(min(period['slat']))) + '°' + shem +\
            ' - ' + str(abs(min(period['elat']))) + '°' + ehem

    return valid_datasets, number_of_subdata, period


def _plot_extratrends(cfg, extratrends, trends, period):
    """Plot trends for ensembles."""
    if 'histmin' in cfg.keys():
        histmin = cfg['histmin']
    else:
       histmin = 0

    if 'histmax' in cfg.keys():
        histmax = cfg['histmax']
    else:
        histmax = 4
    
    if 'ymax' in cfg.keys():
        ymax = cfg['ymax']
    else:
        ymax = 2.5

    res_ar = {'xhist': np.linspace(histmin, histmax, 41),
              'artrend': {},
              'kde1': {}}
    # xhist = np.linspace(0, 4, 41)
    # artrend = {}
    # kde1 = {}
    fig, axx = plt.subplots(figsize=(8, 6))

    valid_datasets = []

    for xtrmdl in cfg['add_model_dist']:
        alias = list(extratrends[xtrmdl].keys())[0]
        valid_datasets.append(select_metadata(cfg['input_data'].values(),
                                              alias=alias)[0]['filename'])
        if alias in trends['cmip6'].keys():
            style = plot.get_dataset_style(xtrmdl, style_file='cmip6')
        elif alias in trends['cmip5'].keys():
            style = plot.get_dataset_style(xtrmdl, style_file='cmip5')
        else:
            style = {'facecolor': (0, 0, 1, 0.2),
                     'color': (0, 0, 1, 1.0)}

        res_ar['artrend'][xtrmdl] = np.fromiter(extratrends[xtrmdl].values(),
                                                dtype=float)
        res_ar['kde1'][xtrmdl] = stats.gaussian_kde(res_ar['artrend'][xtrmdl],
                                                    bw_method="scott")
        axx.hist(res_ar['artrend'][xtrmdl], bins=res_ar['xhist'], density=True,
                 edgecolor=style['color'],
                 facecolor=style['facecolor'])
        axx.plot(res_ar['xhist'], res_ar['kde1'][xtrmdl](res_ar['xhist']),
                 color=style['color'],
                 linewidth=3,
                 label=xtrmdl)

    caption = _plot_obs(trends, axx, ymax)
    caption = _plot_settings(cfg, axx, fig, 'fig2', period) + caption
    plt.close()

    caption = 'Probability density function of the decadal trend in ' + \
        'the Water Vapor Path' + caption

    provenance_record = get_provenance_record(valid_datasets,
                                              caption, ['trend', 'other'],
                                              ['reg'])

    diagnostic_file = get_diagnostic_filename('fig2', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_vars(_write_list_dict(cfg,
                                                 trends,
                                                 res_ar)),
              target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def _plot_obs(trends, axx, maxh):
    """Plot observational or reanalysis data as vertical line."""
    obs_str = ''
    # IPCC colors for obs from
    # https://github.com/IPCC-WG1/colormaps/blob/master/
    # categorical_colors_rgb_0-255/dark_cat.txt
    if trends['obs']:
        obs_str = ' Vertical lines show the trend for'
        for iii, obsname in enumerate(trends['obs'].keys()):
            obscoli = float(iii)
            if iii > 4:
                obscoli = float(iii) - 4.5
            if iii > 8:
                obscoli = float(iii) - 8.75
            if iii > 12:
                obscoli = float(iii) - 12.25
            plotobscol = (1.0 - 0.25 * obscoli, 0.25 * obscoli,
                          0.5 + obscoli * 0.1)
            if obsname == 'RSS':
                plotobscol = (221 / 255.0, 84 / 255.0, 46 / 255.0,)
                # obsnamep = 'RSS'
            if obsname == 'CRU':
                plotobscol = (0.8, 0.4, 0.1, 1)  
            if obsname == 'ERA5':
                # obscol = [(0.8, 0.4, 0.1, 1), (1, 0, 0.2, 1), (0.8, 0.1, 0.4, 1), (1, 0.6, 0, 1)]
                # plotobscol = (1, 0, 0.2, 1)
                plotobscol = (128 / 255.0, 54 / 255.0, 168 / 255.0,)
                # obsnamep = 'ERA5.1'
            if obsname == 'MERRA2':
                plotobscol = (0.8, 0.1, 0.4, 1)            
            if obsname == 'NCEP-NCAR-R1':
                plotobscol = (1, 0.6, 0, 1)
            axx.vlines(trends['obs'][obsname], 0, maxh,
                       colors=plotobscol,
                       linewidth=3,
                       label=obsname)
        obs_str = obs_str + '.'
    return obs_str


def _plot_trends(cfg, trends, valid_datasets, period):
    """Plot probability density function of trends."""
    if 'histmin' in cfg.keys():
        histmin = cfg['histmin']
    else:
        histmin = 1000
        if trends['cmip5']:
            histmin = np.min(np.array([histmin, np.min(np.fromiter(trends['cmip5'].values(),
                                           dtype=float))]))
        if trends['cmip6']:
            histmin = np.min(np.array([histmin, np.min(np.fromiter(trends['cmip6'].values(),
                                           dtype=float))]))
        if trends['obs']:
            histmin = np.min(np.array([histmin, np.min(np.fromiter(trends['obs'].values(),
                                           dtype=float))]))
        
    if 'histmax' in cfg.keys():
        histmax = cfg['histmax']
    else:
        histmax = 4

    if 'ymax' in cfg.keys():
        ymax = cfg['ymax']
    else:
        ymax = 2.5
    res_ar = {'xhist': np.linspace(histmin, histmax, 41)}
    fig, axx = plt.subplots(figsize=(8, 6))

    # IPCC colors for CMIP5 and CMIP6 from
    # https://github.com/IPCC-WG1/colormaps/blob/master/
    # categorical_colors_rgb_0-255/cmip_cat.txt
    # CMIP5
    if trends['cmip5']:
        res_ar['artrend_c5'] = np.fromiter(trends['cmip5'].values(),
                                           dtype=float)
        res_ar['weights_c5'] = np.fromiter(trends['cmip5weights'].values(),
                                           dtype=float)
        # yyy_c5, xxx_c5 = np.histogram(artrend_c5, bins=xhist)
        res_ar['kde1_c5'] = stats.gaussian_kde(res_ar['artrend_c5'],
                                               weights=res_ar['weights_c5'],
                                               bw_method="scott")
        axx.hist(res_ar['artrend_c5'], bins=res_ar['xhist'], density=True,
                 weights=res_ar['weights_c5'],
                 edgecolor=(37 / 250.0, 81 / 255.0, 204 / 255.0, 1.0),
                 facecolor=(37 / 250.0, 81 / 255.0, 204 / 255.0, 0.2))
        axx.plot(res_ar['xhist'], res_ar['kde1_c5'](res_ar['xhist']),
                 color=(37 / 250.0, 81 / 255.0, 204 / 255.0, 1),
                 linewidth=3,
                 label="CMIP5")
        # maxh = np.max(kde1_c5(xhist))

    # CMIP6
    if trends['cmip6']:
        res_ar['artrend_c6'] = np.fromiter(trends['cmip6'].values(),
                                           dtype=float)
        res_ar['weights_c6'] = np.fromiter(trends['cmip6weights'].values(),
                                           dtype=float)
        # yyy_c6, xxx_c6 = np.histogram(artrend_c6, bins=xhist)
        res_ar['kde1_c6'] = stats.gaussian_kde(res_ar['artrend_c6'],
                                               weights=res_ar['weights_c6'],
                                               bw_method="scott")
        axx.hist(res_ar['artrend_c6'], bins=res_ar['xhist'], density=True,
                 weights=res_ar['weights_c6'],
                 edgecolor=(204 / 250.0, 35 / 255.0, 35 / 255.0, 1.0),
                 facecolor=(204 / 250.0, 35 / 255.0, 35 / 255.0, 0.2))
        axx.plot(res_ar['xhist'], res_ar['kde1_c6'](res_ar['xhist']),
                 color=(204 / 250.0, 35 / 255.0, 35 / 255.0, 1),
                 linewidth=3,
                 label="CMIP6")
        # maxh = np.max(kde1_c6(xhist))
    caption = _plot_obs(trends, axx, ymax)
    caption = _plot_settings(cfg, axx, fig, 'fig1', period) + caption
    plt.close()

    caption = 'Probability density function of the decadal trend in ' + \
        'the Water Vapor Path' + caption

    provenance_record = get_provenance_record(valid_datasets,
                                              caption, ['trend', 'other'],
                                              ['reg'])

    diagnostic_file = get_diagnostic_filename('fig1', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = [res_ar['xhist']]
    list_dict["name"] = [{'var_name': 'prw_trends_bins',
                          'long_name': 'Water Vapor Path Trend bins',
                          'units': 'percent'}]
    if trends['cmip5']:
        list_dict["data"].append(res_ar['artrend_c5'])
        list_dict["name"].append({'var_name': 'prw_trends_cmip5',
                                  'long_name': 'Water Vapor Path Trends CMIP5',
                                  'units': 'percent'})
        list_dict["data"].append(res_ar['weights_c5'])
        list_dict["name"].append({'var_name': 'data_set_weights',
                                  'long_name': 'Weights for each data set.',
                                  'units': '1'})
        list_dict["data"].append(res_ar['kde1_c5'](res_ar['xhist']))
        list_dict["name"].append({'var_name': 'prw_trend_distribution_cmip5',
                                  'long_name': 'Water Vapor Path Trends ' +
                                               'distribution CMIP5',
                                  'units': '1'})
    if trends['cmip6']:
        list_dict["data"].append(res_ar['artrend_c6'])
        list_dict["name"].append({'var_name': 'prw_trends_cmip6',
                                  'long_name': 'Water Vapor Path Trends CMIP6',
                                  'units': 'percent'})
        list_dict["data"].append(res_ar['weights_c6'])
        list_dict["name"].append({'var_name': 'data_set_weights',
                                  'long_name': 'Weights for each data set.',
                                  'units': '1'})
        list_dict["data"].append(res_ar['kde1_c6'](res_ar['xhist']))
        list_dict["name"].append({'var_name': 'prw_trend_distribution_cmip6',
                                  'long_name': 'Water Vapor Path Trends ' +
                                               'distribution CMIP6',
                                  'units': '1'})

    if trends['obs']:
        for obsname in trends['obs'].keys():
            list_dict["data"].append(trends['obs'][obsname])
            list_dict["name"].append({'var_name': 'prw_trend_' + obsname,
                                      'long_name': 'Water Vapor Path Trend ' +
                                                   obsname,
                                      'units': 'percent'})

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def _plot_settings(cfg, axx, fig, figname, period):
    """Define Settings for pdf figures."""
    for spine in ["top", "right"]:
        axx.spines[spine].set_visible(False)
    axx.legend(loc=1)
    
    if 'ymax' in cfg.keys():
        ymax = cfg['ymax']
    else:
        ymax = 2.5
    axx.set_ylim([0, ymax])
    axx.set_ylabel('Probability density')
    add = '.'
    if 'common_period' in period.keys():
        add = ' for ' + period['common_period'] + '.'
    elif 'common_span' in period.keys():
        add = ' for ' + period['common_span'] + ' years.'
    if 'common_lats' in period.keys():
        add = ' between ' + period['common_lats'] + add

    vars = list(extract_variables(cfg).keys())
    var = " ".join(vars)
    if var == 'prw':
        var = ' Water Vapour Path'
    else:
        var = ' ' + var
    axx.set_title('Trends in' + var + add)
    axx.set_xlabel('Trend (%/decade)')
    fig.tight_layout()
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)

    return add


def _set_extratrends_dict(cfg):
    """Set dictionary to plot pdf over model ensembles."""
    extratrends = {}
    for extramodel in cfg['add_model_dist']:
        extratrends[extramodel] = OrderedDict()
    return extratrends


def _write_list_dict(cfg, trends, res_ar):
    """Collect data for provenance."""
    list_dict = {}
    list_dict["data"] = [res_ar['xhist']]
    list_dict["name"] = [{'var_name': 'prw_trends_bins',
                          'long_name': 'Water Vapor Path Trend bins',
                          'units': 'percent'}]

    for extramodel in cfg['add_model_dist']:
        list_dict["data"].append(res_ar['artrend'][extramodel])
        list_dict["name"].append({'var_name': 'prw_trends_' + extramodel,
                                  'long_name': 'Water Vapor Path Trends ' +
                                               extramodel,
                                  'units': 'percent'})
        list_dict["data"].append(res_ar['kde1'][extramodel](res_ar['xhist']))
        list_dict["name"].append({'var_name': 'prw_trend_distribution_' +
                                              extramodel,
                                  'long_name': 'Water Vapor Path Trends ' +
                                               'distribution ' + extramodel,
                                  'units': '1'})

    if trends['obs']:
        for obsname in trends['obs'].keys():
            list_dict["data"].append(trends['obs'][obsname])
            list_dict["name"].append({'var_name': 'prw_trend_' + obsname,
                                      'long_name': 'Water Vapor Path Trend ' +
                                                   obsname,
                                      'units': 'percent'})
    return list_dict


def cube_to_save_vars(list_dict):
    """Create cubes to prepare bar plot data for saving to netCDF."""
    # cubes = iris.cube.CubeList()
    for iii, var in enumerate(list_dict["data"]):
        if iii == 0:
            cubes = iris.cube.CubeList([
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units'])])
        else:
            cubes.append(
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units']))

    return cubes


def get_provenance_record(ancestor_files, caption, statistics,
                          domains, plot_type='probability'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': statistics,
        'domains': domains,
        'plot_type': plot_type,
        'realms': ['atmos'],
        'themes': ['atmDyn'],
        'authors': [
            'weigel_katja',
        ],
        'references': [
            'santer20jclim',
        ],
        'ancestors': ancestor_files,
    }
    return record


###############################################################################
# Setup diagnostic
###############################################################################


def main(cfg):
    """Run the diagnostic."""
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Dataset data containers
    input_data = (cfg['input_data'].values())

    if not variables_available(cfg, ['prw']):
        logger.warning("This diagnostic was written and tested only for prw.")

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset and save annual means
    trends = {}
    trends['cmip5'] = OrderedDict()
    trends['cmip6'] = OrderedDict()
    trends['cmip5weights'] = OrderedDict()
    trends['cmip6weights'] = OrderedDict()
    trends['obs'] = {}
    # f_c5 = open(cfg['work_dir'] + "/cmip5_trends.txt", "a")
    # f_c5.write("Model Alias Trend Weight Mean Median " +
    #            "Maximum Mnimum Standard deviation \n")
    # f_c6 = open(cfg['work_dir'] + "/cmip6_trends.txt", "a")
    # f_c6.write("Model Alias Trend Weight Mean Median " +
    #            "Maximum Mnimum Standard deviation \n")

    if 'add_model_dist' in cfg:
        extratrends = _set_extratrends_dict(cfg)

    valid_datasets, number_of_subdata, period = _get_valid_datasets(input_data)

    for dataset_path in input_data:
        project = dataset_path['project']
        cube_load = iris.load(dataset_path['filename'])[0]
        cube = _apply_filter(cfg, cube_load)
        cat.add_month_number(cube, 'time', name='month_number')
        cat.add_year(cube, 'time', name='year')
        alias = dataset_path['alias']
        dataset = dataset_path['dataset']
        # selection = select_metadata(input_data, dataset=dataset)
        # number_of_subdata = float(len(select_metadata(input_data,
        #                                               dataset=dataset)))
        if not _check_full_data(dataset_path, cube):
            continue
        cube_anom = _calculate_anomalies(cube)

        if project == 'CMIP6':
            trend = _calc_trend(cube_anom)
            trends['cmip6'][alias] = trend
            trends['cmip6weights'][alias] = 1. / number_of_subdata[dataset]
            # f_c6.write(dataset + ' ' +
            #            alias + ' ' +
            #            str(round(trend, 4)) + ' ' +
            #            str(round(1. / number_of_subdata[dataset], 4)) + ' ' +
            #            str(round(np.mean(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.median(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.max(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.min(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.std(cube_only_filt.data), 4)) + '\n')
        elif project == 'CMIP5':
            trend = _calc_trend(cube_anom)
            trends['cmip5'][alias] = trend
            trends['cmip5weights'][alias] = 1. / number_of_subdata[dataset]
            # f_c5.write(dataset + ' ' +
            #            alias + ' ' +
            #            str(round(trend, 4)) + ' ' +
            #            str(round(1. / number_of_subdata[dataset], 4)) + ' ' +
            #            str(round(np.mean(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.median(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.max(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.min(cube_only_filt.data), 4)) + ' ' +
            #            str(round(np.std(cube_only_filt.data), 4)) + '\n')
        else:
            trend = _calc_trend(cube_anom)
            trends['obs'][dataset] = trend

        if 'add_model_dist' in cfg:
            if dataset in cfg['add_model_dist']:
                extratrends[dataset][alias] = trend

    # f_c5.close()
    # f_c6.close()
    _plot_trends(cfg, trends, valid_datasets, period)
    if 'add_model_dist' in cfg:
        _plot_extratrends(cfg, extratrends, trends, period)

    ###########################################################################
    # Process data
    ###########################################################################


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
