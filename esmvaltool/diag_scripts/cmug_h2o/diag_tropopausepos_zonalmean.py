#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Interpolate given variable to tropopause height.

###############################################################################
testkw/diag_tropopause.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
ESA-CMUG project
###############################################################################

Description
-----------
    Interpolate given variable to tropopause height.
    Currently only cold point tropopause is used based on ta and plev.

Configuration options
---------------------

###############################################################################

"""


import logging
import os
from pprint import pformat
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import iris
from iris.analysis import Aggregator

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.cmug_h2o.diag_tropopause_zonalmean import (
    read_data_trop_zonal,
    get_sel_lvardata,
    get_range_and_pstring,
    get_provenance_record)

logger = logging.getLogger(os.path.basename(__file__))


def _get_data_for_agg(new_svarcube, new_tacube):
    """Reshape data for use in iris aggregator based on two cubes."""
    dims_to_collapse = set()
    dims_to_collapse.update(new_svarcube.coord_dims('air_pressure'))
    untouched_dims = set(range(new_svarcube.ndim)) -\
        set(dims_to_collapse)
    dims = list(untouched_dims) + list(dims_to_collapse)
    unrolled_data = np.moveaxis(new_tacube.data, dims,
                                range(new_svarcube.ndim))
    return unrolled_data


def _read_data_ta(attributes, svar):
    """Read data for ta instead of variabe from files, average longitudes."""
    print(attributes['filename'])
    print('/' + svar + '/')
    input_file_svar = attributes['filename'].replace('/' + svar + '/', '/ta/')
    input_file_svar = input_file_svar.replace('_' + svar + '_', '_ta_')

    logger.debug("Loading %s", input_file_svar)
    svarcube = iris.load_cube(input_file_svar)
    svarcube = svarcube.collapsed('longitude', iris.analysis.MEAN)

    return svarcube


def find_min(data, data_min, axis):
    """Functions for Iris Aggregator to find the min based on other cube."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    min_ind = np.expand_dims(np.argmin(data_min, axis=axis), axis=axis)
    return_data = np.squeeze(np.take_along_axis(data, min_ind, axis=axis))

    return return_data

def find_min_ind(data, axis):
    """Functions for Iris Aggregator to find the index for min of cube."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    min_ind = np.squeeze(np.expand_dims(np.argmin(data, axis=axis), axis=axis))

    return min_ind


def get_prof_and_plt_data_tropo(cfg, data, available_vars_min_tas):
    """Plot data for singe data sets and get profile for each."""
    profiles = {}
    for svar in available_vars_min_tas:
        profiles[svar] = {}

    # Make an iris aggregator to find value based on minimum different cube.
    min_pos = Aggregator('ag_find_min', find_min,
                         units_func=lambda units: 1)
    
     # Make an iris aggregator to find value based on minimum different cube.
    min_ind = Aggregator('ag_find_min_ind', find_min_ind,
                         units_func=lambda units: 1)

    # interpolate to dense grid, use equal dist points in log(p)
    logpr = np.log(np.array([25000, 2500]))
    sample_points = [('air_pressure', np.exp(np.linspace(logpr[0], logpr[1],
                                                         221)))]
    
    for svar in available_vars_min_tas:
        for attributes in data[svar]:
            dataset = attributes['dataset']
            if dataset == cfg['reference_dataset']:
                att_reference_dataset = attributes
    
    for svar in available_vars_min_tas:
        for attributes in data[svar]:
            dataset = attributes['dataset']
            logger.info("Processing dataset %s", dataset)
            if cfg['ref_ta_all']:
                logger.info("Using ta from dataset %s", cfg['reference_dataset'])
                tacube = _read_data_ta(att_reference_dataset, svar)
            elif (cfg['ref_ta_alternative_dataset'] and dataset == cfg['alternative_dataset']):
                logger.info("Using ta from dataset %s", cfg['reference_dataset'])
                tacube = _read_data_ta(att_reference_dataset, svar)
            else:
                logger.info("Using ta from dataset %s", dataset)
                tacube = _read_data_ta(attributes, svar)
            
            new_tacube = tacube.interpolate(sample_points,
                                            iris.analysis.Linear())
            unrolled_data = _get_data_for_agg(new_tacube, new_tacube)
            plot_zonal_timedev(cfg, new_tacube.collapsed('air_pressure',
                                                         min_pos,
                                                         data_min=unrolled_data),
                               dataset, "Cold point tropopause ",
                               "Air Temperature")

#        for svar in available_vars_min_tas:
            svarcube = read_data_trop_zonal(attributes)

            profiles[svar][dataset] = svarcube.collapsed(['time', 'latitude'],
                                                         iris.analysis.MEAN)

            new_svarcube = svarcube.interpolate(sample_points,
                                                iris.analysis.Linear())

            unrolled_data = _get_data_for_agg(new_svarcube, new_tacube)

            plot_zonal_timedev(cfg,
                               new_svarcube.collapsed('air_pressure',
                                                      min_pos,
                                                      data_min=unrolled_data),
                               dataset, "Cold point tropopause ",
                               new_svarcube.long_name)

            ind_cube = new_tacube.collapsed('air_pressure',
                                                      min_ind)
            
            ind_cube.units = '1'
            newind_cube = deepcopy(ind_cube)
            newind_cube.units = 'Pa'
            # newind_cube.name = 'air_pressure'
            for iii, pii in enumerate(ind_cube.data):
                newind_cube.data[iii] = new_tacube.coord('air_pressure').points[pii]
            
            plot_zonal_mean_trop(cfg,
                                 svarcube.collapsed(['time'],
                                                    iris.analysis.MEAN),
                                 newind_cube.collapsed(['time'],
                                                    iris.analysis.MEAN),
                                dataset, "Zonal mean with Tropopause",
                                svarcube.long_name)
            
            plot_taperec_trop(cfg,
                               svarcube.collapsed('latitude',
                                                      iris.analysis.MEAN),
                               newind_cube.collapsed('latitude',
                                                      iris.analysis.MEAN),
                               dataset, "Tape recorder with Tropopause",
                               svarcube.long_name)
                

    return profiles


def plot_zonal_mean_trop(cfg, mean_cube, mean_plev_cube, dataname, titlestr, variable):
    """Plot zonal mean contour."""
    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(7, 5))

    rps = get_range_and_pstring(variable, mean_cube)

    # draw filled contours
    cnplot = plt.contourf(
        mean_cube.coord('latitude').points,
        mean_cube.coord('air_pressure').points / 100.0,
        mean_cube.data,
        rps['set_range'],
        cmap='jet',
        extend='both')

    axx = plt.gca()
    axx.invert_yaxis()
    axx.set_ylim(250, 1)
    
    axx.plot(mean_plev_cube.coord('latitude').points, mean_plev_cube.data / 100.0, marker='o', color='violet')

    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(rps['print_var'])

    axx.set_ylabel('Pressure [hPa]')
    axx.set_xlabel('Latitude')
    axx.set_title(titlestr + variable)

    fig.tight_layout()

    caption = dataname + " " + titlestr + variable +\
        " between 250 and 1 hPa." +\
        " The diagnostic averages the complete time series."
    figname = 'fig_' + dataname + "_" + titlestr.replace(" ", "_") +\
              variable.replace(" ", "_")
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()

    provenance_record = get_provenance_record(get_sel_lvardata(cfg,
                                                                dataname,
                                                                [variable]),
                                              caption, ['mean'], ['global'])
#
    diagnostic_file = get_diagnostic_filename(figname, cfg)
#
    logger.info("Saving analysis results to %s", diagnostic_file)
#
    iris.save(mean_cube, target=diagnostic_file)
#
    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def plot_zonal_timedev(cfg, mean_cube, dataname, titlestr, variable):
    """Plot zonal mean temporal developement at tropopause."""
    # Plot data
    # create figure and axes instances
    print("mean_cube zonal_timedev")
    print(mean_cube)
    fig, axx = plt.subplots(figsize=(7, 12))

    rps = get_range_and_pstring(variable, mean_cube, tropopause=True)

    iris.coord_categorisation.add_year(mean_cube, 'time', name='year')
    iris.coord_categorisation.add_month_number(mean_cube, 'time',
                                               name='month_number')

    # Adjust (ncdf) time to the format matplotlib expects
    cnplot = plt.contourf(
        mean_cube.coord('latitude').points,
        mean_cube.coord('year').points +
        (mean_cube.coord('month_number').points - 1.0) / 12.0,
        mean_cube.data,
        rps['set_range'],
        cmap='rainbow',
        # cmap='RdBu_r',
        extend='both')

    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(rps['print_var'])

    axx.set_ylabel('Time')
    axx.set_xlabel('Latitude')
    axx.set_title(titlestr + variable)

    fig.tight_layout()
    figname = 'fig_' + dataname + "_" + titlestr.replace(" ", "_") + \
        variable.replace(" ", "_")
    caption = dataname + " " + titlestr + variable
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()

    provenance_record = get_provenance_record(get_sel_lvardata(cfg,
                                                                dataname,
                                                                [variable]),
                                              caption, ['mean'], ['global'])
#
    diagnostic_file = get_diagnostic_filename(figname, cfg)
#
    logger.info("Saving analysis results to %s", diagnostic_file)
#
    iris.save(mean_cube, target=diagnostic_file)
#
    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def plot_taperec_trop(cfg, mean_cube, mean_plev_cube, dataname, titlestr, variable):
    """Plot zonal mean temporal developement at tropopause."""
    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(12, 7))

    rps = get_range_and_pstring(variable, mean_cube, tropopause=True)

    iris.coord_categorisation.add_year(mean_cube, 'time', name='year')
    iris.coord_categorisation.add_month_number(mean_cube, 'time',
                                               name='month_number')

    # Adjust (ncdf) time to the format matplotlib expects
    print("mean_cube.coord('air_pressure').points")
    print(mean_cube.coord('air_pressure').points)
    cnplot = plt.contourf(
        mean_cube.coord('year').points +
        (mean_cube.coord('month_number').points - 1.0) / 12.0,
        mean_cube.coord('air_pressure').points / 100.0,
        np.transpose(mean_cube.data),
        rps['set_range'],
        cmap='jet',
        # cmap='rainbow',
        # cmap='RdBu_r',
        extend='both')

    axx.plot(mean_cube.coord('year').points + (mean_cube.coord('month_number').points - 1.0) / 12.0,
             mean_plev_cube.data / 100.0, marker='o', color='violet')
    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(rps['print_var'])

    axx.set_xlabel('Time')
    axx.set_yscale('log')
    axx.set_ylabel('Pressure [hPa]')
    axx.set_ylim(250, 1)
    axx.set_title(titlestr + variable)

    fig.tight_layout()
    figname = 'fig_' + dataname + "_" + titlestr.replace(" ", "_") + \
        variable.replace(" ", "_")
    caption = dataname + " " + titlestr + variable
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()

    provenance_record = get_provenance_record(get_sel_lvardata(cfg,
                                                                dataname,
                                                                [variable]),
                                              caption, ['mean'], ['global'])
    diagnostic_file = get_diagnostic_filename(figname, cfg)
    logger.info("Saving analysis results to %s", diagnostic_file)
    iris.save(mean_cube, target=diagnostic_file)
    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def main(cfg):
    """Read in data for tropopause calculation."""
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))

    logging.debug("Found variables in recipe:\n%s", available_vars)
    available_vars_min_tas = deepcopy(available_vars)
    available_vars_min_tas.remove('ta')
    # Get input data
    data = {}
    for varname in available_vars:
        data[varname] = select_metadata(cfg['input_data'].values(),
                                        short_name=varname)
        data[varname] = sorted_metadata(data[varname], sort='dataset')

    get_prof_and_plt_data_tropo(cfg, data, available_vars_min_tas)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
