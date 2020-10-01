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

logger = logging.getLogger(os.path.basename(__file__))


def find_min(data, data_min, axis):
    """Functions for Iris Aggregator to find the min based on other cube."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    min_ind = np.expand_dims(np.argmin(data_min, axis=axis), axis=axis)
    return_data = np.squeeze(np.take_along_axis(data, min_ind, axis=axis))

    return return_data


def plot_zonal_mean(cfg, mean_cube, titlestr, variable):
    """Plot zonal mean contour."""
    # if not cfg[n.WRITE_PLOTS]:
    #     return

    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(7, 5))

    if variable == "Air Temperature":
        print_var = "Temperature [K]"
        set_range = np.linspace(180, 230, 21)
    elif variable == "Geopotential Height":
        print_var = "Geopotential Height [m]"
        set_range = np.linspace(10000, 20000, 21)
    elif variable == "Relative Humidity":
        print_var = "Relative Humidity [%]"
        set_range = np.linspace(0, 100, 21)
    elif variable == "Specific Humidity":
        print_var = "Specific Humidity [kg/kg]"
        logval = np.log(np.array([1e-6, 1e-5]))
        set_range = np.exp(np.linspace(logval[0], logval[1], 41))
    else:
        print_var = mean_cube.long_name
        set_range = np.linspace(np.nanmin(mean_cube.data),
                                np.nanmax(mean_cube.data), 21)
    print("mean_cube.units")
    print(mean_cube.units)
    # draw filled contours
    cnplot = plt.contourf(
        mean_cube.coord('latitude').points,
        mean_cube.coord('air_pressure').points / 100.0,
        mean_cube.data,
        set_range,
        cmap='jet',
        # cmap='RdBu_r',
        extend='both')
    
    axx = plt.gca()
    axx.invert_yaxis()
    axx.set_ylim(250, 1)

    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(print_var)

    axx.set_ylabel('Pressure [hPa]')
    axx.set_xlabel('Latitude')
    axx.set_title(titlestr + variable)
    # axx.set_xticks([40, 65, 90, 115])
    # axx.set_xticklabels(['40°E', '65°E', '90°E', '115°E'])
    # axx.set_yticks([-10, 0, 10, 20, 30])
    # axx.set_yticklabels(['10°S', '0°', '10°N', '20°N', '30°N'])

    fig.tight_layout()
    figname = 'fig_' + titlestr.replace(" ", "_") + variable.replace(" ",
                                                                     "_")
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()

#    titlestr = titlestr + ' Box displays the area used to define the ' + \
#        'average ISM (Indian Summer Monsoon) rainfall. Precipitation ' + \
#        'changes are normalized by the corresponding global ' + \
#        'mean SST increase for each model.'
#
#    selection = _get_sel_files_var(cfg, ['pr', 'ts'])
#
#    provenance_record = get_provenance_record(selection,
#                                              titlestr, ['diff'], ['reg'])
#
#    diagnostic_file = get_diagnostic_filename(figname, cfg)
#
#    logger.info("Saving analysis results to %s", diagnostic_file)
#
#    iris.save(cube_to_save_ploted(data, lats, lons, {'var_name': 'd_pr',
#                                                     'long_name': 'Prec' +
#                                                                  'ipita' +
#                                                                  'tion ' +
#                                                                  'Change',
#                                                     'units': 'mm d-1'}),
#              target=diagnostic_file)
#
#    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
#                pformat(provenance_record))
#    with ProvenanceLogger(cfg) as provenance_logger:
#        provenance_logger.log(diagnostic_file, provenance_record)


def plot_zonal_timedev(cfg, mean_cube, titlestr, variable):
    """Plot zonal mean temporal developement at tropopause."""
    # if not cfg[n.WRITE_PLOTS]:
    #     return

    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(7, 12))

    if variable == "Air Temperature":
        print_var = "Temperature [K]"
        set_range = np.linspace(180, 230, 21)
    elif variable == "Geopotential Height":
        print_var = "Geopotential Height [m]"
        set_range = np.linspace(10000, 20000, 21)
    elif variable == "Relative Humidity":
        print_var = "Relative Humidity [%]"
        set_range = np.linspace(0, 100, 21)
    elif variable == "Specific Humidity":
        print_var = "Specific Humidity [kg/kg]"
        set_range = np.linspace(0.1e-5, 2.5e-5, 25)
    else:
        print_var = mean_cube.long_name
        set_range = np.linspace(np.nanmin(mean_cube.data),
                                np.nanmax(mean_cube.data), 21)
    print("mean_cube.units")
    print(mean_cube.units)
    # draw filled contours
    print("mean_cube.data.shape")
    print(mean_cube.data.shape)
    iris.coord_categorisation.add_year(mean_cube, 'time', name='year')
    iris.coord_categorisation.add_month_number(mean_cube, 'time',
                                               name='month_number')
    print("mean_cube.coord('year').points")
    print(mean_cube.coord('year').points)
    # Adjust (ncdf) time to the format matplotlib expects
    cnplot = plt.contourf(
        mean_cube.coord('latitude').points,
        mean_cube.coord('year').points +
        (mean_cube.coord('month_number').points - 1.0)/12.0,
        mean_cube.data,
        set_range,
        cmap='rainbow',
        # cmap='RdBu_r',
        extend='both')

    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(print_var)

    axx.set_ylabel('Time')
    axx.set_xlabel('Latitude')
    axx.set_title(titlestr + variable)
    # axx.set_xticks([40, 65, 90, 115])
    # axx.set_xticklabels(['40°E', '65°E', '90°E', '115°E'])
    # axx.set_yticks([-10, 0, 10, 20, 30])
    # axx.set_yticklabels(['10°S', '0°', '10°N', '20°N', '30°N'])

    fig.tight_layout()
    figname = 'fig_' + titlestr.replace(" ", "_") + variable.replace(" ",
                                                                     "_")
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()


def plot_profiles(cfg, profiles, available_vars_min_tas, available_datasets):
    """Plot zonal mean contour."""
    # if not cfg[n.WRITE_PLOTS]:
    #     return

    # Plot data
    # create figure and axes instances
    for svar in available_vars_min_tas:
        fig, axx = plt.subplots(figsize=(7, 5))
        
        for dataset in available_datasets:
            plt.plot((profiles[svar][dataset]).data,
                     (profiles[svar][dataset]).coord('air_pressure').points /
                     100.0, label = dataset)
    
        axx = plt.gca()
        axx.invert_yaxis()
        
        
        plt.legend(loc='upper right')

        axx.set_ylabel('Pressure [hPa]')
        axx.set_ylim(250, 1)
        axx.set_xlim(0, 1e-4)
        # axx.set_xscale('log')
        if profiles[svar][dataset].long_name ==  "Specific Humidity":
            unitstr = "kg/kg"
        else:
            unitstr = str(profiles[svar][dataset].units)
        
        
        axx.set_xlabel(profiles[svar][dataset].long_name +
                       ' [' + unitstr + ']')
        axx.set_title('Average ' + profiles[svar][dataset].long_name +
                      ' profile')
        fig.tight_layout()
        figname = 'fig_profile_' +\
            profiles[svar][dataset].long_name.replace(" ", "_")
        fig.savefig(get_plot_filename(figname, cfg), dpi=300)
        plt.close()


def main(cfg):
    """Read in data for tropopause calculation."""
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))
    available_datasets = list(group_metadata(cfg['input_data'].values(),
                                         'dataset'))

    logging.debug("Found variables in recipe:\n%s", available_vars)
    available_vars_min_tas = deepcopy(available_vars)
    available_vars_min_tas.remove('ta')
    # Make an iris aggregator to find value based on minimum different cube.
    min_pos = Aggregator('ag_find_min', find_min,
                         units_func=lambda units: 1)
    # Get input data
    data = {}
    for varname in available_vars:
        data[varname] = select_metadata(cfg['input_data'].values(),
                                        short_name=varname)
        data[varname] = sorted_metadata(data[varname], sort='dataset')
    
    
    profiles = {}
    for svar in available_vars_min_tas:
        profiles[svar] = {}

    for attributes in data['ta']:
        logger.info("Processing dataset %s", attributes['dataset'])
        dataset = attributes['dataset']
        for svar in available_vars_min_tas:
            input_file_svar = attributes['filename'].replace('/ta/',
                                                             '/' + svar + '/')
            input_file_svar = input_file_svar.replace('_ta_', '_' + svar + '_')

            logger.debug("Loading %s", input_file_svar)
            svarcube = iris.load_cube(input_file_svar).collapsed('longitude',
                                                        iris.analysis.MEAN)
            profiles[svar][dataset] = svarcube.collapsed(['time', 'latitude'],
                                                    iris.analysis.MEAN)
            logger.debug("Loading %s", attributes['filename'])
            tacube = iris.load_cube(attributes['filename']).collapsed('longitude',
                                                        iris.analysis.MEAN)
            # plev = tacube.coord('air_pressure').points
            plot_zonal_mean(cfg, svarcube.collapsed(['time'],
                                                    iris.analysis.MEAN),
                        dataset + " Zonal mean ",
                        svarcube.long_name)

            # interpolate to dense grid, use equal dist points in log(p)
            logpr = np.log(np.array([25000, 2500]))
            sample_points = [('air_pressure', np.exp(np.linspace(logpr[0],
                                                                 logpr[1],
                                                                 221)))]

            new_tacube = tacube.interpolate(sample_points,
                                            iris.analysis.Linear())
            new_svarcube = svarcube.interpolate(sample_points,
                                                iris.analysis.Linear())

            dims_to_collapse = set()
            dims_to_collapse.update(new_svarcube.coord_dims('air_pressure'))
            untouched_dims = set(range(new_svarcube.ndim)) -\
                set(dims_to_collapse)
            dims = list(untouched_dims) + list(dims_to_collapse)
            unrolled_data = np.moveaxis(new_tacube.data, dims,
                                        range(new_svarcube.ndim))

            trop_svar = new_svarcube.collapsed('air_pressure', min_pos,
                                               data_min=unrolled_data)
            plot_zonal_timedev(cfg, trop_svar,
                        dataset + " Cold point tropopause ",
                        new_svarcube.long_name)

        trop_ta = new_tacube.collapsed('air_pressure', min_pos,
                                       data_min=unrolled_data)

        plot_zonal_timedev(cfg, trop_ta,
                    dataset + " Cold point tropopause ",
                    "Air Temperature")
    print("profiles")
    print(profiles)
    
    plot_profiles(cfg, profiles, available_vars_min_tas, available_datasets)
#        plot_tp_map(cfg, trop_hus.collapsed('time', iris.analysis.MEAN),
#                    dataset + " Cold point tropopause specific humidity",
#                    "Specific humidity [kg/kg]")
#        plot_tp_map(cfg, trop_hur.collapsed('time', iris.analysis.MEAN),
#                    dataset + " Cold point tropopause relative humidity",
#                    "Relative humidity [%]")


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
