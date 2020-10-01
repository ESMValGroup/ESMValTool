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
import cartopy.crs as cart
from cartopy.util import add_cyclic_point
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


def plot_tp_map(cfg, mean_cube, titlestr, variable):
    """Plot contour map."""
    # if not cfg[n.WRITE_PLOTS]:
    #     return

    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(7, 5))
    axx = plt.axes(projection=cart.PlateCarree())
    axx.set_extent([-180, 180, -90, 90], cart.PlateCarree())

    data_c, lon_c = add_cyclic_point(mean_cube.data,
                                     coord=mean_cube.coord('longitude').points)

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
    cnplot = plt.contourf(
        lon_c,
        mean_cube.coord('latitude').points,
        data_c,
        set_range,
        transform=cart.PlateCarree(),
        cmap='rainbow',
        # cmap='RdBu_r',
        extend='both')

    axx.coastlines()

    # add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.8)
    # add colorbar title string
    cbar.set_label(print_var)

    axx.set_xlabel('Longitude')
    axx.set_ylabel('Latitude')
    axx.set_title(titlestr + variable)
    # axx.set_xticks([40, 65, 90, 115])
    # axx.set_xticklabels(['40°E', '65°E', '90°E', '115°E'])
    # axx.set_yticks([-10, 0, 10, 20, 30])
    # axx.set_yticklabels(['10°S', '0°', '10°N', '20°N', '30°N'])

    fig.tight_layout()
    figname = 'fig_' + titlestr.replace(" ", "_") +\
        variable.replace(" ", "_") + '_map'
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


def main(cfg):
    """Read in data for tropopause calculation."""
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))

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

    for attributes in data['ta']:
        logger.info("Processing dataset %s", attributes['dataset'])
        dataset = attributes['dataset']
        for svar in available_vars_min_tas:
            input_file_svar = attributes['filename'].replace('/ta/',
                                                             '/' + svar + '/')
            input_file_svar = input_file_svar.replace('_ta_', '_' + svar + '_')

            logger.debug("Loading %s", input_file_svar)
            svarcube = iris.load_cube(input_file_svar)
            logger.debug("Loading %s", attributes['filename'])
            tacube = iris.load_cube(attributes['filename'])
            # plev = tacube.coord('air_pressure').points

            # interpolate to dense grid, use equal dist points in log(p)
            logpr = np.log(np.array([60000, 2500]))
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
            plot_tp_map(cfg, trop_svar.collapsed('time', iris.analysis.MEAN),
                        dataset + " Cold point tropopause ",
                        new_svarcube.long_name)

        trop_ta = new_tacube.collapsed('air_pressure', min_pos,
                                       data_min=unrolled_data)

        plot_tp_map(cfg, trop_ta.collapsed('time', iris.analysis.MEAN),
                    dataset + " Cold point tropopause ",
                    "Air Temperature")
#        plot_tp_map(cfg, trop_hus.collapsed('time', iris.analysis.MEAN),
#                    dataset + " Cold point tropopause specific humidity",
#                    "Specific humidity [kg/kg]")
#        plot_tp_map(cfg, trop_hur.collapsed('time', iris.analysis.MEAN),
#                    dataset + " Cold point tropopause relative humidity",
#                    "Relative humidity [%]")


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
