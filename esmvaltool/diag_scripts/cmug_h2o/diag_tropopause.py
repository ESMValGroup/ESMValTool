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
from copy import deepcopy
from pprint import pformat

import cartopy.crs as cart
import iris
import matplotlib.pyplot as plt
import numpy as np
from cartopy.util import add_cyclic_point
from iris.analysis import Aggregator

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared._base import (
    get_diagnostic_filename,
    get_plot_filename,
)

logger = logging.getLogger(os.path.basename(__file__))


def _cube_to_save_ploted(var, lats, lons, names):
    """Create cube to prepare plotted data for saving to netCDF."""
    new_cube = iris.cube.Cube(var, var_name=names['var_name'],
                              long_name=names['long_name'],
                              units=names['units'])
    new_cube.add_dim_coord(iris.coords.DimCoord(lats,
                                                var_name='lat',
                                                long_name='latitude',
                                                units='degrees_north'), 0)
    new_cube.add_dim_coord(iris.coords.DimCoord(lons,
                                                var_name='lon',
                                                long_name='longitude',
                                                units='degrees_east'), 1)

    return new_cube


def _get_provenance_record(ancestor_files, caption, statistics,
                           domains, plot_type='geo'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': statistics,
        'domains': domains,
        'plot_type': plot_type,
        'projects': ['cmug'],
        'realms': ['atmos'],
        'themes': ['atmDyn'],
        'authors': [
            'weigel_katja',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def _press_interp_cube(cube, pstart=60000, pend=2500, pgrid=221):
    """Interpolate cube onto dense pressure grid."""
    # interpolate to dense grid, use equal dist points in log(p)
    logpr = np.log(np.array([pstart, pend]))
    sample_points = [('air_pressure', np.exp(np.linspace(logpr[0],
                                                         logpr[1],
                                                         pgrid)))]
    new_cube = cube.interpolate(sample_points, iris.analysis.Linear())

    return new_cube


def _set_new_dims(new_tacube, new_svarcube):
    """Set dimensions to unroll."""
    dims_to_collapse = set()
    dims_to_collapse.update(new_svarcube.coord_dims('air_pressure'))
    untouched_dims = set(range(new_svarcube.ndim)) - set(dims_to_collapse)
    dims = list(untouched_dims) + list(dims_to_collapse)
    unrolled_data = np.moveaxis(new_tacube.data, dims,
                                range(new_svarcube.ndim))
    return unrolled_data


def cube_to_save_ploted_map(var, lats, lons, names):
    """Create cube to prepare plotted data for saving to netCDF."""
    new_cube = iris.cube.Cube(var, var_name=names['var_name'],
                              long_name=names['long_name'],
                              units=names['units'])
    new_cube.add_dim_coord(iris.coords.DimCoord(lats,
                                                var_name='lat',
                                                long_name='latitude',
                                                units='degrees_north'), 0)
    new_cube.add_dim_coord(iris.coords.DimCoord(lons,
                                                var_name='lon',
                                                long_name='longitude',
                                                units='degrees_east'), 1)

    return new_cube


def find_min(data, data_min, axis):
    """Functions for Iris Aggregator to find the min based on other cube."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    min_ind = np.expand_dims(np.argmin(data_min, axis=axis), axis=axis)
    return_data = np.squeeze(np.take_along_axis(data, min_ind, axis=axis))

    return return_data


def plot_tp_map(cfg, mean_cube, titlestr, variable, listdata):
    """Plot contour map."""
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

    fig.tight_layout()
    figname = 'fig_' + titlestr.replace(" ", "_") +\
        variable.replace(" ", "_") + '_map'
    fig.savefig(get_plot_filename(figname, cfg), dpi=300)
    plt.close()

    logger.info("Saving analysis results to %s",
                get_diagnostic_filename(figname, cfg))

    iris.save(_cube_to_save_ploted(data_c,
                                   mean_cube.coord('latitude').points,
                                   lon_c,
                                   {'var_name': mean_cube.var_name,
                                    'long_name': variable,
                                    'units': mean_cube.units}),
              target=get_diagnostic_filename(figname, cfg))

    logger.info("Recording provenance of %s:\n%s",
                get_diagnostic_filename(figname, cfg),
                pformat(_get_provenance_record(listdata,
                                               titlestr + variable,
                                               ['other'], ['global'])))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(get_diagnostic_filename(figname, cfg),
                              _get_provenance_record(listdata,
                                                     titlestr + variable,
                                                     ['other'], ['global']))


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
        logger.debug("Loading %s", attributes['filename'])
        new_tacube = _press_interp_cube(iris.load_cube(attributes['filename']))
        for svar in available_vars_min_tas:
            input_file_svar = attributes['filename'].replace('/ta/',
                                                             '/' + svar + '/')
            input_file_svar = input_file_svar.replace('_ta_', '_' + svar + '_')

            # load, interpolate to dense grid, use equal dist points in log(p)
            logger.debug("Loading %s", input_file_svar)
            new_svarcube = _press_interp_cube(iris.load_cube(input_file_svar))

            # set new dims
            unrolled_data = _set_new_dims(new_tacube, new_svarcube)

            trop_svar = new_svarcube.collapsed('air_pressure', min_pos,
                                               data_min=unrolled_data)
            plot_tp_map(cfg, trop_svar.collapsed('time', iris.analysis.MEAN),
                        dataset + " Cold point tropopause ",
                        new_svarcube.long_name,
                        [data['ta'][0]['filename'], data[svar][0]['filename']])

        trop_ta = new_tacube.collapsed('air_pressure', min_pos,
                                       data_min=unrolled_data)

        plot_tp_map(cfg, trop_ta.collapsed('time', iris.analysis.MEAN),
                    dataset + " Cold point tropopause ",
                    "Air Temperature", [data['ta'][0]['filename']])


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
