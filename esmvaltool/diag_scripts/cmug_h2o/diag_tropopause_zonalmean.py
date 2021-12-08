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


def _get_range_and_pstring(variable, mean_cube, tropopause=False):
    """Get range for color bar and print string."""
    if variable == "Air Temperature":
        print_var = "Temperature [K]"
        set_range = np.linspace(180, 230, 21)
    elif variable == "Geopotential Height":
        print_var = "Geopotential Height [m]"
        set_range = np.linspace(8000, 22000, 25)
    elif variable == "Relative Humidity":
        print_var = "Relative Humidity [%]"
        set_range = np.linspace(0, 100, 21)
    elif variable == "Specific Humidity":
        print_var = "Specific Humidity [kg/kg]"
        if tropopause:
            set_range = np.linspace(0.1e-5, 2.5e-5, 25)
        else:
            logval = np.log(np.array([1e-6, 1e-5]))
            set_range = np.exp(np.linspace(logval[0], logval[1], 41))
    else:
        print_var = mean_cube.long_name
        set_range = np.linspace(np.nanmin(mean_cube.data),
                                np.nanmax(mean_cube.data), 21)

    return {'print_var': print_var, 'set_range': set_range}


def _get_sel_files(cfg, dataname, dim=2):
    """Get filenames from cfg for single models or multi-model mean."""
    selection = []
    if dim == 2:
        for hlp in select_metadata(cfg['input_data'].values(),
                                   dataset=dataname):
            selection.append(hlp['filename'])
    else:
        for hlp in cfg['input_data'].values():
            selection.append(hlp['filename'])

    return selection


def _get_sel_files_var(cfg, varnames):
    """Get filenames from cfg for all model mean and differen variables."""
    selection = []

    for var in varnames:
        for hlp in select_metadata(cfg['input_data'].values(), short_name=var):
            selection.append(hlp['filename'])

    return selection


def _get_sel_lvardata(cfg, dataname, lvarnames):
    """Get filenames from cfg for one model and differen variable(s)."""
    selection = []

    for lvar in lvarnames:
        for hlp in select_metadata(cfg['input_data'].values(), long_name=lvar,
                                   dataset=dataname):
            selection.append(hlp['filename'])

    return selection


def _read_data(attributes, svar):
    """Read data for ta and other variabe from files."""
    input_file_svar = attributes['filename'].replace('/ta/', '/' + svar + '/')
    input_file_svar = input_file_svar.replace('_ta_', '_' + svar + '_')

    logger.debug("Loading %s", input_file_svar)
    svarcube = iris.load_cube(input_file_svar)
    svarcube = svarcube.collapsed('longitude', iris.analysis.MEAN)

    return svarcube


def cube_to_save_profile(var1, var2, names):
    """Create cubes to prepare scatter plot data for saving to netCDF."""
    cubes = iris.cube.CubeList([iris.cube.Cube(var1,
                                               var_name=names['var_name1'],
                                               long_name=names['long_name1'],
                                               units=names['units1'])])
    cubes.append(iris.cube.Cube(var2, var_name=names['var_name2'],
                                long_name=names['long_name2'],
                                units=names['units2']))

    return cubes


def find_min(data, data_min, axis):
    """Functions for Iris Aggregator to find the min based on other cube."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim

    min_ind = np.expand_dims(np.argmin(data_min, axis=axis), axis=axis)
    return_data = np.squeeze(np.take_along_axis(data, min_ind, axis=axis))

    return return_data


def get_provenance_record(ancestor_files, caption, statistics,
                          domains, plot_type='zonal'):
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


def get_prof_and_plt_data(cfg, data, available_vars_min_tas):
    """Plot data for singe data sets and get profile for each."""
    profiles = {}
    for svar in available_vars_min_tas:
        profiles[svar] = {}

    # Make an iris aggregator to find value based on minimum different cube.
    min_pos = Aggregator('ag_find_min', find_min,
                         units_func=lambda units: 1)

    # interpolate to dense grid, use equal dist points in log(p)
    logpr = np.log(np.array([25000, 2500]))
    sample_points = [('air_pressure', np.exp(np.linspace(logpr[0], logpr[1],
                                                         221)))]

    for attributes in data['ta']:
        logger.info("Processing dataset %s", attributes['dataset'])
        dataset = attributes['dataset']

        tacube = _read_data(attributes, 'ta')
        new_tacube = tacube.interpolate(sample_points,
                                        iris.analysis.Linear())
        unrolled_data = _get_data_for_agg(new_tacube, new_tacube)
        plot_zonal_timedev(cfg, new_tacube.collapsed('air_pressure',
                                                     min_pos,
                                                     data_min=unrolled_data),
                           dataset, "Cold point tropopause ",
                           "Air Temperature")

        for svar in available_vars_min_tas:
            svarcube = _read_data(attributes, svar)

            profiles[svar][dataset] = svarcube.collapsed(['time', 'latitude'],
                                                         iris.analysis.MEAN)

            plot_zonal_mean(cfg, svarcube.collapsed(['time'],
                                                    iris.analysis.MEAN),
                            dataset, "Zonal mean ",
                            svarcube.long_name)

            new_svarcube = svarcube.interpolate(sample_points,
                                                iris.analysis.Linear())

            unrolled_data = _get_data_for_agg(new_svarcube, new_tacube)

            plot_zonal_timedev(cfg,
                               new_svarcube.collapsed('air_pressure',
                                                      min_pos,
                                                      data_min=unrolled_data),
                               dataset, "Cold point tropopause ",
                               new_svarcube.long_name)

    return profiles


def plot_zonal_mean(cfg, mean_cube, dataname, titlestr, variable):
    """Plot zonal mean contour."""
    # Plot data
    # create figure and axes instances
    fig, axx = plt.subplots(figsize=(7, 5))

    rps = _get_range_and_pstring(variable, mean_cube)

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

    provenance_record = get_provenance_record(_get_sel_lvardata(cfg,
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
    fig, axx = plt.subplots(figsize=(7, 12))

    rps = _get_range_and_pstring(variable, mean_cube, tropopause=True)

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

    provenance_record = get_provenance_record(_get_sel_lvardata(cfg,
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


def plot_profiles(cfg, profiles, available_vars_min_tas, available_datasets):
    """Plot zonal mean contour."""
    # Plot data
    # create figure and axes instances
    for svar in available_vars_min_tas:
        fig, axx = plt.subplots(figsize=(7, 5))

        for iii, dataset in enumerate(available_datasets):
            plt.plot((profiles[svar][dataset]).data,
                     (profiles[svar][dataset]).coord('air_pressure').points /
                     100.0, label=dataset)
            if iii == 0:
                profiles_save = iris.cube.CubeList([profiles[svar][dataset]])
            else:
                profiles_save.append(profiles[svar][dataset])

        axx = plt.gca()
        axx.invert_yaxis()

        plt.legend(loc='upper right')

        axx.set_ylabel('Pressure [hPa]')
        axx.set_ylim(250, 1)
        onedat = profiles[svar][available_datasets[0]]
        if onedat.long_name == "Specific Humidity":
            unitstr = "kg/kg"
            axx.set_xlim(0, 1e-4)
        else:
            unitstr = str(onedat.units)

        axx.set_xlabel(onedat.long_name +
                       ' [' + unitstr + ']')
        axx.set_title('Average ' +
                      onedat.long_name +
                      ' profile')
        fig.tight_layout()
        figname = 'fig_profile_' +\
            onedat.long_name.replace(" ", "_")
        fig.savefig(get_plot_filename(figname, cfg), dpi=300)
        plt.close()

        provenance_record = get_provenance_record(_get_sel_files_var(cfg,
                                                                     svar),
                                                  'Average ' +
                                                  onedat.long_name +
                                                  ' profile', ['mean'],
                                                  ['global'])

        diagnostic_file = get_diagnostic_filename(figname, cfg)
#
        logger.info("Saving analysis results to %s", diagnostic_file)
#
        iris.save(profiles_save, target=diagnostic_file)
#
        logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                    pformat(provenance_record))
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(diagnostic_file, provenance_record)


def main(cfg):
    """Read in data for tropopause calculation."""
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))
    available_datasets = list(group_metadata(cfg['input_data'].values(),
                                             'dataset'))

    logging.debug("Found variables in recipe:\n%s", available_vars)
    available_vars_min_tas = deepcopy(available_vars)
    available_vars_min_tas.remove('ta')
    # Get input data
    data = {}
    for varname in available_vars:
        data[varname] = select_metadata(cfg['input_data'].values(),
                                        short_name=varname)
        data[varname] = sorted_metadata(data[varname], sort='dataset')

    profiles = get_prof_and_plt_data(cfg, data, available_vars_min_tas)

    plot_profiles(cfg, profiles, available_vars_min_tas, available_datasets)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
