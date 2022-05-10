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
import numpy as np
import matplotlib.pyplot as plt
import iris

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata,
                                            plot)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)

logger = logging.getLogger(os.path.basename(__file__))


def get_range_and_pstring(variable, mean_cube, tropopause=False):
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
            #set_range = np.linspace(0.1e-5, 2.5e-5, 25)
            set_range = np.linspace(0.1e-5, 1e-5, 41)
        else:
            logval = np.log(np.array([1e-6, 1e-5]))
            set_range = np.exp(np.linspace(logval[0], logval[1], 41))
    else:
        print_var = mean_cube.long_name
        set_range = np.linspace(np.nanmin(mean_cube.data),
                                np.nanmax(mean_cube.data), 21)

    return {'print_var': print_var, 'set_range': set_range}


def get_sel_files_var(cfg, varnames):
    """Get filenames from cfg for all model mean and differen variables."""
    selection = []

    for var in varnames:
        for hlp in select_metadata(cfg['input_data'].values(), short_name=var):
            selection.append(hlp['filename'])

    return selection


def get_sel_lvardata(cfg, dataname, lvarnames):
    """Get filenames from cfg for one model and differen variable(s)."""
    selection = []

    for lvar in lvarnames:
        for hlp in select_metadata(cfg['input_data'].values(), long_name=lvar,
                                   dataset=dataname):
            selection.append(hlp['filename'])

    return selection


def read_data_trop_zonal(attributes):
    """Read data from files, average longitudes."""
    logger.debug("Loading %s", attributes['filename'])
    svarcube = iris.load_cube(attributes['filename'])
    svarcube = svarcube.collapsed('longitude', iris.analysis.MEAN)

    return svarcube


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


def get_prof_and_plt_data(cfg, data, available_vars):
    """Plot data for singe data sets and get profile for each."""
    profiles = {}
    projects = {}
    for svar in available_vars:
        profiles[svar] = {}

    
    for svar in available_vars:
        for attributes in data[svar]:
            logger.info("Processing dataset %s", attributes['dataset'])
            dataset = attributes['dataset']
            projects[dataset] = attributes['project']

            svarcube = read_data_trop_zonal(attributes)

            profiles[svar][dataset] = svarcube.collapsed(['time', 'latitude'],
                                                         iris.analysis.MEAN)

            plot_zonal_mean(cfg, svarcube.collapsed(['time'],
                                                    iris.analysis.MEAN),
                            dataset, "Zonal mean ",
                            svarcube.long_name)
                               
            plot_taperec(cfg,
                               svarcube.collapsed('latitude',
                                                      iris.analysis.MEAN),
                               dataset, "Tape recorder ",
                               svarcube.long_name)

    return profiles,projects


def plot_zonal_mean(cfg, mean_cube, dataname, titlestr, variable):
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


def plot_taperec(cfg, mean_cube, dataname, titlestr, variable):
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


def plot_profiles(cfg, profiles, available_vars):
    """Plot zonal mean contour."""
    # Plot data
    # create figure and axes instances
    for svar in available_vars:
        fig, axx = plt.subplots(figsize=(7, 5))

        for iii, dataset in enumerate(list(profiles[svar].keys())):
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
        onedat = profiles[svar][list(profiles[svar].keys())[0]]
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

        provenance_record = get_provenance_record(get_sel_files_var(cfg,
                                                                     svar),
                                                  'Average ' +
                                                  onedat.long_name +
                                                  ' profile', ['mean'],
                                                  ['global'])

        diagnostic_file = get_diagnostic_filename(figname, cfg)
        logger.info("Saving analysis results to %s", diagnostic_file)
        iris.save(profiles_save, target=diagnostic_file)
        logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                    pformat(provenance_record))
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(diagnostic_file, provenance_record)


def plot_logprofiles(cfg, profiles, projects, available_vars):
    """Plot zonal mean contour."""
    # Plot data
    # create figure and axes instances
    for svar in available_vars:
        fig, axx = plt.subplots(figsize=(7, 5))

        for iii, dataset in enumerate(list(profiles[svar].keys())):

            if projects[dataset] == 'CMIP6':
                style = plot.get_dataset_style(dataset, style_file='cmip6')
                lwd = 1
            elif projects[dataset] == 'CMIP5':
                style = plot.get_dataset_style(dataset, style_file='cmip5')
                lwd = 1
            elif dataset == 'ERA5':
                style = {'color': (0, 0, 1, 1.0)}
                lwd = 3
            elif dataset == 'SWOOSH':
                style = {'color': (0, 1, 0, 1.0)}
                lwd = 3
            else:
                style = {'color': (0, 0, 0, 1.0)}
                lwd = 1

            if dataset == 'MultiModelMean':
                style = {'color': (1, 0, 0, 1.0)}
                lwd = 3

            plt.plot((profiles[svar][dataset]).data,
                     (profiles[svar][dataset]).coord('air_pressure').points /
                     100.0, 
                     color=style['color'],
                     linewidth=lwd,
                     label=dataset)
            if iii == 0:
                profiles_save = iris.cube.CubeList([profiles[svar][dataset]])
            else:
                profiles_save.append(profiles[svar][dataset])

        axx = plt.gca()
        axx.invert_yaxis()
        axx.set_yscale('log')
        axx.set_xscale('log')

        plt.legend(loc='upper right')

        axx.set_ylabel('Pressure [hPa]')
        axx.set_ylim(250, 1)
        onedat = profiles[svar][list(profiles[svar].keys())[0]]
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
        figname = 'fig_logprofile_' +\
            onedat.long_name.replace(" ", "_")
        fig.savefig(get_plot_filename(figname, cfg), dpi=300)
        plt.close()

        provenance_record = get_provenance_record(get_sel_files_var(cfg,
                                                                     svar),
                                                  'Average ' +
                                                  onedat.long_name +
                                                  ' profile', ['mean'],
                                                  ['global'])

        diagnostic_file = get_diagnostic_filename(figname, cfg)
        logger.info("Saving analysis results to %s", diagnostic_file)
        iris.save(profiles_save, target=diagnostic_file)
        logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                    pformat(provenance_record))
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(diagnostic_file, provenance_record)


def main(cfg):
    """Read in data for tropopause calculation."""
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))

    logging.debug("Found variables in recipe:\n%s", available_vars)
    # Get input data
    data = {}
    for varname in available_vars:
        data[varname] = select_metadata(cfg['input_data'].values(),
                                        short_name=varname)
        data[varname] = sorted_metadata(data[varname], sort='dataset')

    profiles, projects = get_prof_and_plt_data(cfg, data, available_vars)

    print("projects kw")
    print(projects)

    plot_profiles(cfg, profiles, available_vars)
    plot_logprofiles(cfg, profiles, projects, available_vars)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
