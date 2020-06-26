#!/usr/bin/env  python
# -*- coding: utf-8 -*-


"""Collects SPI or SPEI data comparing historic and future model scenarios.

Applies drought characteristics based on Martin (2018).

###############################################################################
droughtindex/collect_spei.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Collects data produced by diag_save_spi.R or diad_save_spei_all.R
    to plot/process them further.

Configuration options
---------------------
    indexname: "SPI" or "SPEI"
    start_year: year, start of historical time series
    end_year: year, end of future scenario
    comparison_period: should be < (end_year - start_year)/2

###############################################################################

"""
# The imported modules and libraries below allow this script to run accordingly

import os
import glob
import datetime
import iris
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_drought_func import (
    _get_drought_data, _plot_multi_model_maps, _plot_single_maps)


def _get_and_plot_multimodel(cfg, cube, all_drought, input_filenames):
    """Calculate multi-model mean and compare historic and future."""
    all_drought_mean = {}
    for tstype in ['Historic', 'Future']:
        all_drought_mean[tstype] = np.nanmean(all_drought[tstype], axis=-1)

    all_drought_mean['Difference'] = ((all_drought_mean['Future'] -
                                       all_drought_mean['Historic']) /
                                      (all_drought_mean['Future'] +
                                       all_drought_mean['Historic']) * 200)

    # Plot multi model means
    for tstype in ['Historic', 'Future', 'Difference']:
        _plot_multi_model_maps(cfg, all_drought_mean[tstype],
                               [cube.coord('latitude').points,
                                cube.coord('longitude').points],
                               input_filenames,
                               tstype)


def _set_tscube(cfg, cube, time, tstype):
    """Time slice from a cube with start/end given by cfg."""
    if tstype == 'Future':
        # extract time series from rcp model data
        # cfg['end_year'] - cfg['comparison_period'] to cfg['end_year']
        start = datetime.datetime(cfg['end_year'] -
                                  cfg['comparison_period'], 1, 15, 0, 0, 0)
        end = datetime.datetime(cfg['end_year'], 12, 16, 0, 0, 0)
    elif tstype == 'Historic':
        # extract time series from historical model data
        # cfg['start_year'] to cfg['start_year'] + cfg['comparison_period']
        start = datetime.datetime(cfg['start_year'], 1, 15, 0, 0, 0)
        end = datetime.datetime(cfg['start_year'] +
                                cfg['comparison_period'], 12, 16, 0, 0, 0)
    stime = time.nearest_neighbour_index(time.units.date2num(start))
    etime = time.nearest_neighbour_index(time.units.date2num(end))
    tscube = cube[stime:etime, :, :]
    return tscube


def main(cfg):
    """Run the diagnostic.

    Parameters :

    ------------
    cfg : dict
    """
    ######################################################################
    # Read recipe data
    ######################################################################

    # Make an aggregator from the user function.
    # spell_no = Aggregator('spell_count', count_spells,
    #                       units_func=lambda units: 1)

    # Define the parameters of the test.
    first_run = 1

    # Get filenames of input files produced by diag_spei.r
    # input_filenames = (cfg[n.INPUT_FILES])[0] + "/*.nc"
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*_" + \
        (cfg['indexname']).lower() + "_*.nc"

    for iii, spei_file in enumerate(glob.iglob(input_filenames)):
        # Loads the file into a special structure (IRIS cube),
        # which allows us to access the data and additional information
        # with python.
        cube = iris.load(spei_file)[0]
        # lats = cube.coord('latitude').points
        # lons = cube.coord('longitude').points
        time = cube.coord('time')
        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time.
        cube2 = cube.collapsed('time', iris.analysis.MEAN)  # 3D to 2D

        if first_run == 1:
            ncfiles = list(filter(lambda f: f.endswith('.nc'),
                                  os.listdir((cfg[n.INPUT_FILES])[0])))
            all_drought = {}
            all_drought['Historic'] = np.full(cube2.data.shape + (4,) +
                                              (len(ncfiles),), np.nan)
            all_drought['Future'] = np.full(cube2.data.shape + (4,) +
                                            (len(ncfiles),), np.nan)
            first_run = 0
        # Test if time series goes until cfg['end_year']/12
        timecheck = time.units.date2num(datetime.datetime(cfg['end_year'],
                                                          11, 30, 0, 0, 0))

        if time.points[-1] > timecheck:
            for tstype in ['Historic', 'Future']:
                tscube = _set_tscube(cfg, cube, time, tstype)
                drought_show = _get_drought_data(cfg, tscube)
                all_drought[tstype][:, :, :, iii] = drought_show.data
                _plot_single_maps(cfg, cube2, drought_show, tstype, spei_file)

    # Calculating multi model mean and plot it
    _get_and_plot_multimodel(cfg, cube, all_drought,
                             glob.glob(input_filenames))


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
