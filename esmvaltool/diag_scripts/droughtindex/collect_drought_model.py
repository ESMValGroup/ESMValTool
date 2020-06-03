#!/usr/bin/env python
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
from iris.analysis import Aggregator
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_drought_func import (
    _plot_multi_model_maps, _plot_single_maps, count_spells)


def _make_new_cube(tscube, number_drought_charac):
    """Make a new cube with an extra dimension for result of spell count."""
    # make a new cube to increase the size of the data array
    # get two (instead of one) values from the aggregator spell_no
    new_shape = tscube.shape + (number_drought_charac,)
    new_data = iris.util.broadcast_to_shape(
        tscube.data, new_shape, [0, 1, 2])
    new_cube = iris.cube.Cube(new_data)

    new_cube.add_dim_coord(iris.coords.DimCoord(
        tscube.coord('time').points, long_name='time'), 0)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        tscube.coord('latitude').points, long_name='latitude'), 1)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        tscube.coord('longitude').points, long_name='longitude'), 2)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        np.arange(0, number_drought_charac, 1), long_name='z'), 3)
    return new_cube


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
    spell_no = Aggregator('spell_count', count_spells,
                          units_func=lambda units: 1)

    # Define the parameters of the test.
    threshold_spei = -2.0
    number_drought_charac = 4
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
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        time = cube.coord('time')
        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time.
        cube2 = cube.collapsed('time', iris.analysis.MEAN)  # 3D to 2D

        if first_run == 1:
            files = os.listdir((cfg[n.INPUT_FILES])[0])
            ncfiles = list(filter(lambda f: f.endswith('.nc'), files))
            shape_all = cube2.data.shape + (number_drought_charac,) + \
                (len(ncfiles),)
            all_drought_hist = np.full(shape_all, np.nan)
            all_drought_rcp85 = np.full(shape_all, np.nan)
            first_run = 0
        # Test if time series goes until cfg['end_year']/12
        timecheck = time.units.date2num(datetime.datetime(cfg['end_year'],
                                                          11, 30, 0, 0, 0))

        if cube.coord('time').points[-1] > timecheck:
            tscube = _set_tscube(cfg, cube, time, 'Historic')
            new_cube = _make_new_cube(tscube, number_drought_charac)
            # calculate the number of drought events and average duration
            drought_show = new_cube.collapsed('time', spell_no,
                                              threshold=threshold_spei)
            drought_show.rename('Drought characteristics')
            time_len = len(new_cube.coord('time').points) / 12.0
            # Convert number of droughtevents to frequency (per year)
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_hist[:, :, :, iii] = drought_show.data
            _plot_single_maps(cfg, cube2, drought_show, 'Historic')

            tscube = _set_tscube(cfg, cube, time, 'Future')
            new_cube = _make_new_cube(tscube, number_drought_charac)
            # calculate the number of drought events and their average duration
            drought_show = new_cube.collapsed('time', spell_no,
                                              threshold=threshold_spei)
            drought_show.rename('Drought characteristics')
            # length of time series
            time_len = len(new_cube.coord('time').points) / 12.0
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_rcp85[:, :, :, iii] = drought_show.data
            _plot_single_maps(cfg, cube2, drought_show, 'Future')

    # Calculating multi model mean and plot it
    all_drought_hist_mean = np.nanmean(all_drought_hist, axis=-1)
    # to 3D multi model mean
    all_drought_rcp85_mean = np.nanmean(all_drought_rcp85, axis=-1)
    # to 3D multi model mean
    perc_diff = ((all_drought_rcp85_mean - all_drought_hist_mean)
                 / (all_drought_rcp85_mean + all_drought_hist_mean) * 200)

    # Plot multi model means
    _plot_multi_model_maps(cfg, all_drought_hist_mean, lats, lons, 'Historic')
    _plot_multi_model_maps(cfg, all_drought_rcp85_mean, lats, lons, 'Future')
    _plot_multi_model_maps(cfg, perc_diff, lats, lons, 'Difference')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
