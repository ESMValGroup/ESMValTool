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
# from iris.util import rolling_window
from iris.analysis import Aggregator
# from iris.time import PartialDateTime
# import cf_units as unit
import numpy as np
# import cartopy.crs as cart
# import matplotlib.pyplot as plt
# import matplotlib.dates as mda
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_drought_func import (
    count_spells, plot_map_spei_multi, plot_map_spei)


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
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN)  # 3D to 2D

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
        lasttime = cube.coord('time').points[-1]

        if lasttime > timecheck:
            # extract time series from historical model data
            # cfg['start_year'] to cfg['start_year'] + cfg['comparison_period']
            start = datetime.datetime(cfg['start_year'], 1, 15, 0, 0, 0)
            end = datetime.datetime(cfg['start_year'] +
                                    cfg['comparison_period'], 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

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
            drought_numbers_level = np.arange(0, 0.4, 0.05)
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='Historic_No_of_Events_per_year',
                          name='Historic_Number of Events per year')

            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='Historic_Dur_of_Events',
                          name='Historic_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='Historic_Sev_index_of_Events',
                          name='Historic_Severity Index of Events')

            # plot the average spei of drought events
            drought_numbers_level = np.arange(-2.8, -1.8, 0.2)
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='Historic_Avr_' +
                          cfg['indexname'] + '_of_Events',
                          name='Historic_Average ' + cfg['indexname'] +
                          ' of Events')
            # extract time series from rcp model data
            # cfg['end_year'] - cfg['comparison_period'] to cfg['end_year']
            start = datetime.datetime(cfg['end_year'] -
                                      cfg['comparison_period'], 1, 15, 0, 0, 0)
            end = datetime.datetime(cfg['end_year'], 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

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

            # plot the number of drought events
            drought_numbers_level = np.arange(0, 0.4, 0.05)
            # set color levels
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='RCP85_No_of_Events_per_year',
                          name='RCP85_Number of Events per year')

            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='RCP85_Dur_of_Events',
                          name='RCP85_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='RCP85_Sev_index_of_Events',
                          name='RCP85_Severity Index of Events')

            # plot the average spei of drought events
            drought_numbers_level = np.arange(-2.8, -1.8, 0.2)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                          add_to_filename='RCP85_Avr_' + cfg['indexname'] +
                          '_of_Events',
                          name='RCP85_Average ' + cfg['indexname'] +
                          ' of Events')
    # Calculating multi model mean and plot it
    all_drought_hist_mean = np.nanmean(all_drought_hist, axis=-1)
    # to 3D multi model mean
    all_drought_rcp85_mean = np.nanmean(all_drought_rcp85, axis=-1)
    # to 3D multi model mean
    perc_diff = ((all_drought_rcp85_mean - all_drought_hist_mean)
                 / (all_drought_rcp85_mean + all_drought_hist_mean) * 200)
    print("all_drought_rcp85_mean")
    print(all_drought_rcp85_mean)

    # Historic
    data_dict = {}
    data_dict['data'] = all_drought_hist_mean[:, :, 0]
    data_dict['datasetname'] = 'MultiModelMean'
    data_dict['latitude'] = lats
    data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Historic'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'Historic_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'Historic_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Historic_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 3]
    data_dict['drought_char'] = 'Average ' + cfg['indexname'] + ' of Events'
    data_dict['filename'] = 'Historic_Average_' + cfg['indexname'] + \
        '_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # RCP85
    data_dict['data'] = all_drought_rcp85_mean[:, :, 0]
    data_dict['model_kind'] = 'RCP85'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'RCP85_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'RCP85_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'RCP85_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 3]
    data_dict['drought_char'] = 'Average ' + cfg['indexname'] + ' of Events'
    data_dict['filename'] = 'RCP85_Avr_' + cfg['indexname'] + '_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # RCP85 Percentage difference
    data_dict['data'] = perc_diff[:, :, 0]
    data_dict['datasetname'] = 'Percentage'
    # data_dict['latitude'] = lats
    # data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Difference'
    data_dict['drought_char'] = 'Number of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_No_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-100, 110, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-100, 110, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Sev_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 3]
    data_dict['drought_char'] = 'Average ' + cfg['indexname'] + \
        ' of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Avr_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
