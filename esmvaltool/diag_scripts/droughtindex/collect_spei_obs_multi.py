#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Collects data produced by diag_save_spei.r to plot/process them further.

###############################################################################
droughtindex/collect_spei.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Collects data produced by diag_save_spei.r to plot/process them further.

Configuration options
---------------------
    TBD

###############################################################################

"""
import os
import glob
import iris
from iris.analysis import Aggregator
import datetime
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_spei_func import (
    get_latlon_index, count_spells, plot_map_spei_multi,
    plot_map_spei, plot_time_series_spei)


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    #######################################################################
    # Read recipe data
    #######################################################################

    # Get filenames of input files produced by diag_spei.r
    # "cfg[n.INPUT_FILES]" is produced by the ESMValTool and contains
    # information on the SPEI input files produced by diag_spei.r
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*.nc"
    threshold_spei = -2.0
    number_drought_charac = 4
    first_run = 1
    iobs = 0
    # For loop: "glob.iglob" findes all files which match the
    # pattern of "input_filenames".
    # It writes the resulting exact file name onto spei_file
    # and runs the following indented lines for all possibilities
    # for spei_file.
    for iii, spei_file in enumerate(glob.iglob(input_filenames)):
        # Loads the file into a special structure (IRIS cube)
        cube = iris.load(spei_file)[0]
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        time = cube.coord('time')

        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time, i.e. data you need
        # to plot the average over the time series of SPEI on a map
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN)

        if first_run == 1:
            files = os.listdir((cfg[n.INPUT_FILES])[0])
            ncfiles = list(filter(lambda f: f.endswith('.nc'),files))
            shape_all = cube2.data.shape + (number_drought_charac,) + \
                (len(ncfiles) -1 ,)
            all_drought = np.full(shape_all, np.nan)
            first_run = 0

        coords = ('longitude', 'latitude')

        # Bremen
        index_lat = get_latlon_index(lats, 52, 53)
        index_lon = get_latlon_index(lons, 7, 9)
        add_to_filename = 'Bremen'
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        cube_grid_areas = iris.analysis.cartography.area_weights(
            cube[:, index_lat[0]:index_lat[-1] + 1,
                 index_lon[0]:index_lon[-1] + 1])
        cube4 = ((cube[:, index_lat[0]:index_lat[-1] + 1,
                 index_lon[0]:index_lon[-1] + 1]).collapsed(coords,
                                                iris.analysis.MEAN,
                                                weights=cube_grid_areas))

        plot_time_series_spei(cfg, cube4, add_to_filename)

        index_lat = get_latlon_index(lats, 7, 9)
        index_lon = get_latlon_index(lons, 8, 10)
        add_to_filename = 'Nigeria'
        cube_grid_areas = iris.analysis.cartography.area_weights(
            cube[:, index_lat[0]:index_lat[-1] + 1,
                 index_lon[0]:index_lon[-1] + 1])
        cube5 = (cube[:, index_lat[0]:index_lat[-1] + 1,
                 index_lon[0]:index_lon[-1]
                 + 1]).collapsed(coords,
                                 iris.analysis.MEAN, weights=cube_grid_areas)

        start = datetime.datetime(1905, 1, 15, 0, 0, 0)
        end = datetime.datetime(2005, 12, 16, 0, 0, 0)
        time = cube5.coord('time')

        plot_time_series_spei(cfg, cube5, add_to_filename)

        # Make an aggregator from the user function.
        spell_no = Aggregator('spell_count',
                              count_spells,
                              units_func=lambda units: 1)

        # extract time series from 1901-2001 observ and hist. model data
        start = datetime.datetime(1905, 1, 15, 0, 0, 0)
        end = datetime.datetime(2005, 12, 16, 0, 0, 0)
        time = cube.coord('time')

        # make a new cube to increase the size of the data array
        new_shape = cube.shape + (4,)
        new_data = iris.util.broadcast_to_shape(
            cube.data, new_shape, [0, 1, 2])
        new_cube = iris.cube.Cube(new_data)

        new_cube.add_dim_coord(iris.coords.DimCoord(
            cube.coord('time').points, long_name='time'), 0)
        new_cube.add_dim_coord(iris.coords.DimCoord(
            cube.coord('latitude').points, long_name='latitude'), 1)
        new_cube.add_dim_coord(iris.coords.DimCoord(
            cube.coord('longitude').points, long_name='longitude'), 2)
        new_cube.add_dim_coord(iris.coords.DimCoord(
            [0, 1, 2, 3], long_name='z'), 3)

        # calculate the number of drought events and their average duration
        drought_show = new_cube.collapsed('time', spell_no,
                                          threshold=threshold_spei)
        drought_show.rename('Number of 1 month drought')
        # length of time series
        time_length = len(new_cube.coord('time').points) / 12.0
        # Convert number of droughtevents to frequency (per year)
        drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                       0] / time_length
        # dataset_name = cube.metadata.attributes['model_id']
        try:
            dataset_name = cube.metadata.attributes['model_id']
            all_drought[:, :, :, iii - iobs] = drought_show.data
        except KeyError:
            try:
                dataset_name = cube.metadata.attributes['source_id']
                all_drought[:, :, :, iii - iobs] = drought_show.data
            except KeyError:
                dataset_name = 'Observations'
                all_drought_obs = drought_show.data
                iobs = 1

        drought_numbers_level = np.arange(0, 0.4, 0.05)
        # length of time series
        # time_length = len(new_cube.coord('time').points) / 12.0
        cube2.data = drought_show.data[:, :, 0]
        plot_map_spei(cfg, cube2, drought_numbers_level,
                     add_to_filename='No_of_Events',
                     name='Number of Events per year')

        # plot the average duration of drought events
        drought_numbers_level = np.arange(0, 7, 1)
        # use cube2 to get metadata
        cube2.data = drought_show.data[:, :, 1]
        plot_map_spei(cfg, cube2, drought_numbers_level,
                     add_to_filename='Dur_of_Events',
                     name='Duration of Events(month)')

        # plot the average severity index of drought events
        drought_numbers_level = np.arange(0, 9, 1)
        # use cube2 to get metadata
        cube2.data = drought_show.data[:, :, 2]
        plot_map_spei(cfg, cube2, drought_numbers_level,
                     add_to_filename='Sev_index_of_Events',
                     name='Severity Index of Events')

        # plot the average spei of drought events
        drought_numbers_level = np.arange(-2.8, -1.8, 0.2)
        # use cube2 to get metadata
        cube2.data = drought_show.data[:, :, 3]
        plot_map_spei(cfg, cube2, drought_numbers_level,
                     add_to_filename='Avr_spei_of_Events',
                     name='Average spei of Events')

    # Calculating multi model mean and plot it
    all_drought_hist_mean = np.nanmean(all_drought, axis=-1)
    perc_diff = ((all_drought_obs - all_drought_hist_mean)
                 / (all_drought_obs + all_drought_hist_mean) * 200)
    print("all_drought_obs")
    print(all_drought_obs)

    # Historic MultiModelMean
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
    data_dict['drought_numbers_level'] = np.arange(0, 7, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Historic_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'Historic_Avr_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # CRU_OBS MultiModelMean
    data_dict['data'] = all_drought_obs[:, :, 0]
    data_dict['datasetname'] = 'Observations'
    data_dict['model_kind'] = 'Era-Interim'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'Era-Interim_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_obs[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'Era-Interim_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 7, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_obs[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Era-Interim_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_obs[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'Era-Interim_Avr_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # Perc_diff Multimodelmean
    data_dict['data'] = perc_diff[:, :, 0]
    data_dict['datasetname'] = 'Percentage'
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
    data_dict['filename'] = 'Percentage_difference_of_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Avr_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
