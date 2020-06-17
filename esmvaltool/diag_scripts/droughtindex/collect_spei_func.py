#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Collects data produced by diag_save_spei.r to plot/process them further.

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
# import glob
# import datetime
import iris
from iris.util import rolling_window
from iris.analysis import Aggregator
from iris.time import PartialDateTime
import cf_units as unit
import numpy as np
import cartopy.crs as cart
import matplotlib.pyplot as plt
import matplotlib.dates as mda
# import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n


def runs_of_ones_array(bits):
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    run_starts, = np.where(difs > 0)
    run_ends, = np.where(difs < 0)
    return run_ends - run_starts


def runs_of_ones_array_spei(bits, spei):
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], bits, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    run_starts, = np.where(difs > 0)
    run_ends, = np.where(difs < 0)
    spei_sum = np.full(len(run_starts), 0.5)
    for iii, indexs in enumerate(run_starts):
        spei_sum[iii] = np.sum(spei[indexs:run_ends[iii]])

    return [run_ends - run_starts, spei_sum]

# Define a function to perform the custom statistical operation.
# Note: in order to meet the requirements of iris.analysis.Aggregator,
# it must do the calculation over an arbitrary (given) data axis.


def count_spells(data, threshold, axis):

    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim
        # Threshold the data to find the 'significant' points.
    data = data[:, :, 0, :]
    if axis > 2:
        axis = axis - 1

    listshape = []
    inoax = []
    for iii, ishape in enumerate(data.shape):
        if iii != axis:
            listshape.append(ishape)
            inoax.append(iii)
        # else:
        #     listshape.append(1)

    listshape.append(4)
    return_var = np.zeros(tuple(listshape))
    # return_var[:, :, 0] = spell_point_counts_new
    # return_var[:, :, 1] = spell_point_duration_new
    # return_var[:, :, 2] = spell_point_severity_index
    # return_var[:, :, 3] = spell_point_average_spei

    # spell_point_counts_new = np.zeros((listshape[0], listshape[1]))
    # spell_point_duration_new = np.zeros((listshape[0], listshape[1]))
    # spell_point_severity_index = np.zeros((listshape[0], listshape[1]))
    # spell_point_average_spei = np.zeros((listshape[0], listshape[1]))
    for ilat in range(listshape[0]):
        for ilon in range(listshape[1]):
            if axis == 0:
                data_help = (data[:, ilat, ilon])[:, 0]
            if axis == 1:
                data_help = (data[ilat, :, ilon])[:, 0]
            if axis == 2:
                data_help = data[ilat, ilon, :]
            # if np.all(np.isnan(data_help)):

            # print("data_help.count")
            # print(data_help.count)
            if data_help.count() == 0:
                # print("No data")
                return_var[ilat, ilon, 0] = data_help[0]
                return_var[ilat, ilon, 1] = data_help[0]
                return_var[ilat, ilon, 2] = data_help[0]
                return_var[ilat, ilon, 3] = data_help[0]
            else:
                data_hits = data_help < threshold
                [events, spei_sum] = runs_of_ones_array_spei(data_hits,
                                                           data_help)

                return_var[ilat, ilon, 0] = np.count_nonzero(events)
                # np.sum(np.isfinite(events))
                return_var[ilat, ilon, 1] = np.mean(events)
                return_var[ilat, ilon, 2] = np.mean((spei_sum * events) /
                                                    (np.mean(data_help
                                                             [data_hits])
                                                     * np.mean(events)))
                return_var[ilat, ilon, 3] = np.mean(spei_sum / events)

    # listshape.append(4)
    # return_var = np.zeros(tuple(listshape))
    # return_var[:, :, 0] = spell_point_counts_new
    # return_var[:, :, 1] = spell_point_duration_new
    # return_var[:, :, 2] = spell_point_severity_index
    # return_var[:, :, 3] = spell_point_average_spei
    return return_var


def get_latlon_index(coords, lim1, lim2):
    """Get index for given values (1D vector, e.g. lats or lons)
    between two limits"""
    index = (np.where(np.absolute(coords - (lim2 + lim1)
                                  / 2.0) <= (lim2 - lim1)
                      / 2.0))[0]
    return index


def plot_map_spei_multi(cfg, data_dict, colormap='jet'):
    """Plot contour maps for multi model mean."""
    print("data_dict['data']")
    print(data_dict['data'])
    mask = np.isnan(data_dict['data'])
    spei = np.ma.array(data_dict['data'], mask=mask)
    # np.ma.masked_less_equal(spei, 0)
    print("spei")
    print(spei)

    # Get latitudes and longitudes from cube
    lats = data_dict['latitude']
    lons = data_dict['longitude']
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    # mesh_ind = np.ix_(range(360), index)
    spei = spei[np.ix_(range(len(lats)), index)]

    # Get data set name from cube
    dataset_name = data_dict['datasetname']
    # print(dataset_name)

    # Plot data
    # Create figure and axes instances
    fig, axx = plt.subplots(figsize=(8, 4))

    axx = plt.axes(projection=cart.PlateCarree(central_longitude=0.0))
    axx.set_extent([-180.0, 180.0, -90.0, 90.0],
                   cart.PlateCarree(central_longitude=0.0))

    # np.set_printoptions(threshold=np.nan)

    # Draw filled contours
    cnplot = plt.contourf(lons, lats, spei,
                          data_dict['drought_numbers_level'],
                          transform=cart.PlateCarree(central_longitude=0.0),
                          cmap=colormap, extend='both', corner_mask=False)
    # Draw coastlines
    axx.coastlines()

    # Add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation='horizontal')

    # Add colorbar title string
    cbar.set_label(data_dict['model_kind'] + ' ' + data_dict['drought_char'])

    # Set labels and title to each plot
    axx.set_xlabel('Longitude')
    axx.set_ylabel('Latitude')
    axx.set_title(dataset_name + ' ' + data_dict['model_kind'] + ' '
                  + data_dict['drought_char'])

    # Sets number and distance of x ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    # Sets strings for x ticks
    axx.set_xticklabels(['180°W', '120°W', '60°W',
                         '0°', '60°E', '120°E',
                         '180°E'])
    # Sets number and distance of y ticks
    axx.set_yticks(np.linspace(-90, 90, 7))
    # Sets strings for y ticks
    axx.set_yticklabels(['90°S', '60°S', '30°S', '0°',
                         '30°N', '60°N', '90°N'])

    fig.tight_layout()
    fig.savefig(os.path.join(cfg[n.PLOT_DIR],
                             'spei_map' + data_dict['filename'] + '_' +
                             dataset_name + '.' +
                             cfg[n.OUTPUT_FILE_TYPE]))
    plt.close()


def plot_map_spei(cfg, cube, levels, add_to_filename='', name=''):
    """Plot contour map."""
    print("hello map 1")
    # SPEI array to plot
    # spei = cube.data
    print("hello map 2")

    mask = np.isnan(cube.data)
    spei = np.ma.array(cube.data, mask=mask)
    np.ma.masked_less_equal(spei, 0)

    # Get latitudes and longitudes from cube
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    # mesh_ind = np.ix_(range(360), index)
    spei = spei[np.ix_(range(len(lats)), index)]

    # arr1inds = lons.argsort()
    # lons = lons[arr1inds]
    # spei = spei[:,arr1inds,:]

    # Get data set name from cube
    print("cube.metadata.attributes")
    print(cube.metadata.attributes)
    
    # Get data set name from cube    
    try:
        dataset_name = cube.metadata.attributes['model_id']
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes['source_id']
        except KeyError:
            dataset_name = 'Observations'
    # try:
    #     dataset_name = cube.metadata.attributes['model_id']
    # except KeyError:
    #    dataset_name = 'Observations'
    #     # continue
    # dataset_name = cube.metadata.attributes['project_id']
    # if dataset_name in ['CMIP5', 'CMIP6']:
    #     dataset_name = cube.metadata.attributes['model_id']


    print(dataset_name)

    # Plot data
    # Create figure and axes instances
    fig, axx = plt.subplots(figsize=(8, 4))

    axx = plt.axes(projection=cart.PlateCarree(central_longitude=0.0))
    axx.set_extent([-180.0, 180.0, -90.0, 90.0],
                   cart.PlateCarree(central_longitude=0.0))

    # np.set_printoptions(threshold=np.nan)

    # Draw filled contours
    cnplot = plt.contourf(lons, lats, spei,
                          levels,
                          transform=cart.PlateCarree(central_longitude=0.0),
                          cmap='gnuplot', extend='both', corner_mask=False)
    # Draw coastlines
    axx.coastlines()

    # Add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation='horizontal')

    # Add colorbar title string
    cbar.set_label(name)

    # Set labels and title to each plot
    axx.set_xlabel('Longitude')
    axx.set_ylabel('Latitude')
    axx.set_title(dataset_name + ' ' + name)

    # Sets number and distance of x ticks
    axx.set_xticks(np.linspace(-180, 180, 7))
    # Sets strings for x ticks
    axx.set_xticklabels(['180°W', '120°W', '60°W',
                         '0°', '60°E', '120°E',
                         '180°E'])
    # Sets number and distance of y ticks
    axx.set_yticks(np.linspace(-90, 90, 7))
    # Sets strings for y ticks
    axx.set_yticklabels(['90°S', '60°S', '30°S', '0°',
                         '30°N', '60°N', '90°N'])

    fig.tight_layout()
    fig.savefig(os.path.join(cfg[n.PLOT_DIR],
                             'spei_map' + add_to_filename + '_' +
                             dataset_name + '.' +
                             cfg[n.OUTPUT_FILE_TYPE]))
    plt.close()


def plot_time_series_spei(cfg, cube, add_to_filename=''):
    """Plot time series."""

    # SPEI vector to plot
    spei = cube.data
    # Get time from cube
    time = cube.coord('time').points
    # Adjust (ncdf) time to the format matplotlib expects
    add_m_delta = mda.datestr2num('1950-01-01 00:00:00')
    time = time + add_m_delta

    # Get data set name from cube    
    try:
        dataset_name = cube.metadata.attributes['model_id']
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes['source_id']
        except KeyError:
            dataset_name = 'Observations'

    
    # Get data set name from cube
    # try:
    #     dataset_name = cube.metadata.attributes['model_id']
    # except KeyError:
    #     dataset_name = 'Observations'
    # dataset_name = cube.metadata.attributes['project_id']
    # if dataset_name in ['CMIP5', 'CMIP6']:
    #     dataset_name = cube.metadata.attributes['model_id']

    # Plot data
    # Create figure and axes instances
    fig, axx = plt.subplots(figsize=(16, 4))
    axx.plot_date(time, spei, '-', tz=None, xdate=True, ydate=False,
                  color='r', linewidth=4., linestyle='-', alpha=1.,
                  marker='x')
    axx.axhline(y=-2, color='k')

    # Plot labels and title
    axx.set_xlabel('Time')
    axx.set_ylabel('SPEI')
    axx.set_title('Mean SPEI ' + dataset_name + ' ' + add_to_filename)

    # Set limits for y-axis
    axx.set_ylim(-4.0, 4.0)

    # Often improves the layout
    fig.tight_layout()
    # Save plot to file
    fig.savefig(os.path.join(cfg[n.PLOT_DIR],
                             'spei_line' + add_to_filename + '_' + dataset_name
                             + '.' +
                             cfg[n.OUTPUT_FILE_TYPE]))
    plt.close()
