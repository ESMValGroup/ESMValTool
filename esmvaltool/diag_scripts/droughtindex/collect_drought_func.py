#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Drought characteristics and plots based on  Martin (2018).

###############################################################################
droughtindex/collect_drought_obs_multi.py
Author: Katja Weigel, Kemisola Adeniyi (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Functions for:
        collect_drought_obs_multi.py and droughtindex/collect_drought_model.py.

Configuration options
---------------------
    None

###############################################################################

"""


import logging
import os
from pprint import pformat
import numpy as np
import iris
from iris.analysis import Aggregator
import cartopy.crs as cart
import matplotlib.pyplot as plt
import matplotlib.dates as mda
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename)

logger = logging.getLogger(os.path.basename(__file__))


def _get_data_hlp(axis, data, ilat, ilon):
    """Get data_help dependend on axis."""
    if axis == 0:
        data_help = (data[:, ilat, ilon])[:, 0]
    elif axis == 1:
        data_help = (data[ilat, :, ilon])[:, 0]
    elif axis == 2:
        data_help = data[ilat, ilon, :]

    return data_help


def _get_drought_data(cfg, cube):
    """Prepare data and calculate characteristics."""
    # make a new cube to increase the size of the data array
    # Make an aggregator from the user function.
    spell_no = Aggregator('spell_count',
                          count_spells,
                          units_func=lambda units: 1)
    new_cube = _make_new_cube(cube)

    # calculate the number of drought events and their average duration
    drought_show = new_cube.collapsed('time', spell_no,
                                      threshold=cfg['threshold'])
    drought_show.rename('Drought characteristics')
    # length of time series
    time_length = len(new_cube.coord('time').points) / 12.0
    # Convert number of droughtevents to frequency (per year)
    drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                   0] / time_length
    return drought_show


def _provenance_map_spei(cfg, name_dict, spei, dataset_name):
    """Set provenance for plot_map_spei."""
    caption = 'Global map of ' + \
              name_dict['drought_char'] + \
              ' [' + name_dict['unit'] + '] ' + \
              'based on ' + cfg['indexname'] + '.'

    if cfg['indexname'].lower == "spei":
        set_refs = ['martin18grl', 'vicente10jclim', ]
    elif cfg['indexname'].lower == "spi":
        set_refs = ['martin18grl', 'mckee93proc', ]
    else:
        set_refs = ['martin18grl', ]

    provenance_record = get_provenance_record([name_dict['input_filenames']],
                                              caption,
                                              ['global'],
                                              set_refs)

    diagnostic_file = get_diagnostic_filename(cfg['indexname'] + '_map' +
                                              name_dict['add_to_filename'] +
                                              '_' +
                                              dataset_name, cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    cubesave = cube_to_save_ploted(spei, name_dict)
    iris.save(cubesave, target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames):
    """Set provenance for plot_map_spei_multi."""
    caption = 'Global map of the multi-model mean of ' + \
              data_dict['drought_char'] + \
              ' [' + data_dict['unit'] + '] ' + \
              'based on ' + cfg['indexname'] + '.'

    if cfg['indexname'].lower == "spei":
        set_refs = ['martin18grl', 'vicente10jclim', ]
    elif cfg['indexname'].lower == "spi":
        set_refs = ['martin18grl', 'mckee93proc', ]
    else:
        set_refs = ['martin18grl', ]

    provenance_record = get_provenance_record(input_filenames, caption,
                                              ['global'],
                                              set_refs)

    diagnostic_file = get_diagnostic_filename(cfg['indexname'] + '_map' +
                                              data_dict['filename'] + '_' +
                                              data_dict['datasetname'], cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_ploted(spei, data_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def _provenance_time_series_spei(cfg, data_dict):
    """Provenance for time series plots."""
    caption = 'Time series of ' + \
              data_dict['var'] + \
              ' at' + data_dict['area'] + '.'

    if cfg['indexname'].lower == "spei":
        set_refs = ['vicente10jclim', ]
    elif cfg['indexname'].lower == "spi":
        set_refs = ['mckee93proc', ]
    else:
        set_refs = ['martin18grl', ]

    provenance_record = get_provenance_record([data_dict['filename']],
                                              caption,
                                              ['reg'], set_refs,
                                              plot_type='times')

    diagnostic_file = get_diagnostic_filename(cfg['indexname'] +
                                              '_time_series_' +
                                              data_dict['area'] +
                                              '_' +
                                              data_dict['dataset_name'], cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    cubesave = cube_to_save_ploted_ts(data_dict)
    iris.save(cubesave, target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def cube_to_save_ploted(var, data_dict):
    """Create cube to prepare plotted data for saving to netCDF."""
    plot_cube = iris.cube.Cube(var, var_name=data_dict['var'],
                               long_name=data_dict['drought_char'],
                               units=data_dict['unit'])
    plot_cube.add_dim_coord(iris.coords.DimCoord(data_dict['latitude'],
                                                 var_name='lat',
                                                 long_name='latitude',
                                                 units='degrees_north'), 0)
    plot_cube.add_dim_coord(iris.coords.DimCoord(data_dict['longitude'],
                                                 var_name='lon',
                                                 long_name='longitude',
                                                 units='degrees_east'), 1)

    return plot_cube


def cube_to_save_ploted_ts(data_dict):
    """Create cube to prepare plotted time series for saving to netCDF."""
    plot_cube = iris.cube.Cube(data_dict['data'], var_name=data_dict['var'],
                               long_name=data_dict['var'],
                               units=data_dict['unit'])
    plot_cube.add_dim_coord(iris.coords.DimCoord(data_dict['time'],
                                                 var_name='time',
                                                 long_name='Time',
                                                 units='month'), 0)

    return plot_cube


def get_provenance_record(ancestor_files, caption,
                          domains, refs, plot_type='geo'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': domains,
        'plot_type': plot_type,
        'themes': ['phys'],
        'authors': [
            'weigel_katja',
            'adeniyi_kemisola',
        ],
        'references': refs,
        'ancestors': ancestor_files,
    }
    return record


def _make_new_cube(cube):
    """Make a new cube with an extra dimension for result of spell count."""
    new_shape = cube.shape + (4,)
    new_data = iris.util.broadcast_to_shape(cube.data, new_shape, [0, 1, 2])
    new_cube = iris.cube.Cube(new_data)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        cube.coord('time').points, long_name='time'), 0)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        cube.coord('latitude').points, long_name='latitude'), 1)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        cube.coord('longitude').points, long_name='longitude'), 2)
    new_cube.add_dim_coord(iris.coords.DimCoord(
        [0, 1, 2, 3], long_name='z'), 3)
    return new_cube


def _plot_multi_model_maps(cfg, all_drought_mean, lats_lons, input_filenames,
                           tstype):
    """Prepare plots for multi-model mean."""
    data_dict = {'latitude': lats_lons[0],
                 'longitude': lats_lons[1],
                 'model_kind': tstype
                 }
    if tstype == 'Difference':
        # RCP85 Percentage difference
        data_dict.update({'data': all_drought_mean[:, :, 0],
                          'var': 'diffnumber',
                          'datasetname': 'Percentage',
                          'drought_char': 'Number of Events',
                          'unit': '%',
                          'filename': 'Percentage_difference_of_No_of_Events',
                          'drought_numbers_level': np.arange(-100, 110, 10)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='rainbow')

        data_dict.update({'data': all_drought_mean[:, :, 1],
                          'var': 'diffduration',
                          'drought_char': 'Duration of Events',
                          'filename': 'Percentage_difference_of_Dur_of_Events',
                          'drought_numbers_level': np.arange(-100, 110, 10)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='rainbow')

        data_dict.update({'data': all_drought_mean[:, :, 2],
                          'var': 'diffseverity',
                          'drought_char': 'Severity Index of Events [%]',
                          'filename': 'Percentage_difference_of_Sev_of_Events',
                          'drought_numbers_level': np.arange(-50, 60, 10)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='rainbow')

        data_dict.update({'data': all_drought_mean[:, :, 3],
                          'var': 'diff' + (cfg['indexname']).lower(),
                          'drought_char': 'Average ' + cfg['indexname'] +
                                          ' of Events',
                          'filename': 'Percentage_difference_of_Avr_of_Events',
                          'drought_numbers_level': np.arange(-50, 60, 10)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='rainbow')
    else:
        data_dict.update({'data': all_drought_mean[:, :, 0],
                          'var': 'frequency',
                          'unit': 'year-1',
                          'drought_char': 'Number of Events per year',
                          'filename': tstype + '_No_of_Events_per_year',
                          'drought_numbers_level': np.arange(0, 0.4, 0.05)})
        if tstype == 'Observations':
            data_dict['datasetname'] = 'Mean'
        else:
            data_dict['datasetname'] = 'MultiModelMean'
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='gnuplot')

        data_dict.update({'data': all_drought_mean[:, :, 1],
                          'var': 'duration',
                          'unit': 'month',
                          'drought_char': 'Duration of Events [month]',
                          'filename': tstype + '_Dur_of_Events',
                          'drought_numbers_level': np.arange(0, 6, 1)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='gnuplot')

        data_dict.update({'data': all_drought_mean[:, :, 2],
                          'var': 'severity',
                          'unit': '1',
                          'drought_char': 'Severity Index of Events',
                          'filename': tstype + '_Sev_index_of_Events',
                          'drought_numbers_level': np.arange(0, 9, 1)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='gnuplot')
        namehlp = 'Average ' + cfg['indexname'] + ' of Events'
        namehlp2 = tstype + '_Average_' + cfg['indexname'] + '_of_Events'
        data_dict.update({'data': all_drought_mean[:, :, 3],
                          'var': (cfg['indexname']).lower(),
                          'unit': '1',
                          'drought_char': namehlp,
                          'filename': namehlp2,
                          'drought_numbers_level': np.arange(-2.8, -1.8, 0.2)})
        plot_map_spei_multi(cfg, data_dict, input_filenames,
                            colormap='gnuplot')


def _plot_single_maps(cfg, cube2, drought_show, tstype, input_filenames):
    """Plot map of drought characteristics for individual models and times."""
    cube2.data = drought_show.data[:, :, 0]
    name_dict = {'add_to_filename': tstype + '_No_of_Events_per_year',
                 'name': tstype + ' Number of Events per year',
                 'var': 'frequency',
                 'unit': 'year-1',
                 'drought_char': 'Number of Events per year',
                 'input_filenames': input_filenames}
    plot_map_spei(cfg, cube2, np.arange(0, 0.4, 0.05),
                  name_dict)

    # plot the average duration of drought events
    cube2.data = drought_show.data[:, :, 1]
    name_dict.update({'add_to_filename': tstype + '_Dur_of_Events',
                      'name': tstype + ' Duration of Events(month)',
                      'var': 'duration',
                      'unit': 'month',
                      'drought_char': 'Number of Events per year',
                      'input_filenames': input_filenames})
    plot_map_spei(cfg, cube2, np.arange(0, 6, 1), name_dict)

    # plot the average severity index of drought events
    cube2.data = drought_show.data[:, :, 2]
    name_dict.update({'add_to_filename': tstype + '_Sev_index_of_Events',
                      'name': tstype + ' Severity Index of Events',
                      'var': 'severity',
                      'unit': '1',
                      'drought_char': 'Number of Events per year',
                      'input_filenames': input_filenames})
    plot_map_spei(cfg, cube2, np.arange(0, 9, 1), name_dict)

    # plot the average spei of drought events
    cube2.data = drought_show.data[:, :, 3]

    namehlp = tstype + '_Avr_' + cfg['indexname'] + '_of_Events'
    namehlp2 = tstype + '_Average_' + cfg['indexname'] + '_of_Events'
    name_dict.update({'add_to_filename': namehlp,
                      'name': namehlp2,
                      'var': 'severity',
                      'unit': '1',
                      'drought_char': 'Number of Events per year',
                      'input_filenames': input_filenames})
    plot_map_spei(cfg, cube2, np.arange(-2.8, -1.8, 0.2), name_dict)


def runs_of_ones_array_spei(bits, spei):
    """Set 1 at beginning ond -1 at the end of events."""
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


def count_spells(data, threshold, axis):
    """Functions for Iris Aggregator to count spells."""
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim
    data = data[:, :, 0, :]
    if axis > 2:
        axis = axis - 1

    listshape = []
    inoax = []
    for iii, ishape in enumerate(data.shape):
        if iii != axis:
            listshape.append(ishape)
            inoax.append(iii)

    listshape.append(4)
    return_var = np.zeros(tuple(listshape))

    for ilat in range(listshape[0]):
        for ilon in range(listshape[1]):
            data_help = _get_data_hlp(axis, data, ilat, ilon)

            if data_help.count() == 0:
                return_var[ilat, ilon, 0] = data_help[0]
                return_var[ilat, ilon, 1] = data_help[0]
                return_var[ilat, ilon, 2] = data_help[0]
                return_var[ilat, ilon, 3] = data_help[0]
            else:
                data_hits = data_help < threshold
                [events, spei_sum] = runs_of_ones_array_spei(data_hits,
                                                             data_help)

                return_var[ilat, ilon, 0] = np.count_nonzero(events)
                return_var[ilat, ilon, 1] = np.mean(events)
                return_var[ilat, ilon, 2] = np.mean((spei_sum * events) /
                                                    (np.mean(data_help
                                                             [data_hits])
                                                     * np.mean(events)))
                return_var[ilat, ilon, 3] = np.mean(spei_sum / events)

    return return_var


def get_latlon_index(coords, lim1, lim2):
    """Get index for given values between two limits (1D), e.g. lats, lons."""
    index = (np.where(np.absolute(coords - (lim2 + lim1)
                                  / 2.0) <= (lim2 - lim1)
                      / 2.0))[0]
    return index


def plot_map_spei_multi(cfg, data_dict, input_filenames, colormap='jet'):
    """Plot contour maps for multi model mean."""
    print("data_dict['data']")
    print(data_dict['data'])
    # mask = np.isnan(data_dict['data'])
    spei = np.ma.array(data_dict['data'], mask=np.isnan(data_dict['data']))
    # np.ma.masked_less_equal(spei, 0)
    print("spei")
    print(spei)

    # Get latitudes and longitudes from cube
    # lats = data_dict['latitude']
    lons = data_dict['longitude']
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    data_dict.update({'longitude': lons})
    # mesh_ind = np.ix_(range(360), index)
    spei = spei[np.ix_(range(len(data_dict['latitude'])), index)]

    # Plot data
    # Create figure and axes instances
    fig, axx = plt.subplots(figsize=(8, 4))

    axx = plt.axes(projection=cart.PlateCarree(central_longitude=0.0))
    axx.set_extent([-180.0, 180.0, -90.0, 90.0],
                   cart.PlateCarree(central_longitude=0.0))

    # Draw filled contours
    cnplot = plt.contourf(lons, data_dict['latitude'], spei,
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
    axx.set_title(data_dict['datasetname'] + ' ' + data_dict['model_kind'] +
                  ' ' + data_dict['drought_char'])

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
    fig.savefig(get_plot_filename(cfg['indexname'] + '_map' +
                                  data_dict['filename'] + '_' +
                                  data_dict['datasetname'], cfg), dpi=300)
    plt.close()

    _provenance_map_spei_multi(cfg, data_dict, spei, input_filenames)


def plot_map_spei(cfg, cube, levels, name_dict):
    """Plot contour map."""
    print("hello map 1")
    # SPEI array to plot
    # spei = cube.data
    print("hello map 2")

    mask = np.isnan(cube.data)
    spei = np.ma.array(cube.data, mask=mask)
    np.ma.masked_less_equal(spei, 0)

    # Get latitudes and longitudes from cube
    name_dict.update({'latitude': cube.coord('latitude').points})
    lons = cube.coord('longitude').points
    lons = np.where(lons > 180, lons - 360, lons)
    # sort the array
    index = np.argsort(lons)
    lons = lons[index]
    name_dict.update({'longitude': lons})
    # mesh_ind = np.ix_(range(360), index)
    spei = spei[np.ix_(range(len(cube.coord('latitude').points)), index)]

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
    print(dataset_name)

    # Plot data
    # Create figure and axes instances
    fig, axx = plt.subplots(figsize=(8, 4))

    axx = plt.axes(projection=cart.PlateCarree(central_longitude=0.0))
    axx.set_extent([-180.0, 180.0, -90.0, 90.0],
                   cart.PlateCarree(central_longitude=0.0))

    # np.set_printoptions(threshold=np.nan)

    # Draw filled contours
    cnplot = plt.contourf(lons, cube.coord('latitude').points, spei,
                          levels,
                          transform=cart.PlateCarree(central_longitude=0.0),
                          cmap='gnuplot', extend='both', corner_mask=False)
    # Draw coastlines
    axx.coastlines()

    # Add colorbar
    cbar = fig.colorbar(cnplot, ax=axx, shrink=0.6, orientation='horizontal')

    # Add colorbar title string
    cbar.set_label(name_dict['name'])

    # Set labels and title to each plot
    axx.set_xlabel('Longitude')
    axx.set_ylabel('Latitude')
    axx.set_title(dataset_name + ' ' + name_dict['name'])

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

    fig.savefig(get_plot_filename(cfg['indexname'] + '_map' +
                                  name_dict['add_to_filename'] + '_' +
                                  dataset_name, cfg), dpi=300)
    plt.close()

    _provenance_map_spei(cfg, name_dict, spei, dataset_name)


def plot_time_series_spei(cfg, cube, filename, add_to_filename=''):
    """Plot time series."""
    # SPEI vector to plot
    spei = cube.data
    # Get time from cube
    print("cube.coord(time)")
    print(cube.coord('time'))
    time = cube.coord('time').points
    # Adjust (ncdf) time to the format matplotlib expects
    add_m_delta = mda.datestr2num('1850-01-01 00:00:00')
    time = time + add_m_delta

    # Get data set name from cube
    try:
        dataset_name = cube.metadata.attributes['model_id']
    except KeyError:
        try:
            dataset_name = cube.metadata.attributes['source_id']
        except KeyError:
            dataset_name = 'Observations'

    data_dict = {'data': spei,
                 'time': time,
                 'var': cfg['indexname'],
                 'dataset_name': dataset_name,
                 'unit': '1',
                 'filename': filename,
                 'area': add_to_filename}

    fig, axx = plt.subplots(figsize=(16, 4))
    axx.plot_date(time, spei, '-', tz=None, xdate=True, ydate=False,
                  color='r', linewidth=4., linestyle='-', alpha=1.,
                  marker='x')
    axx.axhline(y=-2, color='k')

    # Plot labels and title
    axx.set_xlabel('Time')
    axx.set_ylabel(cfg['indexname'])
    axx.set_title('Mean ' + cfg['indexname'] + ' ' +
                  data_dict['dataset_name'] + ' '
                  + data_dict['area'])

    # Set limits for y-axis
    axx.set_ylim(-4.0, 4.0)

    # Often improves the layout
    fig.tight_layout()
    # Save plot to file
    fig.savefig(get_plot_filename(cfg['indexname'] +
                                  '_time_series_' +
                                  data_dict['area'] +
                                  '_' +
                                  data_dict['dataset_name'], cfg), dpi=300)
    plt.close()

    _provenance_time_series_spei(cfg, data_dict)
