#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Collects SPI or SPEI data comparing models and observations/reanalysis.

Applies drought characteristics based on Martin (2018).

###############################################################################
droughtindex/collect_drought_obs_multi.py
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

###############################################################################

"""
import os
import glob
import iris
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_drought_func import (
    _get_drought_data, _plot_multi_model_maps, _plot_single_maps,
    get_latlon_index, plot_time_series_spei)


def _get_and_plot_obsmodel(cfg, cube, all_drought, all_drought_obs,
                           input_filenames):
    """Calculate multi-model mean and compare it to observations."""
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    all_drought_hist_mean = np.nanmean(all_drought, axis=-1)
    perc_diff = ((all_drought_obs - all_drought_hist_mean)
                 / (all_drought_obs + all_drought_hist_mean) * 200)

    # Plot multi model means
    _plot_multi_model_maps(cfg, all_drought_hist_mean, [lats, lons],
                           input_filenames, 'Historic')
    _plot_multi_model_maps(cfg, all_drought_obs, [lats, lons],
                           input_filenames, 'Observations')
    _plot_multi_model_maps(cfg, perc_diff, [lats, lons],
                           input_filenames, 'Difference')


def ini_time_series_plot(cfg, cube, area, filename):
    """Set up cube for time series plot."""
    coords = ('longitude', 'latitude')
    if area == 'Bremen':
        index_lat = get_latlon_index(cube.coord('latitude').points, 52, 53)
        index_lon = get_latlon_index(cube.coord('longitude').points, 7, 9)
    elif area == 'Nigeria':
        index_lat = get_latlon_index(cube.coord('latitude').points, 7, 9)
        index_lon = get_latlon_index(cube.coord('longitude').points, 8, 10)

    cube_grid_areas = iris.analysis.cartography.area_weights(
        cube[:, index_lat[0]:index_lat[-1] + 1,
             index_lon[0]:index_lon[-1] + 1])
    cube4 = ((cube[:, index_lat[0]:index_lat[-1] + 1,
                   index_lon[0]:index_lon[-1] +
                   1]).collapsed(coords, iris.analysis.MEAN,
                                 weights=cube_grid_areas))

    plot_time_series_spei(cfg, cube4, filename, area)


def main(cfg):
    """Run the diagnostic.

    Parameters :

    ------------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    #######################################################################
    # Read recipe data
    #######################################################################

    # Get filenames of input files produced by diag_spei.r
    # "cfg[n.INPUT_FILES]" is produced by the ESMValTool and contains
    # information on the SPEI input files produced by diag_spei.r
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*_" + \
        (cfg['indexname']).lower() + "_*.nc"
    print(cfg.keys())
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
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        # time = cube.coord('time')

        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time, i.e. data you need
        # to plot the average over the time series of SPEI on a map
        cube2 = cube.collapsed('time', iris.analysis.MEAN)

        # This is only possible because all data must be on the same grid
        if first_run == 1:
            files = os.listdir((cfg[n.INPUT_FILES])[0])
            ncfiles = list(filter(lambda f: f.endswith('.nc'), files))
            shape_all = cube2.data.shape + (4,) + \
                (len(ncfiles) - 1, )
            all_drought = np.full(shape_all, np.nan)
            first_run = 0

        ini_time_series_plot(cfg, cube, 'Bremen', spei_file)
        ini_time_series_plot(cfg, cube, 'Nigeria', spei_file)

        drought_show = _get_drought_data(cfg, cube)

        # Distinguish between model and observations/reanalysis.
        # Collest all model data in one array.
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
        print(dataset_name)
        _plot_single_maps(cfg, cube2, drought_show, 'Historic', spei_file)

    # Calculating multi model mean and plot it
    _get_and_plot_obsmodel(cfg, cube, all_drought, all_drought_obs,
                           glob.glob(input_filenames))


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
