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

Updates:
- changed the filename pattern search to read the metadata.yml produced by
  new spei.R diagnostic.
"""
import os
import glob
import iris
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughts.collect_drought_func import (
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
    """Set up cube for time series plot.
    TODO: This should be configurable in recipe. And maybe find nearest point
    instead of crash if the resolution changes.
    """
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
    # Read input data
    input_data = cfg["input_data"].values()
    proj_groups = e.group_metadata(input_data, 'project')

    all_drought = None
    count = len(proj_groups["OBS"])
    print("Number of datasets: ", count)
    # Loop over all OBS datasets and plot event characteristics
    for iii, meta in enumerate(proj_groups["OBS"]):
        print(meta)
        cube = iris.load_cube(meta["filename"])
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        cube_mean = cube.collapsed('time', iris.analysis.MEAN)
        # we could use one cube(list) instead of np to keep lazyness and meta
        if all_drought is None:
            shape = cube_mean.data.shape + (4, count)
            all_drought = np.full(shape, np.nan)
        # ini_time_series_plot(cfg, cube, 'Bremen', meta["filename"])
        # ini_time_series_plot(cfg, cube, 'Nigeria', meta["filename"])

        drought_show = _get_drought_data(cfg, cube)
        # Distinguish between model and observations/reanalysis.
        # Collest all model data in one array.
        all_drought[:, :, :, iii] = drought_show.data
        _plot_single_maps(cfg, cube_mean, drought_show, 'Historic', meta["filename"])

    # Calculating multi model mean and plot it
    # _get_and_plot_obsmodel(cfg, cube, all_drought, all_drought_obs,
    #                        glob.glob(input_filenames))


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
