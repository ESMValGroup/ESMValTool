#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Compares SPI or SPEI data from models with observations/reanalysis.

###############################################################################
droughs/collect_drought.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Collects data produced by spei.R to plot/process them further.
    Applies drought characteristics based on Martin (2018).

Configuration options
---------------------
    indexname: "SPI" or "SPEI"
    reference_dataset: Dataset name to use for comparison (excluded from MMM)
    threshold: Threshold for binary classifiaction of a drought
    compare_intervals: bool, false
        If true,  begin and end of the time periods are compared instead of
        models and reference.
    comparison_period: should be < (end_year - start_year)/2
    start_year: year, start of historical time series
    end_year: year, end of future scenario
"""

import iris
import numpy as np
import datetime as dt
import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.droughts.collect_drought_func import (
    _get_drought_data, _plot_multi_model_maps, _plot_single_maps)


def _plot_models_vs_obs(cfg, cube, mmm, obs, fnames):
    """Compare drought metrics of multi-model mean to observations."""
    latslons = [cube.coord(i).points for i in ["latitude", "longitude"]]
    perc_diff = ((obs-mmm) / (obs + mmm) * 200)
    _plot_multi_model_maps(cfg, mmm, latslons, fnames, 'Historic')
    _plot_multi_model_maps(cfg, obs, latslons, fnames, 'Observations')
    _plot_multi_model_maps(cfg, perc_diff, latslons, fnames, 'Difference')


def _plot_future_vs_past(cfg, cube, slices, fnames):
    """Compare drought metrics of future and historic time slices."""
    latslons = [cube.coord(i).points for i in ["latitude", "longitude"]]
    slices["Difference"] = ((slices["Future"] - slices["Historic"]) /
        (slices["Future"] + slices["Historic"]) * 200)
    for tstype in ['Historic', 'Future', 'Difference']:
        _plot_multi_model_maps(cfg, slices[tstype], latslons, fnames, tstype)


def _set_tscube(cfg, cube, time, tstype):
    """Time slice from a cube with start/end given by cfg."""
    print("sett time slice")
    if tstype == 'Future':
        print("future")
        start_year = cfg["end_year"] - cfg["comparison_period"]
        start = dt.datetime(start_year , 1, 15, 0, 0, 0)
        end = dt.datetime(cfg['end_year'], 12, 16, 0, 0, 0)
    elif tstype == 'Historic':
        print("historic")
        start = dt.datetime(cfg['start_year'], 1, 15, 0, 0, 0)
        end_year = cfg["start_year"] + cfg["comparison_period"]
        end = dt.datetime(end_year, 12, 16, 0, 0, 0)
    print(start, end)
    stime = time.nearest_neighbour_index(time.units.date2num(start))
    etime = time.nearest_neighbour_index(time.units.date2num(end))
    print(stime, etime)
    print(cube)
    tscube = cube[stime:etime, :, :]
    return tscube


def main(cfg):
    """Run the diagnostic."""
    # Read input data
    input_data = cfg["input_data"].values()
    drought_data = []
    drought_slices = {"Historic": [], "Future": []}
    fnames = []  # why do we need them?
    ref_data = None
    for iii, meta in enumerate(input_data):
        fname = meta["filename"]
        cube = iris.load_cube(fname)
        fnames.append(fname)
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        cube_mean = cube.collapsed('time', iris.analysis.MEAN)
        if cfg.get("compare_intervals", False):
            # calculate and plot metrics per time slice
            for tstype in ['Historic', 'Future']:
                ts_cube = _set_tscube(cfg, cube, cube.coord('time'), tstype)
                drought_show = _get_drought_data(cfg, ts_cube)
                drought_slices[tstype].append(drought_show.data)
                _plot_single_maps(cfg, cube_mean, drought_show, tstype, fname)
        else:
            # calculate and plot metrics per dataset
            drought_show = _get_drought_data(cfg, cube)
            if meta["dataset"] == cfg["reference_dataset"]:
                ref_data = drought_show.data
            else:
                drought_data.append(drought_show.data)
            _plot_single_maps(cfg, cube_mean, drought_show, 'Historic', fname)

    if cfg.get("compare_intervals", False):
        # calculate multi model mean for time slices
        slices = {k: np.array(v) for k, v in drought_slices.items()}
        mean_slices = {k: np.nanmean(v, axis=0) for k, v in slices.items()}
        _plot_future_vs_past(cfg, cube, mean_slices, fnames)
    else:
        # calculate multi model mean and compare with reference dataset
        drought_data = np.array(drought_data)
        mmm = np.nanmean(np.array(drought_data), axis=0)
        _plot_models_vs_obs(cfg, cube, mmm, ref_data, fnames)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
