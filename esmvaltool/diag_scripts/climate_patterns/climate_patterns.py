# (C) Crown Copyright 2022-2024, Met Office.
"""Diagnostic script to build climate patterns from CMIP6 models.

Description
-----------
Builds patterns, anomaly and climatology cubes from CMIP6 models.
This diagnostic needs preprocessed mean monthly cubes, with no
gridding requirements. Default re-grid specification exists to
decrease CPU-load and run-time.

Author
------
Gregory Munday (Met Office, UK)


Configuration options in recipe
-------------------------------
jules_mode: bool, optional (default: false)
    options: true, false
    def: outputs extra data (anomaly, climatology) per variable
         to drive JULES-IMOGEN configuration
parallelise: bool, optional (default: false)
    options: true, false
    def: parallelises code to run N models at once
area: str, optional (default: global)
    options: global, land
    def: area over which to calculate climate patterns
"""

import logging
import os
from pathlib import Path

import iris
import iris.coord_categorisation
import iris.cube
import numpy as np
import sklearn.linear_model
from esmvalcore.preprocessor import (
    area_statistics,
    climate_statistics,
    extract_time,
)

import esmvaltool.diag_scripts.climate_patterns.sub_functions as sf
from esmvaltool.diag_scripts.climate_patterns.plotting import (
    plot_patterns,
    plot_timeseries,
)
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def calculate_climatology(cube, syr=1850, eyr=1889):
    """Handle aggregation to make climatology.

    Parameters
    ----------
    cube : cube
        cube loaded from config dictionary
    syr : int
        set climatology start year
    eyr : int
        set climatology end year

    Returns
    -------
    cube_aggregated : cube
        40 year climatology cube from syr-eyr (default 1850-1889)
    """
    cube_40yr = extract_time(
        cube,
        start_year=syr,
        start_month=1,
        start_day=1,
        end_year=eyr,
        end_month=12,
        end_day=31,
    )
    cube_aggregated = climate_statistics(cube_40yr, "mean", "month")

    return cube_aggregated


def diurnal_temp_range(cubelist):
    """Calculate diurnal range from monthly max and min temperatures.

    Parameters
    ----------
    cubelist : cubelist
        cubelist of tasmin and tasmax

    Returns
    -------
    range_cube : cube
        cube of calculated diurnal range
    """
    range_cube = cubelist[0] - cubelist[1]

    # check in case cubes are wrong way around
    if np.mean(range_cube.data) < 0:
        range_cube = -range_cube

    range_cube.rename("Diurnal Range")
    range_cube.var_name = "range_tl1"

    return range_cube


def calculate_diurnal_range(clim_list, ts_list):
    """Facilitate diurnal range calculation and appending.

    Parameters
    ----------
    clim_list : cubelist
        cubelist of climatology cubes
    ts_list : cubelist
        cubelist of standard timeseries cubes

    Returns
    -------
    clim_list_final : cubelist
        cubelist of climatology cubes including diurnal range
    ts_list_final : cubelist
        cubelist of standard timeseries cubes including diurnal range
    """
    temp_range_list_clim = iris.cube.CubeList([])
    temp_range_list_ts = iris.cube.CubeList([])
    comb_list = [clim_list, ts_list]

    for cube_list in comb_list:
        for cube in cube_list:
            if (cube.var_name in ("tasmax", "tasmin")) and cube in clim_list:
                temp_range_list_clim.append(cube)
            elif (cube.var_name in ("tasmax", "tasmin")) and cube in ts_list:
                temp_range_list_ts.append(cube)
            else:
                pass

    derived_diurnal_clim = diurnal_temp_range(temp_range_list_clim)
    derived_diurnal_ts = diurnal_temp_range(temp_range_list_ts)

    # append diurnal range to lists
    clim_list_final, ts_list_final = append_diurnal_range(
        derived_diurnal_clim, derived_diurnal_ts, clim_list, ts_list
    )

    return clim_list_final, ts_list_final


def append_diurnal_range(
    derived_diurnal_clim, derived_diurnal_ts, clim_list, ts_list
):
    """Append diurnal range to cubelists.

    Parameters
    ----------
    derived_diurnal_clim : cube
        derived diurnal climatology cube
    derived_diurnal_ts : cube
        derived diurnal timeseries cube
    clim_list : cubelist
        existing climatology cubelist, no range
    ts_list : cubelist
        existing timeseries cubelist, no range

    Returns
    -------
    clim_list_final : cubelist
        cubelist of climatology cubes including diurnal range
    ts_list_final : cubelist
        cubelist of standard timeseries cubes including diurnal range
    """
    # creating cube list without tasmax or tasmin
    # (since we just wanted the diurnal range)
    clim_list_final = iris.cube.CubeList([])
    ts_list_final = iris.cube.CubeList([])

    for cube in clim_list:
        if cube.var_name not in ("tasmax", "tasmin"):
            clim_list_final.append(cube)

    for cube in ts_list:
        if cube.var_name not in ("tasmax", "tasmin"):
            ts_list_final.append(cube)

    clim_list_final.append(derived_diurnal_clim)
    ts_list_final.append(derived_diurnal_ts)

    return clim_list_final, ts_list_final


def calculate_anomaly(clim_list, ts_list):
    """Calculate variables as anomalies, and adds diurnal range as variable.

    Parameters
    ----------
    clim_list : cubelist
        cubelist of climatology variables
    ts_list : cubelist
        cubelist of standard variable timeseries

    Returns
    -------
    clim_list_final : cubelist
        cubelist of clim. vars, inc. diurnal range
    anom_list_final : cubelist
        cubelist of anomaly vars, inc. diurnal range
    """
    # calculate diurnal temperature range cube
    clim_list_final, ts_list_final = calculate_diurnal_range(
        clim_list, ts_list
    )

    anom_list_final = ts_list_final.copy()

    # calc the anom by subtracting the monthly climatology from
    # the time series
    for i, _ in enumerate(ts_list_final):
        i_months = (
            anom_list_final[i].coord("month_number").points - 1
        )  # -1 because months are numbered 1..12
        anom_list_final[i].data -= clim_list_final[i][i_months].data

    return clim_list_final, anom_list_final


def regression(tas, cube_data, area, ocean_frac=None, land_frac=None):
    """Calculate coeffs of regression between global surf temp and variable.

    Parameters
    ----------
    tas : cube
        near-surface air temperature
    cube_data : arr
        cube.data array of a variable
    area: str
        area over which to calculate patterns
    ocean_frac: cube
        gridded ocean fraction
    land_frac: cube
        gridded land fraction

    Returns
    -------
    slope_array : arr
        array of grid cells with same shape as initial cube,
        containing the regression slope
    """
    if area == "land":
        # calculate average warming over land
        tas_data = sf.area_avg_landsea(
            tas, ocean_frac, land_frac, land=True, return_cube=False
        )
    else:
        # calculate global average warming
        tas_data = area_statistics(tas, "mean").data

    # Reshape cube for regression
    cube_reshaped = cube_data.reshape(cube_data.shape[0], -1)

    # Perform linear regression on valid values
    model = sklearn.linear_model.LinearRegression(
        fit_intercept=False, copy_X=True
    )
    model.fit(tas_data.reshape(-1, 1), cube_reshaped)

    # Extract regression coefficients
    slopes = model.coef_

    # Reshape the regression coefficients back to the shape of the grid cells
    slope_array = slopes.reshape(cube_data.shape[1:])

    return slope_array


def create_cube(tas_cube, ssp_cube, array, month_number, units=None):
    """Create a new cube from existing metadata, and new aray data.

    Parameters
    ----------
    tas_cube: cube
        near-surface air temperature
    ssp_cube: cube
        cube of a given variable
    array: array
        output array from regression
    month_number: int
        month related to the regression array
    units: str
        units related to the regression variable

    Returns
    -------
    cube: cube
        cube filled with regression array and metadata

    """
    # assigning dim_coords
    coord1 = tas_cube.coord(contains_dimension=1)
    coord2 = tas_cube.coord(contains_dimension=2)
    dim_coords_and_dims = [(coord1, 0), (coord2, 1)]

    # assigning aux_coord
    coord_month = iris.coords.AuxCoord(month_number, var_name="imogen_drive")
    aux_coords_and_dims = [(coord_month, ())]

    cube = sf.rename_variables(ssp_cube, has_orig_vars=False)

    # creating cube
    cube = iris.cube.Cube(
        array,
        units=units,
        dim_coords_and_dims=dim_coords_and_dims,
        aux_coords_and_dims=aux_coords_and_dims,
        var_name=cube.var_name,
        standard_name=cube.standard_name,
    )

    return cube


def calculate_regressions(
    anom_list, area, ocean_frac=None, land_frac=None, yrs=86
):
    """Facilitate the calculation of regression coeffs (climate patterns).

    Also creates of a new cube of patterns per variable.

    Parameters
    ----------
    anom_list : cubelist
        cube list of variables as anomalies
    area: str
        area over which to calculate patterns
    ocean_frac: cube
        gridded ocean fraction
    land_frac: cube
        gridded land fraction
    yrs : int
        int to specify length of scenario

    Returns
    -------
    regr_var_list : cubelist
        cube list of newly created regression slope value cubes, for each var
    """
    regr_var_list = iris.cube.CubeList([])

    for cube in anom_list:
        if cube.var_name == "tl1_anom":
            # convert years to months when selecting
            tas = cube[-yrs * 12 :]

    for cube in anom_list:
        cube = cube[-yrs * 12 :]
        month_list = iris.cube.CubeList([])

        # extracting months, regressing, and merging
        for i in range(1, 13):
            month_cube = cube.extract(iris.Constraint(imogen_drive=i))
            month_tas = tas.extract(iris.Constraint(imogen_drive=i))

            if area == "land":
                regr_array = regression(
                    month_tas,
                    month_cube.data,
                    area=area,
                    ocean_frac=ocean_frac,
                    land_frac=land_frac,
                )
            else:
                regr_array = regression(
                    month_tas,
                    month_cube.data,
                    area=area,
                )

            if cube.var_name in ("swdown_anom", "lwdown_anom"):
                units = "W m-2 K-1"
            else:
                units = cube.units / tas.units

            # create, and append cube of regression values
            month_list.append(
                create_cube(tas, cube.copy(), regr_array, i, units=units)
            )

        month_list = month_list.merge_cube()
        regr_var_list.append(month_list)

    return regr_var_list


def cube_saver(list_of_cubelists, work_path, name_list, jules_mode):
    """Save desired cubelists to work_dir, depending on switch settings.

    Parameters
    ----------
    list_of_cubelists : list
        list containing desired cubelists
    work_path : path
        path to work_dir, to save cubelists
    name_list : list
        list of filename strings for saving
    jules_mode : str
        switch option passed through by ESMValTool config dict

    Returns
    -------
    None
    """
    if jules_mode:
        for i in range(0, 3):
            iris.save(
                list_of_cubelists[i], os.path.join(work_path, name_list[i])
            )
    else:
        for i, cube in enumerate(list_of_cubelists[2]):
            list_of_cubelists[2][i] = sf.rename_variables(
                cube, has_orig_vars=False
            )
        iris.save(list_of_cubelists[2], os.path.join(work_path, name_list[2]))


def save_outputs(cfg, list_of_cubelists, model):
    """Save data and plots to relevant directories.

    Parameters
    ----------
    cfg: dict
        Dictionary passed in by ESMValTool preprocessors
    list_of_cubelists: list
        List of cubelists to save
    model : str
        model name

    Returns
    -------
    None
    """
    work_path, plot_path = sf.make_model_dirs(cfg, model)

    name_list = [
        "climatology_variables.nc",
        "anomaly_variables.nc",
        "patterns.nc",
    ]

    # saving data + plotting
    if cfg["jules_mode"] is True:
        plot_timeseries(
            list_of_cubelists[0],
            plot_path,
            "40 Year Climatologies, 1850-1889",
            "Climatologies",
        )
        plot_timeseries(
            list_of_cubelists[1],
            plot_path,
            "Anomaly Timeseries, 1850-2100",
            "Anomalies",
        )
        plot_patterns(list_of_cubelists[2], plot_path)
        cube_saver(
            list_of_cubelists,
            work_path,
            name_list,
            jules_mode=cfg["jules_mode"],
        )

    else:
        plot_patterns(list_of_cubelists[2], plot_path)
        cube_saver(
            list_of_cubelists,
            work_path,
            name_list,
            jules_mode=cfg["jules_mode"],
        )


def get_provenance_record():
    """Create a provenance record describing the diagnostic data and plot.

    Parameters
    ----------
    None

    Returns
    -------
    record : dict
        provenance record
    """
    record = {
        "caption": ["Generating Climate Patterns from CMIP6 Models"],
        "statistics": ["mean", "other"],
        "domains": ["global"],
        "themes": ["carbon"],
        "realms": ["atmos"],
        "authors": ["munday_gregory"],
    }

    return record


def extract_data_from_cfg(cfg, model):
    """Extract model data from the cfg.

    Parameters
    ----------
    cfg: dict
        Dictionary passed in by ESMValTool preprocessors
    model : str
        model name

    Returns
    -------
    clim_list: cubelist
        cubelist of climatologies
    ts_list: cubelist
        cubelist of spatial timeseries
    sftlf: cube
        land fraction cube
    """
    clim_list = iris.cube.CubeList([])
    ts_list = iris.cube.CubeList([])

    for dataset in cfg["input_data"].values():
        if dataset["dataset"] == model:
            input_file = dataset["filename"]

            # preparing single cube
            cube = sf.load_cube(input_file)

            if dataset["exp"] != "historical-ssp585":
                sftlf = cube
            else:
                # appending to timeseries list
                ts_list.append(cube)

                # making climatology
                clim_cube = calculate_climatology(cube)
                clim_list.append(clim_cube)

    if cfg["area"] == "land":
        return clim_list, ts_list, sftlf

    return clim_list, ts_list, None


def patterns(model, cfg):
    """Driving function for script, taking in model data and saving parameters.

    Parameters
    ----------
    model : str
        model name
    cfg: dict
        Dictionary passed in by ESMValTool preprocessors

    Returns
    -------
    None
    """
    clim_list, ts_list, sftlf = extract_data_from_cfg(cfg, model)

    if cfg["area"] == "land":
        # calculate land/ocean_fracs
        ocean_frac, land_frac = sf.ocean_fraction_calc(sftlf)

    # calculate anomaly over historical + ssp timeseries
    clim_list_final, anom_list_final = calculate_anomaly(clim_list, ts_list)

    for i, cube in enumerate(clim_list_final):
        clim_list_final[i] = sf.rename_variables(
            cube, has_orig_vars=True, new_extension="_clim"
        )
        anom_list_final[i] = sf.rename_variables(
            anom_list_final[i], has_orig_vars=True, new_extension="_anom"
        )

    if cfg["area"] == "land":
        regressions = calculate_regressions(
            anom_list_final,
            cfg["area"],
            ocean_frac=ocean_frac,
            land_frac=land_frac,
        )
    else:
        regressions = calculate_regressions(anom_list_final, cfg["area"])

    list_of_cubelists = [clim_list_final, anom_list_final, regressions]

    save_outputs(cfg, list_of_cubelists, model)

    # Provenance Logging, removed due to sporadic errors. Fix later.

    # model_work_dir, _ = sf.make_model_dirs(
    #     cfg,
    #     model
    # )

    # provenance_record = get_provenance_record()
    # path = os.path.join(model_work_dir, "patterns.nc")
    # with ProvenanceLogger(cfg) as provenance_logger:
    #     provenance_logger.log(path, provenance_record)


def main(cfg):
    """Take in driving data with parallelisation options.

    Parameters
    ----------
    cfg : dict
        the global config dictionary, passed by ESMValTool.

    Returns
    -------
    None
    """
    input_data = cfg["input_data"].values()
    parallelise = cfg["parallelise"]

    models = []
    for mod in input_data:
        model = mod["dataset"]
        if model not in models:
            models.append(model)

    if parallelise is True:
        sf.parallelise(patterns)(models, cfg)
    else:
        for model in models:
            patterns(model, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
