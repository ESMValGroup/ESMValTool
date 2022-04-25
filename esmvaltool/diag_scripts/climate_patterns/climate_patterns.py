import logging
from pathlib import Path

import iris
import iris.coord_categorisation
import iris.cube
import numpy as np
import sklearn.linear_model
import sub_functions as sf
from plotting import plot_cp_timeseries
from rename_variables import (
    rename_anom_variables,
    rename_clim_variables,
    rename_regression_variables,
    rename_variables,
)

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def climatology(cube):
    """Handles aggregation to make climatology."""

    # time extraction

    cube_40yr = cube.extract(
        iris.Constraint(time=lambda t: 1850 <= t.point.year <= 1889,
                        month_number=lambda t: 1 <= t.point <= 12))
    cube_aggregated = make_monthly_climatology(cube_40yr)

    return cube_aggregated


def constrain_latitude(cube):
    """Constraining latitude to meet IMOGEN grid cell specifications."""
    cube_clipped = cube.extract(
        iris.Constraint(latitude=lambda cell: 82.5 >= cell >= -55))
    return cube_clipped


def make_monthly_climatology(cube):
    """Generates a climatology by month_number."""
    if not cube.coords("month_number"):
        iris.coord_categorisation.add_month_number(cube, "time",
                                                   "month_number")
    cube_month_climatol = cube.aggregated_by("month_number",
                                             iris.analysis.MEAN)

    return cube_month_climatol


def diurnal_temp_range(cubelist):
    """Calculates diurnal range from daily max and min temperatures."""
    range_cube = cubelist[0] - cubelist[1]
    range_cube.rename("Diurnal Range")
    range_cube.var_name = ("range_t1p5m")
    range

    return range_cube


def calculate_diurnal_range(clim_list, ts_list):
    """Facilitates diurnal range calculation and appending."""
    temp_range_list_clim = iris.load([])
    temp_range_list_ts = iris.load([])
    comb_list = [clim_list, ts_list]

    for i in range(len(comb_list)):
        for cube in comb_list[i]:
            if (cube.var_name == "tasmax"
                    or cube.var_name == "tasmin") and cube in clim_list:
                temp_range_list_clim.append(cube)
            elif (cube.var_name == "tasmax"
                  or cube.var_name == "tasmin") and cube in ts_list:
                temp_range_list_ts.append(cube)
            else:
                pass

    derived_diurnal_clim = diurnal_temp_range(temp_range_list_clim)
    derived_diurnal_ts = diurnal_temp_range(temp_range_list_ts)

    # append diurnal range to lists
    clim_list_final, ts_list_final = append_diurnal_range(
        derived_diurnal_clim, derived_diurnal_ts, clim_list, ts_list)

    return clim_list_final, ts_list_final


def append_diurnal_range(derived_diurnal_clim, derived_diurnal_ts, clim_list,
                         ts_list):
    """Appends diurnal range to cubelists."""
    # creating cube list without tasmax or tasmin
    # (since we just wanted the diurnal range)
    clim_list_final = iris.load([])
    ts_list_final = iris.load([])

    for cube in clim_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            clim_list_final.append(cube)

    for cube in ts_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            ts_list_final.append(cube)

    clim_list_final.append(derived_diurnal_clim)
    ts_list_final.append(derived_diurnal_ts)

    return clim_list_final, ts_list_final


def calculate_anomaly(clim_list, ts_list):
    """Calculates a climatology."""
    # calculate diurnal temperature range cube
    clim_list_final, ts_list_final = calculate_diurnal_range(
        clim_list, ts_list)

    anom_list_final = ts_list_final.copy()

    # calc the anom by subtracting the monthly climatology from
    # the time series
    for i in range(len(ts_list_final)):
        i_months = anom_list_final[i].coord(
            'month_number').points - 1  # -1 because months are numbered 1..12
        anom_list_final[i].data -= clim_list_final[i][i_months].data

    return clim_list_final, ts_list_final, anom_list_final


def regression(tas, cube_data):
    """Calculates the regression between global surface temp and a variable."""
    slope_array = np.full(tas.data.shape[1:], np.nan)

    # calculate global average
    tas_data = sf.area_avg(tas, return_cube=False)

    for i in range(tas.data.shape[1]):
        for j in range(tas.data.shape[2]):
            if tas.data[0, i, j] is not np.ma.masked:
                model = sklearn.linear_model.LinearRegression(
                    fit_intercept=False, copy_X=True)

                x = tas_data.reshape(-1, 1)
                y = cube_data[:, i, j]

                model.fit(x, y)
                slope_array[i, j] = model.coef_

    return slope_array


def regression_units(tas, cube):
    """Calculates regression coefficient units."""
    print('Cube Units: ', cube.units)
    units = cube.units / tas.units
    print('Regression Units: ', units)

    return units


def calculate_regressions(anom_list):
    """Facilitates the calculation of regression coefficients (climate
    patterns) and the creation of a new cube of patterns per variable."""
    regr_var_list = iris.cube.CubeList([])

    for cube in anom_list:
        if cube.var_name == 't1p5m_anom':

            # getting last 85 years from full timeseries
            tas = cube[-1020:]

    for cube in anom_list:
        cube_ssp = cube[-1020:]
        month_list = iris.cube.CubeList([])

        # exctracting months, regressing, and merging
        for i in range(1, 13):
            month_constraint = iris.Constraint(imogen_drive=i)
            month_cube_ssp = cube_ssp.extract(month_constraint)
            month_tas = tas.extract(month_constraint)

            regr_array = regression(month_tas, month_cube_ssp.data)

            # re-creating cube
            if (cube.var_name == 'swdown_anom'
                    or cube.var_name == 'lwdown_anom'):
                units = 'W m-2 K-1'
            else:
                units = regression_units(tas, cube_ssp)

            # assigning dim_coords
            coord1 = tas.coord(contains_dimension=1)
            coord2 = tas.coord(contains_dimension=2)
            dim_coords_and_dims = [(coord1, 0), (coord2, 1)]

            # assigning aux_coord
            coord_month = iris.coords.AuxCoord(i, var_name='imogen_drive')
            aux_coords_and_dims = [(coord_month, ())]

            cube = rename_regression_variables(cube)

            new_cube = iris.cube.Cube(regr_array,
                                      units=units,
                                      dim_coords_and_dims=dim_coords_and_dims,
                                      aux_coords_and_dims=aux_coords_and_dims,
                                      var_name=cube.var_name,
                                      standard_name=cube.standard_name)
            month_list.append(new_cube)

        conc_cube = month_list.merge_cube()
        regr_var_list.append(conc_cube)

    return regr_var_list


def main(cfg):
    # gets a description of the preprocessed data that we will use as input.
    input_data = cfg["input_data"].values()

    clim_list = iris.load([])
    ts_list = iris.load([])

    for dataset in input_data:
        input_file = dataset["filename"]

        # preparing single cube
        cube_initial = compute_diagnostic(input_file)
        cube_constrained = constrain_latitude(cube_initial)

        # appending to timeseries list
        ts_list.append(cube_constrained)

        # making climatology
        clim_cube = climatology(cube_constrained)
        clim_list.append(clim_cube)

    # calculate anomaly over historical + ssp timeseries
    clim_list_final, ts_list_final, anom_list_final = calculate_anomaly(
        clim_list, ts_list)

    for i in range(len(clim_list_final)):
        clim_list_final[i] = rename_clim_variables(clim_list_final[i])
        ts_list_final[i] = rename_variables(ts_list_final[i])
        anom_list_final[i] = rename_anom_variables(anom_list_final[i])

    regressions = calculate_regressions(anom_list_final.copy())

    # list of variable cube lists
    list_of_cubelists = [
        clim_list_final, ts_list_final, anom_list_final, regressions
    ]
    name_list = [
        "climatology_variables.nc", "ts_variables.nc", "anomaly_variables.nc",
        "regressions.nc"
    ]

    # saving data
    work_path = cfg["work_dir"] + "/"
    for i in range(len(list_of_cubelists)):
        iris.save(list_of_cubelists[i], work_path + name_list[i])

    # saving figures
    plot_path = cfg["plot_dir"] + "/"
    plot_cp_timeseries(list_of_cubelists, plot_path)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
