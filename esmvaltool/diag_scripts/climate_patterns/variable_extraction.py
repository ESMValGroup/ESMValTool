import logging
from pathlib import Path

import gm_functions as gm
import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from rename_variables import (
    rename_anom_variables,
    rename_clim_variables,
    rename_variables,
)
from variable_derivations import diurnal_temp_range

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def plot_timeseries(list_cubelists, plot_path):
    """Plots timeseries of aggregated cubes, across all scenarios."""
    for i in range(len(list_cubelists)):
        cube_list = list_cubelists[i]
        for cube in cube_list:
            avg_cube = gm.area_avg(cube, return_cube=True)
            fig = qplt.plot(avg_cube)
            if i == 0:
                plt.savefig(plot_path + cube.var_name + "_clim")
            if i == 1:
                plt.savefig(plot_path + cube.var_name + "_full_ts")
            if i == 2:
                plt.savefig(plot_path + cube.var_name + "_anom")
            plt.close()

    return fig


def climatology(cube):
    """Handles aggregation to make climatology."""

    # time extraction

    cube_30yr = cube.extract(
        iris.Constraint(time=lambda t: 1851 <= t.point.year <= 1880,
                        month_number=lambda t: 1 <= t.point <= 12))
    cube_aggregated = make_monthly_climatology(cube_30yr)

    return cube_aggregated


def constrain_latitude(cube):
    cube_clipped = cube.extract(
        iris.Constraint(latitude=lambda cell: 82.5 >= cell >= -55))
    return cube_clipped


def make_monthly_climatology(cube):
    if not cube.coords("month_number"):
        iris.coord_categorisation.add_month_number(cube, "time",
                                                   "month_number")
    cube_month_climatol = cube.aggregated_by("month_number",
                                             iris.analysis.MEAN)

    return cube_month_climatol


def calculate_diurnal_range(clim_list, ts_list):
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

    derived_variable_clim = diurnal_temp_range(temp_range_list_clim)
    derived_variable_ts = diurnal_temp_range(temp_range_list_ts)

    return derived_variable_clim, derived_variable_ts


def append_diurnal_range(derived_variable_clim, derived_variable_ts, clim_list,
                         ts_list):
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

    clim_list_final.append(derived_variable_clim)
    ts_list_final.append(derived_variable_ts)

    return clim_list_final, ts_list_final


def calculate_anomaly(clim_list_final, ts_list_final):
    anom_list_final = ts_list_final.copy()

    # calc the anom by subtracting the monthly climatology from
    # the time series
    for i in range(len(ts_list_final)):
        i_months = anom_list_final[i].coord(
            'month_number').points - 1  # -1 because months are numbered 1..12
        anom_list_final[i].data -= clim_list_final[i][i_months].data

    anom_list_renamed = iris.load([])
    for cube in anom_list_final:
        anom_list_renamed.append(rename_anom_variables(cube))

    return anom_list_renamed


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

        # renaming variables, appending to cube list
        ts_renamed = rename_variables(cube_constrained)
        ts_list.append(ts_renamed)

        # making climatology
        clim_cube = climatology(cube_constrained)
        clim_renamed = rename_clim_variables(clim_cube)
        clim_list.append(clim_renamed)

    # calculate diurnal temperature range cube
    (derived_diurnal_clim,
     derived_diurnal_ts) = calculate_diurnal_range(clim_list, ts_list)

    # append diurnal range to lists
    clim_list_final, ts_list_final = append_diurnal_range(
        derived_diurnal_clim, derived_diurnal_ts, clim_list, ts_list)

    # calculate anomaly over historical + ssp timeseries
    anom_list_final = calculate_anomaly(clim_list_final, ts_list_final)

    # list of variable cube lists
    list_of_cubelists = [clim_list_final, ts_list_final, anom_list_final]
    name_list = [
        "climatology_variables.nc", "ts_variables.nc", "anomaly_variables.nc"
    ]

    # saving data
    work_path = cfg["work_dir"] + "/"
    for i in range(len(list_of_cubelists)):
        iris.save(list_of_cubelists[i], work_path + name_list[i])

    # saving figures
    plot_path = cfg["plot_dir"] + "/"
    plot_timeseries(list_of_cubelists, plot_path)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
