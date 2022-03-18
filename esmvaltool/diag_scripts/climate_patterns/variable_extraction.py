import logging
from pathlib import Path

import gm_functions as gm
import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
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
                plt.savefig(plot_path + cube.var_name + "_hist")
            if i == 1:
                plt.savefig(plot_path + cube.var_name + "_scen")
            if i == 2:
                plt.savefig(plot_path + cube.var_name + "_anom")
            plt.close()

    return fig


def prepare_cube(cube, dataset):
    """Handles time extraction, aggregation and latitude extraction."""

    # time extraction
    if dataset["exp"] == "historical":
        cube_30yr = cube.extract(
            iris.Constraint(year=lambda cell: 2010 >= cell > 1980))
    elif dataset["exp"] == "ssp370":
        cube_30yr = cube.extract(
            iris.Constraint(year=lambda cell: 2080 > cell >= 2050))

    # aggregation, latitude clipping
    cube_aggregated = make_monthly_climatology(cube_30yr)

    # anomaly generation
    anomaly = make_variable_anomaly(cube_30yr, cube_aggregated)

    cube_clipped = cube_aggregated.extract(
        iris.Constraint(latitude=lambda cell: 82.5 >= cell >= -55))

    anomaly_clipped = anomaly.extract(
        iris.Constraint(latitude=lambda cell: 82.5 >= cell >= -55))

    return cube_clipped, anomaly_clipped


def make_monthly_climatology(cube):
    if not cube.coords("month"):
        iris.coord_categorisation.add_month(cube, "time", "month")
    cube_month_climatol = cube.aggregated_by(["month"], iris.analysis.MEAN)

    return cube_month_climatol


def make_variable_anomaly(cube_normal, cube_aggregated):
    if not cube_normal.coords("month"):
        iris.coord_categorisation.add_month(cube_normal, "time", "month")
    anomaly = cube_normal.aggregated_by(["month"],
                                        iris.analysis.MAX) - cube_aggregated

    return anomaly


def rename_variables(cube):
    # rename variables to fit in JULES framework
    if cube.var_name == "tas":
        cube.rename("t1p5m_clim")
        cube.var_name = "t1p5m_clim"
    if cube.var_name == "hurs":
        cube.rename("rh1p5m_clim")
        cube.var_name = "rh1p5m_clim"
    if cube.var_name == "huss":
        cube.rename("q1p5m_clim")
        cube.var_name = "q1p5m_clim"
    if cube.var_name == "pr":
        cube.rename("precip_clim")
        cube.var_name = "precip_clim"
    if cube.var_name == "sfcWind":
        cube.rename("wind_clim")
        cube.var_name = "wind_clim"
    if cube.var_name == "ps":
        cube.rename("pstar_clim")
        cube.var_name = "pstar_clim"
    if cube.var_name == "rsds":
        cube.rename("swdown_clim")
        cube.var_name = "swdown_clim"
    if cube.var_name == "rlds":
        cube.rename("lwdown_clim")
        cube.var_name = "lwdown_clim"

    return cube


def main(cfg):
    # Gets a description of the preprocessed data that we will use as input.
    input_data = cfg["input_data"].values()
    hist_list = iris.load([])
    scenario_list = iris.load([])
    anomaly_list = iris.load([])

    for dataset in input_data:
        input_file = dataset["filename"]

        # manipulating single cube
        cube_initial = compute_diagnostic(input_file)
        cube_prepped, anomaly_prepped = prepare_cube(cube_initial, dataset)

        # renaming cube function
        cube_renamed = rename_variables(cube_prepped)
        anomaly_renamed = rename_variables(anomaly_prepped)

        # appending
        if dataset["exp"] == "historical":
            hist_list.append(cube_renamed)
            anomaly_list.append(anomaly_renamed)
        elif dataset["exp"] == "ssp370":
            scenario_list.append(cube_renamed)
        else:
            pass

    # calculate diurnal temperature range cube
    temp_range_list_hist = iris.load([])
    temp_range_list_scen = iris.load([])
    temp_range_list_anom = iris.load([])
    comb_list = [hist_list, scenario_list, anomaly_list]

    for i in range(len(comb_list)):
        for cube in comb_list[i]:
            if (cube.var_name == "tasmax"
                    or cube.var_name == "tasmin") and cube in hist_list:
                temp_range_list_hist.append(cube)
            elif (cube.var_name == "tasmax"
                  or cube.var_name == "tasmin" and cube in scenario_list):
                temp_range_list_scen.append(cube)
            elif (cube.var_name == "tasmax"
                  or cube.var_name == "tasmin" and cube in anomaly_list):
                temp_range_list_anom.append(cube)
            else:
                pass

    derived_variable_hist = diurnal_temp_range(temp_range_list_hist)
    derived_variable_scenario = diurnal_temp_range(temp_range_list_scen)
    derived_variable_anom = diurnal_temp_range(temp_range_list_anom)

    # creating cube list without tasmax or tasmin
    # (since we just wanted the diurnal range)
    hist_list_final = iris.load([])
    scen_list_final = iris.load([])
    anom_list_final = iris.load([])

    for cube in hist_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            hist_list_final.append(cube)

    for cube in scenario_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            scen_list_final.append(cube)

    for cube in anomaly_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            anom_list_final.append(cube)

    hist_list_final.append(derived_variable_hist)
    scen_list_final.append(derived_variable_scenario)
    anom_list_final.append(derived_variable_anom)

    # list of variable cube lists
    list_of_cubelists = [hist_list_final, scen_list_final, anom_list_final]
    name_list = [
        "_historical_variables.nc", "_scenario_variables.nc",
        "_anomaly_variables.nc"
    ]

    # saving data
    work_path = cfg["work_dir"]
    for i in range(len(list_of_cubelists)):
        iris.save(list_of_cubelists[i], work_path + name_list[i])

    # saving figures
    plot_path = cfg["plot_dir"] + "/"
    plot_timeseries(list_of_cubelists, plot_path)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
