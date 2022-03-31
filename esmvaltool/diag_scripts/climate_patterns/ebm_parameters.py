import logging
from pathlib import Path

import gm_functions as gm
import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube


def net_flux_calculation(toa_list):
    for cube in toa_list:
        if cube.var_name == "rsdt":
            rsdt = cube
        if cube.var_name == "rlut":
            rlut = cube
        if cube.var_name == "rsut":
            rsut = cube

    toa_net = rsdt - rlut - rsut
    toa_net.rename("toa")
    toa_net.var_name = "toa"

    return toa_net


def plot_timeseries(list_cubes, plot_path):
    """Plots timeseries of aggregated cubes, across all scenarios."""
    for i, cube in enumerate(list_cubes):
        avg_cube = gm.area_avg(cube, return_cube=True)
        fig = qplt.plot(avg_cube)
        if i == 0:
            plt.savefig(plot_path + cube.var_name + "_ts_global")

        if i == 1:
            plt.savefig(plot_path + cube.var_name + "_ts_global")
        if i == 2:
            plt.savefig(plot_path + cube.var_name + "_ts_land")
        if i == 3:
            plt.savefig(plot_path + cube.var_name + "_ts_land")
        if i == 4:
            plt.savefig(plot_path + cube.var_name + "_ts_ocean")
        if i == 5:
            plt.savefig(plot_path + cube.var_name + "_ts_ocean")
        plt.close()

    for i, cube in enumerate(list_cubes):
        fig = qplt.pcolormesh(cube[0])
        if i == 0:
            plt.savefig(plot_path + cube.var_name + "_mesh_global")
        if i == 1:
            plt.savefig(plot_path + cube.var_name + "_mesh_global")
        if i == 2:
            plt.savefig(plot_path + cube.var_name + "_mesh_land")
        if i == 3:
            plt.savefig(plot_path + cube.var_name + "_mesh_land")
        if i == 4:
            plt.savefig(plot_path + cube.var_name + "_mesh_ocean")
        if i == 5:
            plt.savefig(plot_path + cube.var_name + "_mesh_ocean")
        plt.close()

    return fig


def main(cfg):
    # gets a description of the preprocessed data that we will use as input.
    input_data = cfg["input_data"].values()

    toa_global_list = iris.load([])
    toa_land_list = iris.load([])
    toa_ocean_list = iris.load([])

    for dataset in input_data:
        input_file = dataset["filename"]

        # preparing single cube
        cube_initial = compute_diagnostic(input_file)

        if dataset["preprocessor"] == "global_mean_annual":
            if not cube_initial.var_name == "tas":
                toa_global_list.append(cube_initial)
            else:
                tas_global_cube = cube_initial
        elif dataset["preprocessor"] == "land_mean_annual":
            if not cube_initial.var_name == "tas":
                toa_land_list.append(cube_initial)
            else:
                tas_land_cube = cube_initial
        elif dataset["preprocessor"] == "ocean_mean_annual":
            if not cube_initial.var_name == "tas":
                toa_ocean_list.append(cube_initial)
            else:
                tas_ocean_cube = cube_initial

    toa_global_cube = net_flux_calculation(toa_global_list)
    toa_land_cube = net_flux_calculation(toa_land_list)
    toa_ocean_cube = net_flux_calculation(toa_ocean_list)

    # list of variable cube lists
    list_of_cubes = [
        toa_global_cube, tas_global_cube, toa_land_cube, tas_land_cube,
        toa_ocean_cube, tas_ocean_cube
    ]
    name_list = [
        "TOA_global_timeseries.nc",
        "tas_global_timeseries.nc",
        "TOA_land_timeseries.nc",
        "tas_land_timeseries.nc",
        "TOA_ocean_timeseries.nc",
        "tas_ocean_timeseries.nc",
    ]

    # saving data
    work_path = cfg["work_dir"] + "/"
    for i in range(len(list_of_cubes)):
        iris.save(list_of_cubes[i], work_path + name_list[i])

    # saving figures
    plot_path = cfg["plot_dir"] + "/"
    plot_timeseries(list_of_cubes, plot_path)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
