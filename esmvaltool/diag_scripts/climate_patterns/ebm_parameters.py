import logging
from pathlib import Path

import gm_functions as gm
import iris
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve

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


def ocean_fraction_calc(sftlf):
    sftlf.coord("latitude").coord_system = iris.coord_systems.GeogCS(6371229.0)
    sftlf.coord("longitude").coord_system = iris.coord_systems.GeogCS(
        6371229.0)
    weights = iris.analysis.cartography.area_weights(sftlf)
    lf_area = sftlf.collapsed(
        ['latitude', 'longitude'], iris.analysis.SUM, weights=weights) / 1e12
    sftof = sftlf.copy()
    sftof.data = 100.0 - sftlf.data
    of_area = sftof.collapsed(
        ['latitude', 'longitude'], iris.analysis.SUM, weights=weights) / 1e12

    ocean_frac = of_area.data / (of_area.data + lf_area.data)

    logger.info("Ocean fraction = ", ocean_frac)

    return ocean_frac


def kappa_parameter(f, toa_delta, tas_delta_ocean):
    rms = 10000.0
    kappa = -9999.9
    for i_kappa in range(0, 150):
        kappa_test = 20.0 + 20.0 * float(i_kappa)

        temp_ocean_top = kappa_calc(f, toa_delta, kappa_test)
        rms_local = np.sqrt(
            np.mean((temp_ocean_top - tas_delta_ocean)**2) /
            tas_delta_ocean.shape[0])
        if (rms_local <= rms):
            rms = rms_local
            kappa = kappa_test

    return kappa


def kappa_calc(f, toa_delta, kappa):
    # Set volumetric heat capacity for salt water
    cp = 4.04E6  # (J/K/m3)
    nyr = toa_delta.shape[0]

    n_pde = 20
    dt = (1.0 / float(n_pde)
          ) * 60.0 * 60.0 * 24.0 * 365.0  # Model timestep (seconds)
    temp_ocean_top = np.zeros(nyr)  # Array holds oceanic warming (K)
    # Set number of vertical layers and ocean depth
    n_vert = 254
    depth = 5000.0  # (metres)
    dz = depth / float(n_vert)  # (metres)

    s = (kappa / cp) * (dt / (dz * dz))

    # Now do the pde solving. First set up the arrays (and zero'd)
    t_ocean_old = np.zeros(n_vert + 1)
    t_ocean_new = np.zeros(n_vert + 1)

    C = np.zeros([n_vert + 1, n_vert + 1])
    D = np.zeros([n_vert + 1, n_vert + 1])
    E = np.zeros(n_vert + 1)

    # Conservation check for end of run - add up energy in.
    q_energy = 0.0

    # Now start loop over the different years.
    for j in range(0, nyr):
        factor1 = -toa_delta[j] / (kappa * f)
        factor2 = 0.0

        for k in range(0, n_pde):
            _ = j * n_pde + k

            # First set t_old as t_new ready for next iteration
            t_ocean_old = t_ocean_new

            # Sort out the top points, bottom points and then all for C
            C[0, 0] = s * (1.0 + dz * factor2) + 1
            C[0, 1] = -s
            C[n_vert, n_vert - 1] = -s
            C[n_vert, n_vert] = (1.0 + s)
            for m in range(1, n_vert):
                C[m, m - 1] = -s / 2.0
                C[m, m + 1] = -s / 2.0
                C[m, m] = 1.0 + s

            # Sort out the top points, bottom points and then all for D
            D[0, 0] = -s * (1.0 + dz * factor2) + 1
            D[0, 1] = s
            D[n_vert, n_vert] = (1.0 - s)
            D[n_vert, n_vert - 1] = s
            for m in range(1, n_vert):
                D[m, m - 1] = s / 2.0
                D[m, m + 1] = s / 2.0
                D[m, m] = 1.0 - s

            E[0] = -(kappa / cp) * (dt / dz) * (factor1 + factor1
                                                )  # Assume slow variation in Q

            # Now solve for U_{j+1}.
            # First put bits in to tri-diagonal things for "sparse.spdiags".
            C_sub = np.zeros(n_vert + 1)
            C_main = np.zeros(n_vert + 1)
            C_super = np.zeros(n_vert + 1)
            for m in range(0, n_vert + 1):
                C_main[m] = C[m, m]
            for m in range(0, n_vert):
                C_sub[m] = C[m + 1, m]
                C_super[m] = C[m, m + 1]

            # Now calculate the right-hand side (called b_rhs)
            _ = np.mat(D)
            e_mat = np.mat(E).transpose()
            t_ocean_old_mat = np.mat(t_ocean_old).transpose()
            b_rhs = np.dot(D, t_ocean_old_mat) + e_mat
            b_rhs = np.ravel(b_rhs)

            # Now perform the calculation to update the oceanic temperatures
            t_ocean_new = spsolve(C, b_rhs)

            # Update the cumulative heat entering ocean,
            # and also save to top ocean temperature
            q_energy = q_energy + dt * kappa * (-factor1 -
                                                (t_ocean_new[0] * factor2))

        temp_ocean_top[j] = t_ocean_new[0]

    # Check conservation of energy at end of run.
    # Compare against total heat in the final profile.
    q_energy_derived = 0.0
    for i in range(1, n_vert + 1):
        q_energy_derived = q_energy_derived + cp * 0.5 * (
            t_ocean_new[i] + t_ocean_new[i - 1]) * dz

    conserved = 100.0 * (q_energy_derived / q_energy)
    logger.info("Heat conservation check (%) = ", round(conserved, 2))

    return temp_ocean_top


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

        if cube_initial.var_name == "sftlf":
            sftlf_cube = cube_initial
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

    # calculating global, land and ocean TOA fluxes
    toa_global_cube = net_flux_calculation(toa_global_list)
    toa_land_cube = net_flux_calculation(toa_land_list)
    toa_ocean_cube = net_flux_calculation(toa_ocean_list)

    # calculating ocean fraction
    ocean_frac = ocean_fraction_calc(sftlf_cube)

    # calculating global averages, subtracting mean of first 4 decades
    toa_delta_init = gm.area_avg(toa_global_cube, return_cube=False)
    tas_ocean_init = gm.area_avg(tas_ocean_cube, return_cube=False)
    toa_delta = toa_delta_init - np.mean(toa_delta_init[0:40])
    tas_ocean = tas_ocean_init - np.mean(tas_ocean_init[0:40])

    # calculating EBM parameter, Kappa
    kappa = kappa_parameter(ocean_frac, toa_delta, tas_ocean)
    logger.info("Kappa = ", kappa)

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

    # saving cubes
    work_path = cfg["work_dir"] + "/"
    for i in range(len(list_of_cubes)):
        iris.save(list_of_cubes[i], work_path + name_list[i])

    # saving EBM parameters
    file_frac = open(work_path + 'kappa.dat', 'w')
    data_frac = '{0:10.3f}'.format(kappa)
    file_frac.write(data_frac + '\n')
    file_frac.close()

    # saving figures
    plot_path = cfg["plot_dir"] + "/"
    plot_timeseries(list_of_cubes, plot_path)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
