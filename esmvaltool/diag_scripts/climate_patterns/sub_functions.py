"""Script containing relevant sub-functions for main driving scripts.

Author
------
Gregory Munday (Met Office, UK)
"""

import logging
import os
from pathlib import Path

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import numpy as np
from scipy.sparse.linalg import spsolve

logger = logging.getLogger(Path(__file__).stem)


def area_avg(x, return_cube=None):
    """Calculates the global mean of a variable in a cube, area-weighted.

    Parameters
    ----------
    x (cube): input cube
    return_cube (bool): option to return a cube or array

    Returns
    -------
    x2 (cube): cube with collapsed lat-lons, global mean over time
    x2.data (arr): array with collapsed lat-lons, global mean over time
    """
    if not x.coord('latitude').has_bounds():
        x.coord('latitude').guess_bounds()
    if not x.coord('longitude').has_bounds():
        x.coord('longitude').guess_bounds()
    area = iris.analysis.cartography.area_weights(x, normalize=False)
    x2 = x.collapsed(['latitude', 'longitude'],
                     iris.analysis.MEAN,
                     weights=area)

    if return_cube:
        return x2
    else:
        return x2.data


def area_avg_landsea(x, ocean_frac, land_frac, land=True, return_cube=None):
    """Calculates the global mean of a variable in a cube, with options to be
    land or ocean weighted.

    Parameters
    ----------
    x (cube): input cube
    ocean_frac (cube): ocean fraction cube, found from sftlf
    land_frac (cube): land fraction cube, sftlf
    land (bool): option to weight be land or ocean
    return_cube (bool): option to return a cube or array

    Returns
    -------
    x2 (cube): cube with collapsed lat-lons, global mean over time
    x2.data (arr): array with collapsed lat-lons, global mean over time
    """
    if not x.coord('latitude').has_bounds():
        x.coord('latitude').guess_bounds()
    if not x.coord('longitude').has_bounds():
        x.coord('longitude').guess_bounds()

    global_weights = iris.analysis.cartography.area_weights(x, normalize=False)

    if land is False:
        ocean_frac.data = np.ma.masked_less(ocean_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(ocean_frac,
                                                         normalize=False)
        ocean_area = ocean_frac.collapsed(["latitude", "longitude"],
                                          iris.analysis.SUM,
                                          weights=weights) / 1e12
        print("Ocean area: ", ocean_area.data)
        x2 = x.copy()
        x2.data = x2.data * global_weights * ocean_frac.data

        x2 = x2.collapsed(['latitude', 'longitude'],
                          iris.analysis.SUM) / 1e12 / ocean_area.data

    if land:
        land_frac.data = np.ma.masked_less(land_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(land_frac,
                                                         normalize=False)
        land_area = land_frac.collapsed(["latitude", "longitude"],
                                        iris.analysis.SUM,
                                        weights=weights) / 1e12
        print("Land area: ", land_area.data)
        x2 = x.copy()
        x2.data = x2.data * global_weights * land_frac.data
        x2 = x2.collapsed(['latitude', 'longitude'],
                          iris.analysis.SUM) / 1e12 / land_area.data

    if return_cube:
        return x2
    else:
        return x2.data


def kappa_calc_predict(q, f, kappa, lambda_o, lambda_l, nu):
    """Energy balance model test function, which predicts ocean surface
    temperature given the forcing, kappa, and existing atmospheric parameter
    values.

    Parameters
    ----------
    q (arr): derived effective radiative forcing
    f (float): ocean fraction
    kappa (float): ocean diffusivity parameter (W m-1 K-1)
    lambda_o (float): climate sensitivity of land (W m-2 K-1)
    lambda_l (float): climate sensitivity of ocean (W m-2 K-1)
    nu_ratio (float): land-sea temperature contrast in warming

    Returns
    -------
    temp_ocean_top (arr): ocean surface temp. (Huntingford & Cox, 2000)
    """
    cp = 4.04E6
    nyr = q.shape[0]
    n_pde = 20
    dt = (1.0 / float(n_pde)) * 60.0 * 60.0 * 24.0 * 365.0
    temp_ocean_top = np.zeros(nyr)
    n_vert = 254
    depth = 5000.0
    dz = depth / float(n_vert)

    s = (kappa / cp) * (dt / (dz * dz))

    t_ocean_old = np.zeros(n_vert + 1)
    t_ocean_new = np.zeros(n_vert + 1)

    C = np.zeros([n_vert + 1, n_vert + 1])
    D = np.zeros([n_vert + 1, n_vert + 1])
    E = np.zeros(n_vert + 1)

    q_energy = 0.0

    for j in range(0, nyr):
        factor1 = -q[j] / (kappa * f)
        factor2 = (
            (1.0 - f) * lambda_l * nu) / (f * kappa) + (lambda_o / kappa)

        for k in range(0, n_pde):
            _ = j * n_pde + k
            t_ocean_old = t_ocean_new

            C[0, 0] = s * (1.0 + dz * factor2) + 1
            C[0, 1] = -s
            C[n_vert, n_vert - 1] = -s
            C[n_vert, n_vert] = (1.0 + s)
            for m in range(1, n_vert):
                C[m, m - 1] = -s / 2.0
                C[m, m + 1] = -s / 2.0
                C[m, m] = 1.0 + s

            D[0, 0] = -s * (1.0 + dz * factor2) + 1
            D[0, 1] = s
            D[n_vert, n_vert] = (1.0 - s)
            D[n_vert, n_vert - 1] = s
            for m in range(1, n_vert):
                D[m, m - 1] = s / 2.0
                D[m, m + 1] = s / 2.0
                D[m, m] = 1.0 - s

            E[0] = -(kappa / cp) * (dt / dz) * (factor1 + factor1)

            C_sub = np.zeros(n_vert + 1)
            C_main = np.zeros(n_vert + 1)
            C_super = np.zeros(n_vert + 1)
            for m in range(0, n_vert + 1):
                C_main[m] = C[m, m]
            for m in range(0, n_vert):
                C_sub[m] = C[m + 1, m]
                C_super[m] = C[m, m + 1]

            _ = np.mat(D)
            e_mat = np.mat(E).transpose()
            t_ocean_old_mat = np.mat(t_ocean_old).transpose()
            b_rhs = np.dot(D, t_ocean_old_mat) + e_mat
            b_rhs = np.ravel(b_rhs)

            t_ocean_new = spsolve(C, b_rhs)

            q_energy = q_energy + dt * kappa * (-factor1 -
                                                (t_ocean_new[0] * factor2))

        temp_ocean_top[j] = t_ocean_new[0]

    q_energy_derived = 0.0
    for i in range(1, n_vert + 1):
        q_energy_derived = q_energy_derived + cp * 0.5 * (
            t_ocean_new[i] + t_ocean_new[i - 1]) * dz

    conserved = 100.0 * (q_energy_derived / q_energy)
    logger.info("Heat conservation check (%) = ", round(conserved, 2))

    return temp_ocean_top


def make_model_dirs(cube_initial, work_path, plot_path):
    """Creates directories for each input model for saving.

    Parameters
    ----------
    cube_initial (cube): initial input cube used to retrieve model name
    work_path (path): path to work_dir
    plot_path (path): path to plot_dir

    Returns
    -------
    model_work_dir (path): path to specific model directory in work_dir
    model_plot_dir (path): path to specific plot directory in plot_dir
    """
    w_path = os.path.join(work_path, cube_initial.attributes['source_id'])
    p_path = os.path.join(plot_path, cube_initial.attributes['source_id'])
    os.mkdir(w_path)
    os.mkdir(p_path)

    model_work_dir = w_path + "/"
    model_plot_dir = p_path + "/"

    return model_work_dir, model_plot_dir


def parallelise(f, processes=None):
    """Wrapper to parallelise any function, graciously supplied by George Ford
    (Met Office).

    Parameters
    ----------
    f (func): function to be parallelised
    processes (int): number of threads to be used in parallelisation

    Returns
    -------
    result (any): results of parallelised elements
    """
    import multiprocessing as mp
    if processes is None:
        processes = max(1, mp.cpu_count() - 1)
    if processes <= 0:
        processes = 1

    def easy_parallise(f, sequence):
        pool = mp.Pool(processes=processes)
        result = pool.map_async(f, sequence).get()
        pool.close()
        pool.join()
        return result

    from functools import partial
    return partial(easy_parallise, f)
