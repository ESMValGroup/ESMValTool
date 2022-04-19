import logging
from pathlib import Path

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import numpy as np
from scipy.sparse.linalg import spsolve

logger = logging.getLogger(Path(__file__).stem)


def stats(x):
    if isinstance(x, int):
        print('single number, you numpty. Max/min/mean = ', x)
    elif isinstance(x, str):
        print('cant stat a string!')
    elif isinstance(x, list) or isinstance(x, np.ndarray):
        print('max = ', np.max(x))
        print('min = ', np.min(x))
        print('mean = ', np.mean(x))
        print('st-dev = ', np.std(x))
    else:
        print('what have you given me??')


def months_appended(x):
    mn2 = x.coord("time")
    m2 = mn2.units.num2date(mn2.points)
    m3 = len(m2)
    time = np.linspace(0, m3, m3)

    return time


def load_cubelist_to_cube(filename=None, load_file=True, cube=None):
    if load_file:
        cube = iris.load(filename)
    if isinstance(cube, iris.cube.CubeList):
        for ijk in np.arange(0, int(len(cube))):
            for key in list(cube[ijk].attributes.keys()):
                del cube[ijk].attributes[key]
            if ijk > 0:
                cube[ijk].coord("time").convert_units(
                    cube[0].coord("time").units)
    cube = cube.concatenate_cube()

    return cube


def area_avg(x, return_cube=None):
    # If the cube does not have bounds, add bounds
    if not x.coord('latitude').has_bounds():
        x.coord('latitude').guess_bounds()
    if not x.coord('longitude').has_bounds():
        x.coord('longitude').guess_bounds()
    # Get the area weights using the same cube
    area = iris.analysis.cartography.area_weights(x, normalize=False)
    # Now collapse the lat and lon to find a global mean over time
    x2 = x.collapsed(['latitude', 'longitude'],
                     iris.analysis.MEAN,
                     weights=area)

    if return_cube:
        return x2
    else:
        return x2.data


def area_avg_landsea(x, ocean_frac, land_frac, land=True, return_cube=None):
    # If the cube does not have bounds, add bounds
    if not x.coord('latitude').has_bounds():
        x.coord('latitude').guess_bounds()
    if not x.coord('longitude').has_bounds():
        x.coord('longitude').guess_bounds()
    # Get the area weights using the same cube
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
        # Now collapse the lat and lon to find a global mean over time
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
