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

logger = logging.getLogger(Path(__file__).stem)


def compute_diagnostic(filename):
    """Load cube, remove any dimensions of length: 1.

    Parameters
    ----------
    filename (path): path to load cube file

    Returns
    -------
    cube (cube): a cube
    """
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)
    cube = iris.util.squeeze(cube)

    return cube


def area_avg(x, return_cube=None):
    """Calculate the global mean of a variable in a cube, area-weighted.

    Parameters
    ----------
    x (cube): input cube
    return_cube (bool): option to return a cube or array

    Returns
    -------
    x2 (cube): cube with collapsed lat-lons, global mean over time
    x2.data (arr): array with collapsed lat-lons, global mean over time
    """
    if not x.coord("latitude").has_bounds():
        x.coord("latitude").guess_bounds()
    if not x.coord("longitude").has_bounds():
        x.coord("longitude").guess_bounds()
    area = iris.analysis.cartography.area_weights(x, normalize=False)
    x2 = x.collapsed(["latitude", "longitude"],
                     iris.analysis.MEAN,
                     weights=area)

    if return_cube:
        return x2
    else:
        return x2.data


def area_avg_landsea(x, ocean_frac, land_frac, land=True, return_cube=None):
    """Calculate the global mean of a variable in a cube, with options to be
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
    if not x.coord("latitude").has_bounds():
        x.coord("latitude").guess_bounds()
    if not x.coord("longitude").has_bounds():
        x.coord("longitude").guess_bounds()

    global_weights = iris.analysis.cartography.area_weights(x, normalize=False)

    if land is False:
        ocean_frac.data = np.ma.masked_less(ocean_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(ocean_frac,
                                                         normalize=False)
        ocean_area = (ocean_frac.collapsed(
            ["latitude", "longitude"], iris.analysis.SUM, weights=weights) /
                      1e12)
        print("Ocean area: ", ocean_area.data)
        x2 = x.copy()
        x2.data = x2.data * global_weights * ocean_frac.data

        x2 = (x2.collapsed(["latitude", "longitude"], iris.analysis.SUM) /
              1e12 / ocean_area.data)

    if land:
        land_frac.data = np.ma.masked_less(land_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(land_frac,
                                                         normalize=False)
        land_area = (land_frac.collapsed(
            ["latitude", "longitude"], iris.analysis.SUM, weights=weights) /
                     1e12)
        print("Land area: ", land_area.data)
        x2 = x.copy()
        x2.data = x2.data * global_weights * land_frac.data
        x2 = (x2.collapsed(["latitude", "longitude"], iris.analysis.SUM) /
              1e12 / land_area.data)

    if return_cube:
        return x2
    else:
        return x2.data


def make_model_dirs(cube_initial, work_path, plot_path):
    """Create directories for each input model for saving.

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
    w_path = os.path.join(work_path, cube_initial.attributes["source_id"])
    p_path = os.path.join(plot_path, cube_initial.attributes["source_id"])
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
