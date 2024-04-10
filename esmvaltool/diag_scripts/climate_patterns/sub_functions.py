# (C) Crown Copyright 2022-2024, Met Office.
"""Script containing relevant sub-functions for driving scripts.

Author
------
Gregory Munday (Met Office, UK)
"""

import logging
import multiprocessing as mp
import os
from functools import partial
from pathlib import Path

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import numpy as np

logger = logging.getLogger(Path(__file__).stem)


def load_cube(filename):
    """Load cube, remove any dimensions of length: 1.

    Parameters
    ----------
    filename : path
        path to load cube file

    Returns
    -------
    cube : cube
        a cube
    """
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)
    cube = iris.util.squeeze(cube)

    return cube


def ocean_fraction_calc(sftlf):
    """Calculate gridded land and ocean fractions.

    Parameters
    ----------
    sftlf: cube
        land-fraction cube from piControl experiment

    Returns
    -------
    ocean_frac: cube
        ocean_fraction cube for area-weights
    land_frac: cube
        land_fraction cube for area-weights
    of: float
        ocean fraction float for EBM calculations
    """
    sftlf.coord("latitude").coord_system = iris.coord_systems.GeogCS(
        6371229.0
    )
    sftlf.coord("longitude").coord_system = iris.coord_systems.GeogCS(
        6371229.0
    )
    sftof = 100 - sftlf

    ocean_frac = sftof / 100
    land_frac = sftlf / 100

    return ocean_frac, land_frac


def area_avg_landsea(cube, ocean_frac, land_frac, land=True, return_cube=None):
    """Calculate the global mean of a variable in a cube.

    Parameters
    ----------
    cube : cube
        input cube
    ocean_frac : cube
        ocean fraction cube, found from sftlf
    land_frac : cube
        land fraction cube, sftlf
    land : bool
        option to weight be land or ocean
    return_cube : bool
        option to return a cube or array

    Returns
    -------
    cube2 : cube
        cube with collapsed lat-lons, global mean over time
    cube2.data : arr
        array with collapsed lat-lons, global mean over time
    """
    if not cube.coord("latitude").has_bounds():
        cube.coord("latitude").guess_bounds()
    if not cube.coord("longitude").has_bounds():
        cube.coord("longitude").guess_bounds()

    global_weights = iris.analysis.cartography.area_weights(
        cube,
        normalize=False
    )

    if land is False:
        ocean_frac.data = np.ma.masked_less(ocean_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(
            ocean_frac,
            normalize=False
        )
        ocean_area = (
            ocean_frac.collapsed(
                ["latitude", "longitude"], iris.analysis.SUM, weights=weights
            )
            / 1e12
        )
        cube2 = cube * global_weights * ocean_frac

        cube2 = (
            cube2.collapsed(["latitude", "longitude"], iris.analysis.SUM)
            / 1e12
            / ocean_area
        )

    if land:
        land_frac.data = np.ma.masked_less(land_frac.data, 0.01)
        weights = iris.analysis.cartography.area_weights(
            land_frac,
            normalize=False
        )
        land_area = (
            land_frac.collapsed(
                ["latitude", "longitude"], iris.analysis.SUM, weights=weights
            )
            / 1e12
        )
        cube2 = cube * global_weights * land_frac
        cube2 = (
            cube2.collapsed(["latitude", "longitude"], iris.analysis.SUM)
            / 1e12
            / land_area
        )

    if return_cube:
        return cube2

    return cube2.data


def make_model_dirs(cfg, model):
    """Create directories for each input model for saving.

    Parameters
    ----------
    cfg: dict
        Dictionary passed in by ESMValTool preprocessors
    model : str
        model name

    Returns
    -------
    model_work_dir : path
        path to specific model directory in work_dir
    model_plot_dir : path
        path to specific plot directory in plot_dir
    """
    work_path = cfg["work_dir"] + "/"
    plot_path = cfg["plot_dir"] + "/"
    w_path = os.path.join(work_path, model)
    p_path = os.path.join(plot_path, model)

    if not os.path.exists(w_path):
        os.mkdir(w_path)
    if not os.path.exists(p_path):
        os.mkdir(p_path)

    model_work_dir = w_path + "/"
    model_plot_dir = p_path + "/"

    return model_work_dir, model_plot_dir


def parallelise(function, processes=None):
    """Parallelise any function, by George Ford, Met Office.

    Parameters
    ----------
    function : function
        function to be parallelised
    processes : int
        number of threads to be used in parallelisation

    Returns
    -------
    result : any
        results of parallelised elements
    """
    if processes is None:
        processes = max(1, mp.cpu_count() - 1)
    if processes <= 0:
        processes = 1

    def easy_parallise(func, sequence, cfg):
        with mp.Pool(processes=processes) as pool:
            config_wrapper = partial(func, cfg=cfg)
            result = pool.map_async(config_wrapper, sequence).get()
            pool.close()
            pool.join()
            return result

    return partial(easy_parallise, function)
