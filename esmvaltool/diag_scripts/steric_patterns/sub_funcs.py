# (C) Crown Copyright 2022-2025, Met Office.
"""Script containing sub-functions for main scripts.

Author
------
Gregory Munday (Met Office, UK)
"""
from __future__ import annotations

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import iris
import iris.analysis.cartography
import iris.coord_categorisation

logger = logging.getLogger(Path(__file__).stem)


def load_cube(filename: str) -> iris.cube.Cube:
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
    return iris.util.squeeze(cube)


def make_model_dirs(cfg: dict, model: str) -> tuple[Path, Path]:
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
    work_path = cfg["work_dir"]
    plot_path = cfg["plot_dir"]
    model_work_dir = Path(work_path) / model
    model_plot_dir = Path(plot_path) / model

    if not Path.exists(model_work_dir):
        Path.mkdir(model_work_dir)
    if not Path.exists(model_plot_dir):
        Path.mkdir(model_plot_dir)
    return model_work_dir, model_plot_dir


def parallelise(function: callable, processes: int | None = None) -> callable:
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

    def easy_parallise(func: callable, sequence: list, cfg: dict) -> list:
        with mp.Pool(processes=processes) as pool:
            config_wrapper = partial(func, cfg=cfg)
            result = pool.map_async(config_wrapper, sequence).get()
            pool.close()
            pool.join()
            return result
    return partial(easy_parallise, function)
