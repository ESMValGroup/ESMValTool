# (C) Crown Copyright 2022-2025, Met Office.
"""Diagnostic script to calculate patterns of sterodynamic sea-level change from CMIP6 models.

Description
-----------
Calculates regressions between the global thermal expansion (zostoga) and the
dynamic sea-level change (zos) for the CMIP6 models. This gives you patterns
of steric sea-levelchange for the different scenarios.

Author
------
Gregory Munday (Met Office, UK)
"""

import logging
from pathlib import Path

import iris
import iris.coord_categorisation
import iris.cube
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.steric_patterns import sub_funcs as sf

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record() -> dict:
    """Create a provenance record describing the diagnostic data and plot.

    Parameters
    ----------
    None

    Returns
    -------
    record : dict
        provenance record
    """
    record = {
        "caption": ["Calculating sterodynamic sea-level patterns"],
        "statistics": ["mean", "other"],
        "domains": ["global"],
        "themes": ["ocean"],
        "realms": ["ocean"],
        "authors": ["munday_gregory"],
    }
    return record


def calculate_drift(cube: iris.cube.Cube) -> float:
    """Calculate zostoga model drift over the PiControl experiment.

    Parameters
    ----------
    cube : iris.cube.Cube
        cube of global thermal expansion

    Returns
    -------
    drift : float
        model drift (slope of regression)
    """
    model = LinearRegression(fit_intercept=True)
    model.fit(cube.coord("time").points.reshape(-1, 1), cube.data)
    return model.coef_


def detrend_zostoga(
    zostoga: iris.cube.Cube, drift_coef: float,
    plot_path: Path) -> iris.cube.Cube:
    """Detrend the zostoga cube.

    Parameters
    ----------
    zostoga: iris.cube.Cube
        raw zostoga cube
    drift_coef: float
        drift coefficient to detrend zostoga
    plot_path: Path
        path to plot folder

    Returns
    -------
    zostoga_dedrended: iris.cube.Cube
        detrended zostoga cube
    """
    trend = drift_coef * zostoga.coord("time").points
    trend -= trend[0]
    zostoga_detrended = zostoga - trend

    # Check how the detrending looks
    fig = plt.figure(figsize=(10, 6))

    ax = fig.add_subplot(111)
    ax.scatter(
        np.arange(len(zostoga.data)), zostoga.data, s=2,
        alpha=0.8, color="navy", label="Original")
    ax.scatter(
        np.arange(len(zostoga.data)), zostoga_detrended.data, s=2,
        alpha=0.8, color="darkorchid", label="Detrended")
    ax.set_xlabel("Time from PiControl start (yrs)")
    ax.set_ylabel("Global thermal expansion (m)")
    ax.legend(loc="upper left", frameon=False)

    fig.savefig(Path(plot_path) / f'detrended_{zostoga.attributes["source_id"]}.png')
    return zostoga_detrended


def dyn_steric_regression(
    zostoga: iris.cube.Cube, zos: iris.cube.Cube) -> tuple[np.array]:
    """Calculate the zostoga/zos regression.

    Parameters
    ----------
    zostoga: iris.cube.Cube
        zostoga cube
    zos: iris.cube.Cube
        zos cube

    Returns
    -------
    slopes: np.array
        regression slopes over the grid
    mask: np.array
        model-specific mask
    """
    time = zos.coord("time")
    dates = time.units.num2date(time.points)
    yrs = np.array([date.year for date in dates])
    start_yr = 2005
    end_yr = 2100

    # Calculate regression coefficients for period 2005-2100
    idx = ((yrs >= start_yr) & (yrs <= end_yr))

    zostoga.data = zostoga.data - zostoga.data[0:10].mean()
    zos.data = zos.data - zos.data[0:10].mean()

    # Calculate the slope and intercepts of linear fits of the
    # global and local sea level projections
    regr = LinearRegression()
    regr.fit(
        zostoga.data[idx].reshape(-1, 1),
        zos.data[idx].reshape(zos.data[idx].shape[0], -1))
    slopes = regr.coef_.reshape(180, 360)

    # Deal with dodgy land mask, which won't really matter in the end anyway
    if zos.attributes["source_id"] == "GISS-E2-1-H":
        mask = np.full(zos.data.shape, np.nan)
    else:
        mask = zos.data.mask[0]

    return slopes, mask


def save_data(
        slopes: np.array, mask: np.array,
        work_path: Path, model: str, scenario: str) -> None:
    """Save the zostoga/zos regression slopes and model mask.

    Parameters
    ----------
    slopes: np.array
        zostoga/zos regression slopes
    mask: np.array
        model mask
    work_path: Path
        path to save data to

    Returns
    -------
    None
    """
    np.save(Path(work_path) / f"zos_regression_{scenario}_{model}.npy", slopes)
    np.save(Path(work_path) / f"zos_mask_{scenario}_{model}.npy", mask)


def evaluate_patterns(
        zostoga: list, zos: list, slopes: list,
        plot_path: Path, model: str) -> None:
    """Evaluate the patterns.

    Parameters
    ----------
    zostoga: list
        list of zostoga scenarios
    zos: list
        list of zos scenarios
    slopes: list
        list of slopes
    plot_path: Path
        path to save plots to
    model: str
        model name

    Returns
    -------
    None
    """


def extract_data_from_cfg(model: str, cfg: dict) -> tuple[list]:
    """Extract model data from the cfg.

    Parameters
    ----------
    model : str
        model name
    cfg: dict
        dictionary passed in by ESMValTool preprocessors

    Returns
    -------
    zostoga: list
        list of zostoga scenarios
    zos: list
        list of zos scenarios
    """
    for dataset in cfg["input_data"].values():
        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zostoga_piControl")):
            input_file = dataset["filename"]
            zostoga_picontrol = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zostoga_245")):
            input_file = dataset["filename"]
            zostoga_245 = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zostoga_370")):
            input_file = dataset["filename"]
            zostoga_370 = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zostoga_585")):
            input_file = dataset["filename"]
            zostoga_585 = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zos_245")):
            input_file = dataset["filename"]
            zos_245 = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zos_370")):
            input_file = dataset["filename"]
            zos_370 = sf.load_cube(input_file)

        if ((dataset["dataset"] == model) and
            (dataset["variable_group"] == "zos_585")):
            input_file = dataset["filename"]
            zos_585 = sf.load_cube(input_file)


    zostoga = [zostoga_picontrol, zostoga_245, zostoga_370, zostoga_585]
    zos = [zos_245, zos_370, zos_585]
    return zostoga, zos


def patterns(model: str, cfg: dict) -> None:
    """Run pattern calculations.

    Parameters
    ----------
    model : str
        mame of the model to calculate patterns for
    cfg : dict
        global config dictionary

    Returns
    -------
    None
    """
    work_path, plot_path = sf.make_model_dirs(cfg, model)
    zostoga_list, zos_list = extract_data_from_cfg(model, cfg)

    # Calculate drift from PiControl zostoga
    zostoga_drift = calculate_drift(zostoga_list[0])

    # Detrend the scenario zostogas
    zostoga_detrended = [
        detrend_zostoga(z, zostoga_drift, plot_path) for z in zostoga_list[1:]]

    # Calculate regression between zostoga and zos
    slopes_masks = [
        dyn_steric_regression(z_dtr, zos) for (z_dtr, zos)
        in zip(zostoga_detrended, zos_list)]
    slopes, masks = zip(*slopes_masks)

    scenarios = ["ssp245", "ssp370", "ssp585"]
    for i, (s, m) in enumerate(zip(slopes, masks)):
        save_data(s, m, work_path, model, scenarios[i])

    # Test the patterns
    evaluate_patterns(slopes, masks, plot_path, model, scenarios)


def main(cfg: dict) -> None:
    """Take in driving data with parallelisation options.

    Parameters
    ----------
    cfg : dict
        the global config dictionary, passed by ESMValTool.

    Returns
    -------
    None
    """
    input_data = cfg["input_data"].values()

    models = []
    for mod in input_data:
        model = mod["dataset"]
        if model not in models:
            models.append(model)

    sf.parallelise(patterns)(models, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
