# (C) Crown Copyright 2025, Met Office.
"""Diagnostic script to calculate patterns of sterodynamic sea-level change from CMIP6 models.

Description
-----------
Calculates regressions between the global thermal expansion (zostoga) and the
dynamic sea-level change (zos) for the CMIP6 models. This gives you patterns
of steric sea-level change for the different scenarios.

Author
------
Gregory Munday (Met Office, UK)
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import iris.coord_categorisation
import iris.cube
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    run_diagnostic,
    save_figure,
)
from esmvaltool.diag_scripts.steric_patterns import sub_funcs as sf
import iris.quickplot as qplt

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
    return {
        "caption": ["Calculating sterodynamic sea-level patterns"],
        "statistics": ["mean", "other"],
        "domains": ["global"],
        "themes": ["phys"],
        "realms": ["ocean"],
        "authors": ["munday_gregory"],
        "plot_types": ["scatter", "map"],
        "references": ["palmer2020", "perks2023"]}


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

    fig.savefig(
        Path(plot_path)
        / f'detrended_{zostoga.attributes["source_id"]}.png', dpi=150)
    return zostoga_detrended


def dyn_steric_regression(
        zostoga: iris.cube.Cube, zos: iris.cube.Cube,
        plot_path: Path, scenario: str) -> tuple[np.array]:
    """Calculate the zostoga/zos regression.

    Parameters
    ----------
    zostoga: iris.cube.Cube
        zostoga cube
    zos: iris.cube.Cube
        zos cube
    plot_path: Path
        path to save plots to
    scenario: str
        scenario name

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
    index = (yrs >= start_yr) & (yrs <= end_yr)
    zostoga = zostoga[index]
    zos = zos[index]

    zostoga.data = zostoga.data - np.mean(zostoga.data[0:10])
    zos.data = zos.data - np.mean(zos.data[0:10], axis=0)

    # Calculate the slope and intercepts of linear fits of the
    # global and local sea level projections
    regr = LinearRegression()
    regr.fit(
        zostoga.data.reshape(-1, 1),
        zostoga.data.reshape(-1, 1) + zos.data.reshape(zos.data.shape[0], -1))
    slopes = regr.coef_.reshape(180, 360)

    # Deal with dodgy land mask, which won't really matter in the end anyway
    if zos.attributes["source_id"] == "GISS-E2-1-H":
        mask = np.full(zos.data.shape, np.nan)
    else:
        mask = zos.data.mask[0]

    fig = evaluate_regression(zostoga.data, zos.data, slopes, zos)
    fig.savefig(
        Path(plot_path)
        / f"regression_{zostoga.attributes['source_id']}_{scenario}.png",
        dpi=150)

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
    model: str
        model name
    scenario: str
        scenario name

    Returns
    -------
    None
    """
    np.save(
        Path(work_path)
        / f"zos_regression_{scenario}_{model}.npy", slopes)
    np.save(Path(work_path) / f"zos_mask_{scenario}_{model}.npy", mask)


def evaluate_regression(
        zostoga_data: np.array, zos_data: np.array,
        slopes: np.array, zos_cube: iris.cube.Cube) -> plt.figure:
    """Evaluate the regression.

    Parameters
    ----------
    zostoga_data: np.array
        global thermal expansion data
    zos_data: np.array
        dynamic sea-level change data
    slopes: np.array
        regression slopes

    Returns
    -------
    fig: plt.figure
        figure of the regression
    """
    x_vals = np.linspace(
        zostoga_data.min(),
        zostoga_data.max(),
        len(zostoga_data))

    # Plot the regression at a few grid points
    fig = plt.figure(figsize=(12, 6), layout="constrained")

    lat_idxs = [150, 26, 170]
    lon_idxs = [110, 30, 340]
    lats = zos_cube.coord("latitude").points[lat_idxs]
    lons = zos_cube.coord("longitude").points[lon_idxs]
    for index in range(3):
        ax = fig.add_subplot(1, 3, index+1)
        ax.scatter(
            zostoga_data,
            zostoga_data + zos_data[:, lat_idxs[index], lon_idxs[index]],
            s=2,
            alpha=0.8, color="navy", label="Model")
        ax.plot(
            x_vals, slopes[lat_idxs[index], lon_idxs[index]] * x_vals,
            color="darkorchid", label="Regression")
        ax.text(
            0.62,
            0.95,
            f"({float(lats[index]):.1f}$\degree$, "
            f"{float(lons[index]):.1f}$\degree$)",
            transform=ax.transAxes)
        ax.set_ylim([-0.2, 0.8])
        ax.set_xlabel("Global thermal expansion (m)")
        ax.set_ylabel("Dynamic sea level (m)")
        ax.legend(loc="upper left", frameon=False)

    return fig


def plot_evals(
        diff_list: list, zos: list, mse_list: list) -> plt.figure:
    """Plot the evaluation of the thermosteric patterns.

    Parameters
    ----------
    diff_list: list
        list of difference maps
    zos: list
        list of zos scenarios
    mse_list: list
        list of mean squared errors

    Returns
    -------
    fig: plt.figure
        figure of the evaluation
    """
    # Plot the mse for each scenario
    fig = plt.figure(figsize=(12, 5), layout="constrained")
    vmin = -0.5
    vmax = 0.5
    titles = ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
    for index in range(3):
        ax = fig.add_subplot(2, 3, index+1, projection=ccrs.PlateCarree())
        ax.pcolormesh(
            zos[index].coord("longitude").points,
            zos[index].coord("latitude").points,
            diff_list[index], transform=ccrs.PlateCarree(),
            vmin=vmin, vmax=vmax,
            cmap="RdBu_r")
        ax.set_title(titles[index])
        cbar = plt.colorbar(
            ax.collections[0], ax=ax, orientation="horizontal")
        cbar.set_label("Prediction - ESM (m)")

    ax = fig.add_subplot(2, 3, (4, 6))
    time = np.linspace(2015, 2100, 1032)
    ax.plot(time, mse_list[0][:(86 * 12)], label="SSP245", color="navy")
    ax.plot(time, mse_list[1][:(86 * 12)], label="SSP370", color="orange")
    ax.plot(time, mse_list[2][:(86 * 12)], label="SSP585", color="darkorchid")
    ax.set_xlabel("Year")
    ax.set_ylabel("Global mean squared error (m)")
    ax.legend(loc="upper left", frameon=False)
    return fig


def evaluate_patterns(
        zostoga_list: list, zos_list: list, slopes: list,
        plot_path: Path, model: str) -> None:
    """Evaluate the patterns.

    Parameters
    ----------
    zostoga_list: list
        list of zostoga scenarios
    zos_list: list
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
    # Scale patterns for each scenario
    mse_list = []
    diff_list = []
    for index, (zostoga, slope) in enumerate(zip(
                                                zostoga_list,
                                                slopes,
                                                strict=True)):
        zostoga.data = zostoga.data - np.mean(zostoga.data[0:10])
        zos_list[index].data = (zos_list[index].data - 
                                np.mean(zos_list[index].data[0:10], axis=0))

        p_scaled = (zostoga.data[:, np.newaxis, np.newaxis] *
                    slope[np.newaxis, :, :])

        # Diff maps for end of century (20 yr mean)
        scenario_yrs = 86
        months = 12
        end_idx = (scenario_yrs * months) - 1
        start_idx = end_idx - (20 * months)
        diff_list.append(
            np.mean(
                (zos_list[index][start_idx:end_idx].data +
                zostoga.data[start_idx:end_idx, np.newaxis, np.newaxis]),
                axis=0)
            - np.mean(p_scaled[start_idx:end_idx], axis=0))

        # Calculate mean squared error
        mse = np.nanmean(
            ((zos_list[index].data +
            zostoga.data[:, np.newaxis, np.newaxis]) - p_scaled) ** 2,
            axis=(1, 2))
        mse_list.append(mse)

    fig = plot_evals(diff_list, zos_list, mse_list)
    fig.savefig(Path(plot_path) / f"mse_{model}.png", dpi=150)


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
    zostoga_names = [
        "zostoga_piControl", "zostoga_245",
        "zostoga_370", "zostoga_585"]
    zos_names = ["zos_245", "zos_370","zos_585"]
    zostoga_cubes = []
    zos_cubes = []
    for dataset in cfg["input_data"].values():
        if dataset["dataset"] == model:
            for zostoga_name in zostoga_names:
                if dataset["variable_group"] == zostoga_name:
                    input_file = dataset["filename"]
                    zostoga_cubes.append(sf.load_cube(input_file))

            for zos_name in zos_names:
                if dataset["variable_group"] == zos_name:
                    input_file = dataset["filename"]
                    zos_cubes.append(sf.load_cube(input_file))

    return zostoga_cubes, zos_cubes


def prepare_zostoga(
        zostoga_list: list, plot_path: Path) -> list:
    """Prepare the zostoga cube for regression.

    Parameters
    ----------
    zostoga_list: list
        list of raw zostoga cubes
    plot_path: Path
        path to plot folder

    Returns
    -------
    zostoga_list: list
        list of zostoga cubes
    """
    # Calculate drift from PiControl zostoga
    zostoga_drift = calculate_drift(zostoga_list[0])

    # Detrend the scenario zostogas
    return [
        detrend_zostoga(zostoga, zostoga_drift, plot_path)
        for zostoga in zostoga_list[1:]]  # [1:] to skip PiControl


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

    # Prepare zostoga for regression
    zostoga_list = prepare_zostoga(zostoga_list, plot_path)

    # Calculate regression between zostoga and zos
    scenarios = ["ssp245", "ssp370", "ssp585"]
    slopes, masks = [], []
    for index, (z_dtr, zos) in enumerate(zip(
                                            zostoga_list,
                                            zos_list,
                                            strict=True)):
        slopes_arr, masks_arr = dyn_steric_regression(
            z_dtr, zos, plot_path, scenarios[index])
        save_data(slopes_arr, masks_arr, work_path, model, scenarios[index])
        slopes.append(slopes_arr)
        masks.append(masks_arr)

    # Test the patterns
    evaluate_patterns(zostoga_list, zos_list, slopes, plot_path, model)


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

    sf.parallelise(patterns, processes=8)(models, cfg)

    # Log provenance
    model_work_dir = Path(cfg["work_dir"])
    provenance_record = get_provenance_record()
    path = Path(model_work_dir / "patterns.nc")
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(path, provenance_record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
