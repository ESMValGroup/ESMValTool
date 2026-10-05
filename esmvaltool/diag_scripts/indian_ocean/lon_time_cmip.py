import logging
import os
from pathlib import Path

import cmocean
import iris  # type: ignore
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from basic_functions import (
    get_provenance_record,
    iso_depth_4d,
    load_and_update_dict,
)
from matplotlib.patches import Rectangle

from esmvaltool.diag_scripts.shared import (  # type: ignore
    group_metadata,
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

COASTAL_MASK_LAT_MAX = -4
COASTAL_MASK_LON_MIN = 54
COASTAL_MASK_LON_MAX = 56
OBS_DATASETS = {"NCEP", "HadISST", "EN4"}


def _replace_fill_values(cube, fill_value=1e20):
    """Replace values >= fill_value with NaN in a cube's data array."""
    data = np.asarray(cube.data, dtype=float)
    data[data >= fill_value] = np.nan
    cube.data = ma.masked_invalid(data)
    return cube


def _collapse_latitude(cube, lat_band=None):
    """Collapse latitude by mean while ignoring masked/invalid values.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube with latitude coordinate.
    lat_band : float, optional
        If specified, only collapse over latitudes within 
        ±lat_band from the equator.
        E.g., lat_band=2.0 collapses over [-2, 2] latitude.

    Returns
    -------
    collapsed_cube : iris.cube.Cube
        Cube with latitude collapsed by mean.
    """
    try:
        cube.coord("latitude")
        if lat_band is not None:
            # Extract only the specified latitude band
            cube = cube.extract(
                iris.Constraint(latitude=lambda x: -lat_band <= x <= lat_band)
            )
            if cube is None:
                logger.warning(
                    f"Could not extract latitude band ±{lat_band}"
                )
                return None
        masked = ma.masked_invalid(np.asarray(cube.data, dtype=float))
        return cube.copy(data=masked).collapsed("latitude", iris.analysis.MEAN)
    except iris.exceptions.CoordinateNotFoundError:
        return cube


def _compute_multimodel_bias_and_std(model_cubes, obs_data):
    """
    Compute multi-model median bias and inter-model std dev against obs_data.

    Parameters
    ----------
    model_cubes : list of iris.cube.Cube
        Model cubes with fill values already replaced by NaN.
    obs_data : np.ndarray
        Observation array with the same spatial/temporal shape as each model
        cube.

    Returns
    -------
    mm_median_bias : np.ndarray
        Multi-model median of (model - obs).
    mm_std : np.ndarray
        Inter-model standard deviation of model values.
    """
    bias_list = []
    model_list = []
    for mc in model_cubes:
        model_data = np.ma.filled(np.ma.asarray(mc.data, dtype=float), np.nan)
        bias_list.append(model_data - obs_data)
        model_list.append(model_data)
    bias_stack = np.array(bias_list)  # (n_models, ...)
    model_stack = np.array(model_list)  # (n_models, ...)
    mm_median_bias = np.nanmedian(bias_stack, axis=0)
    mm_std = np.nanstd(model_stack, axis=0)
    return mm_median_bias, mm_std


def _build_lontime_panel_spec(
    data,
    cmap,
    colorbar_label,
    panel_title,
    vmin=None,
    vmax=None,
):
    """Bundle panel rendering settings into a single object."""
    return {
        "data": data,
        "cmap": cmap,
        "colorbar_label": colorbar_label,
        "panel_title": panel_title,
        "vmin": vmin,
        "vmax": vmax,
    }


def _build_monthly_map_spec(
    cmap,
    title,
    output_basename,
    input_filenames,
    vmin=None,
    vmax=None,
    lat_band=2.0,
):
    """Bundle monthly map rendering settings into a single object."""
    return {
        "cmap": cmap,
        "title": title,
        "output_basename": output_basename,
        "input_filenames": input_filenames,
        "vmin": vmin,
        "vmax": vmax,
        "lat_band": lat_band,
    }


def _separate_obs_and_models(plot_dict):
    """
    Split a plot_dict into
    (obs_name, obs_cube, obs_dataset, model_cubes, input_filenames).
    """
    obs_entry = None
    model_cubes = []
    input_filenames = set()
    for dataset, info in plot_dict.items():
        file = info["filename"]
        print(f"Processing dataset {dataset} with file(s): {file}")
        input_filenames.update(file if isinstance(file, list) else [file])
        cube = info["cube"]
        if os.path.basename(file).startswith("OBS"):
            obs_entry = (dataset, cube)
        else:
            model_cubes.append(cube)
    return obs_entry, model_cubes, input_filenames


def _plot_lontime_panel(ax, fig, lon_points, month_points, panel_spec):
    """
    Render a single lon-time pcolormesh panel with a horizontal colorbar.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw onto.
    fig : matplotlib.figure.Figure
        Parent figure (needed for colorbar).
    lon_points, month_points : array-like
        Longitude and month coordinate values taken directly from the cube.
    panel_spec : dict
        Panel metadata and data to render.
    """
    data = panel_spec["data"]
    cmap = panel_spec["cmap"]
    colorbar_label = panel_spec["colorbar_label"]
    panel_title = panel_spec["panel_title"]
    vmin = panel_spec["vmin"]
    vmax = panel_spec["vmax"]
    month_labels = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    n_months = data.shape[0]

    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)

    im = ax.pcolormesh(
        lon_points,
        month_points,
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="nearest",
    )
    fig.colorbar(
        im, ax=ax, orientation="horizontal", pad=0.15, label=colorbar_label
    )
    ax.set_title(panel_title, fontsize=12)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Month")
    ax.set_yticks(list(range(1, n_months + 1))[:n_months])
    ax.set_yticklabels(month_labels[:n_months], fontsize=9)
    return im


def _prepare_lontime_inputs(plot_dict, lat_band):
    """Split obs and model cubes and collapse latitude for lon-time plots."""
    obs_entry, model_cubes, input_filenames = (
        _separate_obs_and_models(plot_dict)
    )
    if obs_entry is None:
        return None, None, None, None, None

    obs_dataset, raw_obs_cube = obs_entry
    obs_cube = _collapse_latitude(
        _replace_fill_values(raw_obs_cube.copy()), lat_band=lat_band
    )
    if obs_cube is None:
        return None, None, None, None, None

    model_cubes_processed = []
    for mc in model_cubes:
        processed = _collapse_latitude(
            _replace_fill_values(mc.copy()), lat_band=lat_band
        )
        if processed is not None:
            model_cubes_processed.append(processed)

    if not model_cubes_processed:
        return None, None, None, None, None

    obs_data = np.ma.filled(np.ma.asarray(obs_cube.data, dtype=float), np.nan)
    return (
        obs_dataset,
        obs_cube,
        model_cubes_processed,
        obs_data,
        input_filenames,
    )


def _collect_highlight_biases(
    plot_dict, highlight_datasets, obs_dataset, obs_data, lat_band=None
):
    """Compute bias arrays for highlighted datasets."""
    highlight_biases = {}
    for dataset, info in plot_dict.items():
        if dataset in highlight_datasets and dataset != obs_dataset:
            cube = info["cube"]
            cube_processed = _collapse_latitude(
                _replace_fill_values(cube.copy()), lat_band=lat_band
            )
            if cube_processed is None:
                continue
            model_data = np.ma.filled(
                np.ma.asarray(cube_processed.data, dtype=float), np.nan
            )
            highlight_biases[dataset] = model_data - obs_data
    return highlight_biases


def _derive_lontime_limits(
        variable, obs_data, mm_median_bias, highlight_biases):
    """Derive colour limits for the lon-time panels."""
    if variable == "tos":
        obs_limits = (25.0, 30.0)
    elif variable in ("ua", "va", "wind"):
        vabs = np.nanmax(np.abs(obs_data))
        obs_limits = (-vabs, vabs)
    else:
        obs_limits = (np.nanmin(obs_data), np.nanmax(obs_data))

    max_bias = np.nanmax(np.abs(mm_median_bias))
    for highlight_bias in highlight_biases.values():
        max_bias = max(max_bias, np.nanmax(np.abs(highlight_bias)))

    return (
        obs_limits[0],
        obs_limits[1],
        round(max_bias, 1),
        np.nanmax(np.abs(mm_median_bias)),
    )


def _render_lontime_panels(
    panel_context,
):
    """Render and save the lon-time panel figure."""
    cfg = panel_context["cfg"]
    obs_dataset = panel_context["obs_dataset"]
    obs_data = panel_context["obs_data"]
    mm_median_bias = panel_context["mm_median_bias"]
    mm_std = panel_context["mm_std"]
    highlight_biases = panel_context["highlight_biases"]
    lon_points = panel_context["lon_points"]
    month_points = panel_context["month_points"]
    cmap_list = panel_context["cmap_list"]
    title = panel_context["title"]
    output_basename = panel_context["output_basename"]
    input_filenames = panel_context["input_filenames"]
    obs_vmin = panel_context["obs_vmin"]
    obs_vmax = panel_context["obs_vmax"]
    bias_vlim = panel_context["bias_vlim"]
    std_vmax = panel_context["std_vmax"]

    n_rows = 2 if highlight_biases else 1
    fig, axes_raw = plt.subplots(
        n_rows, 3, figsize=(18, 5 * n_rows), constrained_layout=True
    )
    axes = axes_raw.reshape(1, -1) if n_rows == 1 else axes_raw

    _plot_lontime_panel(
        axes[0, 0],
        fig,
        lon_points,
        month_points,
        _build_lontime_panel_spec(
            obs_data,
            cmap_list[0],
            obs_dataset,
            f"Obs ({obs_dataset})",
            vmin=obs_vmin,
            vmax=obs_vmax,
        ),
    )
    _plot_lontime_panel(
        axes[0, 1],
        fig,
        lon_points,
        month_points,
        _build_lontime_panel_spec(
            mm_median_bias,
            cmap_list[1],
            "Bias (model - obs)",
            "MM-median bias",
            vmin=-bias_vlim,
            vmax=bias_vlim,
        ),
    )
    _plot_lontime_panel(
        axes[0, 2],
        fig,
        lon_points,
        month_points,
        _build_lontime_panel_spec(
            mm_std,
            cmap_list[2],
            "Std dev (models)",
            "Inter-model std dev",
            vmin=0,
            vmax=std_vmax,
        ),
    )

    for i, dataset in enumerate(list(highlight_biases.keys())[:3]):
        highlight_bias = highlight_biases[dataset]
        _plot_lontime_panel(
            axes[1, i],
            fig,
            lon_points,
            month_points,
            _build_lontime_panel_spec(
                highlight_bias,
                cmap_list[1],
                f"Bias {dataset}",
                f"{dataset} bias",
                vmin=-bias_vlim,
                vmax=bias_vlim,
            ),
        )

    if highlight_biases:
        if len(highlight_biases) == 2:
            bias_values = list(highlight_biases.values())
            bias_diff = bias_values[0] - bias_values[1]
            dataset_names = list(highlight_biases.keys())
            _plot_lontime_panel(
                axes[1, 2],
                fig,
                lon_points,
                month_points,
                _build_lontime_panel_spec(
                    bias_diff,
                    cmap_list[1],
                    "Difference",
                    f"{dataset_names[0]} - {dataset_names[1]}",
                    vmin=-bias_vlim,
                    vmax=bias_vlim,
                ),
            )
        else:
            for i in range(len(highlight_biases), 3):
                axes[1, i].axis("off")

    fig.suptitle(title, fontsize=14)
    provenance_record = get_provenance_record(
        output_basename, sorted(list(input_filenames))
    )
    save_figure(output_basename, provenance_record, cfg, bbox_inches="tight")
    logger.info("Lon-time plot saved: %s", output_basename)
    plt.close(fig)


def plot_lon_time_multimodel(
    cfg,
    plot_dict,
    plot_spec,
):
    """
    Plot lon-time climatology: obs, MM-median bias, and inter-model std dev.

    Creates a three-panel figure:
      Panel 1 — Observed climatology (lon x month).
      Panel 2 — Multi-model median bias (model - obs).
      Panel 3 — Inter-model standard deviation.
      Panel 4-6 - Bias or comparison of up to three highlighted 
                datasets (if specified in cfg).

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    plot_spec : dict
        Plot settings including colormaps, titles, limits, and latitude band.
    """
    cmap_list = plot_spec["cmap_list"]
    title = plot_spec["title"]
    output_basename = plot_spec["output_basename"]
    variable = plot_spec.get("variable")
    obs_vmin = plot_spec.get("obs_vmin")
    obs_vmax = plot_spec.get("obs_vmax")
    bias_vlim = plot_spec.get("bias_vlim")
    std_vmax = plot_spec.get("std_vmax")
    lat_band = plot_spec.get("lat_band")

    logger.info("Plotting lon-time plots: %s", output_basename)

    prepared = _prepare_lontime_inputs(plot_dict, lat_band)
    if prepared[0] is None:
        logger.warning(
            "Unable to prepare lon-time plot data for %s, skipping.",
            output_basename,
        )
        return

    obs_dataset, obs_cube, model_cubes, obs_data, input_filenames = prepared
    mm_median_bias, mm_std = _compute_multimodel_bias_and_std(
        model_cubes, obs_data
    )

    highlight_datasets = cfg.get("highlight_datasets", [])[:3]
    highlight_biases = _collect_highlight_biases(
        plot_dict,
        highlight_datasets,
        obs_dataset,
        obs_data,
        lat_band=lat_band,
    )

    try:
        lon_centres = obs_cube.coord("longitude").points
    except iris.exceptions.CoordinateNotFoundError:
        lon_centres = np.arange(obs_data.shape[1])
    try:
        month_centres = obs_cube.coord("month_number").points
    except iris.exceptions.CoordinateNotFoundError:
        month_centres = np.arange(1, obs_data.shape[0] + 1)

    obs_vmin, obs_vmax, derived_bias_vlim, derived_std_vmax = (
        _derive_lontime_limits(variable, obs_data, 
                               mm_median_bias, highlight_biases)
    )
    bias_vlim = bias_vlim if bias_vlim is not None else derived_bias_vlim*0.8
    std_vmax = std_vmax if std_vmax is not None else derived_std_vmax*0.8

    _render_lontime_panels(
        {
            "cfg": cfg,
            "obs_dataset": obs_dataset,
            "obs_data": obs_data,
            "mm_median_bias": mm_median_bias,
            "mm_std": mm_std,
            "highlight_biases": highlight_biases,
            "lon_points": lon_centres,
            "month_points": month_centres,
            "cmap_list": cmap_list,
            "title": title,
            "output_basename": output_basename,
            "input_filenames": input_filenames,
            "obs_vmin": obs_vmin,
            "obs_vmax": obs_vmax,
            "bias_vlim": bias_vlim,
            "std_vmax": std_vmax,
        }
    )


def plot_monthly_maps(
    cfg,
    data,
    lon_centres,
    lat_centres,
    plot_spec,
):
    """
    Plot 12 monthly map panels (3 rows x 4 columns) with a shared colorbar.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    data : np.ndarray
        3-D array shaped (month, latitude, longitude).
    lon_centres, lat_centres : array-like
        Coordinate centre points.
    plot_spec : dict
        Plot settings including colormap, title, output name, limits, and
        provenance inputs.
    """
    cmap = plot_spec["cmap"]
    title = plot_spec["title"]
    output_basename = plot_spec["output_basename"]
    input_filenames = plot_spec["input_filenames"]
    vmin = plot_spec.get("vmin")
    vmax = plot_spec.get("vmax")
    lat_band = plot_spec.get("lat_band", 2.0)

    data = np.asarray(data)
    is_3d = data.ndim == 3
    if not is_3d:
        logger.warning(
            "Expected 3D data (month, lat, lon) for %s, got shape %s",
            output_basename,
            data.shape,
        )
        return

    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)

    n_months = min(12, data.shape[0])
    month_labels = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    fig, axes = plt.subplots(3, 4, figsize=(16, 10), constrained_layout=True)
    axes = axes.flatten()

    im = None
    for m in range(12):
        ax = axes[m]
        if m < n_months:
            im = ax.pcolormesh(
                lon_centres,
                lat_centres,
                data[m, :, :],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                shading="nearest",
            )
            ax.set_title(month_labels[m], fontsize=10)
            ax.set_xlabel("Longitude (°E)")
            ax.set_ylabel("Latitude (°N)")
            # Add black box highlighting equatorial band
            rect = Rectangle(
                (lon_centres[0], -lat_band),
                lon_centres[-1] - lon_centres[0],
                2 * lat_band,
                linewidth=2,
                edgecolor="black",
                facecolor="none",
            )
            ax.add_patch(rect)
        else:
            ax.axis("off")

    if im is not None:
        cbar = fig.colorbar(
            im,
            ax=axes.tolist(),
            orientation="horizontal",
            pad=0.06,
            shrink=0.9,
        )
        cbar.set_label(title)

    fig.suptitle(title, fontsize=14)
    provenance_record = get_provenance_record(
        output_basename, sorted(list(input_filenames))
    )
    save_figure(output_basename, provenance_record, cfg, bbox_inches="tight")
    logger.info("Monthly map saved: %s", output_basename)
    plt.close(fig)


def plot_map_multimodel(
    cfg, plot_dict, cmap_list, title, output_basename, lat_band=2.0
):
    """
    Plot monthly maps for obs, MM-median bias, and inter-model std dev.

    Calls plot_monthly_maps three times — once per panel type — saving
    separate figures for obs, bias, and std dev.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    cmap_list : list of str
        Three colourmap names for [obs, bias, std dev].
    title : str
        Base title string appended with the panel type for each figure.
    output_basename : str
        Base stem for output filenames; suffixes '_obs', '_bias', '_stddev' are
        appended.
    lat_band : float, optional
        Latitude band (symmetric about equator) to highlight with a box
        (default 2.0 for ±2°).
    """
    logger.info("Plotting monthly map plots: %s", output_basename)

    obs_entry, model_cubes, input_filenames = _separate_obs_and_models(
        plot_dict
    )
    if obs_entry is None:
        logger.warning(
            "No obs found for map plot %s, skipping.", output_basename
        )
        return
    if not model_cubes:
        logger.warning(
            "No model data found for map plot %s, skipping.", output_basename
        )
        return

    obs_name, obs_cube = obs_entry
    obs_cube = _replace_fill_values(obs_cube.copy())
    model_cubes = [_replace_fill_values(mc.copy()) for mc in model_cubes]
    obs_data = np.ma.filled(np.ma.asarray(obs_cube.data, dtype=float), np.nan)

    mm_median_bias, mm_std = _compute_multimodel_bias_and_std(
        model_cubes, obs_data
    )

    try:
        lon_centres = obs_cube.coord("longitude").points
    except iris.exceptions.CoordinateNotFoundError:
        lon_centres = np.arange(obs_data.shape[-1])
    try:
        lat_centres = obs_cube.coord("latitude").points
    except iris.exceptions.CoordinateNotFoundError:
        lat_centres = np.arange(obs_data.shape[-2])

    plot_monthly_maps(
        cfg,
        obs_data,
        lon_centres,
        lat_centres,
        _build_monthly_map_spec(
            cmap_list[0],
            f"Obs ({obs_name}) - {title}",
            output_basename + "_clim_obs",
            input_filenames,
            vmin=np.nanmin(obs_data),
            vmax=np.nanmax(obs_data),
            lat_band=lat_band,
        ),
    )
    bias_vlim = np.nanmax(np.abs(mm_median_bias))
    plot_monthly_maps(
        cfg,
        mm_median_bias,
        lon_centres,
        lat_centres,
        _build_monthly_map_spec(
            cmap_list[1],
            f"MM-median bias - {title}",
            output_basename + "_model_bias",
            input_filenames,
            vmin=-bias_vlim,
            vmax=bias_vlim,
            lat_band=lat_band,
        ),
    )
    plot_monthly_maps(
        cfg,
        mm_std,
        lon_centres,
        lat_centres,
        _build_monthly_map_spec(
            cmap_list[2],
            f"Inter-model std dev - {title}",
            output_basename + "_stddev",
            input_filenames,
            vmin=0,
            vmax=np.nanmax(mm_std),
            lat_band=lat_band,
        ),
    )


def _extract_pressure_level(plot_dict, pressure_level):
    """
    Extract a specific pressure level from all cubes in plot_dict.

    Parameters
    ----------
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    pressure_level : float
        Pressure level to extract (e.g., 100000. for surface, 85000., 20000.).

    Returns
    -------
    new_plot_dict : dict
        New plot_dict with cubes extracted at the specified pressure level.
    """
    new_plot_dict = {}
    for dataset, info in plot_dict.items():
        cube = info["cube"]
        try:
            # Extract the specific pressure level
            extracted_cube = cube.extract(
                iris.Constraint(air_pressure=pressure_level)
            )
            if extracted_cube is None:
                logger.warning(
                    f"Could not extract pressure level "
                    f"{pressure_level} from {dataset}"
                )
                continue
            new_plot_dict[dataset] = {
                "cube": extracted_cube,
                "filename": info["filename"],
            }
        except Exception as e:
            logger.warning(
                f"Error extracting pressure level "
                f"{pressure_level} from {dataset}: {e}"
            )
            continue
    return new_plot_dict


def _compute_wind_shear(plot_dict, pressure_top, pressure_bottom):
    """
    Compute wind shear (top - bottom) between two pressure levels.

    Parameters
    ----------
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    pressure_top : float
        Upper pressure level (e.g., 20000. for 200 hPa).
    pressure_bottom : float
        Lower pressure level (e.g., 85000. for 850 hPa).

    Returns
    -------
    shear_plot_dict : dict
        New plot_dict with wind shear (u_top - u_bottom) cubes.
    """
    shear_plot_dict = {}
    for dataset, info in plot_dict.items():
        cube = info["cube"]
        try:
            # Extract both pressure levels
            cube_top = cube.extract(iris.Constraint(air_pressure=pressure_top))
            cube_bottom = cube.extract(
                iris.Constraint(air_pressure=pressure_bottom)
            )

            if cube_top is None or cube_bottom is None:
                logger.warning(
                    f"Could not extract both pressure levels from {dataset}"
                )
                continue

            # Compute shear
            shear_cube = cube_top - cube_bottom
            shear_plot_dict[dataset] = {
                "cube": shear_cube,
                "filename": info["filename"],
            }
        except Exception as e:
            logger.warning(f"Error computing wind shear for {dataset}: {e}")
            continue
    return shear_plot_dict


def _create_iso_depth_dict(
    cfg,
    plot_dict,
    iso_level=20.0,
    time_measure="month_number",
    apply_coastal_mask=True,
):
    """Create a plot_dict with 4D cubes converted to isotherm depth.

    Optionally applies a robust local-neighborhood coastal mask to remove
    isolated shallow outliers (Seychelles) that can skew 
    zonal means and multi-model stats.
    """
    new_plot_dict = {}
    for dataset, info in plot_dict.items():
        cube = info["cube"]
        is_4d = cube.ndim == 4
        if is_4d:
            iso_cube = iso_depth_4d(cube, iso_level, time_measure=time_measure)
            iso_cube = _replace_fill_values(iso_cube)
            if apply_coastal_mask:
                lat_all = iso_cube.coord("latitude").points
                lon_all = iso_cube.coord("longitude").points
                lat_mask = lat_all <= COASTAL_MASK_LAT_MAX
                lon_mask = (lon_all >= COASTAL_MASK_LON_MIN) \
                            & (lon_all <= COASTAL_MASK_LON_MAX)

                box_mask = lat_mask[:, None] & lon_mask[None, :]  # (lat, lon)

                # Mask out the box
                iso_cube.data[:, box_mask] = np.nan

            new_plot_dict[dataset] = {
                "cube": iso_cube,
                "filename": info["filename"],
            }
        else:
            new_plot_dict[dataset] = info
    return new_plot_dict


def main(cfg):
    """
    Plot monthly climatologies for multiple datasets and observations.
    """
    LAT_BAND_AVG = cfg.get(
        "lat_band_avg", 2.0
    )  # Default to ±2° if not specified
    input_data = cfg["input_data"].values()
    grouped_data = group_metadata(input_data, "dataset")
    eio_wind_monthly, eio_sst_monthly, eio_theta_monthly, eio_pr_monthly = (
        {},
        {},
        {},
        {},
    )
    for _group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, "eio_wind_monthly", eio_wind_monthly)
        load_and_update_dict(group_md, "eio_sst_monthly", eio_sst_monthly)
        load_and_update_dict(group_md, "eio_theta_monthly", eio_theta_monthly)
        load_and_update_dict(group_md, "eio_pr_monthly", eio_pr_monthly)

    logger.info("Data loaded, now plotting.")

    eio_t20d_monthly = _create_iso_depth_dict(
        cfg,
        eio_theta_monthly,
        iso_level=20.0,
        time_measure="month_number",
        apply_coastal_mask=True,
    )

    # Extract surface wind (100000 Pa) 
    # Compute wind shear (200 hPa - 850 hPa)
    eio_wind_surface = _extract_pressure_level(eio_wind_monthly, 100000.0)
    eio_wind_850 = _extract_pressure_level(eio_wind_monthly, 85000.0)
    eio_wind_200 = _extract_pressure_level(eio_wind_monthly, 20000.0)
    eio_wind_shear = _compute_wind_shear(eio_wind_monthly, 20000.0, 85000.0)

    plot_lon_time_multimodel(
        cfg,
        eio_wind_surface,
        {
            "cmap_list": ["cmo.delta", "BrBG", "RdPu"],
            "title": (
                "Indian Ocean equatorial zonal wind (surface, 1000 hPa) "
                "— monthly climatology"
            ),
            "output_basename": "lon_time_eio_wind_surface",
            "variable": "ua",
            "lat_band": LAT_BAND_AVG,
        },
    )
    plot_lon_time_multimodel(
        cfg,
        eio_wind_shear,
        {
            "cmap_list": ["cmo.delta", "BrBG", "RdPu"],
            "title": (
                "Indian Ocean equatorial zonal wind shear "
                "(200 hPa - 850 hPa) — monthly climatology"
            ),
            "output_basename": "lon_time_eio_wind_shear",
            "variable": "ua",
            "lat_band": LAT_BAND_AVG,
        },
    )
    plot_lon_time_multimodel(
        cfg,
        eio_sst_monthly,
        {
            "cmap_list": ["RdYlBu_r", "RdBu_r", "RdPu"],
            "title": ("Indian Ocean equatorial SST " "— monthly climatology"),
            "output_basename": "lon_time_eio_sst",
            "variable": "tos",
            "lat_band": LAT_BAND_AVG,
        },
    )
    plot_lon_time_multimodel(
        cfg,
        eio_t20d_monthly,
        {
            "cmap_list": ["cmo.deep", "cmo.tarn", "RdPu"],
            "title": (
                "Indian Ocean equatorial 20°C isotherm depth "
                "— monthly climatology"
            ),
            "output_basename": "lon_time_eio_t20d",
            "variable": "t20d",
            "lat_band": LAT_BAND_AVG,
        },
    )
    plot_lon_time_multimodel(
        cfg,
        eio_pr_monthly,
        {
            "cmap_list": ["cmo.rain", "BrBG", "RdPu"],
            "title": (
                "Indian Ocean equatorial precipitation "
                "— monthly climatology"
            ),
            "output_basename": "lon_time_eio_pr",
            "variable": "pr",
            "bias_vlim": 0.00008,
            "lat_band": LAT_BAND_AVG,
        },
    )

    plot_map_multimodel(
        cfg,
        eio_wind_surface,
        cmap_list=["BrBG", "BrBG", "RdPu"],
        title="Indian Ocean equatorial zonal wind — monthly climatology",
        output_basename="map_eio_wind_surface",
        lat_band=LAT_BAND_AVG,
    )
    plot_map_multimodel(
        cfg,
        eio_sst_monthly,
        cmap_list=["RdYlBu_r", "RdBu_r", "RdPu"],
        title="Indian Ocean equatorial SST — monthly climatology",
        output_basename="map_eio_sst",
        lat_band=LAT_BAND_AVG,
    )
    plot_map_multimodel(
        cfg,
        eio_t20d_monthly,
        cmap_list=["cmo.deep", "cmo.tarn", "RdPu"],
        title="Indian Ocean equatorial 20°C " \
        "isotherm depth — monthly climatology",
        output_basename="map_eio_t20d",
        lat_band=LAT_BAND_AVG,
    )
    plot_map_multimodel(
        cfg,
        eio_pr_monthly,
        cmap_list=["cmo.rain", "BrBG", "RdPu"],
        title="Indian Ocean equatorial precipitation — monthly climatology",
        output_basename="map_eio_pr",
        lat_band=LAT_BAND_AVG,
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
