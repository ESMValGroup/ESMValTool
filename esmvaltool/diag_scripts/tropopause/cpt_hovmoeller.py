"""Create tropical CPT Hovmoeller plots with contour overlays.

This diagnostic identifies the tropical cold-point tropopause (CPT) by
searching for the minimum temperature between configurable pressure bounds and
creates, per dataset, one standalone Hovmoeller plot per variable for each of
two views:

- longitude-time (mean over the tropical latitude band, full longitude range)
- latitude-time (restricted to the tropical latitude band, mean over the
  full longitude range)

In both views:

- filled colors (all plots): CPT pressure
- line contours: CPT temperature, CPT specific humidity, or CPT RH over ice
  (one variable per figure)

Configuration options
---------------------
- temperature_var: Short name for temperature variable (default: ``ta``).
- humidity_var: Short name for specific humidity variable (default: ``hus``).
- geopotential_height_var: Optional short name for geopotential height
  (default: ``zg``). Used for output only.
- pressure_bounds_hpa: Two values [top, bottom] in hPa for CPT search
  (default: [10.0, 500.0]).
- tropical_lat_bounds: Two values [south, north] in degrees (default:
  [-20.0, 20.0]).
- group_by: Metadata facet used to group inputs (default: ``alias``).
- lon_resolution: Longitude bin width in degrees (default: 2.5).
- lat_resolution: Latitude bin width in degrees  (default: 2.5).
- plot_filename: Prefix for output file names
  (default: ``tropical_cpt_hovmoeller``).
- dpi: Figure DPI (default: 300).
"""

from __future__ import annotations

import calendar
import logging
import re
from pathlib import Path

import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
from iris.cube import Cube

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)


def _safe_name(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z_.-]+", "_", text)


def _get_provenance(caption: str, ancestors: list[str]) -> dict:
    return {
        "caption": caption,
        "statistics": ["mean", "other"],
        "domains": ["trop"],
        "plot_types": ["zonal"],
        "authors": ["mcquaid_aoife"],
        "ancestors": ancestors,
    }


def _get_variable_attrs(group: list[dict], short_name: str) -> dict:
    for attrs in group:
        if attrs.get("short_name") == short_name:
            return attrs
    available = sorted({attrs.get("short_name") for attrs in group})
    msg = (
        f"Required variable '{short_name}' not found in grouped input; "
        f"available short_names: {available}"
    )
    raise ValueError(msg)


def _get_pressure_coord(cube: Cube):
    for candidate in ("air_pressure", "plev", "pressure"):
        if cube.coords(candidate):
            return cube.coord(candidate)
    for coord in cube.coords(axis="Z"):
        if coord.units.is_convertible("Pa"):
            return coord
    msg = f"Could not identify pressure coordinate for cube '{cube.summary(shorten=True)}'"
    raise ValueError(msg)


def _build_3d_cube_from_level_selection(
    source_cube: Cube,
    z_coord_name: str,
    data_3d: np.ndarray,
    *,
    var_name: str,
    long_name: str,
    units: str,
) -> Cube:
    z_dim = source_cube.coord_dims(z_coord_name)[0]
    index = [slice(None)] * source_cube.ndim
    index[z_dim] = 0
    out_cube = source_cube[tuple(index)].copy(data=data_3d)
    if out_cube.coords(z_coord_name):
        out_cube.remove_coord(z_coord_name)
    out_cube.var_name = var_name
    out_cube.long_name = long_name
    out_cube.units = units
    return out_cube


def _select_at_cpt(
    ta_cube: Cube,
    target_cube: Cube,
    p_bounds_hpa: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    p_coord = _get_pressure_coord(ta_cube)
    ta_z_dim = ta_cube.coord_dims(p_coord.name())[0]

    p_hpa = p_coord.copy()
    p_hpa.convert_units("hPa")
    p_points = p_hpa.points
    pmin, pmax = min(p_bounds_hpa), max(p_bounds_hpa)
    level_mask = (p_points >= pmin) & (p_points <= pmax)
    if not np.any(level_mask):
        msg = (
            f"No pressure levels found between {pmin} and {pmax} hPa "
            f"for coordinate '{p_coord.name()}'"
        )
        raise ValueError(msg)

    ta_data = np.asarray(ta_cube.core_data())
    target_data = np.asarray(target_cube.core_data())
    selector_shape = [1] * ta_cube.ndim
    selector_shape[ta_z_dim] = level_mask.size
    level_selector = level_mask.reshape(selector_shape)

    ta_masked = np.where(
        level_selector & np.isfinite(ta_data),
        ta_data,
        np.inf,
    )
    invalid = np.all(np.isinf(ta_masked), axis=ta_z_dim)
    cpt_idx = np.argmin(ta_masked, axis=ta_z_dim)
    target_z_dim = target_cube.coord_dims(
        _get_pressure_coord(target_cube).name()
    )[0]
    expanded_idx = np.expand_dims(cpt_idx, axis=target_z_dim)
    selected = np.squeeze(
        np.take_along_axis(target_data, expanded_idx, axis=target_z_dim),
        axis=target_z_dim,
    )
    selected = selected.astype(float)
    selected[invalid] = np.nan

    cpt_p_hpa = p_points[cpt_idx].astype(float)
    cpt_p_hpa[invalid] = np.nan

    return selected, cpt_p_hpa


def _compute_rhi_percent(
    specific_humidity: np.ndarray,
    temperature_k: np.ndarray,
    pressure_hpa: np.ndarray,
) -> np.ndarray:
    eps = 0.622
    pressure_pa = pressure_hpa * 100.0
    q = np.asarray(specific_humidity, dtype=float)
    t = np.asarray(temperature_k, dtype=float)
    p = np.asarray(pressure_pa, dtype=float)

    vapor_pressure = q * p / (eps + (1.0 - eps) * q)
    log_es_ice = (
        9.550426 - (5723.265 / t) + (3.53068 * np.log(t)) - (0.00728332 * t)
    )
    es_ice = np.exp(log_es_ice)
    rhi = 100.0 * vapor_pressure / es_ice
    rhi[~np.isfinite(rhi)] = np.nan
    return rhi


def _monthly_tropical_mean(
    cube: Cube,
    lat_bounds: tuple[float, float],
) -> Cube:
    work = cube.copy()
    if not work.coords("year"):
        iris.coord_categorisation.add_year(work, "time", name="year")
    if not work.coords("month_number"):
        iris.coord_categorisation.add_month_number(
            work,
            "time",
            name="month_number",
        )
    monthly = work.aggregated_by(["year", "month_number"], iris.analysis.MEAN)

    lat_name = monthly.coord(axis="Y").name()
    lat_min, lat_max = min(lat_bounds), max(lat_bounds)
    tropical = monthly.extract(
        iris.Constraint(**{lat_name: lambda lat: lat_min <= lat <= lat_max}),
    )
    if tropical is None:
        msg = f"No latitude data found inside tropical band [{lat_min}, {lat_max}]"
        raise ValueError(msg)
    return tropical.collapsed(lat_name, iris.analysis.MEAN)


def _monthly_tropical_mean_lat(
    cube: Cube,
    lat_bounds: tuple[float, float],
) -> Cube:
    """Time-vs-latitude counterpart of ``_monthly_tropical_mean``.

    Restricts latitude to the tropical band (kept as the plotted axis)
    and collapses over the full longitude range with a plain (unweighted)
    mean, mirroring the unweighted latitude collapse used on the
    longitude-time path.
    """
    work = cube.copy()
    if not work.coords("year"):
        iris.coord_categorisation.add_year(work, "time", name="year")
    if not work.coords("month_number"):
        iris.coord_categorisation.add_month_number(
            work,
            "time",
            name="month_number",
        )
    monthly = work.aggregated_by(["year", "month_number"], iris.analysis.MEAN)

    lat_name = monthly.coord(axis="Y").name()
    lon_name = monthly.coord(axis="X").name()
    lat_min, lat_max = min(lat_bounds), max(lat_bounds)
    tropical = monthly.extract(
        iris.Constraint(**{lat_name: lambda lat: lat_min <= lat <= lat_max}),
    )
    if tropical is None:
        msg = f"No latitude data found inside tropical band [{lat_min}, {lat_max}]"
        raise ValueError(msg)
    return tropical.collapsed(lon_name, iris.analysis.MEAN)


def _as_time_axis(cube: Cube, axis: str) -> Cube:
    """Ensure ``cube`` has dims ordered (time, <axis>), axis in {'X', 'Y'}."""
    time_name = cube.coord(axis="T").name()
    other_name = cube.coord(axis=axis).name()
    time_dim = cube.coord_dims(time_name)[0]
    other_dim = cube.coord_dims(other_name)[0]
    if (time_dim, other_dim) == (0, 1):
        return cube
    out_cube = cube.copy()
    out_cube.transpose((time_dim, other_dim))
    return out_cube


def _finite_levels(data: np.ndarray, n_levels: int) -> np.ndarray:
    finite = np.asarray(data, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.linspace(0.0, 1.0, n_levels)
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if np.isclose(vmin, vmax):
        delta = 1.0 if np.isclose(vmin, 0.0) else 0.05 * abs(vmin)
        vmin -= delta
        vmax += delta
    return np.linspace(vmin, vmax, n_levels)


def _plot_hovmoeller_set(
    cfg: dict,
    dataset_label: str,
    cpt_p_cube: Cube,
    cpt_t_cube: Cube,
    cpt_q_cube: Cube,
    cpt_rhi_cube: Cube,
    ancestors: list[str],
    *,
    axis: str,
) -> None:
    """Save one standalone Hovmoeller figure per variable for this dataset.

    Each figure shows CPT pressure as filled colors (shared across all
    variables) with a single variable's CPT contours overlaid, so every
    output PNG is a complete, self-contained plot rather than a subplot.

    ``axis`` selects the plotted view: ``"X"`` for time-vs-longitude
    (``cpt_*_cube`` must already be longitude-resolved, i.e. reduced over
    latitude), or ``"Y"`` for time-vs-latitude (``cpt_*_cube`` must already
    be latitude-resolved, i.e. reduced over longitude).
    """
    cpt_p = _as_time_axis(cpt_p_cube, axis)
    cpt_t = _as_time_axis(cpt_t_cube, axis)
    cpt_q = _as_time_axis(cpt_q_cube, axis)
    cpt_rhi = _as_time_axis(cpt_rhi_cube, axis)

    coord = cpt_p.coord(axis=axis)
    coord_points = coord.points
    coord_name = coord.name()
    coord_units = coord.units

    y = np.arange(cpt_p.shape[0])
    p_data = np.asarray(cpt_p.core_data())

    year = cpt_p.coord("year").points.astype(int)
    month = cpt_p.coord("month_number").points.astype(int)

    y_ticks = np.where((month == 1) | (month == 7))[0]
    if y_ticks.size == 0:
        y_ticks = np.arange(0, len(y), 8)

    month_abbrs = [calendar.month_abbr[m].upper() for m in month]
    y_tick_labels = [f"{month_abbrs[i]}{year[i]}" for i in y_ticks]

    p_levels = _finite_levels(p_data, 21)

    # (filename token, panel title, data cube, clabel format)
    panels = [
        ("temperature", "CPT temperature contours [K]", cpt_t, "%.0f"),
        (
            "specific_humidity",
            "CPT specific humidity contours [x10\u2076 kg kg-1]",
            cpt_q,
            lambda x: f"{x * 1e6:.2f}",
        ),
        ("rhi", "CPT RH over ice contours [%]", cpt_rhi, "%.0f"),
    ]

    suffix = "lon" if axis == "X" else "lat"
    view_desc = "longitude-time" if axis == "X" else "latitude-time"
    basename = f"{cfg.get('plot_filename', 'tropical_cpt_hovmoeller')}_{_safe_name(dataset_label)}"

    for token, title, contour_cube, clabel_fmt in panels:
        contour_data = np.asarray(contour_cube.core_data())
        contour_levels = _finite_levels(contour_data, 8)

        fig, ax = plt.subplots(figsize=(8, 10), constrained_layout=True)
        cf = ax.contourf(
            coord_points,
            y,
            p_data,
            levels=p_levels,
            cmap="Reds_r",
            extend="both",
        )
        cs = ax.contour(
            coord_points,
            y,
            contour_data,
            levels=contour_levels,
            colors="k",
            linewidths=0.8,
        )
        ax.clabel(cs, inline=True, fontsize=12, fmt=clabel_fmt)

        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
        ax.set_title(title, fontsize=16)
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        ax.set_ylabel("Time", fontsize=16)
        ax.set_xlabel(f"{coord_name} [{coord_units}]", fontsize=16)
        fig.suptitle(f"Tropical CPT Hovmoeller ({dataset_label})", fontsize=18)

        cbar = fig.colorbar(cf, ax=ax, orientation="vertical", shrink=0.95)
        cbar.set_label("CPT pressure [hPa]", fontsize=16)
        cbar.set_ticks([80, 90, 100, 110, 120])
        cbar.ax.tick_params(labelsize=12)

        caption = (
            f"Tropical CPT {view_desc} Hovmoeller for {dataset_label}: "
            f"colors show CPT pressure; contours show {title.lower()}."
        )
        provenance = _get_provenance(caption, ancestors)
        save_figure(
            f"{basename}_{token}_{suffix}",
            provenance,
            cfg,
            figure=fig,
            dpi=cfg.get("dpi", 300),
        )


def _dataset_label(group_key: str, group_items: list[dict]) -> str:
    if group_key not in (None, "None"):
        return str(group_key)
    return str(group_items[0].get("dataset", "dataset"))


def main(cfg: dict) -> None:
    cfg.setdefault("temperature_var", "ta")
    cfg.setdefault("humidity_var", "hus")
    cfg.setdefault("geopotential_height_var", "zg")
    cfg.setdefault("pressure_bounds_hpa", [10.0, 500.0])
    cfg.setdefault("tropical_lat_bounds", [-20.0, 20.0])
    cfg.setdefault("group_by", "alias")
    cfg.setdefault("plot_filename", "tropical_cpt_hovmoeller")
    cfg.setdefault("lon_resolution", 2.5)
    cfg.setdefault("lat_resolution", 2.5)

    input_data = list(cfg.get("input_data", {}).values())
    if not input_data:
        raise ValueError(
            "Provide preprocessed ESMValTool inputs",
        )

    grouped = group_metadata(input_data, cfg["group_by"], sort="short_name")

    for group_key, group_items in grouped.items():
        label = _dataset_label(group_key, group_items)
        logger.info("Processing %s", label)
        ta_attrs = _get_variable_attrs(group_items, cfg["temperature_var"])
        hus_attrs = _get_variable_attrs(group_items, cfg["humidity_var"])
        zg_attrs = None
        try:
            zg_attrs = _get_variable_attrs(
                group_items, cfg["geopotential_height_var"]
            )
        except ValueError:
            logger.info(
                "Optional geopotential height variable '%s' not found for %s",
                cfg["geopotential_height_var"],
                label,
            )

        ta_cube = iris.load_cube(ta_attrs["filename"])
        hus_cube = iris.load_cube(hus_attrs["filename"])
        p_bounds_hpa = tuple(cfg["pressure_bounds_hpa"])
        lat_bounds = tuple(cfg["tropical_lat_bounds"])

        cpt_t_data, cpt_p_hpa = _select_at_cpt(ta_cube, ta_cube, p_bounds_hpa)
        cpt_q_data, _ = _select_at_cpt(ta_cube, hus_cube, p_bounds_hpa)

        p_coord_name = _get_pressure_coord(ta_cube).name()
        cpt_t = _build_3d_cube_from_level_selection(
            ta_cube,
            p_coord_name,
            cpt_t_data,
            var_name="cpt_ta",
            long_name="Cold point tropopause temperature",
            units=str(ta_cube.units),
        )
        cpt_p = _build_3d_cube_from_level_selection(
            ta_cube,
            p_coord_name,
            cpt_p_hpa,
            var_name="cpt_p",
            long_name="Cold point tropopause pressure",
            units="hPa",
        )
        cpt_q = _build_3d_cube_from_level_selection(
            hus_cube,
            _get_pressure_coord(hus_cube).name(),
            cpt_q_data,
            var_name="cpt_hus",
            long_name="Cold point tropopause specific humidity",
            units=str(hus_cube.units),
        )
        cpt_rhi_data = _compute_rhi_percent(cpt_q_data, cpt_t_data, cpt_p_hpa)
        cpt_rhi = _build_3d_cube_from_level_selection(
            ta_cube,
            p_coord_name,
            cpt_rhi_data,
            var_name="cpt_rhi",
            long_name="Cold point tropopause relative humidity over ice",
            units="%",
        )

        cpt_z_monthly = None
        if zg_attrs is not None:
            zg_cube = iris.load_cube(zg_attrs["filename"])
            cpt_z_data, _ = _select_at_cpt(ta_cube, zg_cube, p_bounds_hpa)
            cpt_z = _build_3d_cube_from_level_selection(
                zg_cube,
                _get_pressure_coord(zg_cube).name(),
                cpt_z_data,
                var_name="cpt_zg",
                long_name="Cold point tropopause geopotential height",
                units=str(zg_cube.units),
            )
            cpt_z_monthly = _monthly_tropical_mean(cpt_z, lat_bounds)

        cpt_p_monthly_lon = _monthly_tropical_mean(cpt_p, lat_bounds)
        cpt_t_monthly_lon = _monthly_tropical_mean(cpt_t, lat_bounds)
        cpt_q_monthly_lon = _monthly_tropical_mean(cpt_q, lat_bounds)
        cpt_rhi_monthly_lon = _monthly_tropical_mean(cpt_rhi, lat_bounds)

        cpt_p_monthly_lat = _monthly_tropical_mean_lat(cpt_p, lat_bounds)
        cpt_t_monthly_lat = _monthly_tropical_mean_lat(cpt_t, lat_bounds)
        cpt_q_monthly_lat = _monthly_tropical_mean_lat(cpt_q, lat_bounds)
        cpt_rhi_monthly_lat = _monthly_tropical_mean_lat(cpt_rhi, lat_bounds)

        ancestors = [ta_attrs["filename"], hus_attrs["filename"]]
        if zg_attrs is not None:
            ancestors.append(zg_attrs["filename"])

        base = f"{cfg['plot_filename']}_{_safe_name(label)}"
        provenance = _get_provenance(
            "Monthly tropical CPT diagnostics at longitude-time resolution.",
            ancestors,
        )
        save_data(f"{base}_pressure", provenance, cfg, cpt_p_monthly_lon)
        save_data(f"{base}_temperature", provenance, cfg, cpt_t_monthly_lon)
        save_data(
            f"{base}_specific_humidity", provenance, cfg, cpt_q_monthly_lon
        )
        save_data(f"{base}_rhi", provenance, cfg, cpt_rhi_monthly_lon)
        if cpt_z_monthly is not None:
            save_data(f"{base}_height", provenance, cfg, cpt_z_monthly)

        _plot_hovmoeller_set(
            cfg,
            label,
            cpt_p_monthly_lon,
            cpt_t_monthly_lon,
            cpt_q_monthly_lon,
            cpt_rhi_monthly_lon,
            ancestors,
            axis="X",
        )
        _plot_hovmoeller_set(
            cfg,
            label,
            cpt_p_monthly_lat,
            cpt_t_monthly_lat,
            cpt_q_monthly_lat,
            cpt_rhi_monthly_lat,
            ancestors,
            axis="Y",
        )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
