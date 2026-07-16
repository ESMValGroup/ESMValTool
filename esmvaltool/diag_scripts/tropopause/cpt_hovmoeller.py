"""Create tropical CPT Hovmoeller plots with contour overlays.

This diagnostic identifies the tropical cold-point tropopause (CPT) by
searching for the minimum temperature between configurable pressure bounds and
creates a longitude-time Hovmoeller plot (20S-20N mean) with:

- filled colors: CPT pressure
- line contours: CPT temperature, CPT specific humidity, CPT RH over ice

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
- plot_filename: Prefix for output file names
  (default: ``tropical_cpt_hovmoeller``).
- dpi: Figure DPI (default: 300).
- icon: Optional block for raw ICON input (same pattern as the
  wheeler_kiladis diagnostic): file patterns, grid file, variable names,
  and time metadata.
"""

from __future__ import annotations

import glob
import logging
import re
from pathlib import Path
import calendar

import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from iris.coords import AuxCoord, DimCoord
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
        "authors": ["weigel_katja", "schlund_manuel"],
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
    target_z_dim = target_cube.coord_dims(_get_pressure_coord(target_cube).name())[0]
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
        9.550426
        - (5723.265 / t)
        + (3.53068 * np.log(t))
        - (0.00728332 * t)
    )
    es_ice = np.exp(log_es_ice)
    rhi = 100.0 * vapor_pressure / es_ice
    rhi[~np.isfinite(rhi)] = np.nan
    return rhi


def _ensure_datetime_time_axis(
    data: xr.DataArray,
    *,
    start_date: str,
    time_step_hours: int | None = None,
) -> xr.DataArray:
    if np.issubdtype(data["time"].dtype, np.datetime64):
        return data

    raw_time = np.asarray(data["time"].values)
    if np.issubdtype(raw_time.dtype, np.number):
        rounded = np.rint(raw_time).astype(np.int64)
        if np.all((rounded >= 10000101) & (rounded <= 99991231)):
            parsed = pd.to_datetime(
                rounded.astype(str),
                format="%Y%m%d",
                errors="coerce",
            )
            if np.all(~pd.isna(parsed)):
                return data.assign_coords(time=parsed.to_numpy(dtype="datetime64[ns]"))
    elif raw_time.dtype.kind in {"U", "S", "O"}:
        parsed = pd.to_datetime(raw_time, errors="coerce")
        if np.all(~pd.isna(parsed)):
            return data.assign_coords(time=parsed.to_numpy(dtype="datetime64[ns]"))

    try:
        tmp_name = data.name if data.name is not None else "var"
        decoded = xr.decode_cf(data.to_dataset(name=tmp_name))
        if np.issubdtype(decoded["time"].dtype, np.datetime64):
            return decoded[tmp_name]
    except Exception:
        pass

    nt = data.sizes["time"]
    if time_step_hours is None:
        new_time = pd.date_range(start=start_date, periods=nt, freq="1D")
    else:
        new_time = pd.date_range(
            start=start_date,
            periods=nt,
            freq=pd.to_timedelta(time_step_hours, unit="h"),
        )
    return data.assign_coords(time=new_time)


def _open_icon_dataarray(
    *,
    file_pattern: str,
    variable: str,
    start_date: str,
    time_step_hours: int,
    chunks_time: int,
) -> xr.DataArray:
    files = sorted(glob.glob(file_pattern))
    if not files:
        raise FileNotFoundError(f"No ICON files found for pattern: {file_pattern}")

    dataset = xr.open_mfdataset(
        files,
        combine="by_coords",
        chunks={"time": chunks_time},
        data_vars="minimal",
        coords="minimal",
        compat="override",
    )
    if variable not in dataset:
        available = list(dataset.data_vars)
        raise KeyError(
            f"{variable!r} not found in ICON files. Available variables: {available}",
        )
    data = dataset[variable]
    if "time" not in data.dims:
        raise ValueError(f"ICON variable {variable!r} has no 'time' dimension.")
    return _ensure_datetime_time_axis(
        data,
        start_date=start_date,
        time_step_hours=time_step_hours,
    )


def _icon_grid_lat_lon(
    data: xr.DataArray,
    grid_file: str,
) -> tuple[xr.DataArray, xr.DataArray, str]:
    grid = xr.open_dataset(grid_file)
    non_time_dims = [dim for dim in data.dims if dim != "time"]
    if not non_time_dims:
        raise ValueError(f"Could not infer ICON dimensions from {data.dims}")

    horiz_dim = None
    grid_dim = None
    for dim in non_time_dims:
        size = data.sizes[dim]
        for gdim in grid.dims:
            if grid.sizes[gdim] == size:
                horiz_dim = dim
                grid_dim = gdim
                break
        if horiz_dim is not None:
            break

    if horiz_dim is None:
        raise ValueError(
            "Could not match an ICON horizontal dimension "
            f"for data dims {data.dims} and grid dims {dict(grid.sizes)}",
        )

    if "clat" not in grid or "clon" not in grid:
        raise KeyError("ICON grid file must contain 'clat' and 'clon'.")

    clat = grid["clat"]
    clon = grid["clon"]
    if grid_dim != horiz_dim:
        clat = clat.rename({grid_dim: horiz_dim})
        clon = clon.rename({grid_dim: horiz_dim})

    lat_values = clat.values
    lon_values = clon.values
    if np.nanmax(np.abs(lat_values)) <= np.pi / 2.0 + 0.1:
        lat_values = np.rad2deg(lat_values)
    if np.nanmax(np.abs(lon_values)) <= 2.0 * np.pi + 0.1:
        lon_values = np.rad2deg(lon_values)
    lon_values = np.mod(lon_values, 360.0)

    lat_deg = xr.DataArray(
        lat_values,
        dims=(horiz_dim,),
        coords={horiz_dim: data[horiz_dim]},
        name="lat",
    )
    lon_deg = xr.DataArray(
        lon_values,
        dims=(horiz_dim,),
        coords={horiz_dim: data[horiz_dim]},
        name="lon",
    )
    return lat_deg, lon_deg, horiz_dim


def _get_vertical_dim(data: xr.DataArray, horiz_dim: str) -> str:
    vertical_candidates = [
        "plev",
        "lev",
        "level",
        "height",
        "height_2",
        "altitude",
        "model_level",
    ]
    for name in vertical_candidates:
        if name in data.dims and name not in ("time", horiz_dim):
            return name
    dims = [dim for dim in data.dims if dim not in ("time", horiz_dim)]
    if len(dims) != 1:
        raise ValueError(
            "Could not infer ICON vertical dimension. "
            f"Expected one non-time/non-horizontal dim, got {dims}",
        )
    return dims[0]


def _pressure_hpa_from_icon(
    ta_data: xr.DataArray,
    pressure_data: xr.DataArray | None,
    vertical_dim: str,
) -> xr.DataArray:
    if pressure_data is not None:
        p = pressure_data
    elif vertical_dim in ta_data.coords:
        p = ta_data[vertical_dim]
    else:
        raise ValueError(
            "No pressure information found. Set icon.pressure_variable "
            "or provide a pressure coordinate on the temperature field.",
        )

    units = str(p.attrs.get("units", "")).lower()
    if "hpa" in units or "mbar" in units:
        return p.astype(float)
    if "pa" in units:
        return p.astype(float) / 100.0

    p_values = p.astype(float)
    if float(np.nanmax(p_values)) > 2000.0:
        return p_values / 100.0
    return p_values


def _select_at_cpt_icon(
    ta_data: xr.DataArray,
    target_data: xr.DataArray,
    pressure_hpa: xr.DataArray,
    vertical_dim: str,
    p_bounds_hpa: tuple[float, float],
) -> tuple[xr.DataArray, xr.DataArray]:
    pmin, pmax = min(p_bounds_hpa), max(p_bounds_hpa)
    p_level = pressure_hpa.broadcast_like(ta_data)
    valid_level = (p_level >= pmin) & (p_level <= pmax)
    ta_valid = ta_data.where(valid_level)
    # xarray cannot use a dask-chunked DataArray as vectorized indexer.
    # Compute indexer/mask eagerly, while keeping target arrays lazy.
    idx = (
        ta_valid.fillna(np.inf)
        .argmin(dim=vertical_dim)
        .astype(np.int64)
        .compute()
    )
    has_valid = valid_level.any(dim=vertical_dim).compute()

    selected = target_data.isel({vertical_dim: idx})
    selected = selected.where(has_valid)
    cpt_pressure_hpa = p_level.isel({vertical_dim: idx}).where(has_valid)
    return selected, cpt_pressure_hpa


def _icon_tropical_hovmoeller(
    data: xr.DataArray,
    lat_deg: xr.DataArray,
    lon_deg: xr.DataArray,
    horiz_dim: str,
    lat_bounds: tuple[float, float],
    lon_resolution: float,
) -> xr.DataArray:
    lat_min, lat_max = min(lat_bounds), max(lat_bounds)
    tropical_mask = (lat_deg >= lat_min) & (lat_deg <= lat_max)
    tropical_data = data.where(tropical_mask, drop=True)
    tropical_lat = lat_deg.where(tropical_mask, drop=True)
    tropical_lon = lon_deg.where(tropical_mask, drop=True)

    nlon = int(round(360.0 / lon_resolution))
    lon_centers = np.arange(nlon) * lon_resolution + 0.5 * lon_resolution
    lon_bin_values = np.floor(tropical_lon.values / lon_resolution).astype(int)
    lon_bin_values = np.clip(lon_bin_values, 0, nlon - 1)

    lon_bin = xr.DataArray(
        lon_bin_values,
        dims=(horiz_dim,),
        coords={horiz_dim: tropical_data[horiz_dim]},
        name="lon_bin",
    )
    weights = xr.DataArray(
        np.cos(np.deg2rad(tropical_lat.values)),
        dims=(horiz_dim,),
        coords={horiz_dim: tropical_data[horiz_dim]},
        name="weights",
    )

    numerator = (tropical_data * weights).groupby(lon_bin).sum(dim=horiz_dim)
    denominator = weights.groupby(lon_bin).sum(dim=horiz_dim)
    hov = numerator / denominator
    hov = hov.reindex(lon_bin=np.arange(nlon))
    hov = hov.assign_coords(longitude=("lon_bin", lon_centers))
    hov = hov.swap_dims({"lon_bin": "longitude"}).drop_vars("lon_bin")
    hov = hov.transpose("time", "longitude")
    return hov


def _time_lon_cube_from_xarray(
    data: xr.DataArray,
    *,
    var_name: str,
    long_name: str,
    units: str,
) -> Cube:
    values = np.asarray(data.values, dtype=float)
    lon = np.asarray(data["longitude"].values, dtype=float)
    times = pd.to_datetime(data["time"].values)
    time_coord = DimCoord(
        np.arange(values.shape[0], dtype=float),
        standard_name="time",
        units="days since 1970-01-01",
    )
    lon_coord = DimCoord(
        lon,
        standard_name="longitude",
        units="degrees_east",
    )
    cube = Cube(
        values,
        var_name=var_name,
        long_name=long_name,
        units=units,
        dim_coords_and_dims=[(time_coord, 0), (lon_coord, 1)],
    )
    cube.add_aux_coord(AuxCoord(times.year, long_name="year"), 0)
    cube.add_aux_coord(AuxCoord(times.month, long_name="month_number"), 0)
    return cube


def _process_icon_dataset(cfg: dict) -> None:
    icon_cfg = cfg["icon"]
    label = cfg.get("dataset_label", icon_cfg.get("label", "ICON"))
    p_bounds_hpa = tuple(cfg["pressure_bounds_hpa"])
    lat_bounds = tuple(cfg["tropical_lat_bounds"])
    lon_resolution = float(cfg.get("lon_resolution", 2.5))
    start_date = icon_cfg.get("start_date", "1979-01-01")
    time_step_hours = int(icon_cfg.get("time_step_hours", 6))
    chunks_time = int(icon_cfg.get("chunks_time", 240))

    ta_da = _open_icon_dataarray(
        file_pattern=icon_cfg["ta_file_pattern"],
        variable=icon_cfg.get("ta_variable", "ta"),
        start_date=start_date,
        time_step_hours=time_step_hours,
        chunks_time=chunks_time,
    )
    hus_da = _open_icon_dataarray(
        file_pattern=icon_cfg["hus_file_pattern"],
        variable=icon_cfg.get("hus_variable", "hus"),
        start_date=start_date,
        time_step_hours=time_step_hours,
        chunks_time=chunks_time,
    )
    zg_da = None
    if icon_cfg.get("zg_file_pattern") and icon_cfg.get("zg_variable"):
        zg_da = _open_icon_dataarray(
            file_pattern=icon_cfg["zg_file_pattern"],
            variable=icon_cfg["zg_variable"],
            start_date=start_date,
            time_step_hours=time_step_hours,
            chunks_time=chunks_time,
        )

    pressure_da = None
    if icon_cfg.get("pressure_file_pattern") and icon_cfg.get("pressure_variable"):
        pressure_da = _open_icon_dataarray(
            file_pattern=icon_cfg["pressure_file_pattern"],
            variable=icon_cfg["pressure_variable"],
            start_date=start_date,
            time_step_hours=time_step_hours,
            chunks_time=chunks_time,
        )

    ta_da, hus_da = xr.align(ta_da, hus_da, join="inner")
    if zg_da is not None:
        ta_da, zg_da = xr.align(ta_da, zg_da, join="inner")
    if pressure_da is not None:
        ta_da, pressure_da = xr.align(ta_da, pressure_da, join="inner")

    lat_deg, lon_deg, horiz_dim = _icon_grid_lat_lon(ta_da, icon_cfg["grid_file"])
    vertical_dim = _get_vertical_dim(ta_da, horiz_dim)
    pressure_hpa = _pressure_hpa_from_icon(ta_da, pressure_da, vertical_dim)

    cpt_t, cpt_p = _select_at_cpt_icon(
        ta_da,
        ta_da,
        pressure_hpa,
        vertical_dim,
        p_bounds_hpa,
    )
    cpt_q, _ = _select_at_cpt_icon(
        ta_da,
        hus_da,
        pressure_hpa,
        vertical_dim,
        p_bounds_hpa,
    )
    cpt_rhi = xr.DataArray(
        _compute_rhi_percent(cpt_q.values, cpt_t.values, cpt_p.values),
        dims=cpt_t.dims,
        coords=cpt_t.coords,
        name="cpt_rhi",
    )

    cpt_z = None
    if zg_da is not None:
        cpt_z, _ = _select_at_cpt_icon(
            ta_da,
            zg_da,
            pressure_hpa,
            vertical_dim,
            p_bounds_hpa,
        )

    monthly = {
        "pressure": cpt_p.resample(time="1MS").mean(),
        "temperature": cpt_t.resample(time="1MS").mean(),
        "specific_humidity": cpt_q.resample(time="1MS").mean(),
        "rhi": cpt_rhi.resample(time="1MS").mean(),
    }
    if cpt_z is not None:
        monthly["height"] = cpt_z.resample(time="1MS").mean()

    hov = {
        name: _icon_tropical_hovmoeller(
            field,
            lat_deg,
            lon_deg,
            horiz_dim,
            lat_bounds,
            lon_resolution,
        )
        for name, field in monthly.items()
    }

    cpt_p_monthly = _time_lon_cube_from_xarray(
        hov["pressure"],
        var_name="cpt_p",
        long_name="Cold point tropopause pressure",
        units="hPa",
    )
    cpt_t_monthly = _time_lon_cube_from_xarray(
        hov["temperature"],
        var_name="cpt_ta",
        long_name="Cold point tropopause temperature",
        units=str(ta_da.attrs.get("units", "K")),
    )
    cpt_q_monthly = _time_lon_cube_from_xarray(
        hov["specific_humidity"],
        var_name="cpt_hus",
        long_name="Cold point tropopause specific humidity",
        units=str(hus_da.attrs.get("units", "1")),
    )
    cpt_rhi_monthly = _time_lon_cube_from_xarray(
        hov["rhi"],
        var_name="cpt_rhi",
        long_name="Cold point tropopause relative humidity over ice",
        units="%",
    )
    cpt_z_monthly = None
    if "height" in hov:
        cpt_z_monthly = _time_lon_cube_from_xarray(
            hov["height"],
            var_name="cpt_zg",
            long_name="Cold point tropopause geopotential height",
            units=str(zg_da.attrs.get("units", "m")),
        )

    ancestors = [icon_cfg["ta_file_pattern"], icon_cfg["hus_file_pattern"]]
    if icon_cfg.get("zg_file_pattern"):
        ancestors.append(icon_cfg["zg_file_pattern"])
    if icon_cfg.get("pressure_file_pattern"):
        ancestors.append(icon_cfg["pressure_file_pattern"])
    ancestors.append(icon_cfg["grid_file"])

    base = f"{cfg['plot_filename']}_{_safe_name(label)}"
    provenance = _get_provenance(
        "Monthly tropical CPT diagnostics at longitude-time resolution.",
        ancestors,
    )
    save_data(f"{base}_pressure", provenance, cfg, cpt_p_monthly)
    save_data(f"{base}_temperature", provenance, cfg, cpt_t_monthly)
    save_data(f"{base}_specific_humidity", provenance, cfg, cpt_q_monthly)
    save_data(f"{base}_rhi", provenance, cfg, cpt_rhi_monthly)
    if cpt_z_monthly is not None:
        save_data(f"{base}_height", provenance, cfg, cpt_z_monthly)

    _plot_one_dataset(
        cfg,
        label,
        cpt_p_monthly,
        cpt_t_monthly,
        cpt_q_monthly,
        cpt_rhi_monthly,
        ancestors,
    )


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


def _as_time_lon(cube: Cube) -> Cube:
    time_name = cube.coord(axis="T").name()
    lon_name = cube.coord(axis="X").name()
    time_dim = cube.coord_dims(time_name)[0]
    lon_dim = cube.coord_dims(lon_name)[0]
    if (time_dim, lon_dim) == (0, 1):
        return cube
    out_cube = cube.copy()
    out_cube.transpose((time_dim, lon_dim))
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


def _plot_one_dataset(
    cfg: dict,
    dataset_label: str,
    cpt_p_cube: Cube,
    cpt_t_cube: Cube,
    cpt_q_cube: Cube,
    cpt_rhi_cube: Cube,
    ancestors: list[str],
) -> None:
    cpt_p = _as_time_lon(cpt_p_cube)
    cpt_t = _as_time_lon(cpt_t_cube)
    cpt_q = _as_time_lon(cpt_q_cube)
    cpt_rhi = _as_time_lon(cpt_rhi_cube)

    lon = cpt_p.coord(axis="X").points
    y = np.arange(cpt_p.shape[0])
    p_data = np.asarray(cpt_p.core_data())

    year = cpt_p.coord("year").points.astype(int)
    month = cpt_p.coord("month_number").points.astype(int)

    y_ticks = np.where((month == 1) | (month == 7))[0]
    if y_ticks.size == 0:
        y_ticks = np.arange(0, len(y), 8)

    month_abbrs = [calendar.month_abbr[m].upper() for m in month]
    y_tick_labels = [f"{month_abbrs[i]}{year[i]}" for i in y_ticks]

    panels = [
        ("CPT temperature contours [K]", cpt_t, "black"),
        ("CPT specific humidity contours [x10\u2076 kg kg-1]", cpt_q, "darkgreen"),
        ("CPT RH over ice contours [%]", cpt_rhi, "darkred"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 12), constrained_layout=True)
    p_levels = _finite_levels(p_data, 21)
    contour_sets = []
    label_sets = []
    for ax, (title, contour_cube, color) in zip(axes, panels, strict=True):
        contour_data = np.asarray(contour_cube.core_data())
        contour_levels = _finite_levels(contour_data, 8)
        cf = ax.contourf(lon, y, p_data, levels=p_levels, cmap="Reds_r", extend="both")
        cs = ax.contour(
            lon,
            y,
            contour_data,
            levels=contour_levels,
            colors='k',
            linewidths=0.8,
        )
        contour_sets.append(cf)
        label_sets.append(cs)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
        ax.set_title(title, fontsize=16)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

    # axes[0] (temperature) and axes[2] (RH): plain whole-number labels
    for ax, cs in zip((axes[0], axes[2]), (label_sets[0], label_sets[2]), strict=True):
        ax.clabel(cs, inline=True, fontsize=12, fmt="%.0f")

    # axes[1] (specific humidity): scaled by x10^6, 1 decimal place
    axes[1].clabel(
        label_sets[1],
        inline=True,
        fontsize=12,
        fmt=lambda x: f"{x * 1e6:.2f}",
    )

    axes[0].set_ylabel("Time", fontsize=16)
    axes[1].set_xlabel(f"{cpt_p.coord(axis='X').name()} [{cpt_p.coord(axis='X').units}]", fontsize=16)
    fig.suptitle(
        f"Tropical CPT Hovmoeller ({dataset_label})", fontsize=18
    )
    cbar = fig.colorbar(contour_sets[0], ax=axes, orientation="vertical", shrink=0.95)
    cbar.set_label("CPT pressure [hPa]", fontsize=16)
    cbar.set_ticks([80, 90, 100, 110, 120])
    cbar.ax.tick_params(labelsize=12)

    basename = f"{cfg.get('plot_filename', 'tropical_cpt_hovmoeller')}_{_safe_name(dataset_label)}"
    caption = (
        "Tropical (20S-20N) CPT longitude-time Hovmoeller: colors show CPT pressure; "
        "contours show CPT temperature, CPT specific humidity, and CPT RH over ice."
    )
    provenance = _get_provenance(caption, ancestors)
    save_figure(basename, provenance, cfg, figure=fig, dpi=cfg.get("dpi", 300))


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

    if cfg.get("icon"):
        logger.info("Processing raw ICON input from script configuration")
        _process_icon_dataset(cfg)
        return

    input_data = list(cfg["input_data"].values())
    if not input_data:
        raise ValueError(
            "No input_data provided and no icon block configured. "
            "Provide preprocessed ESMValTool inputs or script.icon settings.",
        )
    grouped = group_metadata(input_data, cfg["group_by"], sort="short_name")

    for group_key, group_items in grouped.items():
        label = _dataset_label(group_key, group_items)
        logger.info("Processing %s", label)
        ta_attrs = _get_variable_attrs(group_items, cfg["temperature_var"])
        hus_attrs = _get_variable_attrs(group_items, cfg["humidity_var"])
        zg_attrs = None
        try:
            zg_attrs = _get_variable_attrs(group_items, cfg["geopotential_height_var"])
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

        cpt_p_monthly = _monthly_tropical_mean(cpt_p, lat_bounds)
        cpt_t_monthly = _monthly_tropical_mean(cpt_t, lat_bounds)
        cpt_q_monthly = _monthly_tropical_mean(cpt_q, lat_bounds)
        cpt_rhi_monthly = _monthly_tropical_mean(cpt_rhi, lat_bounds)

        ancestors = [ta_attrs["filename"], hus_attrs["filename"]]
        if zg_attrs is not None:
            ancestors.append(zg_attrs["filename"])

        base = f"{cfg['plot_filename']}_{_safe_name(label)}"
        provenance = _get_provenance(
            "Monthly tropical CPT diagnostics at longitude-time resolution.",
            ancestors,
        )
        save_data(f"{base}_pressure", provenance, cfg, cpt_p_monthly)
        save_data(f"{base}_temperature", provenance, cfg, cpt_t_monthly)
        save_data(f"{base}_specific_humidity", provenance, cfg, cpt_q_monthly)
        save_data(f"{base}_rhi", provenance, cfg, cpt_rhi_monthly)
        if cpt_z_monthly is not None:
            save_data(f"{base}_height", provenance, cfg, cpt_z_monthly)

        _plot_one_dataset(
            cfg,
            label,
            cpt_p_monthly,
            cpt_t_monthly,
            cpt_q_monthly,
            cpt_rhi_monthly,
            ancestors,
        )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
