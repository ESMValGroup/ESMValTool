"""Area-weighted remapping of regular lat/lon cubes onto a common 1° grid."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr

GridKind = Literal["era5", "centers"]
_LAT_NAMES = ("lat", "latitude", "nav_lat")
_LON_NAMES = ("lon", "longitude", "nav_lon")
_TIME_NAMES = ("time", "valid_time")


def one_degree_coordinates(kind: GridKind = "era5") -> tuple[np.ndarray, np.ndarray]:
    """Latitude and longitude centres of a global 1° grid, in degrees.

    ``era5`` matches the monthly ERA5 files already saved for this project:
    latitude 90, 89, …, −90 and longitude 0, 1, …, 359. ``centers`` uses the
    middle of each 1° cell: latitude −89.5 … 89.5 and longitude 0.5 … 359.5.
    """
    if kind == "era5":
        lat = np.arange(90.0, -90.0 - 1e-9, -1.0)
        lon = np.arange(0.0, 360.0, 1.0)
    elif kind == "centers":
        lat = np.arange(89.5, -90.0, -1.0)
        lon = np.arange(0.5, 360.0, 1.0)
    else:
        raise ValueError("kind must be 'era5' or 'centers'")
    return lat, lon


def _axis_name(obj: xr.Dataset | xr.DataArray, candidates: tuple[str, ...], role: str) -> str:
    for name in candidates:
        if name not in obj.dims and name not in getattr(obj, "coords", {}):
            continue
        values = np.asarray(obj[name].values)
        if values.ndim != 1:
            raise ValueError(f"{role} axis {name!r} is not a regular 1D coordinate")
        return name
    raise ValueError(f"No regular {role} axis among {candidates}. Dims={list(obj.dims)}")


def _ascending(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        raise ValueError("An axis needs at least two coordinates")
    order = np.argsort(values)
    sorted_values = values[order]
    if np.any(np.diff(sorted_values) <= 0):
        raise ValueError("Axis coordinates must be strictly monotonic")
    flipped = not np.array_equal(order, np.arange(values.size))
    return sorted_values, order, flipped


def _edges(centers: np.ndarray, *, limit: tuple[float, float] | None) -> np.ndarray:
    mid = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - (mid[0] - centers[0])
    last = centers[-1] + (centers[-1] - mid[-1])
    edges = np.concatenate([[first], mid, [last]])
    if limit is not None:
        edges = np.clip(edges, limit[0], limit[1])
    if np.any(np.diff(edges) < 0):
        raise ValueError("Axis edges are not increasing")
    return edges


def _overlap(src_edges: np.ndarray, tgt_edges: np.ndarray, *, latitude: bool) -> np.ndarray:
    src0, src1 = src_edges[:-1], src_edges[1:]
    tgt0, tgt1 = tgt_edges[:-1], tgt_edges[1:]

    def pair(shift: float) -> np.ndarray:
        left = np.maximum(src0[:, None] + shift, tgt0[None, :])
        right = np.minimum(src1[:, None] + shift, tgt1[None, :])
        overlap = np.clip(right - left, 0.0, None)
        if not latitude:
            return np.deg2rad(overlap)
        lo = np.deg2rad(left)
        hi = np.deg2rad(right)
        return np.where(overlap > 0, np.sin(hi) - np.sin(lo), 0.0)

    weight = pair(0.0)
    if not latitude:
        weight = weight + pair(-360.0) + pair(360.0)
    return weight


def _remap(values: np.ndarray, lat_weight: np.ndarray, lon_weight: np.ndarray) -> np.ndarray:
    squeeze = values.ndim == 2
    if squeeze:
        values = values[None, ...]
    out = np.empty((values.shape[0], lat_weight.shape[1], lon_weight.shape[1]), dtype=np.float64)
    for start in range(0, values.shape[0], 24):
        block = np.asarray(values[start : start + 24], dtype=np.float64)
        mask = np.isfinite(block)
        filled = np.where(mask, block, 0.0)
        numerator = np.matmul(lat_weight.T, np.matmul(filled, lon_weight))
        denominator = np.matmul(lat_weight.T, np.matmul(mask.astype(np.float64), lon_weight))
        with np.errstate(invalid="ignore", divide="ignore"):
            remapped = numerator / denominator
        out[start : start + 24] = np.where(denominator > 0, remapped, np.nan)
    return out[0] if squeeze else out


def _prepare_field(data: xr.DataArray, lat_name: str, lon_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    lat, lat_order, _ = _ascending(np.asarray(data[lat_name].values, dtype=float))
    lon_raw = np.mod(np.asarray(data[lon_name].values, dtype=float), 360.0)
    lon, lon_order, _ = _ascending(lon_raw)
    field = data.transpose(..., lat_name, lon_name).values
    field = np.take(field, lat_order, axis=-2)
    field = np.take(field, lon_order, axis=-1)
    other = [name for name in data.dims if name not in {lat_name, lon_name}]
    if len(other) > 1:
        raise ValueError(f"{data.name or 'field'} has extra dimensions {other}")
    return field, lat, lon, other


def regrid_to_1deg(
    data: xr.Dataset | xr.DataArray,
    *,
    grid: GridKind = "era5",
) -> xr.Dataset | xr.DataArray:
    """Remap a regular lat/lon cube onto a global 1° grid.

    Overlapping source and target cells are combined with a first-order
    area weight (``sin`` latitude by longitude span). The input is not modified.
    Coordinates in the result are named ``lat`` and ``lon``. A time axis is
    kept and named ``time``.
    """
    return_array = isinstance(data, xr.DataArray)
    dataset = data.to_dataset() if return_array else data
    lat_name = _axis_name(dataset, _LAT_NAMES, "latitude")
    lon_name = _axis_name(dataset, _LON_NAMES, "longitude")
    target_lat, target_lon = one_degree_coordinates(grid)
    lat_weight = _overlap(
        _edges(np.sort(np.asarray(dataset[lat_name].values, dtype=float)), limit=(-90.0, 90.0)),
        _edges(np.sort(target_lat), limit=(-90.0, 90.0)),
        latitude=True,
    )
    lon_weight = _overlap(
        _edges(np.sort(np.mod(np.asarray(dataset[lon_name].values, dtype=float), 360.0)), limit=None),
        _edges(np.sort(np.mod(target_lon, 360.0)), limit=None),
        latitude=False,
    )
    time_name = next((name for name in _TIME_NAMES if name in dataset.coords), None)

    data_vars = {}
    for name, variable in dataset.data_vars.items():
        if lat_name not in variable.dims or lon_name not in variable.dims:
            continue
        field, _, _, other = _prepare_field(variable, lat_name, lon_name)
        remapped = _restore_axis(_remap(field, lat_weight, lon_weight), target_lat, target_lon)
        dims = ("time", "lat", "lon") if other else ("lat", "lon")
        attrs = dict(variable.attrs)
        attrs["regrid"] = "area-weighted 1 degree"
        attrs["regrid_grid"] = grid
        data_vars[name] = (dims, remapped, attrs)
    if not data_vars:
        raise ValueError("No regular lat/lon data variables to regrid")

    coords: dict = {
        "lat": ("lat", target_lat, {"units": "degrees_north", "long_name": "latitude"}),
        "lon": ("lon", target_lon, {"units": "degrees_east", "long_name": "longitude"}),
    }
    needs_time = any(len(spec[0]) == 3 for spec in data_vars.values())
    if needs_time:
        if time_name is None:
            raise ValueError("Spatial fields have a time dimension but no time coordinate")
        coords["time"] = ("time", np.asarray(dataset[time_name].values), dict(dataset[time_name].attrs))
    result = xr.Dataset(data_vars, coords=coords, attrs=dict(dataset.attrs))
    result.attrs["regrid"] = "area-weighted 1 degree"
    result.attrs["regrid_grid"] = grid
    if return_array:
        array_name = data.name or next(iter(result.data_vars))
        return result[array_name]
    return result


def _lonlat_to_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat_r = np.deg2rad(np.asarray(lat, dtype=np.float64))
    lon_r = np.deg2rad(np.mod(np.asarray(lon, dtype=np.float64), 360.0))
    cos_lat = np.cos(lat_r)
    return np.stack((cos_lat * np.cos(lon_r), cos_lat * np.sin(lon_r), np.sin(lat_r)), axis=-1)


def _curvilinear_latlon(dataset: xr.Dataset) -> tuple[str, str]:
    lat_name = next((name for name in _LAT_NAMES if name in dataset.coords or name in dataset.variables), None)
    lon_name = next((name for name in _LON_NAMES if name in dataset.coords or name in dataset.variables), None)
    if lat_name is None or lon_name is None:
        raise ValueError(f"Need 2D latitude/longitude coordinates, got {list(dataset.coords)}")
    lat = dataset[lat_name]
    lon = dataset[lon_name]
    if lat.ndim != 2 or lon.ndim != 2:
        raise ValueError(
            f"{lat_name}/{lon_name} must be 2-dimensional for nearest curvilinear regrid "
            f"(got {lat.ndim}D / {lon.ndim}D)"
        )
    return lat_name, lon_name


def regrid_curvilinear_nearest_to_1deg(
    data: xr.Dataset | xr.DataArray,
    *,
    grid: GridKind = "era5",
    max_distance_deg: float = 1.5,
) -> xr.Dataset | xr.DataArray:
    """Nearest-neighbour remap of a 2D-lat/lon (curvilinear) cube onto 1°.

    Used for native CMIP ocean grids such as ACCESS ``(j, i)``. Target cells
    farther than ``max_distance_deg`` from any source cell stay missing.
    """
    from scipy.spatial import cKDTree

    return_array = isinstance(data, xr.DataArray)
    dataset = data.to_dataset(name=data.name or "z") if return_array else data
    lat_name, lon_name = _curvilinear_latlon(dataset)
    src_lat = np.asarray(dataset[lat_name].values, dtype=np.float64)
    src_lon = np.asarray(dataset[lon_name].values, dtype=np.float64)
    src_xyz = _lonlat_to_xyz(src_lat.ravel(), src_lon.ravel())
    finite_src = np.isfinite(src_xyz).all(axis=1)
    tree = cKDTree(src_xyz[finite_src])
    src_index = np.flatnonzero(finite_src)

    target_lat, target_lon = one_degree_coordinates(grid)
    tgt_lat, tgt_lon = np.meshgrid(target_lat, target_lon, indexing="ij")
    tgt_xyz = _lonlat_to_xyz(tgt_lat.ravel(), tgt_lon.ravel())
    distance, neighbours = tree.query(tgt_xyz, k=1)
    max_chord = 2.0 * np.sin(np.deg2rad(max_distance_deg) / 2.0)
    take = src_index[neighbours]
    valid = distance <= max_chord

    time_name = next((name for name in _TIME_NAMES if name in dataset.coords), None)
    data_vars = {}
    spatial = tuple(dataset[lat_name].dims)
    for name, variable in dataset.data_vars.items():
        if any(dim not in variable.dims for dim in spatial):
            continue
        values = variable.transpose(..., *spatial).values
        flat = values.reshape(*values.shape[: values.ndim - 2], -1)
        remapped = flat[..., take]
        remapped[..., ~valid] = np.nan
        remapped = remapped.reshape(*values.shape[: values.ndim - 2], target_lat.size, target_lon.size)
        other = [dim for dim in variable.dims if dim not in spatial]
        if len(other) > 1:
            raise ValueError(f"{name} has extra dimensions {other}")
        dims = ("time", "lat", "lon") if other else ("lat", "lon")
        attrs = dict(variable.attrs)
        attrs["regrid"] = "nearest curvilinear to 1 degree"
        attrs["regrid_grid"] = grid
        data_vars[name] = (dims, remapped, attrs)
    if not data_vars:
        raise ValueError("No curvilinear data variables to regrid")

    coords: dict = {
        "lat": ("lat", target_lat, {"units": "degrees_north", "long_name": "latitude"}),
        "lon": ("lon", target_lon, {"units": "degrees_east", "long_name": "longitude"}),
    }
    if any(len(spec[0]) == 3 for spec in data_vars.values()):
        if time_name is None:
            raise ValueError("Spatial fields have a time dimension but no time coordinate")
        coords["time"] = ("time", np.asarray(dataset[time_name].values), dict(dataset[time_name].attrs))
    result = xr.Dataset(data_vars, coords=coords, attrs=dict(dataset.attrs))
    result.attrs["regrid"] = "nearest curvilinear to 1 degree"
    result.attrs["regrid_grid"] = grid
    if return_array:
        array_name = data.name or next(iter(result.data_vars))
        return result[array_name]
    return result


def _restore_axis(remapped: np.ndarray, target_lat: np.ndarray, target_lon: np.ndarray) -> np.ndarray:
    """Move ascending remap axes back to the requested lat/lon order."""
    lat_ascending_index = np.argsort(target_lat)
    lon_ascending_index = np.argsort(np.mod(target_lon, 360.0))
    public = np.empty_like(remapped)
    public[..., lat_ascending_index, :] = remapped
    ordered = np.empty_like(public)
    ordered[..., lon_ascending_index] = public
    return ordered


def regrid_file_to_1deg(
    source: str | Path,
    destination: str | Path,
    *,
    grid: GridKind = "era5",
    overwrite: bool = False,
) -> Path:
    """Write a 1° copy of ``source`` to ``destination``.

    The source path is never replaced. An existing destination is left in
    place unless ``overwrite`` is true.
    """
    source_path = Path(source).expanduser().resolve()
    destination_path = Path(destination).expanduser().resolve()
    if source_path == destination_path:
        raise ValueError(f"Refusing to overwrite the original file: {source_path}")
    if destination_path.exists() and not overwrite:
        raise FileExistsError(f"Destination already exists: {destination_path}")
    with xr.open_dataset(source_path) as dataset:
        remapped = regrid_to_1deg(dataset, grid=grid)
    remapped.attrs["regrid_source"] = str(source_path)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination_path.with_suffix(destination_path.suffix + ".part")
    encoding = {
        name: {"zlib": True, "complevel": 4, "dtype": "float32"}
        for name, variable in remapped.data_vars.items()
        if np.issubdtype(variable.dtype, np.floating)
    }
    remapped.to_netcdf(temporary, encoding=encoding)
    temporary.replace(destination_path)
    return destination_path
