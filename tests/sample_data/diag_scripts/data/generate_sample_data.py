#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import iris
import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube

iris.FUTURE.save_split_attrs = True

ROOT = Path(__file__).parent
VAR_METADATA: dict[str, dict[str, str]] = {
    "tas": {
        "var_name": "tas",
        "standard_name": "air_temperature",
        "long_name": "Near-Surface Air Temperature",
        "units": "K",
    },
}


def data(arr: list | np.ndarray) -> np.ndarray:
    """Get float32 data (all values in [-0.5, 0.5] will be masked)."""
    return np.ma.masked_inside(np.array(arr, dtype=np.float32), -0.5, 0.5)


def save_cube(
    filename: str,
    var_name: str,
    coords: list[DimCoord],
    cube_data: list | np.ndarray,
) -> None:
    """Save cube."""
    dim_coord_spec = [(c, i) for (i, c) in enumerate(coords)]
    kwargs = VAR_METADATA[var_name]
    cube = Cube(data(cube_data), dim_coords_and_dims=dim_coord_spec, **kwargs)
    path = ROOT / filename
    iris.save(cube, path)
    print("Saved", path)


def main() -> None:
    """Generate sample data."""
    # Define coordinates
    alt16_coord = DimCoord(
        [2.0, 4.0, 6.0],
        bounds=[[1.0, 3.0], [3.0, 5.0], [5.0, 7.0]],
        var_name="alt16",
        standard_name="altitude",
        long_name="altitude",
        units="m",
        attributes={"positive": "up"},
    )
    hour_coord = DimCoord(
        [1, 2, 4],
        var_name="hour",
        long_name="hour",
        units="1",
    )
    lat_coord = DimCoord(
        [10.0, 30.0, 50.0],
        bounds=[[0.0, 20.0], [20.0, 40.0], [40.0, 60.0]],
        var_name="lat",
        standard_name="latitude",
        long_name="Latitude",
        units="degrees_north",
    )
    lon_coord = DimCoord(
        [20.0, 60.0, 100.0],
        bounds=[[0.0, 40.0], [40.0, 80.0], [80.0, 120.0]],
        var_name="lon",
        standard_name="longitude",
        long_name="Longitude",
        units="degrees_east",
    )
    month_number_coord = DimCoord(
        [1, 2, 4],
        var_name="month_number",
        long_name="month_number",
        units="1",
    )
    plev_coord = DimCoord(
        [1000.0, 500.0, 200.0],
        var_name="plev",
        standard_name="air_pressure",
        long_name="pressure",
        units="Pa",
        attributes={"positive": "down"},
    )
    time_coord = DimCoord(
        [15.0, 46.0, 75.0],
        bounds=[[0.0, 31.0], [31.0, 60.0], [60.0, 91.0]],
        var_name="time",
        standard_name="time",
        long_name="time",
        units="days since 2000-01-01",
    )

    # 1D hour
    save_cube("tas_amon_1d_hour_0.nc", "tas", [hour_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_hour_1.nc", "tas", [hour_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_hour_2.nc", "tas", [hour_coord], [-1.0, -3.0, 1.5])

    # 1D month_number
    save_cube(
        "tas_amon_1d_month_number_0.nc",
        "tas",
        [month_number_coord],
        [2.0, 3.0, np.nan],
    )
    save_cube(
        "tas_amon_1d_month_number_1.nc",
        "tas",
        [month_number_coord],
        [1.0, -1.0, 2.0],
    )
    save_cube(
        "tas_amon_1d_month_number_2.nc",
        "tas",
        [month_number_coord],
        [-1.0, -3.0, 1.5],
    )

    # 1D time
    save_cube("tas_amon_1d_time_0.nc", "tas", [time_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_time_1.nc", "tas", [time_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_time_2.nc", "tas", [time_coord], [-1.0, -3.0, 1.5])

    # 2D alt16/lat
    save_cube(
        "tas_amon_2d_alt16_lat_0.nc",
        "tas",
        [alt16_coord, lat_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_alt16_lat_1.nc",
        "tas",
        [alt16_coord, lat_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_alt16_lat_2.nc",
        "tas",
        [alt16_coord, lat_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D lat/lon
    save_cube(
        "tas_amon_2d_lat_lon_0.nc",
        "tas",
        [lat_coord, lon_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_lat_lon_1.nc",
        "tas",
        [lat_coord, lon_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_lat_lon_2.nc",
        "tas",
        [lat_coord, lon_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D plev/lat
    save_cube(
        "tas_amon_2d_plev_lat_0.nc",
        "tas",
        [plev_coord, lat_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_plev_lat_1.nc",
        "tas",
        [plev_coord, lat_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_plev_lat_2.nc",
        "tas",
        [plev_coord, lat_coord],
        np.arange(9).reshape(3, 3) * -1,
    )


if __name__ == "__main__":
    main()
