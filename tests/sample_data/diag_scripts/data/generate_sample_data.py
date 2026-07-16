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
    "emd_tas": {
        "var_name": "emd_tas",
        "standard_name": None,
        "long_name": "EMD of Near-Surface Air Temperature",
        "units": "K",
    },
    "pr": {
        "var_name": "pr",
        "standard_name": "precipitation_flux",
        "long_name": "Precipitation",
        "units": "kg m-2 s-1",
    },
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

    # 0D
    save_cube("pr_amon_0d_0.nc", "pr", [], 6.0)
    save_cube("pr_amon_0d_1.nc", "pr", [], 10.0)
    save_cube("pr_amon_0d_2.nc", "pr", [], 4.0)
    save_cube("tas_amon_0d_0.nc", "tas", [], 3.0)
    save_cube("tas_amon_0d_1.nc", "tas", [], 5.0)
    save_cube("tas_amon_0d_2.nc", "tas", [], 2.0)

    # 1D alt16
    save_cube("tas_amon_1d_alt16_0.nc", "tas", [alt16_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_alt16_1.nc", "tas", [alt16_coord], [1.0, -1.0, 2.0])
    save_cube(
        "tas_amon_1d_alt16_2.nc",
        "tas",
        [alt16_coord],
        [-1.0, -3.0, 1.5],
    )

    # 1D hour
    save_cube("tas_amon_1d_hour_0.nc", "tas", [hour_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_hour_1.nc", "tas", [hour_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_hour_2.nc", "tas", [hour_coord], [-1.0, -3.0, 1.5])

    # 1D lat
    save_cube("tas_amon_1d_lat_0.nc", "tas", [lat_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_lat_1.nc", "tas", [lat_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_lat_2.nc", "tas", [lat_coord], [-1.0, -3.0, 1.5])

    # 1D month_number
    save_cube(
        "emd_tas_amon_1d_month_number_0.nc",
        "emd_tas",
        [month_number_coord],
        [2.0, 3.0, np.nan],
    )
    save_cube(
        "emd_tas_amon_1d_month_number_1.nc",
        "emd_tas",
        [month_number_coord],
        [2.5, 4.5, 4.0],
    )
    save_cube(
        "emd_tas_amon_1d_month_number_2.nc",
        "emd_tas",
        [month_number_coord],
        [2.0, 1.5, 2.5],
    )
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

    # 1D plev
    save_cube("tas_amon_1d_plev_0.nc", "tas", [plev_coord], [2.0, 3.0, 0.0])
    save_cube("tas_amon_1d_plev_1.nc", "tas", [plev_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_plev_2.nc", "tas", [plev_coord], [-1.0, -3.0, 1.5])

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

    # 2D month_number/lat
    save_cube(
        "tas_amon_2d_month_number_lat_0.nc",
        "tas",
        [month_number_coord, lat_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_month_number_lat_1.nc",
        "tas",
        [month_number_coord, lat_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_month_number_lat_2.nc",
        "tas",
        [month_number_coord, lat_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D month_number/lon
    save_cube(
        "tas_amon_2d_month_number_lon_0.nc",
        "tas",
        [month_number_coord, lon_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_month_number_lon_1.nc",
        "tas",
        [month_number_coord, lon_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_month_number_lon_2.nc",
        "tas",
        [month_number_coord, lon_coord],
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

    # 2D time/alt16
    save_cube(
        "tas_amon_2d_time_alt16_0.nc",
        "tas",
        [time_coord, alt16_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_time_alt16_1.nc",
        "tas",
        [time_coord, alt16_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_time_alt16_2.nc",
        "tas",
        [time_coord, alt16_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D time/lat
    save_cube(
        "tas_amon_2d_time_lat_0.nc",
        "tas",
        [time_coord, lat_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_time_lat_1.nc",
        "tas",
        [time_coord, lat_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_time_lat_2.nc",
        "tas",
        [time_coord, lat_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D time/lon
    save_cube(
        "tas_amon_2d_time_lon_0.nc",
        "tas",
        [time_coord, lon_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_time_lon_1.nc",
        "tas",
        [time_coord, lon_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_time_lon_2.nc",
        "tas",
        [time_coord, lon_coord],
        np.arange(9).reshape(3, 3) * -1,
    )

    # 2D time/plev
    save_cube(
        "tas_amon_2d_time_plev_0.nc",
        "tas",
        [time_coord, plev_coord],
        np.arange(9).reshape(3, 3),
    )
    save_cube(
        "tas_amon_2d_time_plev_1.nc",
        "tas",
        [time_coord, plev_coord],
        np.arange(9).reshape(3, 3) * 2 + 2,
    )
    save_cube(
        "tas_amon_2d_time_plev_2.nc",
        "tas",
        [time_coord, plev_coord],
        np.arange(9).reshape(3, 3) * -1,
    )


if __name__ == "__main__":
    main()
