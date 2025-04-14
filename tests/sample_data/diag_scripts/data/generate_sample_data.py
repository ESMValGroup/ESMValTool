#!/usr/bin/env python
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


def data(arr: list) -> np.ndarray:
    """Get masked float32 data."""
    return np.ma.masked_invalid(np.array(arr, dtype=np.float32))


def save_cube(
    filename: str,
    var_name: str,
    coords: list[DimCoord],
    cube_data: list,
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
    time_coord = DimCoord(
        [15.0, 46.0, 75.0],
        bounds=[[0.0, 31.0], [31.0, 60.0], [60.0, 91.0]],
        var_name="time",
        standard_name="time",
        long_name="time",
        units="days since 2000-01-01",
    )
    month_number_coord = DimCoord(
        [1, 2, 4],
        var_name="month_number",
        long_name="month_number",
        units="1",
    )
    hour_coord = DimCoord(
        [1, 2, 4], var_name="hour", long_name="hour", units="1"
    )

    # 1D time
    save_cube("tas_amon_1d_time_0.nc", "tas", [time_coord], [2.0, 3.0, np.nan])
    save_cube("tas_amon_1d_time_1.nc", "tas", [time_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_time_2.nc", "tas", [time_coord], [-1.0, -3.0, 1.5])

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

    # 1D hour
    save_cube("tas_amon_1d_hour_0.nc", "tas", [hour_coord], [2.0, 3.0, np.nan])
    save_cube("tas_amon_1d_hour_1.nc", "tas", [hour_coord], [1.0, -1.0, 2.0])
    save_cube("tas_amon_1d_hour_2.nc", "tas", [hour_coord], [-1.0, -3.0, 1.5])


if __name__ == "__main__":
    main()
