"""Tests for Russell et al. (2018) diagnostic calculations."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import iris
import iris.coords
import iris.cube
import numpy as np

from esmvaltool.diag_scripts.russell18jgr.russell_common import (
    add_latitude_coord,
    southern_ocean_flux_sum,
)
from esmvaltool.diag_scripts.russell18jgr.russell_fig6_shared import (
    prepare_section,
)


def _volume_cube(data, lat, lon, *, monthly=False):
    """Build a fixed or monthly cell-volume cube for transport tests."""
    depth = iris.coords.DimCoord(
        [5.0, 15.0],
        standard_name="depth",
        units="m",
        attributes={"positive": "down"},
    )
    latitude = iris.coords.DimCoord(
        lat, standard_name="latitude", units="degrees"
    )
    longitude = iris.coords.DimCoord(
        lon, standard_name="longitude", units="degrees"
    )
    coords = [(depth, 0), (latitude, 1), (longitude, 2)]
    if monthly:
        time = iris.coords.DimCoord(
            [15.0, 45.0],
            standard_name="time",
            units="days since 2000-01-01",
            bounds=[[0.0, 30.0], [30.0, 60.0]],
        )
        coords = [(time, 0)] + [(coord, dim + 1) for coord, dim in coords]
    return iris.cube.Cube(data, dim_coords_and_dims=coords)


def test_figure4_monthly_cell_volume_matches_mean_fixed_volume():
    """Drake transport uses the period mean when volume varies monthly."""
    script = (
        Path(__file__).resolve().parents[3]
        / "esmvaltool/diag_scripts/russell18jgr/russell18jgr-fig4.py"
    )
    spec = importlib.util.spec_from_file_location("russell_fig4", script)
    module = importlib.util.module_from_spec(spec)
    shared = ModuleType("esmvaltool.diag_scripts.shared")
    for name in (
        "ProvenanceLogger",
        "get_diagnostic_filename",
        "get_plot_filename",
        "run_diagnostic",
    ):
        setattr(shared, name, None)
    with patch.dict(sys.modules, {"esmvaltool.diag_scripts.shared": shared}):
        spec.loader.exec_module(module)

    lat = np.array([-76.0, -60.0, -48.0])
    lon = np.array([290.0, 291.0, 292.0])
    static = _volume_cube(np.full((2, 3, 3), 2.0), lat, lon)
    monthly = _volume_cube(
        np.stack([np.full((2, 3, 3), 1.0), np.full((2, 3, 3), 3.0)]),
        lat,
        lon,
        monthly=True,
    )
    section = np.ones((2, 3))
    expected = module.total_transport(
        section, static, lat, lon, 1, [5.0, 15.0]
    )
    result = module.total_transport(section, monthly, lat, lon, 1, [5.0, 15.0])

    np.testing.assert_allclose(result, expected)


def test_figure6_monthly_cell_volume_matches_mean_fixed_volume():
    """Density-layer sections use a 3D period-mean volume field."""
    lat = np.array([-40.0, -30.0, -20.0])
    lon = np.array([0.0, 1.0, 2.0])
    static = _volume_cube(np.full((2, 3, 3), 2.0), lat, lon)
    monthly = _volume_cube(
        np.stack([np.full((2, 3, 3), 1.0), np.full((2, 3, 3), 3.0)]),
        lat,
        lon,
        monthly=True,
    )
    theta = _volume_cube(np.full((2, 3, 3), 10.0), lat, lon)
    salt = _volume_cube(np.full((2, 3, 3), 35.0), lat, lon)
    velocity = _volume_cube(np.full((2, 3, 3), 1.0), lat, lon)

    expected = prepare_section(theta, salt, velocity, static)
    result = prepare_section(theta, salt, velocity, monthly)

    np.testing.assert_allclose(result["area"], expected["area"])


def test_southern_ocean_flux_sum_preserves_input_masks():
    """Masked land cells must not contribute to integrated ocean flux."""
    flux = np.ma.array([[1.0, 2.0], [3.0, 4.0]], mask=[[0, 1], [0, 0]])
    area = np.ma.array([[10.0, 10.0], [10.0, 10.0]], mask=[[0, 0], [1, 0]])

    result = southern_ocean_flux_sum(flux, area, [-60.0, -30.0], 1.0)

    assert result == 50.0


def test_add_latitude_coord_accepts_nonmonotonic_rows():
    """Curvilinear row latitudes remain attached to their data rows."""
    cube = iris.cube.Cube([1.0, 2.0, 3.0])
    lat = [-75.0, -60.0, -65.0]

    add_latitude_coord(cube, lat)

    assert cube.coord("latitude").points.tolist() == lat
    assert cube.coord_dims("latitude") == (0,)
    assert cube.coord("latitude") in cube.coords(dim_coords=False)
