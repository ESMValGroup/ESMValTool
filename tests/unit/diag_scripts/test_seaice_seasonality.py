"""Tests for daily Antarctic sea ice seasonality."""

import iris
import numpy as np
import pytest

from esmvaltool.diag_scripts.seaice import seaice_seasonality as seasonality


def make_cube(data):
    """Build a daily siconc cube for the leap-year 2000/01 ice year."""
    cube = iris.cube.Cube(data, var_name="siconc", units="%")
    cube.add_dim_coord(
        iris.coords.DimCoord(
            np.arange(data.shape[0]) + 0.5,
            standard_name="time",
            units="days since 2000-02-15",
        ),
        0,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            np.linspace(-80, -60, data.shape[1]),
            standard_name="latitude",
            units="degrees_north",
        ),
        1,
    )
    cube.add_dim_coord(
        iris.coords.DimCoord(
            np.linspace(0, 270, data.shape[2]),
            standard_name="longitude",
            units="degrees_east",
        ),
        2,
    )
    return cube


def test_advance_requires_consecutive_days_and_complete_data():
    """Scattered ice and partial coverage cannot produce event dates."""
    data = np.ma.zeros((365, 1, 4))
    data[100:300, 0, 0] = 0.8
    data[[10, 20, 30, 40, 50], 0, 1] = 0.8
    data[:, 0, 2] = 0.9
    data[100:300, 0, 3] = 0.8
    data[150, 0, 3] = np.ma.masked

    advance, retreat, duration = seasonality.seasonality_fields(data)

    assert advance[0, 0] == 101
    assert retreat[0, 0] == 301
    assert duration[0, 0] == 200
    assert all(field.mask[0, 1] for field in (advance, retreat, duration))
    assert (advance[0, 2], retreat[0, 2], duration[0, 2]) == (
        1,
        365,
        364,
    )
    assert all(field.mask[0, 3] for field in (advance, retreat, duration))


def test_validate_daily_ice_year_rejects_missing_days():
    """The diagnostic must reject irregular time sampling."""
    cube = make_cube(np.zeros((366, 2, 2)))
    assert seasonality.validate_daily_ice_year(cube) == 366
    points = cube.coord("time").points.copy()
    points[10] += 0.25
    cube.coord("time").points = points
    with pytest.raises(ValueError, match="every day"):
        seasonality.validate_daily_ice_year(cube)


def test_main_writes_multimodel_outputs(tmp_path, monkeypatch):
    """Two model inputs yield two NetCDF files and a comparison map."""
    records = []

    class ProvenanceLogger:
        def __init__(self, _cfg):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def log(self, path, record):
            records.append((path, record))

    monkeypatch.setattr(seasonality, "ProvenanceLogger", ProvenanceLogger)
    monkeypatch.setattr(
        seasonality,
        "get_diagnostic_filename",
        lambda name, _cfg: str(tmp_path / f"{name}.nc"),
    )
    monkeypatch.setattr(
        seasonality,
        "get_plot_filename",
        lambda name, _cfg: str(tmp_path / f"{name}.png"),
    )
    inputs = {}
    for model, first_day in (("MODEL-A", 100), ("MODEL-B", 120)):
        data = np.ma.zeros((366, 3, 4))
        data[first_day:300, :, :] = 80
        cube = make_cube(data)
        path = tmp_path / f"{model}.nc"
        iris.save(cube, path)
        inputs[model] = {
            "dataset": model,
            "ensemble": "r1i1p1f1",
            "filename": str(path),
        }

    seasonality.main({"input_data": inputs})

    for model in inputs:
        path = tmp_path / f"seaice_seasonality_{model}_r1i1p1f1.nc"
        assert path.is_file()
        assert {cube.var_name for cube in iris.load(path)} == {
            "siseason_advance",
            "siseason_retreat",
            "siseason_duration",
        }
    assert (tmp_path / "seaice_seasonality_comparison.png").is_file()
    assert len(records) == 3
