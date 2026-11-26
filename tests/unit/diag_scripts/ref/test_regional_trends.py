"""Check that sequential region processing preserves trend values."""

import iris
import numpy as np
from esmvalcore import preprocessor as pp
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.diag_scripts.ref.regional_trends import _regional_trend


def test_regional_trend_matches_decomposed_preprocessor():
    """Compare the original and sequential paths, including missing cells."""
    time = DimCoord(
        np.arange(36) * 30,
        standard_name="time",
        units="days since 2000-01-01",
    )
    latitude = DimCoord(
        np.arange(-85.0, 86.0, 10.0),
        standard_name="latitude",
        units="degrees_north",
    )
    longitude = DimCoord(
        np.arange(5.0, 360.0, 10.0),
        standard_name="longitude",
        units="degrees_east",
    )
    latitude.guess_bounds()
    longitude.guess_bounds()
    rng = np.random.default_rng(1)
    data = rng.normal(size=(36, 18, 36)) + np.arange(36)[:, None, None]
    mask = np.zeros(data.shape, dtype=bool)
    mask[:, 13, 24] = True  # 45°N, 115°W in WNA.
    mask[:6, 12, 11] = True  # 35°N, 115°E in EAS.
    cube = Cube(
        np.ma.array(data, mask=mask),
        standard_name="air_temperature",
        units="K",
        dim_coords_and_dims=[(time, 0), (latitude, 1), (longitude, 2)],
    )
    regions = ["WNA", "EAS", "SAH", "NAU"]
    unmasked_cube = cube.copy(data=data)
    for region, missing_count in (("WNA", 36), ("EAS", 6)):
        kwargs = {"decomposed": True, "ids": {"Acronym": [region]}}
        selected = pp.extract_shape(cube, "ar6", **kwargs)
        unmasked = pp.extract_shape(unmasked_cube, "ar6", **kwargs)
        assert (
            np.ma.count_masked(selected.data)
            - np.ma.count_masked(
                unmasked.data,
            )
            == missing_count
        )

    original = pp.extract_shape(
        cube,
        "ar6",
        decomposed=True,
        ids={"Acronym": regions},
    )
    original = pp.area_statistics(original, operator="mean")
    original = pp.linear_trend(original)
    original = pp.convert_units(original, "K yr-1")

    sequential = _regional_trend(cube, regions, "K yr-1")
    assert set(original.coord("shape_id").points) == set(
        sequential.coord("shape_id").points,
    )
    assert original.units == sequential.units
    for region in regions:
        constraint = iris.Constraint(shape_id=region)
        expected = original.extract(constraint).data
        actual = sequential.extract(constraint).data
        np.testing.assert_array_equal(
            np.ma.getmaskarray(expected),
            np.ma.getmaskarray(actual),
        )
        np.testing.assert_allclose(expected, actual, rtol=1e-10)
