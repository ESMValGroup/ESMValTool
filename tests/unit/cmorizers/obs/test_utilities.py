"""Tests for the module :mod:`esmvaltool.cmorizers.obs.utilities`."""

from unittest.mock import Mock

import dask.array as da
import iris
import numpy as np
import pytest
from cf_units import Unit

import esmvaltool.cmorizers.obs.utilities as utils


def np_to_da(array, lazy):
    """Convert numpy array to dask array."""
    if not lazy:
        return array
    if array is None:
        return array
    return da.from_array(array)


def is_lazy(cube):
    """Check if data is lazy."""
    if not cube.has_lazy_data():
        return False
    for coord in cube.coords(dim_coords=False):
        if not coord.has_lazy_points():
            return False
        if coord.has_bounds():
            if not coord.has_lazy_bounds():
                return False
    return True


def cubes_generator(lazy=True):
    """Generate a list of cubes via test parametrization."""
    cube_datas = [
        np.array([[0, 1], [-1, 0]], dtype=np.int),
        np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float32),
        np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float64),
        np.ma.masked_equal([[0, 1], [2, 3]], 3).astype(np.int),
        np.ma.masked_values([[0.0, 1.0], [2.0, 3.0]], 3.0).astype(np.float32),
        np.ma.masked_values([[0.0, 1.0], [2.0, 3.0]], 3.0).astype(np.float64),
    ]
    x_coords = [
        (np.array([1, 3], dtype=np.int), None),
        (np.array([1, 3],
                  dtype=np.int), np.array([[0, 2], [2, 4]], dtype=np.int)),
        (np.array([1.0, 3.0], dtype=np.float32),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float32)),
        (np.array([1.0, 3.0], dtype=np.float64), None),
        (np.array([1.0, 3.0], dtype=np.float64),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float64)),
    ]
    y_coords = [
        (np.array([1, 3], dtype=np.int),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float32)),
        (np.array([1.0, 3.0], dtype=np.float32),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float64)),
        (np.array([1.0, 3.0],
                  dtype=np.float64), np.array([[0, 2], [2, 4]], dtype=np.int)),
    ]
    for cube_data in cube_datas:
        cube_data = np_to_da(cube_data, lazy)
        for x_val in x_coords:
            x_val = (np_to_da(x_val[0], lazy), np_to_da(x_val[1], lazy))
            x_coord = iris.coords.DimCoord(x_val[0],
                                           bounds=x_val[1],
                                           var_name='x')
            for y_val in y_coords:
                y_val = (np_to_da(y_val[0], lazy), np_to_da(y_val[1], lazy))
                y_coord = iris.coords.DimCoord(y_val[0],
                                               bounds=y_val[1],
                                               var_name='y')
                aux_coord = iris.coords.AuxCoord(y_val[0],
                                                 bounds=y_val[1],
                                                 var_name='aux')
                cube = iris.cube.Cube(
                    cube_data,
                    var_name='test_var',
                    dim_coords_and_dims=[(x_coord, 0), (y_coord, 1)],
                    aux_coords_and_dims=[(aux_coord, 0)],
                )
                yield cube


@pytest.mark.parametrize('cube', cubes_generator(lazy=True))
def test_fix_dtype_lazy(cube):
    """Test fix for lazy data."""
    assert is_lazy(cube)
    utils._fix_dtype(cube)
    assert cube.dtype == np.float32
    for coord in cube.coords():
        assert coord.dtype == np.float64
        if coord.has_bounds():
            assert coord.bounds_dtype == np.float64
    assert is_lazy(cube)


@pytest.mark.parametrize('cube', cubes_generator(lazy=False))
def test_fix_dtype_not_lazy(cube):
    """Test fix for realized data."""
    assert not is_lazy(cube)
    utils._fix_dtype(cube)
    assert cube.dtype == np.float32
    for coord in cube.coords():
        assert coord.dtype == np.float64
        if coord.has_bounds():
            assert coord.bounds_dtype == np.float64
    assert not is_lazy(cube)


def mock_var_info(var_dict):
    mock_dict = Mock()
    mock_dict.__dict__ = var_dict
    return mock_dict


def _create_sample_cube():
    """Create a quick CMOR-compliant sample cube."""
    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube_data = np.ones((2, 3, 2, 2))
    cube_data[1, 1, 1, 1] = 22.
    time = iris.coords.DimCoord([15, 45],
                                standard_name='time',
                                bounds=[[1., 30.], [30., 60.]],
                                units=Unit('days since 1950-01-01',
                                           calendar='gregorian'))
    zcoord = iris.coords.DimCoord([0.5, 5., 50.],
                                  var_name='depth',
                                  standard_name='depth',
                                  bounds=[[0., 2.5], [2.5, 25.], [25., 250.]],
                                  units='m',
                                  attributes={'positive': 'down'})
    lons = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='longitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='degrees_east',
                                coord_system=coord_sys)
    lats = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='latitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='degrees_north',
                                coord_system=coord_sys)
    coords_spec = [(time, 0), (zcoord, 1), (lats, 2), (lons, 3)]
    cube = iris.cube.Cube(cube_data, dim_coords_and_dims=coords_spec)
    return cube


def test_add_scalar_height_coord():
    """Test add height aux coord."""
    cube = _create_sample_cube()
    utils.add_scalar_height_coord(cube, height=10.)
    assert cube.coord("height").points[0] == 10.
    assert "positive" in cube.coord("height").attributes
    assert cube.coord("height").attributes["positive"] == "up"


@pytest.mark.parametrize('time_units',
                         ['months since 1950-01-01 00:00:00',
                          'days since 0000-01-01 00:00:00',
                          'days since 1950-1-1',
                          'days since 1950-1-1 00:00:00'])
def test_convert_time_units(time_units):
    """Test convert time units functionlity."""
    cube = _create_sample_cube()
    cube.coord("time").units = time_units
    utils.convert_timeunits(cube, "1950")
    converted_units = cube.coord("time").units
    if time_units == 'months since 1950-01-01 00:00:00':
        assert converted_units == 'months since 1950-01-01 00:00:00'
    else:
        assert converted_units == 'days since 1950-01-01 00:00:00'


def test_fix_coords():
    """Test fix coordinates."""
    cube = _create_sample_cube()
    cube.coord("time").bounds = None
    cube.coord('time').convert_units(
        Unit('days since 1850-1-1 00:00:00', calendar='gregorian'))
    cube.coord("longitude").bounds = None
    cube.coord("latitude").bounds = None
    cube.coord("depth").bounds = None
    cube.coord("longitude").points = cube.coord("longitude").points - 3.
    cube.coord("time").var_name = "cows"
    cube.coord("longitude").var_name = "cows"
    cube.coord("latitude").var_name = "cows"
    cube.coord("longitude").units = "m"
    cube.coord("latitude").units = "K"
    utils.fix_coords(cube)
    assert cube.coord("time").var_name == "time"
    assert cube.coord("longitude").var_name == "lon"
    assert cube.coord("latitude").var_name == "lat"
    assert cube.coord("longitude").standard_name == "longitude"
    assert cube.coord("latitude").standard_name == "latitude"
    assert cube.coord("longitude").long_name == "longitude coordinate"
    assert cube.coord("latitude").long_name == "latitude coordinate"
    assert cube.coord("longitude").units == "degrees"
    assert cube.coord("latitude").units == "degrees"
    assert cube.coord("depth").var_name == "lev"
    assert cube.coord("depth").attributes['positive'] == "down"
    assert cube.coord("time").has_bounds()
    assert cube.coord("time").bounds[0][1] == 30.
    assert cube.coord("time").units == 'days since 1950-1-1 00:00:00'
    assert cube.coord("time").units.calendar == "gregorian"
    assert cube.coord("longitude").points[0] == 178.5
    assert cube.coord("longitude").points[1] == 179.5
    assert cube.coord("longitude").has_bounds()
    assert cube.coord("longitude").bounds[1][1] == 180.
    assert cube.data[1, 1, 1, 0] == 22.
    assert cube.coord("latitude").has_bounds()
    assert cube.coord("depth").has_bounds()
    assert cube.coord('latitude').coord_system is None
    assert cube.coord('longitude').coord_system is None


def test_fix_var_metadata():
    """Test fixing the variable metadata."""
    cube = _create_sample_cube()
    cube.var_name = "cows"
    cube.long_name = "flying cows"
    cube.units = "m"
    var_info = {
        "short_name": "tas",
        "frequency": "mon",
        "modeling_realm": "atmos",
        "standard_name": "air_temperature",
        "units": "K",
        "cell_methods": "area: time: mean",
        "cell_measures": "area: areacella",
        "long_name": "Near-Surface Air Temperature",
        "comment": "near-surface (usually, 2 meter) air temperature",
        "dimensions": "longitude latitude time height2m",
        "out_name": "tas",
        "type": "real",
        "positive": "",
        "valid_min": "",
        "valid_max": "",
        "ok_min_mean_abs": "",
        "ok_max_mean_abs": ""
    }
    var_info = mock_var_info(var_info)
    utils.fix_var_metadata(cube, var_info)
    assert cube.var_name == "tas"
    assert cube.long_name == "Near-Surface Air Temperature"
    assert cube.units == "K"
    assert cube.standard_name == "air_temperature"


def test_set_global_atts_correct():
    """Test set global attributes."""
    cube = _create_sample_cube()
    global_attrs = {
        'dataset_id': '1',
        'version': '2',
        'tier': '3',
        'source': '4',
        'reference': 'acknow_author',
        'comment': '6',
        'project_id': '7',
    }
    utils.set_global_atts(cube, global_attrs)
    attrs = cube.attributes
    assert '1 ' in attrs['title']
    assert attrs['version'] == '2'
    assert attrs['tier'] == '3'
    assert attrs['source'] == '4'
    assert attrs['reference'] == 'doi not found'
    assert attrs['comment'] == '6'
    assert attrs['project_id'] == '7'


def test_set_global_atts_incorrect():
    """Test set global attributes."""
    cube = _create_sample_cube()
    global_attrs = {
        'version': '2',
        'tier': '3',
        'source': '4',
        'reference': 'acknow_author',
        'comment': '6',
        'project_id': '7',
    }
    msg = \
        "".join(["All CMORized datasets need the ",
                 "global attributes 'dataset_id', ",
                 "'version', 'tier', 'source', 'reference', 'comment' and ",
                 "'project_id' specified in the configuration file"])
    with pytest.raises(KeyError) as key_err:
        utils.set_global_atts(cube, global_attrs)
        assert msg in key_err


def test_flip_dim_coord():
    """Test flip dimensional coordinate."""
    cube = _create_sample_cube()
    assert cube.data[1, 1, 1, 1] == 22.
    utils.flip_dim_coord(cube, "latitude")
    assert cube.data[1, 1, 0, 1] == 22.


def test_read_cmor_config():
    """Test the cmor table reading."""
    cfg = utils.read_cmor_config("WOA")
    assert cfg['attributes']['dataset_id'] == 'WOA'
    assert 'thetao' in cfg['variables']
    assert 'Omon' in cfg['cmor_table'].tables
    assert 'thetao' in cfg['cmor_table'].tables['Omon']
