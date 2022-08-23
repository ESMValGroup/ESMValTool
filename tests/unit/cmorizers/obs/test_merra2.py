import iris
import iris.coord_systems
import iris.fileformats
import numpy as np
import pytest
from cf_units import Unit

from esmvaltool.cmorizers.data.formatters.datasets.merra2 import (
    _load_cube
)


def _create_sample_cube():
    """Create a quick CMOR-compliant sample cube."""
    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube_data = np.ones((1, 3, 2, 2))
    cube_data[0, 1, 1, 1] = 22.
    time = iris.coords.DimCoord([
        1.0,
    ],
                                standard_name='time',
                                bounds=[[0.5, 1.5]],
                                units=Unit('days since 0000-01-01',
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
                                units='K',
                                coord_system=coord_sys)
    lats = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='latitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='K',
                                coord_system=coord_sys)
    coords_spec = [(time, 0), (zcoord, 1), (lats, 2), (lons, 3)]
    cube = iris.cube.Cube(cube_data, dim_coords_and_dims=coords_spec)
    drop_attrs = [
        'History', 'Filename', 'Comment', 'RangeBeginningDate',
        'RangeEndingDate', 'GranuleID', 'ProductionDateTime', 'Source'
    ]
    for attr in drop_attrs:
        cube.attributes[attr] = "cow"
    drop_time_attrs = [
        'begin_date', 'begin_time', 'time_increment', 'valid_range', 'vmax',
        'vmin'
    ]
    for attr in drop_time_attrs:
        cube.coord('time').attributes[attr] = "1982"
    cube.coord('time').attributes["valid_range"] = [1982, 1989]

    return cube


def test_load_cube_pairwise_vars(tmp_path):
    """Test loading MERRA2 cubes."""
    path_cubes = tmp_path / "cubes.nc"
    cube_1 = _create_sample_cube()
    cube_1.var_name = "SWTDN"
    cube_2 = _create_sample_cube()
    cube_2.var_name = "SWTNT"
    cubes = iris.cube.CubeList([cube_1, cube_2])
    iris.save(cubes, str(path_cubes))
    var = {
        'short_name': 'rsut',
        'mip': 'Amon', 'raw': 'SWTDN-SWTNT',
        'file': 'MERRA2_???.tavgM_2d_rad_Nx.{year}??.nc4'
    }
    in_files = str(tmp_path / "cubes.nc")
    selection = _load_cube(in_files, var)
    assert np.mean(selection.data) == 0.


def test_load_cube_threewise_vars(tmp_path):
    """Test loading MERRA2 cubes."""
    path_cubes = tmp_path / "cubes.nc"
    cube_1 = _create_sample_cube()
    cube_1.var_name = "SWTDN"
    cube_2 = _create_sample_cube()
    cube_2.var_name = "SWTNT"
    cube_3 = _create_sample_cube()
    cube_3.var_name = "COW"
    cubes = iris.cube.CubeList([cube_1, cube_2, cube_3])
    iris.save(cubes, str(path_cubes))
    var = {
        'short_name': 'rsut',
        'mip': 'Amon', 'raw': 'SWTDN-SWTNT-COW',
        'file': 'MERRA2_???.tavgM_2d_rad_Nx.{year}??.nc4'
    }
    in_files = str(tmp_path / "cubes.nc")
    with pytest.raises(NotImplementedError) as exc:
        _load_cube(in_files, var)
    print(exc)
