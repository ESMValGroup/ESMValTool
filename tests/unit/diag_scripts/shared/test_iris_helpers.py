"""Tests for the module :mod:`esmvaltool.diag_scripts.shared.iris_helpers`."""
from unittest import mock

import iris
import iris.coords
import iris.cube
import iris.exceptions
import numpy as np
import pytest
from cf_units import Unit

from esmvaltool import ESMValToolDeprecationWarning
from esmvaltool.diag_scripts.shared import iris_helpers as ih

LONG_NAME = 'x'
DIM_COORD_1 = iris.coords.DimCoord(np.arange(3.0) - 1.0, long_name=LONG_NAME)
AUX_COORD_1 = iris.coords.AuxCoord(np.arange(3.0) - 1.0, long_name=LONG_NAME)
AUX_COORD_2 = iris.coords.AuxCoord([10.0, 20.0, 30.0], long_name='longer')
SMALL_COORD = iris.coords.DimCoord([0.0], long_name=LONG_NAME)
LONG_COORD_1 = iris.coords.AuxCoord([-1.0, 0.0, 1.0, 1.], long_name=LONG_NAME)
LONG_COORD_2 = iris.coords.DimCoord([-1.0, -0.5, 0.0, 1.0],
                                    long_name=LONG_NAME)
WRONG_COORD = iris.coords.DimCoord([-200.0, +200.0], long_name=LONG_NAME)
SCALAR_COORD = iris.coords.AuxCoord(2.71, long_name='e')
DUP_COORD = iris.coords.AuxCoord([-1.0, 0.0, 1.0, 1.0], long_name=LONG_NAME)
CUBE_1 = iris.cube.Cube(
    np.ma.masked_invalid([-1.0, np.nan, 2.0]),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(DIM_COORD_1, 0)],
    aux_coords_and_dims=[(SCALAR_COORD, []), (AUX_COORD_2, 0)])
CUBE_2 = iris.cube.Cube(
    np.ma.masked_invalid([-1.0, np.nan, 2.0]),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(DIM_COORD_1, 0)],
    aux_coords_and_dims=[(SCALAR_COORD, [])])
CUBE_3 = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, 3.14, np.nan]),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(DIM_COORD_1, 0)])
CUBE_4 = iris.cube.Cube(
    np.ma.masked_invalid([1.0, 2.0, 3.0, 3.0]),
    var_name='a',
    attributes={'1': '2'},
    aux_coords_and_dims=[(SCALAR_COORD, []), (LONG_COORD_1, 0)])
CUBE_5 = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, 3.14, np.nan, np.nan]),
    var_name='a',
    attributes={'1': '2'},
    aux_coords_and_dims=[(LONG_COORD_1, 0)])
CUBE_SMALL = iris.cube.Cube([3.14],
                            var_name='a',
                            attributes={'1': '2'},
                            dim_coords_and_dims=[(SMALL_COORD, 0)])
CUBE_LONG = iris.cube.Cube(
    np.ma.masked_invalid([-1.0, np.nan, np.nan, 2.0]),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(LONG_COORD_2, 0)],
    aux_coords_and_dims=[(SCALAR_COORD, [])])
CUBE_SMALL_LONG = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, np.nan, 3.14, np.nan]),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(LONG_COORD_2, 0)])
CUBE_WRONG = iris.cube.Cube(
    np.arange(2.0),
    var_name='a',
    attributes={'1': '2'},
    dim_coords_and_dims=[(WRONG_COORD, 0)])
CUBE_DUP = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, 3.14, 2.71, 6.28]),
    var_name='a',
    attributes={'1': '2'},
    aux_coords_and_dims=[(DUP_COORD, 0)])
TEST_TRANSFORM_COORD_TO_REF = [
    (DIM_COORD_1, [CUBE_1, CUBE_1], [CUBE_2, CUBE_2]),
    (DIM_COORD_1, [CUBE_SMALL, CUBE_1], [CUBE_3, CUBE_2]),
    (DIM_COORD_1, [CUBE_WRONG, CUBE_1], ValueError),
    (DIM_COORD_1, [CUBE_DUP, CUBE_1], ValueError),
    (AUX_COORD_1, [CUBE_1, CUBE_1], [CUBE_2, CUBE_2]),
    (AUX_COORD_1, [CUBE_SMALL, CUBE_1], [CUBE_3, CUBE_2]),
    (AUX_COORD_1, [CUBE_WRONG, CUBE_1], ValueError),
    (AUX_COORD_1, [CUBE_DUP, CUBE_1], ValueError),
    (LONG_COORD_1, [CUBE_1, CUBE_1], ValueError),
    (LONG_COORD_1, [CUBE_SMALL, CUBE_1], ValueError),
    (LONG_COORD_1, [CUBE_WRONG, CUBE_1], ValueError),
    (LONG_COORD_1, [CUBE_DUP, CUBE_1], ValueError),
    (LONG_COORD_2, [CUBE_1, CUBE_1], [CUBE_LONG, CUBE_LONG]),
    (LONG_COORD_2, [CUBE_SMALL, CUBE_1], [CUBE_SMALL_LONG, CUBE_LONG]),
    (LONG_COORD_2, [CUBE_WRONG, CUBE_1], ValueError),
    (LONG_COORD_2, [CUBE_DUP, CUBE_1], ValueError),
    (DIM_COORD_1, [CUBE_1], [CUBE_2]),
    (DIM_COORD_1, [CUBE_SMALL], [CUBE_3]),
    (DIM_COORD_1, [CUBE_WRONG], ValueError),
    (DIM_COORD_1, [CUBE_DUP], ValueError),
    (AUX_COORD_1, [CUBE_1], [CUBE_2]),
    (AUX_COORD_1, [CUBE_SMALL], [CUBE_3]),
    (AUX_COORD_1, [CUBE_WRONG], ValueError),
    (AUX_COORD_1, [CUBE_DUP], ValueError),
    (LONG_COORD_1, [CUBE_1], ValueError),
    (LONG_COORD_1, [CUBE_SMALL], ValueError),
    (LONG_COORD_1, [CUBE_WRONG], ValueError),
    (LONG_COORD_1, [CUBE_DUP], ValueError),
    (LONG_COORD_2, [CUBE_1], [CUBE_LONG]),
    (LONG_COORD_2, [CUBE_SMALL], [CUBE_SMALL_LONG]),
    (LONG_COORD_2, [CUBE_WRONG], ValueError),
    (LONG_COORD_2, [CUBE_DUP], ValueError),
]


@pytest.mark.parametrize('ref_coord,cubes,output', TEST_TRANSFORM_COORD_TO_REF)
def test_transform_coord_to_ref(ref_coord, cubes, output):
    """Test transforming coordinate to reference."""
    # ValueErrors
    if isinstance(output, type):
        with pytest.raises(output):
            new_cubes = ih._transform_coord_to_ref(cubes, ref_coord)
        return

    # Working examples
    cubes = iris.cube.CubeList(cubes)
    output = iris.cube.CubeList(output)
    new_cubes = ih._transform_coord_to_ref(cubes, ref_coord)
    assert new_cubes == output


DIM_COORD_2 = iris.coords.DimCoord(np.arange(3.0) - 1.0, long_name='aaa')
DIM_COORD_3 = iris.coords.DimCoord(np.arange(3.0) + 1.0, long_name=LONG_NAME)
CUBE_6 = iris.cube.Cube(
    np.ma.arange(3.0) + 100.0,
    var_name='a',
    dim_coords_and_dims=[(DIM_COORD_2, 0)])
CUBE_7 = iris.cube.Cube(
    np.ma.arange(3.0) - 100.0,
    var_name='a',
    dim_coords_and_dims=[(DIM_COORD_3, 0)])
TEST_CHECK_COORDINATE = [
    ([CUBE_1, CUBE_1, CUBE_1], DIM_COORD_1.points),
    ([CUBE_1], DIM_COORD_1.points),
    ([CUBE_1, CUBE_6], iris.exceptions.CoordinateNotFoundError),
    ([CUBE_1, CUBE_7], ValueError),
]


@pytest.mark.parametrize('cubes,output', TEST_CHECK_COORDINATE)
def test_check_coordinate(cubes, output):
    """Test checking of coordinates."""
    if isinstance(output, type):
        with pytest.raises(output):
            out = ih.check_coordinate(cubes, LONG_NAME)
    else:
        out = ih.check_coordinate(cubes, LONG_NAME)
        assert np.array_equal(out, output)


DICT_1 = {'a': 'b', 'c': 'd'}
DICT_2 = {'short_name': 'x'}
DICT_3 = {'var_name': 'x'}
TEST_CONVERT_TO_IRIS = [
    (DICT_1, DICT_1),
    (DICT_2, DICT_3),
    (DICT_3, DICT_3),
    ({
        **DICT_1,
        **DICT_2,
    }, {
        **DICT_1,
        **DICT_3,
    }),
    ({
        **DICT_1,
        **DICT_3,
    }, {
        **DICT_1,
        **DICT_3,
    }),
    ({
        **DICT_1,
        **DICT_2,
        'var_name': ':(',
    }, {
        **DICT_1,
        **DICT_3,
    }),
]


@pytest.mark.parametrize('dict_in,dict_out', TEST_CONVERT_TO_IRIS)
def test_convert_to_iris(dict_in, dict_out):
    """Test converting metadata dictionary checking of coordinates."""
    if 'short_name' in dict_in and 'var_name' in dict_in:
        with pytest.raises(KeyError):
            ih.convert_to_iris(dict_in)
        return
    new_dict = ih.convert_to_iris(dict_in)
    assert new_dict == dict_out
    assert new_dict is not dict_in


@mock.patch('esmvaltool.diag_scripts.shared.iris_helpers.iris.load_cube',
            autospec=True)
def test_get_mean_cube(mock_load_cube):
    """Test calculation of mean cubes."""
    datasets = [
        {'test': 'x', 'filename': 'a/b.nc'},
        {'test': 'y', 'filename': 'a/b/c.nc'},
        {'test': 'z', 'filename': 'c/d.nc'},
    ]
    cube = CUBE_1.copy([-4.0, 2.0, -4.0])
    cube.coord(DIM_COORD_1).attributes = {'test': 1}
    cube.cell_methods = [iris.coords.CellMethod('mean', coords=LONG_NAME)]
    cubes = [CUBE_1, CUBE_2, cube]
    mock_load_cube.side_effect = cubes
    cube_out = iris.cube.Cube(
        [-2.0, 2.0, 0.0],
        var_name='a',
        dim_coords_and_dims=[(DIM_COORD_1, 0)],
        cell_methods=[iris.coords.CellMethod('mean', coords='cube_label')],
    )
    result = ih.get_mean_cube(datasets)
    assert result == cube_out


TEST_IRIS_PROJECT_CONSTRAINT = [
    (['ONE'], False, [2.0, 6.0], ['a', 'e']),
    (['ONE'], True, [3.0, 4.0, 5.0], ['b', 'c', 'd']),
    (['ONE', 'THREE'], False, [2.0, 4.0, 6.0], ['a', 'c', 'e']),
    (['ONE', 'THREE'], True, [3.0, 5.0], ['b', 'd']),
]


@pytest.mark.parametrize('constr,negate,data,points',
                         TEST_IRIS_PROJECT_CONSTRAINT)
def test_iris_project_constraint(constr, negate, data, points):
    """Test iris constraint for projects."""
    input_data = [{
        'project': 'ONE',
        'dataset': 'a',
    }, {
        'project': 'TWO',
        'dataset': 'b',
    }, {
        'project': 'THREE',
        'dataset': 'c',
    }, {
        'project': 'ONE',
        'dataset': 'e',
    }]
    dataset_coord = iris.coords.AuxCoord(['a', 'b', 'c', 'd', 'e'],
                                         long_name='dataset')
    cube = iris.cube.Cube(
        np.arange(5.0) + 2.0, aux_coords_and_dims=[(dataset_coord, 0)])
    new_cube = iris.cube.Cube(
        data,
        aux_coords_and_dims=[(iris.coords.AuxCoord(
            points, long_name='dataset'), 0)])
    constraint = ih.iris_project_constraint(constr, input_data, negate=negate)
    assert cube.extract(constraint) == new_cube


ATTRS = [
    {
        'test': 1,
        'oh': 'yeah',
    },
    {
        'a2': 'c2',
    },
]
VAR_ATTRS = [
    {
        'var_name': 'var',
        'long_name': 'LOOONG NAME',
    },
    {
        'standard_name': 'air_temperature',
        'units': 'K',
    },
]
DATSET_COORD_1 = iris.coords.AuxCoord(['x', 'b', 'c', 'a', 'y', 'z'],
                                      long_name='dataset')
DATSET_COORD_1_SORTED = iris.coords.AuxCoord(['a', 'b', 'c', 'x', 'y', 'z'],
                                             long_name='dataset')
DATSET_COORD_2 = iris.coords.AuxCoord(['t', 'w', 'z', 'b', 'x'],
                                      long_name='dataset')
DATSET_COORD_3 = iris.coords.AuxCoord(['r', 's'], long_name='dataset')
DATSET_COORD_4 = iris.coords.AuxCoord(['c', 'c', 'b', 'a'],
                                      long_name='dataset')
DATSET_COORD_5 = iris.coords.AuxCoord(['b', 'x', 'z'], long_name='dataset')
CUBE_DAT_1 = iris.cube.Cube(
    np.arange(6.0) - 2.0,
    aux_coords_and_dims=[(DATSET_COORD_1, 0)],
    attributes=ATTRS[0],
    **VAR_ATTRS[0])
CUBE_DAT_1_SORTED = iris.cube.Cube([1.0, -1.0, 0.0, -2.0, 2.0, 3.0],
                                   aux_coords_and_dims=[(DATSET_COORD_1_SORTED,
                                                         0)],
                                   attributes=ATTRS[0],
                                   **VAR_ATTRS[0])
CUBE_DAT_1_OUT = iris.cube.Cube([-1.0, -2.0, 3.0],
                                aux_coords_and_dims=[(DATSET_COORD_5, 0)],
                                attributes=ATTRS[0],
                                **VAR_ATTRS[0])
CUBE_DAT_2 = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, 0.0, np.nan, 3.14, 2.71]),
    aux_coords_and_dims=[(DATSET_COORD_2, 0)],
    attributes=ATTRS[1],
    **VAR_ATTRS[1])
CUBE_DAT_2_OUT = iris.cube.Cube(
    np.ma.masked_invalid([3.14, 2.71, np.nan]),
    aux_coords_and_dims=[(DATSET_COORD_5, 0)],
    attributes=ATTRS[1],
    **VAR_ATTRS[1])
CUBE_DAT_3 = iris.cube.Cube(
    np.arange(2.0),
    aux_coords_and_dims=[(DATSET_COORD_3, 0)],
    attributes=ATTRS[0],
    **VAR_ATTRS[0])
CUBE_DAT_4 = iris.cube.Cube(
    np.ma.masked_invalid([np.nan, 2.0, 3.0, 42.0]),
    aux_coords_and_dims=[(DATSET_COORD_4, 0)],
    attributes=ATTRS[1],
    **VAR_ATTRS[1])
TEST_INTERSECT_DATASET_COORDS = [
    ([CUBE_DAT_1, CUBE_1], iris.exceptions.CoordinateNotFoundError),
    ([CUBE_DAT_1, CUBE_DAT_4], ValueError),
    ([CUBE_DAT_1, CUBE_DAT_3], ValueError),
    ([CUBE_DAT_1], [CUBE_DAT_1_SORTED]),
    ([CUBE_DAT_1, CUBE_DAT_1], [CUBE_DAT_1_SORTED, CUBE_DAT_1_SORTED]),
    ([CUBE_DAT_1, CUBE_DAT_2], [CUBE_DAT_1_OUT, CUBE_DAT_2_OUT]),
    ([CUBE_DAT_2, CUBE_DAT_1], [CUBE_DAT_2_OUT, CUBE_DAT_1_OUT]),
]


@pytest.mark.parametrize('cubes,output', TEST_INTERSECT_DATASET_COORDS)
def test_intersect_dataset_coords(cubes, output):
    """Test intersecting dataset coordinates."""
    # ValueErrors
    if isinstance(output, type):
        with pytest.raises(output):
            new_cubes = ih.intersect_dataset_coordinates(cubes)
        return

    # Working examples
    cubes = iris.cube.CubeList(cubes)
    output = iris.cube.CubeList(output)
    new_cubes = ih.intersect_dataset_coordinates(cubes)
    assert new_cubes == output


def test_prepare_cube_for_merging():
    """Test preprocessing cubes before merging."""
    label = 'abcde'
    aux_coord = iris.coords.AuxCoord(label,
                                     var_name='cube_label',
                                     long_name='cube_label')
    cube_in = CUBE_1.copy()
    cube_in.coord(DIM_COORD_1).attributes = {'test_attr': 1}
    cube_in.cell_methods = [iris.coords.CellMethod('mean', coords=LONG_NAME)]
    cube_out = iris.cube.Cube(
        np.ma.masked_invalid([-1.0, np.nan, 2.0]),
        var_name='a',
        dim_coords_and_dims=[(DIM_COORD_1, 0)],
        aux_coords_and_dims=[(aux_coord, [])],
    )
    assert cube_in != cube_out
    ih.prepare_cube_for_merging(cube_in, label)
    assert cube_in == cube_out


DIM_COORD_4 = DIM_COORD_1.copy([100.0, 150.0, 160.0])
DIM_COORD_4.rename('time')
DIM_COORD_LONGEST = DIM_COORD_1.copy([-200.0, -1.0, 0.0, 1.0, 2.0, 3.0, 200.0])
CUBE_8 = CUBE_1.copy()
CUBE_8.coord(LONG_NAME).points = np.array([100.0, 150.0, 160.0])
CUBE_8.coord(LONG_NAME).rename('time')
CUBE_WRONG_COORD = CUBE_WRONG.copy()
CUBE_WRONG_COORD.coord(LONG_NAME).rename('wrooong')
TEST_UNIFY_1D_CUBES = [
    ([CUBE_1, iris.cube.Cube([[1.0]])], LONG_NAME, ValueError),
    ([CUBE_1, iris.cube.Cube(0.0)], LONG_NAME, ValueError),
    (
        [iris.cube.Cube([0.0])],
        LONG_NAME,
        iris.exceptions.CoordinateNotFoundError,
    ),
    (
        [CUBE_1, CUBE_WRONG_COORD, CUBE_3],
        LONG_NAME,
        iris.exceptions.CoordinateNotFoundError,
    ),
    ([CUBE_1, CUBE_4, CUBE_3], LONG_NAME, ValueError),
    ([CUBE_7, CUBE_1, CUBE_WRONG], LONG_NAME, DIM_COORD_LONGEST),
    ([CUBE_8], 'time', DIM_COORD_4),
]


@pytest.mark.parametrize('cubes,coord_name,output', TEST_UNIFY_1D_CUBES)
@mock.patch.object(ih, '_transform_coord_to_ref', autospec=True)
@mock.patch(
    'esmvaltool.diag_scripts.shared.io.iris.util.unify_time_units',
    autospec=True)
def test_unify_1d_cubes(mock_unify_time, mock_transform, cubes, coord_name,
                        output):
    """Test unifying 1D cubes."""
    # ValueErrors
    if isinstance(output, type):
        with pytest.raises(output):
            ih.unify_1d_cubes(cubes, coord_name)
        return

    # Working examples
    cubes = iris.cube.CubeList(cubes)
    ih.unify_1d_cubes(cubes, coord_name)
    assert mock_transform.call_args_list == [mock.call(cubes, output)]
    mock_transform.reset_mock()
    if coord_name == 'time':
        assert mock_unify_time.call_count == 1
    else:
        assert not mock_unify_time.called


@pytest.fixture
def cube_with_time():
    """Cube that includes time coordinate."""
    time_coord = iris.coords.DimCoord(
        [1, 3],
        bounds=[[0, 2], [2, 4]],
        standard_name='time',
        units='days since 1850-01-03',
    )
    cube = iris.cube.Cube(
        [1, 2],
        var_name='x',
        dim_coords_and_dims=[(time_coord, 0)],
    )
    return cube


def test_unify_time_coord_str(cube_with_time):
    """Test ``unify_time_coord``."""
    ih.unify_time_coord(cube_with_time)

    expected_units = Unit('days since 1850-01-01 00:00:00',
                          calendar='standard')
    time_coord = cube_with_time.coord('time')

    assert time_coord.var_name == 'time'
    assert time_coord.standard_name == 'time'
    assert time_coord.long_name == 'time'
    assert time_coord.units == expected_units
    assert time_coord.attributes == {}

    np.testing.assert_array_equal(time_coord.points, [3, 5])
    np.testing.assert_array_equal(time_coord.bounds, [[2, 4], [4, 6]])


def test_unify_time_coord_unit(cube_with_time):
    """Test ``unify_time_coord``."""
    target_units = Unit('days since 1850-01-02 00:00:00', calendar='gregorian')
    ih.unify_time_coord(cube_with_time, target_units=target_units)

    expected_units = Unit('days since 1850-01-02 00:00:00',
                          calendar='gregorian')
    time_coord = cube_with_time.coord('time')

    assert time_coord.var_name == 'time'
    assert time_coord.standard_name == 'time'
    assert time_coord.long_name == 'time'
    assert time_coord.units == expected_units
    assert time_coord.attributes == {}
    assert time_coord.units == expected_units

    np.testing.assert_array_equal(time_coord.points, [2, 4])
    np.testing.assert_array_equal(time_coord.bounds, [[1, 3], [3, 5]])


def test_unify_time_coord_no_bounds(cube_with_time):
    """Test ``unify_time_coord``."""
    cube_with_time.coord('time').bounds = None
    ih.unify_time_coord(cube_with_time, 'days since 1850-01-04')

    expected_units = Unit('days since 1850-01-04 00:00:00')
    time_coord = cube_with_time.coord('time')

    assert time_coord.var_name == 'time'
    assert time_coord.standard_name == 'time'
    assert time_coord.long_name == 'time'
    assert time_coord.units == expected_units
    assert time_coord.attributes == {}

    np.testing.assert_array_equal(time_coord.points, [0, 2])
    assert time_coord.bounds is None


def test_unify_time_coord_no_time(cube_with_time):
    """Test ``unify_time_coord``."""
    cube_with_time.remove_coord('time')
    with pytest.raises(iris.exceptions.CoordinateNotFoundError):
        ih.unify_time_coord(cube_with_time)


def test_var_name_constraint():
    """Test var_name constraint."""
    cubes_in = iris.cube.CubeList([
        iris.cube.Cube(0, var_name='a', long_name='aaa'),
        iris.cube.Cube(1, var_name='a', long_name='bbb'),
        iris.cube.Cube(2, var_name='b', long_name='a'),
        iris.cube.Cube(3, var_name='c', long_name='aaa'),
    ])
    cubes_out = cubes_in[:2].copy()
    with pytest.warns(ESMValToolDeprecationWarning):
        constraint = ih.var_name_constraint('a')
    assert cubes_in is not cubes_out
    result = cubes_in.extract(constraint)
    assert cubes_in is not result
    assert result == cubes_out
