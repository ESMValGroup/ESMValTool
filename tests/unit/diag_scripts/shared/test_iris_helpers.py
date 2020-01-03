"""Tests for the module :mod:`esmvaltool.diag_scripts.shared.iris_helpers`."""

import iris
import mock
import numpy as np
import pytest

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
CUBES_TO_TRANSFORM = [
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


@pytest.mark.parametrize('ref_coord,cubes,output', CUBES_TO_TRANSFORM)
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
CUBES_TO_CHECK_COORD = [
    ([CUBE_1, CUBE_1, CUBE_1], DIM_COORD_1.points),
    ([CUBE_1], DIM_COORD_1.points),
    ([CUBE_1, CUBE_6], iris.exceptions.CoordinateNotFoundError),
    ([CUBE_1, CUBE_7], ValueError),
]


@pytest.mark.parametrize('cubes,output', CUBES_TO_CHECK_COORD)
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
DICTS_TO_CONVERT = [
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


@pytest.mark.parametrize('dict_in,dict_out', DICTS_TO_CONVERT)
@mock.patch.object(ih, 'logger', autospec=True)
def test_convert_to_iris(mock_logger, dict_in, dict_out):
    """Test converting metadata dictionary checking of coordinates."""
    new_dict = ih.convert_to_iris(dict_in)
    assert new_dict == dict_out
    assert new_dict is not dict_in
    if 'short_name' in dict_in and 'var_name' in dict_in:
        mock_logger.warning.assert_called()
    else:
        mock_logger.warning.assert_not_called()


PROJECT_CONSTRAINTS = [
    (['ONE'], False, [2.0, 6.0], ['a', 'e']),
    (['ONE'], True, [3.0, 4.0, 5.0], ['b', 'c', 'd']),
    (['ONE', 'THREE'], False, [2.0, 4.0, 6.0], ['a', 'c', 'e']),
    (['ONE', 'THREE'], True, [3.0, 5.0], ['b', 'd']),
]


@pytest.mark.parametrize('constr,negate,data,points', PROJECT_CONSTRAINTS)
@mock.patch.object(ih, 'logger', autospec=True)
def test_iris_project_constraint(mock_logger, constr, negate, data, points):
    """Test iris constraint for projects."""
    cfg = {
        'input_data': {
            'p1': {
                'project': 'ONE',
                'dataset': 'a',
            },
            'p2': {
                'project': 'TWO',
                'dataset': 'b',
            },
            'p3': {
                'project': 'THREE',
                'dataset': 'c',
            },
            'p4': {
                'project': 'ONE',
                'dataset': 'e',
            },
        },
        'does_not_matter': 'oh no',
    }
    dataset_coord = iris.coords.AuxCoord(['a', 'b', 'c', 'd', 'e'],
                                         long_name='dataset')
    cube = iris.cube.Cube(
        np.arange(5.0) + 2.0, aux_coords_and_dims=[(dataset_coord, 0)])
    new_cube = iris.cube.Cube(
        data,
        aux_coords_and_dims=[(iris.coords.AuxCoord(
            points, long_name='dataset'), 0)])
    constraint = ih.iris_project_constraint(constr, cfg, negate=negate)
    assert cube.extract(constraint) == new_cube
    mock_logger.warning.assert_not_called()
    mock_logger.reset_mock()
    cfg['input_data']['p5'] = {'project': 'ONE', 'ohhh': 1}
    constraint = ih.iris_project_constraint(constr, cfg, negate=negate)
    assert cube.extract(constraint) == new_cube
    mock_logger.warning.assert_called_once()


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
CUBES_TO_INTERSECT = [
    ([CUBE_DAT_1, CUBE_1], iris.exceptions.CoordinateNotFoundError),
    ([CUBE_DAT_1, CUBE_DAT_4], ValueError),
    ([CUBE_DAT_1, CUBE_DAT_3], ValueError),
    ([CUBE_DAT_1], [CUBE_DAT_1_SORTED]),
    ([CUBE_DAT_1, CUBE_DAT_1], [CUBE_DAT_1_SORTED, CUBE_DAT_1_SORTED]),
    ([CUBE_DAT_1, CUBE_DAT_2], [CUBE_DAT_1_OUT, CUBE_DAT_2_OUT]),
    ([CUBE_DAT_2, CUBE_DAT_1], [CUBE_DAT_2_OUT, CUBE_DAT_1_OUT]),
]


@pytest.mark.parametrize('cubes,output', CUBES_TO_INTERSECT)
def test_intersect_dataset_coords(cubes, output):
    """Test unifying 1D cubes."""
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


DIM_COORD_4 = DIM_COORD_1.copy([100.0, 150.0, 160.0])
DIM_COORD_4.rename('time')
DIM_COORD_LONGEST = DIM_COORD_1.copy([-200.0, -1.0, 0.0, 1.0, 2.0, 3.0, 200.0])
CUBE_8 = CUBE_1.copy()
CUBE_8.coord(LONG_NAME).points = np.array([100.0, 150.0, 160.0])
CUBE_8.coord(LONG_NAME).rename('time')
CUBE_WRONG_COORD = CUBE_WRONG.copy()
CUBE_WRONG_COORD.coord(LONG_NAME).rename('wrooong')
CUBES_TO_UNIFY = [
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


@pytest.mark.parametrize('cubes,coord_name,output', CUBES_TO_UNIFY)
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
