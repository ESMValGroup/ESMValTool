"""Tests for the module :mod:`esmvaltool.diag_scripts.shared.io`."""
import os
from collections import OrderedDict
from copy import deepcopy
from unittest import mock

import iris
import numpy as np
import pytest
import yaml

from esmvaltool.diag_scripts.shared import io

with open(os.path.join(os.path.dirname(__file__), 'configs',
                       'test_io.yml')) as file_:
    CONFIG = yaml.safe_load(file_)


@pytest.mark.parametrize('data', CONFIG['_has_necessary_attributes'])
@mock.patch.object(io, 'logger', autospec=True)
def test_has_necessary_attributes(mock_logger, data):
    """Test attribute checks."""
    for log_level in ('debug', 'info', 'warning', 'error', 'exception'):
        metadata = data['input']
        kwargs = data.get('kwargs', {})
        has_atts = io._has_necessary_attributes(metadata,
                                                log_level=log_level,
                                                **kwargs)
        assert has_atts == data['output']
        logger_func = getattr(mock_logger, log_level)
        if has_atts:
            logger_func.assert_not_called()
        else:
            logger_func.assert_called()
        assert logger_func.call_count == data['n_logger']
        mock_logger.reset_mock()


CFG = {
    'input_files': [
        'metadata.yml',
        'test_metadata.yml',
        'valid/dir/1',
        'valid/dir/2',
    ],
    'other_attr':
    'I am not used!',
}
ROOT_DIR = '/root/to/something'
TEST_GET_ALL_ANCESTOR_FILES = [
    (None, [
        os.path.join(ROOT_DIR, 'egg.yml'),
        os.path.join(ROOT_DIR, 'root2', 'x.nc'),
        os.path.join(ROOT_DIR, 'root2', 'y.png'),
        os.path.join(ROOT_DIR, 'root3', 'egg.nc'),
        os.path.join(ROOT_DIR, 'root4', 'egg.nc'),
        os.path.join(ROOT_DIR, 'test.nc'),
        os.path.join(ROOT_DIR, 'test_1.nc'),
        os.path.join(ROOT_DIR, 'test_2.yml'),
    ]),
    ('*', [
        os.path.join(ROOT_DIR, 'egg.yml'),
        os.path.join(ROOT_DIR, 'root2', 'x.nc'),
        os.path.join(ROOT_DIR, 'root2', 'y.png'),
        os.path.join(ROOT_DIR, 'root3', 'egg.nc'),
        os.path.join(ROOT_DIR, 'root4', 'egg.nc'),
        os.path.join(ROOT_DIR, 'test.nc'),
        os.path.join(ROOT_DIR, 'test_1.nc'),
        os.path.join(ROOT_DIR, 'test_2.yml'),
    ]),
    ('*.nc', [
        os.path.join(ROOT_DIR, 'root2', 'x.nc'),
        os.path.join(ROOT_DIR, 'root3', 'egg.nc'),
        os.path.join(ROOT_DIR, 'root4', 'egg.nc'),
        os.path.join(ROOT_DIR, 'test.nc'),
        os.path.join(ROOT_DIR, 'test_1.nc'),
    ]),
    ('test*', [
        os.path.join(ROOT_DIR, 'test.nc'),
        os.path.join(ROOT_DIR, 'test_1.nc'),
        os.path.join(ROOT_DIR, 'test_2.yml'),
    ]),
    ('*.yml', [
        os.path.join(ROOT_DIR, 'egg.yml'),
        os.path.join(ROOT_DIR, 'test_2.yml'),
    ]),
    ('egg.nc*', [
        os.path.join(ROOT_DIR, 'root3', 'egg.nc'),
        os.path.join(ROOT_DIR, 'root4', 'egg.nc'),
    ]),
]


@pytest.mark.parametrize('pattern,output', TEST_GET_ALL_ANCESTOR_FILES)
@mock.patch('esmvaltool.diag_scripts.shared.io.os.walk', autospec=True)
def test_get_all_ancestor_files(mock_walk, pattern, output):
    """Test retrieving of ancestor files."""
    input_dirs = [
        [
            (ROOT_DIR, ['dir', '__pycache__'], ['test.nc', 'egg.yml']),
            (os.path.join(ROOT_DIR, 'root2'), ['d'], ['x.nc', 'y.png']),
            (os.path.join(ROOT_DIR, 'root3'), [], ['egg.nc']),
        ],
        [
            (ROOT_DIR, ['dir', '__pycache__'], ['test_1.nc', 'test_2.yml']),
            (os.path.join(ROOT_DIR, 'root4'), ['d2'], ['egg.nc']),
        ],
    ]
    mock_walk.side_effect = input_dirs
    files = io.get_all_ancestor_files(CFG, pattern=pattern)
    assert files == output


TEST_GET_ANCESTOR_FILE = [
    ([], ValueError),
    (['I/am/a/cool/file.nc'], 'I/am/a/cool/file.nc'),
    (['I/am/a/cool/file.nc', 'oh/no/file_2.nc'], ValueError),
]


@pytest.mark.parametrize('files,output', TEST_GET_ANCESTOR_FILE)
@mock.patch.object(io, 'get_all_ancestor_files', autospec=True)
def test_get_ancestor_file(mock_get_all_ancestors, files, output):
    """Test retrieving of single ancestor file."""
    mock_get_all_ancestors.return_value = files
    if isinstance(output, type):
        with pytest.raises(output):
            io.get_ancestor_file(CFG, pattern='*')
        return
    returned_file = io.get_ancestor_file(CFG, pattern='*')
    assert returned_file == output


INVALID_STANDARD_NAME = 'I_am_an_invalid_standard_name'
LONG_NAME = 'Loooong name'
SHORT_NAME = 'var'
STANDARD_NAME = 'air_temperature'
UNITS = 'K'

A_1 = {
    'dataset': 'model',
    'filename': 'r/a.nc',
    'project': 'CMIP42',
}
V_1 = {
    'long_name': LONG_NAME,
    'var_name': SHORT_NAME,
    'units': UNITS,
}
C_1 = iris.cube.Cube(0, **V_1, attributes=A_1)
A_2 = {
    'dataset': 'model',
    'filename': 'r1/b.ps',
}
V_2 = {
    'long_name': LONG_NAME,
    'var_name': SHORT_NAME,
    'units': UNITS,
}
C_2 = iris.cube.Cube(0, **V_2, attributes=A_2)
A_3 = {
    'filename': 'r/a.nc',
}
V_3 = {
    'long_name': LONG_NAME,
    'var_name': SHORT_NAME,
    'units': UNITS,
}
C_3 = iris.cube.Cube(0, **V_3, attributes=A_3)
A_4 = {
    'dataset': 'model',
    'filename': 'r1/b.nc',
    'project': 'CMIP42',
}
V_4 = {
    'long_name': LONG_NAME,
    'var_name': SHORT_NAME,
    'standard_name': STANDARD_NAME,
    'units': UNITS,
}
C_4 = iris.cube.Cube(0, **V_4, attributes=A_4)
A_5 = {
    'dataset': 'model',
    'filename': 'r/a.nc',
    'project': 'CMIP42',
}
V_5 = {
    'long_name': LONG_NAME,
    'var_name': SHORT_NAME,
    'standard_name': None,
    'units': UNITS,
}
C_5 = iris.cube.Cube(0, **V_5, attributes=A_5)


W_1 = [('r', [], ['a.nc'])]
W_2 = [('r', [], ['a.nc']), ('r1', ['d1'], ['b.nc'])]
W_2_X = [('r1', [], ['b.nc', 'b.ps'])]

TEST_NETCDF_TO_METADATA = [
    ([C_1], W_1, None, [{**A_1, **V_1}], 0),
    ([C_1], W_1, '*', [{**A_1, **V_1}], 0),
    ([C_1, C_4], W_2, None, [{**A_1, **V_1}, {**A_4, **V_4}], 0),
    ([C_1, C_4], W_2, '*', [{**A_1, **V_1}, {**A_4, **V_4}], 0),
    ([C_5, C_4], W_2, None, [{**A_5, **V_5}, {**A_4, **V_4}], 0),
    ([C_5, C_4], W_2, '*', [{**A_5, **V_5}, {**A_4, **V_4}], 0),
    ([C_4], W_2_X, None, [{**A_4, **V_4}], 0),
    ([C_4], W_2_X, '*', [{**A_4, **V_4}], 0),
    ([C_2], W_1, None, ValueError, 1),
    ([C_2], W_1, '*', ValueError, 1),
    ([C_3], W_1, None, ValueError, 2),
    ([C_3], W_1, '*', ValueError, 2),
    ([C_2, C_3], W_2, None, ValueError, 3),
    ([C_2, C_3], W_2, '*', ValueError, 3),
    ([C_1, C_3], W_2, None, ValueError, 2),
    ([C_1, C_3], W_2, '*', ValueError, 2),
]


@pytest.mark.parametrize('cubes,walk_out,root,output,n_logger',
                         TEST_NETCDF_TO_METADATA)
@mock.patch.object(io, 'get_all_ancestor_files', autospec=True)
@mock.patch.object(io, 'logger', autospec=True)
@mock.patch('esmvaltool.diag_scripts.shared.io.iris.load_cube', autospec=True)
@mock.patch('esmvaltool.diag_scripts.shared.io.os.walk', autospec=True)
def test_netcdf_to_metadata(mock_walk, mock_load_cube, mock_logger,
                            mock_get_all_ancestors, cubes, walk_out, root,
                            output, n_logger):
    """Test cube to metadata."""
    ancestors = []
    for (files_root, _, files) in walk_out:
        new_files = [os.path.join(files_root, f) for f in files]
        ancestors.extend(new_files)
    mock_get_all_ancestors.return_value = ancestors
    mock_walk.return_value = walk_out
    mock_load_cube.side_effect = cubes
    if isinstance(output, type):
        with pytest.raises(output):
            io.netcdf_to_metadata({}, pattern=root, root=root)
    else:
        for dataset in output:
            dataset['short_name'] = dataset.pop('var_name')
            dataset.setdefault('standard_name', None)
        metadata = io.netcdf_to_metadata({}, pattern=root, root=root)
        assert metadata == output
    assert mock_logger.error.call_count == n_logger


ATTRS_IN = [
    {
        'dataset': 'a',
        'filename': 'path/to/model1.nc',
        'project': 'CMIP42',
        'bool': True,
    },
    {
        'dataset': 'b',
        'filename': 'path/to/model2.nc',
        'project': 'CMIP42',
    },
    {
        'dataset': 'c',
        'filename': 'path/to/model3.nc',
    },
    {
        'dataset': 'd',
        'filename': 'path/to/model4.nc',
        'project': 'CMIP42',
    },
]
ATTRS_OUT = [
    {
        'dataset': 'a',
        'filename': 'path/to/model1.nc',
        'project': 'CMIP42',
        'bool': 'True',
        'invalid_standard_name': INVALID_STANDARD_NAME,
        'attr': 'test',
    },
    {
        'dataset': 'b',
        'filename': 'path/to/model2.nc',
        'project': 'CMIP42',
        'attr': 'test',
    },
    {},
    {
        'dataset': 'd',
        'filename': 'path/to/model4.nc',
        'project': 'CMIP42',
        'attr': 'test',
    },
]
VAR_ATTRS_IN = [
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'units': UNITS,
    },
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'units': UNITS,
    },
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
    },
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'standard_name': STANDARD_NAME,
        'units': UNITS,
    },
]
VAR_ATTRS_OUT = [
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'standard_name': None,
        'units': UNITS,
    },
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'standard_name': None,
        'units': UNITS,
    },
    {},
    {
        'long_name': LONG_NAME,
        'var_name': SHORT_NAME,
        'standard_name': STANDARD_NAME,
        'units': UNITS,
    },
]
ADD_ATTRS = {'project': 'PROJECT', 'attr': 'test'}
ADD_VAR_ATTRS = {'standard_name': STANDARD_NAME, 'var_name': 'test'}
CUBES_IN = [
    iris.cube.Cube(0, attributes=ADD_ATTRS, **ADD_VAR_ATTRS) for _ in range(4)
]
OUTPUT = [
    iris.cube.Cube(0, attributes=ATTRS_OUT[idx], **VAR_ATTRS_OUT[idx]) for idx
    in range(4)
]
OUTPUT[2] = ValueError
for var_attr in VAR_ATTRS_IN:
    var_attr['short_name'] = var_attr.pop('var_name')
ATTRS_IN[0]['standard_name'] = INVALID_STANDARD_NAME
METADATA = [{**a, **VAR_ATTRS_IN[idx]} for (idx, a) in enumerate(ATTRS_IN)]
TEST_METADATA_TO_NETDCF = zip(METADATA, CUBES_IN, OUTPUT)


@pytest.mark.parametrize('metadata,cube,output', TEST_METADATA_TO_NETDCF)
@mock.patch.object(io, 'iris_save', autospec=True)
@mock.patch.object(io, 'logger', autospec=True)
def test_metadata_to_netcdf(mock_logger, mock_save, metadata, cube, output):
    """Test metadata to cube."""
    if isinstance(output, type):
        with pytest.raises(output):
            io.metadata_to_netcdf(cube, metadata)
        assert not mock_save.called
        return
    io.metadata_to_netcdf(cube, metadata)
    if metadata.get('standard_name') == INVALID_STANDARD_NAME:
        mock_logger.warning.assert_called()
        assert 'invalid_standard_name' in output.attributes
    else:
        mock_logger.warning.assert_not_called()
        assert 'invalid_standard_name' not in output.attributes
    save_args = (output, metadata['filename'])
    assert mock_save.call_args_list == [mock.call(*save_args)]


PATH = 'path/to/super/cube'
VAR_ATTRS_NEW = [
    {
        'long_name': 'I do not have units :(',
        'short_name': 'sad',
    },
    {
        'long_name': 'Long name',
        'short_name': 'var',
        'units': '1',
    },
    {
        'short_name': SHORT_NAME,
        'long_name': LONG_NAME,
        'standard_name': STANDARD_NAME,
        'units': UNITS,
    },
]
ATTRS_NEW = [
    {},
    {},
    {
        'test': '123',
        'answer': 42,
    },
]
TEST_SAVE_1D_DATA = zip(VAR_ATTRS_NEW, ATTRS_NEW)


@pytest.mark.parametrize('var_attrs,attrs', TEST_SAVE_1D_DATA)
@mock.patch.object(io, 'iris_save', autospec=True)
@mock.patch.object(io, 'logger', autospec=True)
def test_save_1d_data(mock_logger, mock_save, var_attrs, attrs):
    """Test saving of 1 dimensional data."""
    coord_name = 'inclination'
    data = [
        np.ma.masked_invalid([1.0, np.nan, -1.0]),
        np.arange(2.0) + 100.0,
        np.ma.masked_invalid([33.0, 22.0, np.nan, np.nan, -77.0]),
    ]
    coords = [
        iris.coords.DimCoord(np.arange(3.0) - 3.0, long_name=coord_name),
        iris.coords.DimCoord(np.arange(2.0) + 2.0, long_name=coord_name),
        iris.coords.DimCoord(np.array([-7.0, -3.0, -2.71, 3.0, 314.15]),
                             long_name=coord_name),
    ]
    cubes = OrderedDict([
        ('model1',
         iris.cube.Cube(data[0],
                        var_name='xy',
                        units='kg',
                        attributes={'hi': '!'},
                        dim_coords_and_dims=[(coords[0], 0)])),
        ('model2',
         iris.cube.Cube(data[1],
                        var_name='zr',
                        units='1',
                        attributes={},
                        dim_coords_and_dims=[(coords[1], 0)])),
        ('model3',
         iris.cube.Cube(data[2],
                        var_name='wa',
                        units='unknown',
                        attributes={'very': 'long cube'},
                        dim_coords_and_dims=[(coords[2], 0)])),
    ])
    dataset_dim = iris.coords.AuxCoord(list(cubes.keys()), long_name='dataset')
    dim_1 = coords[0].copy([-7.0, -3.0, -2.71, -2.0, -1.0, 2.0, 3.0, 314.15])
    output_data = np.ma.masked_invalid(
        [[np.nan, 1.0, np.nan, np.nan, -1.0, np.nan, np.nan, np.nan],
         [np.nan, np.nan, np.nan, np.nan, np.nan, 100.0, 101.0, np.nan],
         [33.0, 22.0, np.nan, np.nan, np.nan, np.nan, np.nan, -77.0]])
    output_dims = [(dataset_dim, 0), (dim_1, 1)]

    # Without cubes
    with pytest.raises(ValueError):
        io.save_1d_data({}, PATH, coord_name, var_attrs, attrs)
    mock_logger.error.assert_not_called()
    assert not mock_save.called
    mock_logger.reset_mock()
    mock_save.reset_mock()

    # With cubes
    if 'units' not in var_attrs:
        with pytest.raises(ValueError):
            io.save_1d_data(cubes, PATH, coord_name, var_attrs, attrs)
        mock_logger.error.assert_called_once()
        assert not mock_save.called
        return
    io.save_1d_data(cubes, PATH, coord_name, var_attrs, attrs)
    iris_var_attrs = deepcopy(var_attrs)
    iris_var_attrs['var_name'] = iris_var_attrs.pop('short_name')
    new_cube = iris.cube.Cube(output_data,
                              aux_coords_and_dims=output_dims,
                              attributes=attrs,
                              **iris_var_attrs)
    mock_logger.error.assert_not_called()
    assert mock_save.call_args_list == [mock.call(new_cube, PATH)]


CUBELIST = [
    iris.cube.Cube(1),
    iris.cube.Cube(2, attributes={
        'filename': 'a',
        'x': 'y',
    }),
]
CUBELIST_OUT = [
    iris.cube.Cube(1, attributes={'filename': PATH}),
    iris.cube.Cube(2, attributes={
        'filename': PATH,
        'x': 'y',
    }),
]
CUBES_TO_SAVE = [
    (iris.cube.Cube(0), iris.cube.Cube(0, attributes={'filename': PATH})),
    (CUBELIST, CUBELIST_OUT),
    (iris.cube.CubeList(CUBELIST), iris.cube.CubeList(CUBELIST_OUT)),
]


@pytest.mark.parametrize('source,output', CUBES_TO_SAVE)
@mock.patch('esmvaltool.diag_scripts.shared.io.iris.save', autospec=True)
@mock.patch.object(io, 'logger', autospec=True)
def test_iris_save(mock_logger, mock_save, source, output):
    """Test iris save function."""
    io.iris_save(source, PATH)
    assert mock_save.call_args_list == [mock.call(output, PATH)]
    mock_logger.info.assert_called_once()


AUX_COORDS = [
    None,
    None,
    iris.coords.AuxCoord([2, 3, 5], long_name='Primes!'),
]
TEST_SAVE_SCALAR_DATA = zip(VAR_ATTRS_NEW, ATTRS_NEW, AUX_COORDS)


@pytest.mark.parametrize('var_attrs,attrs,aux_coord', TEST_SAVE_SCALAR_DATA)
@mock.patch.object(io, 'iris_save', autospec=True)
@mock.patch.object(io, 'logger', autospec=True)
def test_save_scalar_data(mock_logger, mock_save, var_attrs, attrs, aux_coord):
    """Test saving of scalar data."""
    data = OrderedDict([
        ('model1', np.nan),
        ('model2', 1.0),
        ('model3', 3.14),
    ])
    dataset_dim = iris.coords.AuxCoord(list(data.keys()), long_name='dataset')
    output_data = np.ma.masked_invalid([np.nan, 1.0, 3.14])

    # Without data
    with pytest.raises(ValueError):
        io.save_scalar_data({}, PATH, var_attrs)
    mock_logger.error.assert_not_called()
    assert not mock_save.called
    mock_logger.reset_mock()
    mock_save.reset_mock()

    # With data
    if 'units' not in var_attrs:
        with pytest.raises(ValueError):
            io.save_scalar_data(data, PATH, var_attrs, aux_coord, attrs)
        mock_logger.error.assert_called_once()
        assert not mock_save.called
        return
    io.save_scalar_data(data, PATH, var_attrs, aux_coord, attrs)
    iris_var_attrs = deepcopy(var_attrs)
    iris_var_attrs['var_name'] = iris_var_attrs.pop('short_name')
    new_cube = iris.cube.Cube(output_data,
                              aux_coords_and_dims=[(dataset_dim, 0)],
                              attributes=attrs,
                              **iris_var_attrs)
    if aux_coord is not None:
        new_cube.add_aux_coord(aux_coord, 0)
    mock_logger.error.assert_not_called()
    assert mock_save.call_args_list == [mock.call(new_cube, PATH)]
