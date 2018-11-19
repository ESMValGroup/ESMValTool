from __future__ import absolute_import, division, print_function

import cf_units
import iris
import mock
import numpy as np

import tests
from esmvaltool.preprocessor._mapping import (get_empty_data,
                                              map_slices,
                                              ref_to_dims_index)


class TestHelpers(tests.Test):
    def setUp(self):
        self.coord_system = mock.Mock(return_value=None)
        self.scalar_coord = mock.sentinel.scalar_coord
        self.scalar_coord.name = lambda: 'scalar_coord'
        self.coord = mock.sentinel.coord
        self.coords = mock.Mock(return_value=[self.scalar_coord, self.coord])

        def coord(name_or_coord):
            if name_or_coord == 'coord':
                return self.coord
            elif name_or_coord == 'scalar_coord':
                return self.scalar_coord
            else:
                raise iris.exceptions.CoordinateNotFoundError('')

        def coord_dims(coord):
            if coord == self.coord:
                return [0]
            elif coord == self.scalar_coord:
                return []
            else:
                raise iris.exceptions.CoordinateNotFoundError('')
        self.cube = mock.Mock(
            spec=iris.cube.Cube,
            coord_system=self.coord_system,
            coords=self.coords,
            coord=coord,
            coord_dims=coord_dims,
            ndim=4,
        )

    def test_get_empty_data(self):
        shape = (3, 3)
        data = get_empty_data(shape)
        self.assertIsInstance(data, np.ma.MaskedArray)
        self.assertEqual(data.shape, shape)

    def test_ref_to_dims_index__int(self):
        ref_to_dims_index(self.cube, 0)

    def test_ref_to_dims_index__invalid_int(self):
        self.assertRaises(ValueError,
                          ref_to_dims_index, self.cube, -1)
        self.assertRaises(ValueError,
                          ref_to_dims_index, self.cube, 100)

    def test_ref_to_dims_index__scalar_coord(self):
        self.assertRaises(ValueError,
                          ref_to_dims_index, self.cube, 'scalar_coord')

    def test_ref_to_dims_index__valid_coordinate_name(self):
        ref_to_dims_index(self.cube, 'coord')

    def test_ref_to_dims_index__invalid_coordinate_name(self):
        self.assertRaises(iris.exceptions.CoordinateNotFoundError,
                          ref_to_dims_index, self.cube, 'test')

    def test_ref_to_dims_index__invalid_type(self):
        self.assertRaises(ValueError,
                          ref_to_dims_index,
                          self.cube, mock.sentinel.something)


class Test(tests.Test):
    def setUp(self):
        self.coord_system = mock.Mock(return_value=None)
        self.scalar_coord = mock.sentinel.scalar_coord
        self.scalar_coord.name = lambda: 'scalar_coord'
        self.time = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='time',
            long_name='time',
            shape=(3,),
        )
        self.z = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='height',
            long_name='height',
            shape=(4,),
        )
        self.src_latitude = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='latitude',
            long_name='latitude',
            shape=(5,),
            points=np.array([1.1, 2.2, 3.3, 4.4, 5.5]),
        )
        self.src_longitude = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='longitude',
            long_name='longitude',
            shape=(6,),
            points=np.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6]),
        )
        self.dst_latitude = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='latitude',
            long_name='latitude',
            shape=(2,),
            points=np.array([1.1, 2.2]),
        )
        self.dst_longitude = mock.Mock(
            spec=iris.coords.DimCoord,
            standard_name='longitude',
            long_name='longitude',
            shape=(2,),
            points=np.array([1.1, 2.2]),
        )

        def src_coord(name_or_coord):
            if name_or_coord in ['latitude', self.src_latitude]:
                return self.src_latitude
            elif name_or_coord in ['longitude', self.src_longitude]:
                return self.src_longitude
            elif name_or_coord == 'scalar_coord':
                return self.scalar_coord
            else:
                raise iris.exceptions.CoordinateNotFoundError('')

        def coord_dims(coord):
            if coord in [self.time, self.dst_latitude]:
                return [0]
            elif coord in [self.z, self.dst_longitude]:
                return [1]
            elif coord in [self.src_latitude]:
                return [2]
            elif coord in [self.src_longitude]:
                return [3]
            elif coord == self.scalar_coord:
                return []
            else:
                raise iris.exceptions.CoordinateNotFoundError('')

        def src_coords(*args, contains_dimension=None, dim_coords=None):
            dim_coords_list = [self.time,
                               self.z,
                               self.src_latitude,
                               self.src_longitude]
            if contains_dimension is not None:
                return [dim_coords_list[contains_dimension]]
            if dim_coords:
                return dim_coords_list
            return [self.scalar_coord] + dim_coords_list

        def src_repr_coords(*args, **kwargs):
            dim_coords = [self.src_latitude, self.src_longitude]
            if kwargs.get('dim_coords', False):
                return dim_coords
            if 'contains_dimension' in kwargs:
                return dim_coords
            return [self.scalar_coord] + dim_coords

        def dst_repr_coords(*args, **kwargs):
            dim_coords = [self.dst_latitude, self.dst_longitude]
            if kwargs.get('dim_coords', False):
                return dim_coords
            return [self.scalar_coord] + dim_coords
        self.src_cube = mock.Mock(
            spec=iris.cube.Cube,
            coord_system=self.coord_system,
            coords=src_coords,
            coord=src_coord,
            coord_dims=coord_dims,
            ndim=4,
            shape=(3, 4, 5, 6),
            standard_name='sea_surface_temperature',
            long_name='Sea surface temperature',
            var_name='tos',
            units=cf_units.Unit('K'),
            attributes={},
            cell_methods={},
            __getitem__=lambda a, b: mock.sentinel.src_data,
        )
        self.src_repr = mock.Mock(
            spec=iris.cube.Cube,
            coords=src_repr_coords,
            ndim=2,
        )
        self.dst_repr = mock.Mock(
            spec=iris.cube.Cube,
            coords=dst_repr_coords,
            shape=(2, 2),
        )

    @mock.patch('esmvaltool.preprocessor._mapping.get_empty_data',
                return_value=mock.sentinel.empty_data)
    @mock.patch('iris.cube.Cube')
    def test_map_slices(self, mock_cube, mock_get_empty_data):
        dst = map_slices(self.src_cube, lambda s: np.ones((2, 2)),
                         self.src_repr, self.dst_repr)
        self.assertEqual(dst, mock_cube.return_value)
        dim_coords = self.src_cube.coords(dim_coords=True)[:2] \
            + self.dst_repr.coords(dim_coords=True)
        dim_coords_and_dims = [(c, i) for i, c in enumerate(dim_coords)]
        mock_cube.assert_called_once_with(
            data=mock.sentinel.empty_data,
            standard_name=self.src_cube.standard_name,
            long_name=self.src_cube.long_name,
            var_name=self.src_cube.var_name,
            units=self.src_cube.units,
            attributes=self.src_cube.attributes,
            cell_methods=self.src_cube.cell_methods,
            dim_coords_and_dims=dim_coords_and_dims,
        )
