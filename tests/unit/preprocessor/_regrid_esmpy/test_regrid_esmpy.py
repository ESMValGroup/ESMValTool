"""Unit tests for the esmvaltool.preprocessor._regrid_esmpy module."""
import cf_units
import iris
import mock
import numpy as np
from iris.exceptions import CoordinateNotFoundError

import tests
from esmvaltool.preprocessor._regrid_esmpy import (
    build_regridder, build_regridder_2d, coords_iris_to_esmpy,
    cube_to_empty_field, get_grid, get_grid_representant,
    get_grid_representants, get_representant, is_lon_circular, regrid)


def identity(*args, **kwargs):
    """Return args, acting as identity for mocking functions."""
    # pylint: disable=unused-argument
    # Here, kwargs will be ignored.
    if len(args) == 1:
        return args[0]
    return args


def mock_cube_to_empty_field(cube):
    """Return associated field for mock cube."""
    return cube.field


class MockGrid(mock.MagicMock):
    """Mock ESMF grid."""

    get_coords = mock.Mock(return_value=mock.MagicMock())
    add_coords = mock.Mock()
    add_item = mock.Mock()
    get_item = mock.Mock(return_value=mock.MagicMock())


class MockGridItem(mock.Mock):
    """Mock ESMF enum for grid items."""

    MASK = mock.sentinel.gi_mask


class MockRegridMethod(mock.Mock):
    """Mock ESMF enum for regridding methods."""

    BILINEAR = mock.sentinel.rm_bilinear
    CONSERVE = mock.sentinel.rm_conserve
    NEAREST_STOD = mock.sentinel.rm_nearest_stod


class MockStaggerLoc(mock.Mock):
    """Mock ESMF enum for stagger locations."""

    CENTER = mock.sentinel.sl_center
    CORNER = mock.sentinel.sl_corner


class MockUnmappedAction(mock.Mock):
    """Mock ESMF enum for unmapped actions."""

    IGNORE = mock.sentinel.ua_ignore


ESMF_REGRID_METHODS = {
    'linear': MockRegridMethod.BILINEAR,
    'area_weighted': MockRegridMethod.CONSERVE,
    'nearest': MockRegridMethod.NEAREST_STOD,
}

MASK_REGRIDDING_MASK_VALUE = {
    mock.sentinel.rm_bilinear: np.array([1]),
    mock.sentinel.rm_conserve: np.array([1]),
    mock.sentinel.rm_nearest_stod: np.array([]),
}


@mock.patch('esmvaltool.preprocessor._regrid_esmpy.MASK_REGRIDDING_MASK_VALUE',
            MASK_REGRIDDING_MASK_VALUE)
@mock.patch('esmvaltool.preprocessor._regrid_esmpy.ESMF_REGRID_METHODS',
            ESMF_REGRID_METHODS)
@mock.patch('ESMF.Manager', mock.Mock)
@mock.patch('ESMF.GridItem', MockGridItem)
@mock.patch('ESMF.RegridMethod', MockRegridMethod)
@mock.patch('ESMF.StaggerLoc', MockStaggerLoc)
@mock.patch('ESMF.UnmappedAction', MockUnmappedAction)
class TestHelpers(tests.Test):
    """Unit tests for helper functions."""

    # pylint: disable=too-many-instance-attributes, too-many-public-methods
    def setUp(self):
        """Set up fixtures."""
        # pylint: disable=too-many-locals
        lat_1d_pre_bounds = np.linspace(-90, 90, 5)
        lat_1d_bounds = np.stack(
            [lat_1d_pre_bounds[:-1], lat_1d_pre_bounds[1:]], axis=1)
        lat_1d_points = lat_1d_bounds.mean(axis=1)
        lon_1d_pre_bounds = np.linspace(0, 360, 5)
        lon_1d_bounds = np.stack(
            [lon_1d_pre_bounds[:-1], lon_1d_pre_bounds[1:]], axis=1)
        lon_1d_points = lon_1d_bounds.mean(axis=1)
        lon_2d_points, lat_2d_points = np.meshgrid(lon_1d_points,
                                                   lat_1d_points)
        (lon_2d_pre_bounds, lat_2d_pre_bounds) = np.meshgrid(
            lon_1d_pre_bounds, lat_1d_pre_bounds)
        lat_2d_bounds = np.stack([
            lat_2d_pre_bounds[:-1, :-1], lat_2d_pre_bounds[:-1, 1:],
            lat_2d_pre_bounds[1:, 1:], lat_2d_pre_bounds[1:, :-1]
        ],
                                 axis=2)
        lon_2d_bounds = np.stack([
            lon_2d_pre_bounds[:-1, :-1], lon_2d_pre_bounds[:-1, 1:],
            lon_2d_pre_bounds[1:, 1:], lon_2d_pre_bounds[1:, :-1]
        ],
                                 axis=2)
        self.lat_1d = mock.Mock(
            iris.coords.DimCoord,
            standard_name='latitude',
            long_name='latitude',
            ndim=1,
            points=lat_1d_points,
            bounds=lat_1d_bounds,
            has_bounds=mock.Mock(return_value=True))
        self.lat_1d_no_bounds = mock.Mock(
            iris.coords.DimCoord,
            standard_name='latitude',
            ndim=1,
            points=lat_1d_points,
            has_bounds=mock.Mock(return_value=False),
            bounds=lat_1d_bounds,
            guess_bounds=mock.Mock())
        self.lon_1d = mock.Mock(
            iris.coords.DimCoord,
            standard_name='longitude',
            long_name='longitude',
            ndim=1,
            points=lon_1d_points,
            bounds=lon_1d_bounds,
            has_bounds=mock.Mock(return_value=True),
            circular=True)
        self.lon_1d_aux = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='longitude',
            long_name='longitude',
            ndim=1,
            shape=lon_1d_points.shape,
            points=lon_1d_points,
            bounds=lon_1d_bounds,
            has_bounds=mock.Mock(return_value=True))
        self.lat_2d = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='latitude',
            long_name='latitude',
            ndim=2,
            points=lat_2d_points,
            bounds=lat_2d_bounds,
            has_bounds=mock.Mock(return_value=True))
        self.lon_2d = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='longitude',
            long_name='longitude',
            ndim=2,
            points=lon_2d_points,
            bounds=lon_2d_bounds,
            has_bounds=mock.Mock(return_value=True))
        self.lon_2d_non_circular = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='longitude',
            ndim=2,
            points=lon_2d_points[:, 1:-1],
            bounds=lon_2d_bounds[:, 1:-1],
            has_bounds=mock.Mock(return_value=True))
        self.lat_3d = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='latitude',
            long_name='latitude',
            ndim=3)
        self.lon_3d = mock.Mock(
            iris.coords.AuxCoord,
            standard_name='longitude',
            long_name='longitude',
            ndim=3)
        depth_pre_bounds = np.linspace(0, 5000, 5)
        depth_bounds = np.stack([depth_pre_bounds[:-1], depth_pre_bounds[1:]],
                                axis=1)
        depth_points = depth_bounds.mean(axis=1)
        self.depth = mock.Mock(
            iris.coords.DimCoord,
            standard_name='depth',
            long_name='depth',
            ndim=1,
            shape=depth_points.shape,
            points=depth_points,
            bounds=depth_bounds,
            has_bounds=mock.Mock(return_value=True))
        data_shape = lon_2d_points.shape
        raw_data = np.arange(np.prod(data_shape)).reshape(data_shape)
        mask = np.zeros(data_shape)
        mask[:data_shape[0] // 2] = True
        self.data = np.ma.masked_array(raw_data, mask)
        self.data_3d = np.repeat(
            self.data[..., np.newaxis], depth_points.shape[0], axis=-1)
        self.expected_esmpy_lat = np.array([[-67.5, -22.5, 22.5, 67.5],
                                            [-67.5, -22.5, 22.5, 67.5],
                                            [-67.5, -22.5, 22.5, 67.5],
                                            [-67.5, -22.5, 22.5, 67.5]])
        self.expected_esmpy_lon = np.array([[45., 45., 45., 45.],
                                            [135., 135., 135., 135.],
                                            [225., 225., 225., 225.],
                                            [315., 315., 315., 315.]])
        self.expected_esmpy_lat_corners = np.array([[-90., -45., 0., 45., 90.],
                                                    [-90., -45., 0., 45., 90.],
                                                    [-90., -45., 0., 45., 90.],
                                                    [-90., -45., 0., 45., 90.],
                                                    [-90., -45., 0., 45.,
                                                     90.]])
        self.expected_esmpy_lon_corners = np.array(
            [[0., 0., 0., 0., 0.], [90., 90., 90., 90., 90.],
             [180., 180., 180., 180., 180.], [270., 270., 270., 270., 270.],
             [360., 360., 360., 360., 360.]])
        self.coords = {
            'latitude': self.lat_2d,
            'longitude': self.lon_2d,
            'depth': self.depth
        }
        self.coord_dims = {
            'latitude': (0, 1),
            'longitude': (0, 1),
            self.lat_2d: (0, 1),
            self.lon_2d: (0, 1),
        }

        def coord(name=None, axis=None):
            """Return selected coordinate for mock cube."""
            if axis == 'Z':
                raise CoordinateNotFoundError()
            return self.coords[name]

        def coords(dim_coords=None):
            """Return coordinates for mock cube."""
            if dim_coords:
                return []
            return list(self.coords.values())

        self.cube = mock.Mock(
            spec=iris.cube.Cube,
            dtype=np.float32,
            long_name='longname',
            ndim=2,
            shape=self.data.shape,
            data=self.data,
            coord=coord,
            coord_dims=lambda name: self.coord_dims[name],
            coords=coords,
        )
        self.cube.__getitem__ = mock.Mock(return_value=self.cube)
        self.unmasked_cube = mock.Mock(
            spec=iris.cube.Cube,
            dtype=np.float32,
            long_name='longname',
        )
        self.coord_dims_3d = {
            'latitude': (1, 2),
            'longitude': (1, 2),
            self.lat_2d: (1, 2),
            self.lon_2d: (1, 2),
            'depth': (0, ),
            self.depth: (0, ),
        }

        def coord_3d(name=None, dimensions=None, dim_coords=None, axis=None):
            """Return coord for 3d mock cube."""
            # pylint: disable=unused-argument
            if axis == 'Z' or dimensions == [0]:
                return self.coords['depth']
            return self.coords[name]

        self.cube_3d = mock.Mock(
            spec=iris.cube.Cube,
            dtype=np.float32,
            standard_name=None,
            long_name='longname',
            var_name='ln',
            units=cf_units.Unit('1'),
            attributes={},
            cell_methods=[],
            ndim=3,
            shape=self.data_3d.shape,
            data=self.data_3d,
            coord=coord_3d,
            coord_dims=lambda name: self.coord_dims_3d[name],
        )
        self.cube.__getitem__ = mock.Mock(return_value=self.cube)

    def test_coords_iris_to_esmpy_mismatched_dimensions(self):
        """Test coord conversion with mismatched dimensions."""
        self.assertRaises(ValueError, coords_iris_to_esmpy, self.lat_1d,
                          self.lon_2d, True)

    def test_coords_iris_to_esmpy_invalid_dimensions(self):
        """Test coord conversion with invalid dimensions."""
        self.assertRaises(NotImplementedError, coords_iris_to_esmpy,
                          self.lat_3d, self.lon_3d, True)

    def test_coords_iris_to_esmpy_call_guess_bounds(self):
        """Test coord conversion with missing bounds."""
        coords_iris_to_esmpy(self.lat_1d_no_bounds, self.lon_1d, True)
        self.lat_1d_no_bounds.guess_bounds.assert_called_once()

    def test_coords_iris_to_esmpy_1d_circular(self):
        """Test coord conversion with 1d coords and circular longitudes."""
        (esmpy_lat, esmpy_lon,
         esmpy_lat_corners, esmpy_lon_corners) = coords_iris_to_esmpy(
             self.lat_1d, self.lon_1d, True)
        self.assertArrayEqual(esmpy_lat, self.expected_esmpy_lat)
        self.assertArrayEqual(esmpy_lon, self.expected_esmpy_lon)
        self.assertArrayEqual(esmpy_lat_corners,
                              self.expected_esmpy_lat_corners[:-1])
        self.assertArrayEqual(esmpy_lon_corners,
                              self.expected_esmpy_lon_corners[:-1])

    def test_coords_iris_to_esmpy_1d_non_circular(self):
        """Test coord conversion with 1d coords and non circular longitudes."""
        (esmpy_lat, esmpy_lon,
         esmpy_lat_corners, esmpy_lon_corners) = coords_iris_to_esmpy(
             self.lat_1d, self.lon_1d, False)
        self.assertArrayEqual(esmpy_lat, self.expected_esmpy_lat)
        self.assertArrayEqual(esmpy_lon, self.expected_esmpy_lon)
        self.assertArrayEqual(esmpy_lat_corners,
                              self.expected_esmpy_lat_corners)
        self.assertArrayEqual(esmpy_lon_corners,
                              self.expected_esmpy_lon_corners)

    def test_coords_iris_to_esmpy_2d_circular(self):
        """Test coord conversion with 2d coords and circular longitudes."""
        (esmpy_lat, esmpy_lon,
         esmpy_lat_corners, esmpy_lon_corners) = coords_iris_to_esmpy(
             self.lat_2d, self.lon_2d, True)
        self.assertArrayEqual(esmpy_lat, self.expected_esmpy_lat)
        self.assertArrayEqual(esmpy_lon, self.expected_esmpy_lon)
        self.assertArrayEqual(esmpy_lat_corners,
                              self.expected_esmpy_lat_corners[:-1])
        self.assertArrayEqual(esmpy_lon_corners,
                              self.expected_esmpy_lon_corners[:-1])

    def test_coords_iris_to_esmpy_2d_non_circular(self):
        """Test coord conversion with 2d coords and non circular longitudes."""
        (esmpy_lat, esmpy_lon,
         esmpy_lat_corners, esmpy_lon_corners) = coords_iris_to_esmpy(
             self.lat_2d, self.lon_2d, False)
        self.assertArrayEqual(esmpy_lat, self.expected_esmpy_lat)
        self.assertArrayEqual(esmpy_lon, self.expected_esmpy_lon)
        self.assertArrayEqual(esmpy_lat_corners,
                              self.expected_esmpy_lat_corners)
        self.assertArrayEqual(esmpy_lon_corners,
                              self.expected_esmpy_lon_corners)

    def test_get_grid_circular(self):
        """Test building of ESMF grid from iris cube circular longitude."""
        expected_get_coords_calls = [
            mock.call(0),
            mock.call(1),
            mock.call(0, staggerloc=mock.sentinel.sl_corner),
            mock.call(1, staggerloc=mock.sentinel.sl_corner),
        ]
        with mock.patch('ESMF.Grid', MockGrid) as mg:
            mg.get_coords.reset_mock()
            mg.add_coords.reset_mock()
            mg.add_item.reset_mock()
            get_grid(self.expected_esmpy_lat, self.expected_esmpy_lon,
                     self.expected_esmpy_lat_corners[:-1],
                     self.expected_esmpy_lon_corners[:-1], True)
            mg.get_coords.assert_has_calls(expected_get_coords_calls)
            mg.add_coords.assert_called_once_with([mock.sentinel.sl_corner])
            mg.add_item.assert_called_once_with(mock.sentinel.gi_mask,
                                                mock.sentinel.sl_center)

    def test_get_grid_non_circular(self):
        """Test building of ESMF grid from iris cube non circular longitude."""
        expected_get_coords_calls = [
            mock.call(0),
            mock.call(1),
            mock.call(0, staggerloc=mock.sentinel.sl_corner),
            mock.call(1, staggerloc=mock.sentinel.sl_corner),
        ]
        with mock.patch('ESMF.Grid', MockGrid) as mg:
            mg.get_coords.reset_mock()
            mg.add_coords.reset_mock()
            mg.add_item.reset_mock()
            get_grid(self.expected_esmpy_lat, self.expected_esmpy_lon,
                     self.expected_esmpy_lat_corners,
                     self.expected_esmpy_lon_corners, False)
            mg.get_coords.assert_has_calls(expected_get_coords_calls)
            mg.add_coords.assert_called_once_with([mock.sentinel.sl_corner])
            mg.add_item.assert_called_once_with(mock.sentinel.gi_mask,
                                                mock.sentinel.sl_center)

    def test_is_lon_circular_dim_coords_true(self):
        """Test detection of circular longitudes 1d dim coords."""
        is_circ = is_lon_circular(self.lon_1d)
        self.assertTrue(is_circ)

    def test_is_lon_circular_dim_coords_false(self):
        """Test detection of non circular longitudes 1d dim coords."""
        self.lon_1d.circular = False
        is_circ = is_lon_circular(self.lon_1d)
        self.assertFalse(is_circ)

    def test_is_lon_circular_1d_aux_coords(self):
        """Test detection of circular longitudes 1d aux coords."""
        is_circ = is_lon_circular(self.lon_1d_aux)
        self.assertTrue(is_circ)

    def test_is_lon_circular_invalid_dimension(self):
        """Test detection of circular longitudes, invalid coordinates."""
        self.assertRaises(NotImplementedError, is_lon_circular, self.lon_3d)

    def test_is_lon_circular_invalid_argument(self):
        """Test detection of circular longitudes, invalid argument."""
        self.assertRaises(ValueError, is_lon_circular, None)

    def test_is_lon_circular_2d_aux_coords(self):
        """Test detection of circular longitudes 2d aux coords."""
        is_circ = is_lon_circular(self.lon_2d)
        self.assertTrue(is_circ)

    def test_is_lon_circular_2d_aux_coords_non_circ(self):
        """Test detection of non circular longitudes 2d aux coords."""
        is_circ = is_lon_circular(self.lon_2d_non_circular)
        self.assertFalse(is_circ)

    @mock.patch('ESMF.Grid', MockGrid)
    @mock.patch('ESMF.Field')
    def test_cube_to_empty_field(self, mock_field):
        """Test building of empty field from iris cube."""
        field = cube_to_empty_field(self.cube)
        self.assertEqual(mock_field.return_value, field)
        mock_field.assert_called_once()
        ckwargs = mock_field.call_args[1]
        self.assertEqual('longname', ckwargs['name'])
        self.assertEqual(mock.sentinel.sl_center, ckwargs['staggerloc'])

    def test_get_representant(self):
        """Test extraction of horizontal representant from iris cube."""
        horizontal_slice = ['latitude', 'longitude']
        get_representant(self.cube, horizontal_slice)
        self.cube.__getitem__.assert_called_once_with((slice(None, None, None),
                                                       slice(None, None,
                                                             None)))

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.cube_to_empty_field',
                mock_cube_to_empty_field)
    @mock.patch('ESMF.Regrid')
    def test_build_regridder_2d_unmasked_data(self, mock_regrid):
        """Test building of 2d regridder for unmasked data."""
        self.cube.data = self.cube.data.data
        self.cube.field = mock.Mock()
        mock.sentinel.dst_rep.field = mock.Mock()
        build_regridder_2d(self.cube, mock.sentinel.dst_rep,
                           mock.sentinel.regrid_method, .99)
        expected_kwargs = {
            'src_mask_values': np.array([1]),
            'dst_mask_values': np.array([1]),
            'regrid_method': mock.sentinel.regrid_method,
            'srcfield': self.cube.field,
            'dstfield': mock.sentinel.dst_rep.field,
            'unmapped_action': mock.sentinel.ua_ignore,
            'ignore_degenerate': True,
        }
        mock_regrid.assert_called_once_with(**expected_kwargs)

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.cube_to_empty_field',
                mock_cube_to_empty_field)
    @mock.patch('ESMF.Regrid')
    def test_build_regridder_2d_masked_data(self, mock_regrid):
        """Test building of 2d regridder for masked data."""
        mock_regrid.return_value = mock.Mock(
            return_value=mock.Mock(data=self.data.T))
        regrid_method = mock.sentinel.rm_bilinear
        src_rep = mock.MagicMock(data=self.data)
        dst_rep = mock.MagicMock()
        src_rep.field = mock.MagicMock(data=self.data.copy())
        dst_rep.field = mock.MagicMock()
        build_regridder_2d(src_rep, dst_rep, regrid_method, .99)
        expected_calls = [
            mock.call(
                src_mask_values=np.array([]),
                dst_mask_values=np.array([]),
                srcfield=src_rep.field,
                dstfield=dst_rep.field,
                unmapped_action=mock.sentinel.ua_ignore,
                ignore_degenerate=True,
                regrid_method=regrid_method),
            mock.call(
                src_mask_values=np.array([1]),
                dst_mask_values=np.array([1]),
                regrid_method=regrid_method,
                srcfield=src_rep.field,
                dstfield=dst_rep.field,
                unmapped_action=mock.sentinel.ua_ignore,
                ignore_degenerate=True),
        ]
        kwargs = mock_regrid.call_args_list[0][-1]
        expected_kwargs = expected_calls[0][-1]
        self.assertEqual(expected_kwargs.keys(), kwargs.keys())
        array_keys = set(['src_mask_values', 'dst_mask_values'])
        for key in kwargs.keys():
            if key in array_keys:
                self.assertTrue((expected_kwargs[key] == kwargs[key]).all())
            else:
                self.assertEqual(expected_kwargs[key], kwargs[key])
        self.assertTrue(mock_regrid.call_args_list[1] == expected_calls[1])

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.cube_to_empty_field',
                mock_cube_to_empty_field)
    @mock.patch('ESMF.Regrid')
    def test_regridder_2d_unmasked_data(self, mock_regrid):
        """Test regridder for unmasked 2d data."""
        field_regridder = mock.Mock(return_value=mock.Mock(data=self.data.T))
        mock_regrid.return_value = field_regridder
        regrid_method = mock.sentinel.rm_bilinear
        src_rep = mock.MagicMock(data=self.data, dtype=np.float32)
        dst_rep = mock.MagicMock(shape=(4, 4))
        regridder = build_regridder_2d(src_rep, dst_rep, regrid_method, .99)
        field_regridder.reset_mock()
        regridder(src_rep)
        field_regridder.assert_called_once_with(src_rep.field, dst_rep.field)

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.cube_to_empty_field',
                mock_cube_to_empty_field)
    @mock.patch('ESMF.Regrid')
    def test_regridder_2d_masked_data(self, mock_regrid):
        """Test regridder for masked 2d data."""
        field_regridder = mock.Mock(return_value=mock.Mock(data=self.data.T))
        mock_regrid.return_value = field_regridder
        regrid_method = mock.sentinel.rm_bilinear
        src_rep = mock.MagicMock(data=self.data)
        dst_rep = mock.MagicMock(shape=(4, 4))
        regridder = build_regridder_2d(src_rep, dst_rep, regrid_method, .99)
        field_regridder.reset_mock()
        regridder(self.cube)
        field_regridder.assert_called_once_with(src_rep.field, dst_rep.field)

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.build_regridder_3d')
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.build_regridder_2d')
    def test_build_regridder_2(self, mock_regridder_2d, mock_regridder_3d):
        """Test build regridder for 2d data."""
        # pylint: disable=no-self-use
        src_rep = mock.Mock(ndim=2)
        dst_rep = mock.Mock(ndim=2)
        build_regridder(src_rep, dst_rep, 'nearest')
        mock_regridder_2d.assert_called_once_with(
            src_rep, dst_rep, mock.sentinel.rm_nearest_stod, .99)
        mock_regridder_3d.assert_not_called()

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.build_regridder_3d')
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.build_regridder_2d')
    def test_build_regridder_3(self, mock_regridder_2d, mock_regridder_3d):
        """Test build regridder for 3d data."""
        # pylint: disable=no-self-use
        src_rep = mock.Mock(ndim=3)
        dst_rep = mock.Mock(ndim=3)
        build_regridder(src_rep, dst_rep, 'nearest')
        mock_regridder_3d.assert_called_once_with(
            src_rep, dst_rep, mock.sentinel.rm_nearest_stod, .99)
        mock_regridder_2d.assert_not_called()

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_representant')
    def test_get_grid_representant_2d(self, mock_get_representant):
        """Test extraction of 2d grid representant from 2 spatial d cube."""
        mock_get_representant.return_value = mock.sentinel.ret
        ret = get_grid_representant(self.cube)
        self.assertEqual(mock.sentinel.ret, ret)
        mock_get_representant.assert_called_once_with(
            self.cube, ['latitude', 'longitude'])

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_representant')
    def test_get_grid_representant_2d_horiz_only(self, mock_get_representant):
        """Test extraction of forced 2d grid representant from 2d cube."""
        mock_get_representant.return_value = mock.sentinel.ret
        ret = get_grid_representant(self.cube, True)
        self.assertEqual(mock.sentinel.ret, ret)
        mock_get_representant.assert_called_once_with(
            self.cube, ['latitude', 'longitude'])

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_representant')
    def test_get_grid_representant_3d(self, mock_get_representant):
        """Test extraction of 3d grid representant from 3 spatial d cube."""
        mock_get_representant.return_value = mock.sentinel.ret
        ret = get_grid_representant(self.cube_3d)
        self.assertEqual(mock.sentinel.ret, ret)
        mock_get_representant.assert_called_once_with(
            self.cube_3d, [self.depth, 'latitude', 'longitude'])

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_representant')
    def test_get_grid_representant_3d_horiz_only(self, mock_get_representant):
        """Test extraction of 2d grid representant from 3 spatial d cube."""
        mock_get_representant.return_value = mock.sentinel.ret
        ret = get_grid_representant(self.cube_3d, True)
        self.assertEqual(mock.sentinel.ret, ret)
        mock_get_representant.assert_called_once_with(
            self.cube_3d, ['latitude', 'longitude'])

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_grid_representant',
                mock.Mock(side_effect=identity))
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_empty_data')
    @mock.patch('iris.cube.Cube')
    def test_get_grid_representants_3d_src(self, mock_cube,
                                           mock_get_empty_data):
        """Test extraction of grid representants from 3 spatial d cube."""
        src = self.cube_3d
        mock_get_empty_data.return_value = mock.sentinel.empty_data
        src_rep = get_grid_representants(src, self.cube)[0]
        self.assertEqual(src, src_rep)
        mock_cube.assert_called_once_with(
            data=mock.sentinel.empty_data,
            standard_name=src.standard_name,
            long_name=src.long_name,
            var_name=src.var_name,
            units=src.units,
            attributes=src.attributes,
            cell_methods=src.cell_methods,
            dim_coords_and_dims=[(self.depth, 0)],
        )

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_grid_representant',
                mock.Mock(side_effect=identity))
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_empty_data')
    @mock.patch('iris.cube.Cube')
    def test_get_grid_representants_2d_src(self, mock_cube,
                                           mock_get_empty_data):
        """Test extraction of grid representants from 2 spatial d cube."""
        src = self.cube
        mock_get_empty_data.return_value = mock.sentinel.empty_data
        src_rep = get_grid_representants(src, self.cube)[0]
        self.assertEqual(src, src_rep)
        mock_cube.assert_called_once_with(
            data=mock.sentinel.empty_data,
            standard_name=src.standard_name,
            long_name=src.long_name,
            var_name=src.var_name,
            units=src.units,
            attributes=src.attributes,
            cell_methods=src.cell_methods,
            dim_coords_and_dims=[],
        )

    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.map_slices')
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.build_regridder')
    @mock.patch('esmvaltool.preprocessor._regrid_esmpy.get_grid_representants',
                mock.Mock(side_effect=identity))
    def test_regrid(self, mock_build_regridder, mock_map_slices):
        """Test full regrid method."""
        mock_build_regridder.return_value = mock.sentinel.regridder
        mock_map_slices.return_value = mock.sentinel.regridded
        regrid(self.cube_3d, self.cube)
        mock_build_regridder.assert_called_once_with(self.cube_3d, self.cube,
                                                     'linear')
        mock_map_slices.assert_called_once_with(
            self.cube_3d, mock.sentinel.regridder, self.cube_3d, self.cube)
