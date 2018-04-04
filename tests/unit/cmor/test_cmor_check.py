"""Unit tests for the CMORCheck class."""

import sys
import unittest

import iris
import iris.coord_categorisation
import iris.coords
import iris.util
import numpy
from cf_units import Unit

from esmvaltool.cmor.check import CMORCheck, CMORCheckError


class VariableInfoMock:
    """Mock for the variables defintion"""

    def __init__(self):
        self.short_name = 'short_name'
        self.standard_name = 'age_of_sea_ice'  # Iris don't accept fakes ...
        self.long_name = 'Long Name'
        self.units = 'years'  # ... nor in the units
        self.valid_min = '0'
        self.valid_max = '100'
        self.frequency = 'day'
        self.positive = 'up'

        generic_level = CoordinateInfoMock('depth')
        generic_level.generic_level = True
        generic_level.axis = 'Z'

        requested = CoordinateInfoMock('air_pressure')
        requested.requested = [str(number) for number in range(20)]

        self.coordinates = {
            'time': CoordinateInfoMock('time'),
            'lat': CoordinateInfoMock('lat'),
            'lon': CoordinateInfoMock('lon'),
            'air_pressure': requested,
            'depth': generic_level
        }


class CoordinateInfoMock:
    """Mock for the coordinates info"""

    def __init__(self, name):
        self.name = name
        self.generic_level = False

        self.axis = ""
        self.value = ""
        standard_names = {'lat': 'latitude', 'lon': 'longitude'}
        if name in standard_names:
            self.standard_name = standard_names[name]
        else:
            self.standard_name = name
        self.long_name = "Long name"
        self.out_name = self.name
        self.var_name = self.name

        units = {
            'lat': 'degrees_north',
            'lon': 'degrees_east',
            'time': 'days since 1950-01-01 00:00:00'
        }
        if name in units:
            self.units = units[name]
        else:
            self.units = "units"

        self.stored_direction = "increasing"
        self.requested = []

        valid_limits = {'lat': ('-90', '90'), 'lon': ('0', '360')}
        if name in valid_limits:
            self.valid_min = valid_limits[name][0]
            self.valid_max = valid_limits[name][1]
        else:
            self.valid_min = ""
            self.valid_max = ""


class TestCMORCheck(unittest.TestCase):
    """Test CMORCheck class"""

    def setUp(self):
        self.var_info = VariableInfoMock()
        self.cube = self.get_cube(self.var_info)

    def test_report_error(self):
        """Test report error function"""
        checker = CMORCheck(self.cube, self.var_info)
        self.assertFalse(checker.has_errors())
        checker.report_error('New error: {}', 'something failed')
        self.assertTrue(checker.has_errors())

    def test_fail_on_error(self):
        """Test exception is raised if fail_on_error is activated"""
        checker = CMORCheck(self.cube, self.var_info, fail_on_error=True)
        with self.assertRaises(CMORCheckError):
            checker.report_error('New error: {}', 'something failed')

    def test_report_warning(self):
        """Test report warning function"""
        checker = CMORCheck(self.cube, self.var_info)
        self.assertFalse(checker.has_errors())
        checker.report_warning('New error: {}', 'something failed')
        self.assertTrue(checker.has_warnings())

    def test_warning_fail_on_error(self):
        """Test report warning function with fail_on_error"""
        if sys.version_info[0] == 2:
            from StringIO import StringIO
        else:
            from io import StringIO
        checker = CMORCheck(self.cube, self.var_info, fail_on_error=True)
        stdout = sys.stdout
        sys.stdout = StringIO()
        checker.report_warning('New error: {}', 'something failed')
        output = sys.stdout.getvalue().strip()
        sys.stdout = stdout
        self.assertEqual(output, 'WARNING: New error: something failed')

    def test_check(self):
        """Test checks succeeds for a good cube"""
        self._check_cube()

    def _check_cube(self, automatic_fixes=False):
        checker = CMORCheck(
            self.cube, self.var_info, automatic_fixes=automatic_fixes)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_month_number(self):
        """Test checks succeeds for a good cube with month number"""
        iris.coord_categorisation.add_month_number(self.cube, 'time')
        self._check_cube()

    def test_check_with_day_of_month(self):
        """Test checks succeeds for a good cube with day of month"""
        iris.coord_categorisation.add_day_of_month(self.cube, 'time')
        self._check_cube()

    def test_check_with_day_of_year(self):
        """Test checks succeeds for a good cube with day of year"""
        iris.coord_categorisation.add_day_of_year(self.cube, 'time')
        self._check_cube()

    def test_check_with_year(self):
        """Test checks succeeds for a good cube with year"""
        iris.coord_categorisation.add_year(self.cube, 'time')
        self._check_cube()

    def test_check_with_unit_conversion(self):
        """Test check succeds for a good cube requiring unit converision"""
        self.cube.units = 'days'
        self._check_cube()

    def test_check_with_positive(self):
        self.var_info.positive = 'up'
        self.cube = self.get_cube(self.var_info)
        self._check_cube()

    def test_invalid_rank(self):
        """Test check fails in metadata step when rank is not correct"""
        self.cube.remove_coord('latitude')
        self._check_fails_in_metadata()

    def test_rank_with_aux_coords(self):
        """Check succeeds even if a required coordinate is an aux coord"""
        iris.util.demote_dim_coord_to_aux_coord(self.cube, 'latitude')
        self._check_cube()

    def _check_fails_in_metadata(self, automatic_fixes=False, frequency=None):
        checker = CMORCheck(
            self.cube,
            self.var_info,
            automatic_fixes=automatic_fixes,
            frequency=frequency)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_non_requested(self):
        """
        Warning if requested values are not present

        Check issue a warning if a values requested
        for a coordinate are not correct in the metadata step
        """
        coord = self.cube.coord('air_pressure')
        values = numpy.linspace(0, 40, len(coord.points))
        self._update_coordinate_values(self.cube, coord, values)
        checker = CMORCheck(self.cube, self.var_info)
        checker.check_metadata()
        self.assertTrue(checker.has_warnings())

    def test_non_increasing(self):
        """Fail in metadata if increasing coordinate is decreasing"""
        coord = self.cube.coord('latitude')
        values = numpy.linspace(coord.points[-1], coord.points[0],
                                len(coord.points))
        self._update_coordinate_values(self.cube, coord, values)
        self._check_fails_in_metadata()

    def test_non_decreasing(self):
        """Fail in metadata if decreasing coordinate is increasing"""
        self.var_info.coordinates['lat'].stored_direction = 'decreasing'
        self._check_fails_in_metadata()

    def test_non_decreasing_fix(self):
        """Check automatic fix for non decreasing coordinate"""
        self.cube.data[0, 0, 0, 0, 0] = 70
        self.var_info.coordinates['lat'].stored_direction = 'decreasing'
        self._check_cube(automatic_fixes=True)
        self._check_cube()
        index = [0, 0, 0, 0, 0]
        index[self.cube.coord_dims('latitude')[0]] = -1
        self.assertEqual(self.cube.data.item(tuple(index)), 70)
        self.assertEqual(self.cube.data[0, 0, 0, 0, 0], 50)
        cube_points = self.cube.coord('latitude').points
        reference = numpy.linspace(90, -90, 20, endpoint=True)
        for index in range(20):
            self.assertTrue(iris.util.approx_equal(cube_points[index],
                                                   reference[index]))

    def test_not_correct_lons(self):
        """Fail if longitudes are not correct in metadata step"""
        self.cube = self.cube.intersection(longitude=(-180., 180.))
        self._check_fails_in_metadata()

    def test_lons_automatic_fix(self):
        """Test automatic fixes for bad longitudes"""
        self.cube = self.cube.intersection(longitude=(-180., 180.))
        self._check_cube(automatic_fixes=True)

    def test_high_lons_automatic_fix(self):
        """Test automatic fixes for high longitudes"""
        self.cube = self.cube.intersection(longitude=(180., 520.))
        self._check_cube(automatic_fixes=True)

    def test_not_valid_min(self):
        """Fail if coordinate values below valid_min"""
        coord = self.cube.coord('latitude')
        values = numpy.linspace(coord.points[0] - 1, coord.points[-1],
                                len(coord.points))
        self._update_coordinate_values(self.cube, coord, values)
        self._check_fails_in_metadata()

    def test_not_valid_max(self):
        """Fail if coordinate values above valid_max"""
        coord = self.cube.coord('latitude')
        values = numpy.linspace(coord.points[0], coord.points[-1] + 1,
                                len(coord.points))
        self._update_coordinate_values(self.cube, coord, values)
        self._check_fails_in_metadata()

    @staticmethod
    def _update_coordinate_values(cube, coord, values):
        [dimension] = cube.coord_dims(coord)
        cube.remove_coord(coord)
        new_coord = iris.coords.DimCoord(
            values,
            standard_name=coord.standard_name,
            long_name=coord.long_name,
            var_name=coord.var_name,
            units=coord.units)
        cube.add_dim_coord(new_coord, dimension)

    def test_bad_units(self):
        """Fail if coordinates have bad units"""
        self.cube.coord('latitude').units = 'degrees_n'
        self._check_fails_in_metadata()

    def test_units_automatic_fix(self):
        """Test automatic fix for bad coordinate units"""
        self.cube.coord('latitude').units = 'degrees_n'
        self._check_cube(automatic_fixes=True)

    def test_units_automatic_fix_failed(self):
        """Test automatic fix fail for incompatible coordinate units"""
        self.cube.coord('latitude').units = 'degC'
        self._check_fails_in_metadata(automatic_fixes=True)

    def test_bad_time(self):
        """Fail if time have bad units"""
        self.cube.coord('time').units = 'days'
        self._check_fails_in_metadata()

    def test_time_automatic_fix(self):
        """Test automatic fix for time units"""
        self.cube.coord('time').units = 'days since 1950-1-1 00:00:00'
        self._check_cube(automatic_fixes=True)

    def test_time_automatic_fix_failed(self):
        """Test automatic fix fail for incompatible time units"""
        self.cube.coord('time').units = 'K'
        self._check_fails_in_metadata(automatic_fixes=True)

    def test_bad_standard_name(self):
        """Fail if coordinates have bad standard names at metadata step"""
        self.cube.coord('time').standard_name = 'region'
        self._check_fails_in_metadata()

    def test_bad_data_units(self):
        """Fail if data has bad units at metadata step"""
        self.cube.units = 'hPa'
        self._check_fails_in_metadata()

    def test_bad_data_standard_name(self):
        """Fail if data have bad standard_name at metadata step"""
        self.cube.standard_name = 'wind_speed'
        self._check_fails_in_metadata()

    def test_bad_positive(self):
        """Fail if positive value is incorrect at metadata step"""
        self.cube.attributes['positive'] = 'up'
        self.var_info.positive = 'down'
        self._check_fails_in_metadata()

    def test_bad_standard_name_genlevel(self):
        """Fail if a generic level has a bad standard name at metadata step"""
        self.cube.coord('depth').standard_name = None
        self._check_fails_in_metadata()

    def test_bad_frequency_day(self):
        """Fail at metadata if frequency (day) not matches data frequency"""
        self.cube = self.get_cube(self.var_info, frequency='mon')
        self._check_fails_in_metadata(frequency='day')

    def test_bad_frequency_subhr(self):
        """Fail at metadata if frequency (subhr) not matches data frequency"""
        self._check_fails_in_metadata(frequency='subhr')

    def test_bad_frequency_dec(self):
        """Fail at metadata if frequency (dec) not matches data frequency"""
        self._check_fails_in_metadata(frequency='d')

    def test_bad_frequency_yr(self):
        """Fail at metadata if frequency (yr) not matches data frequency"""
        self._check_fails_in_metadata(frequency='yr')

    def test_bad_frequency_mon(self):
        """Fail at metadata if frequency (mon) not matches data frequency"""
        self._check_fails_in_metadata(frequency='mon')

    def test_bad_frequency_hourly(self):
        """Fail at metadata if frequency (3hr) not matches data frequency"""
        self._check_fails_in_metadata(frequency='3hr')

    def test_frequency_not_supported(self):
        """Fail at metadata if frequency is not supported"""
        self._check_fails_in_metadata(frequency='wrong_freq')

    # For the moment, we don't have a variable definition with these values
    # to test

    def test_data_not_valid_max(self):
        """Warning if data is above valid_max in data step"""
        self.var_info.valid_max = '10000'
        self.cube.data[0] = 100000000000
        self._check_warnings_on_data()

    def test_data_not_valid_min(self):
        """Warning if data is below valid_min in data step"""
        self.var_info.valid_min = '-100'
        self.cube.data[0] = -100000000000
        self._check_warnings_on_data()

    def _check_fails_on_data(self):
        checker = CMORCheck(self.cube, self.var_info)
        checker.check_metadata()
        with self.assertRaises(CMORCheckError):
            checker.check_data()

    def _check_warnings_on_data(self):
        checker = CMORCheck(self.cube, self.var_info)
        checker.check_metadata()
        checker.check_data()
        self.assertTrue(checker.has_warnings())

    def get_cube(self,
                 var_info,
                 set_time_units="days since 1850-1-1 00:00:00",
                 frequency=None):
        """
        Create a cube based on a specification

        Parameters
        ----------
        var_info:
            variable specification
        set_time_units: str
            units for the time coordinate
        frequency: None or str
            frequency of the generated data

        Returns
        -------
        iris.cube.Cube

        """
        coords = []
        scalar_coords = []
        index = 0
        if not frequency:
            frequency = var_info.frequency
        for dim_spec in var_info.coordinates.values():
            coord = self._create_coord_from_spec(dim_spec, set_time_units,
                                                 frequency)
            if isinstance(coord, iris.coords.DimCoord):
                coords.append((coord, index))
                index += 1
            elif isinstance(coord, iris.coords.AuxCoord):
                scalar_coords.append(coord)

        if var_info.valid_min:
            valid_min = float(var_info.valid_min)
        else:
            valid_min = 0

        if var_info.valid_max:
            valid_max = float(var_info.valid_max)
        else:
            valid_max = valid_min + 100

        var_data = (numpy.ones(len(coords) * [20], 'f') *
                    (valid_min + (valid_max - valid_min) / 2))
        cube = iris.cube.Cube(
            var_data,
            standard_name=var_info.standard_name,
            long_name=var_info.long_name,
            var_name=var_info.short_name,
            units=var_info.units,
            attributes=None,
        )
        if var_info.positive:
            cube.attributes['positive'] = var_info.positive

        for coord, i in coords:
            cube.add_dim_coord(coord, i)

        for coord in scalar_coords:
            cube.add_aux_coord(coord)

        return cube

    @staticmethod
    def _construct_scalar_coord(coord_spec):
        return iris.coords.AuxCoord(
            coord_spec.value,
            standard_name=coord_spec.standard_name,
            long_name=coord_spec.long_name,
            var_name=coord_spec.out_name,
            units=coord_spec.units,
            attributes=None)

    def _create_coord_from_spec(self, coord_spec, set_time_units, frequency):
        if coord_spec.units.startswith("days since "):
            coord_spec.units = set_time_units
        coord_spec.frequency = frequency

        if coord_spec.value:
            return self._construct_scalar_coord(coord_spec)

        return self._construct_array_coord(coord_spec)

    def _construct_array_coord(self, dim_spec):
        if dim_spec.units.startswith("days since "):
            values = self._get_time_values(dim_spec)
            unit = Unit(dim_spec.units, calendar='360_day')
        else:
            values = self._get_values(dim_spec)
            unit = Unit(dim_spec.units)
        # Set up attributes dictionary
        coord_atts = {'stored_direction': dim_spec.stored_direction}
        coord = iris.coords.DimCoord(
            values,
            standard_name=dim_spec.standard_name,
            long_name=dim_spec.long_name,
            var_name=dim_spec.out_name,
            attributes=coord_atts,
            units=unit,
        )
        return coord

    @staticmethod
    def _get_values(dim_spec):
        valid_min = dim_spec.valid_min
        if valid_min:
            valid_min = float(valid_min)
        else:
            valid_min = 0.0
        valid_max = dim_spec.valid_max
        if valid_max:
            valid_max = float(valid_max)
        else:
            valid_max = 100.0
        decreasing = dim_spec.stored_direction == 'decreasing'
        endpoint = not dim_spec.standard_name == 'longitude'
        if decreasing:
            values = numpy.linspace(
                valid_max, valid_min, 20, endpoint=endpoint)
        else:
            values = numpy.linspace(
                valid_min, valid_max, 20, endpoint=endpoint)
        values = numpy.array(values)
        if dim_spec.requested:
            requested = [float(val) for val in dim_spec.requested]
            requested.sort(reverse=decreasing)
            for j, request in enumerate(requested):
                values[j] = request
            if decreasing:
                extra_values = numpy.linspace(
                    len(requested), valid_min, 20 - len(requested))
            else:
                extra_values = numpy.linspace(
                    len(requested), valid_max, 20 - len(requested))

            for j in range(len(requested), 20):
                values[j] = extra_values[j - len(requested)]

        return values

    @staticmethod
    def _get_time_values(dim_spec):
        frequency = dim_spec.frequency
        if frequency == 'mon':
            delta = 30
        elif frequency == 'day':
            delta = 1
        elif frequency.ends_with('hr'):
            delta = float(frequency[:-2]) / 24
        else:
            raise Exception('Frequency {} not supported'.format(frequency))
        start = 0
        end = start + delta * 20
        return numpy.arange(start, end, step=delta)


if __name__ == "__main__":
    unittest.main()
