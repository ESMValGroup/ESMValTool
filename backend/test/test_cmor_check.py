"""
test_cmor_check.py
==================

Unit tests for the CMORCheck class.

"""

# Standard library imports
import os
import unittest
import json
import sys

# Third-party imports
import numpy
import iris
import iris.coords
import iris.coord_categorisation
from cf_units import Unit
from backend.variable_definition import VariablesInfo

# Local imports
from backend.cmor_check import CMORCheck, CMORCheckError
MIP_TABLE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             '..', 'cmip6-cmor-tables', 'Tables')


class TestCMORCheckErrorReporting(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestCMORCheckErrorReporting, cls).setUpClass()
        cls.variables_info = VariablesInfo()

    def setUp(self):
        self.varid = "tas"
        self.table = 'Amon'
        self.cube = CubeCreator().get_cube(self.table, self.varid)

    def test_report_error(self):
        checker = CMORCheck(self.cube, self.table, self.variables_info,)
        self.assertFalse(checker.has_errors())
        checker.report_error('New error: {}', 'something failed')
        self.assertTrue(checker.has_errors())

    def test_fail_on_error(self):
        checker = CMORCheck(self.cube, self.table, self.variables_info, fail_on_error=True)
        with self.assertRaises(CMORCheckError):
            checker.report_error('New error: {}', 'something failed')

    def test_report_warning(self):
        checker = CMORCheck(self.cube, self.table, self.variables_info)
        self.assertFalse(checker.has_errors())
        checker.report_warning('New error: {}', 'something failed')
        self.assertTrue(checker.has_warnings())

    def test_report_warning_with_fail_error(self):
        checker = CMORCheck(self.cube, self.table, self.variables_info, fail_on_error=True)
        sys.stdout.truncate(0)
        checker.report_warning('New error: {}', 'something failed')
        output = sys.stdout.getvalue().strip()  # because stdout is an StringIO instance
        self.assertEquals(output, 'WARNING: New error: something failed')


class TestCMORCheckGoodCube(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestCMORCheckGoodCube, cls).setUpClass()
        cls.variables_info = VariablesInfo()

    def setUp(self):
        self.varid = "tas"
        self.table = "Amon"
        self.cube_creator = CubeCreator()

    def test_check(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_month_number(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        iris.coord_categorisation.add_month_number(cube, 'time')
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_day_of_month(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        iris.coord_categorisation.add_day_of_month(cube, 'time')
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_day_of_year(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        iris.coord_categorisation.add_day_of_year(cube, 'time')
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_year(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        iris.coord_categorisation.add_year(cube, 'time')
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_unit_conversion(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.units = 'degC'
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_ocean_levels(self):
        cube = self.cube_creator.get_cube('Omon', 'uo')
        checker = CMORCheck(cube, 'Omon', self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_band(self):
        cube = self.cube_creator.get_cube('Amon', 'cl')
        checker = CMORCheck(cube, 'Amon', self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_postive_attributte(self):
        cube = self.cube_creator.get_cube(self.table, 'tauu')
        checker = CMORCheck(cube, self.table, self.variables_info)
        checker.check_metadata()
        checker.check_data()

    def test_check_with_var_and_table_correction(self):
        cube = self.cube_creator.get_cube('AERmon', 'o3')
        cube.var_name = 'tro3'
        checker = CMORCheck(cube, 'AERmon', self.variables_info)
        checker.check_metadata()
        checker.check_data()


class TestCMORCheckBadCube(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestCMORCheckBadCube, cls).setUpClass()
        cls.variables_info = VariablesInfo()

    def setUp(self):
        self.varid = "ta"
        self.table = "Amon"
        self.cube_creator = CubeCreator()

    def test_invalid_rank(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info)
        cube.remove_coord('latitude')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_non_increasing(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[-1], coord.points[0],
                                len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_non_decreasing(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info)
        coord = cube.coord('air_pressure')
        values = numpy.linspace(coord.points[-1], coord.points[0],
                                len(coord.points))
        self._update_coordinate_values(cube, coord, values, 2)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_not_correct_lons(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube = cube.intersection(longitude=(-180., 180.))
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_not_correct_lons_automatic_fix(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube = cube.intersection(longitude=(-180., 180.))
        checker = CMORCheck(cube, self.table, self.variables_info, automatic_fixes=True)
        checker.check_metadata()

    def test_not_valid_min(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0]-1, coord.points[-1],
                                len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_not_valid_max(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0], coord.points[-1]+1,
                                len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def _update_coordinate_values(self, cube, coord, values, position):
        cube.remove_coord(coord)
        new_coord = iris.coords.DimCoord(values,
                                         standard_name=coord.standard_name,
                                         long_name=coord.long_name,
                                         var_name=coord.var_name,
                                         units=coord.units)
        cube.add_dim_coord(new_coord, position)

    def test_bad_units(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.coord('latitude').units = 'degrees_n'
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_time(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info)
        cube.coord('time').units = 'days'
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_time_automatic_fix(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, automatic_fixes=True)
        cube.coord('time').units = 'days since 1950-1-1 00:00:00'
        checker.check_metadata()

    def test_bad_time_automatic_fix_failed(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, automatic_fixes=True)
        cube.coord('time').units = 'K'
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_standard_name(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.coord('time').standard_name = 'region'
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_data_units(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.units = 'hPa'
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_data_standard_name(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.standard_name = 'wind_speed'
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_data_var_name(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        cube.var_name = 'wind_speed'
        with self.assertRaises(CMORCheckError):
            checker = CMORCheck(cube, self.table, self.variables_info)
            checker.check_metadata()

    def test_bad_postive_attributte(self):
        cube = self.cube_creator.get_cube(self.table, 'tauu')
        cube.attributes['positive'] = 'up'
        checker = CMORCheck(cube, self.table, self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_standard_name_generic_level_attributte(self):
        cube = self.cube_creator.get_cube('Omon', 'uo')
        cube.coord('depth').standard_name = None
        checker = CMORCheck(cube, 'Omon', self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_axis_generic_level_attribute(self):
        cube = self.cube_creator.get_cube('Omon', 'uo')
        cube.coord('depth').attributes['positive'] = ''
        checker = CMORCheck(cube, 'Omon', self.variables_info)
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_day(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='day')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_subhr(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='subhr')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_dec(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='dec')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_yr(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='yr')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_mon(self):
        cube = self.cube_creator.get_cube('day', 'clt')
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='mon')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_hourly(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='3hr')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    def test_bad_frequency_not_supported(self):
        cube = self.cube_creator.get_cube(self.table, self.varid)
        checker = CMORCheck(cube, self.table, self.variables_info, frequency='wrong_freq')
        with self.assertRaises(CMORCheckError):
            checker.check_metadata()

    # For the moment, we don't have a variable definition with these values
    # to test

    # def test_data_not_valid_max(self):
    #     cube = _create_good_cube(self.varid)
    #     checker = CMORCheck(cube, self.table, self.variables_info)
    #     cube.data[0] = 100000000000
    #     checker.check_metadata()
    #     with self.assertRaises(CMORCheckError):
    #         checker.check_data()
    #
    # def test_data_not_valid_min(self):
    #     cube = _create_good_cube(self.varid)
    #     checker = CMORCheck(cube, self.table, self.variables_info)
    #     cube.data[0] = -100000000000
    #     checker.check_metadata()
    #     with self.assertRaises(CMORCheckError):
    #         checker.check_data()


class CubeCreator(object):
    def __init__(self):
        self.coord_specs = self._parse_mip_table("coordinate")
        self.vars_spec = {}

    def get_cube(self, table, var_name,
                 set_time_units="days since 1850-1-1 00:00:00",
                 frequency=None):
        """
        Creates a cube based on a specification

        :param table: variable's table
        :param var_name: variable's name
        :param set_time_units: time units to use
        :return:
        """
        # Get a specification and build a cube from it.
        if var_name not in self.vars_spec:
            self.vars_spec[var_name] = self._parse_mip_table(table,
                                                             get_var=var_name)

        spec = self.vars_spec[var_name]
        if frequency is None:
            frequency = spec['frequency']

        coords = []
        scalar_coords = []
        index = 0
        for i, dim in enumerate(spec["dimensions"].split()):
            coord = self._create_coord_from_spec(table, dim, set_time_units,
                                                 frequency)
            if isinstance(coord, iris.coords.DimCoord):
                coords.append((coord, index))
                index += 1
            elif isinstance(coord, iris.coords.AuxCoord):
                scalar_coords.append(coord)

        if spec['valid_min']:
            valid_min = float(spec['valid_min'])
        else:
            valid_min = 0

        if spec['valid_max']:
            valid_max = float(spec['valid_max'])
        else:
            valid_max = valid_min + 100

        var_data = (numpy.ones(len(coords) * [20], 'f') *
                    (valid_min + (valid_max - valid_min) / 2))
        cube = iris.cube.Cube(var_data,
                              standard_name=spec["standard_name"],
                              long_name=spec["long_name"],
                              var_name=spec["out_name"],
                              units=spec["units"],
                              attributes=None,
                              cell_methods=spec["cell_methods"],
                              )
        if spec['positive']:
            cube.attributes['positive'] = spec['positive']

        for coord, i in coords:
            cube.add_dim_coord(coord, i)

        for coord in scalar_coords:
            cube.add_aux_coord(coord)

        return cube

    def _construct_scalar_coord(self, coord_spec):
        return iris.coords.AuxCoord(float(coord_spec['value']),
                                    standard_name=coord_spec["standard_name"],
                                    long_name=coord_spec["long_name"],
                                    var_name=coord_spec["out_name"],
                                    units=coord_spec["units"],
                                    attributes=None)
#                                   {'axis': coord_spec["axis"]})

    def _get_coord_spec(self, table, dim, set_time_units):
        try:
            dim_spec = self.coord_specs["axis_entry"][dim]
        except KeyError:
            table_spec = self._parse_mip_table(table)
            if dim in table_spec['Header']['generic_levels']:
                if dim == 'olevel':
                    dim_spec = {'standard_name': 'depth', 'units': '',
                                'value': '', 'valid_min': '', 'valid_max': '',
                                'stored_direction': '', 'requested': '',
                                'long_name': 'ocean_depth_coordinate',
                                'out_name': 'lev', 'axis': 'Z',
                                'positive': 'down'}
                elif dim == 'alevel':
                    dim_spec = {'standard_name': 'height', 'units': '',
                                'value': '', 'valid_min': '', 'valid_max': '',
                                'stored_direction': '', 'requested': '',
                                'long_name': 'ocean_depth_coordinate',
                                'out_name': 'lev', 'axis': 'Z',
                                'positive': 'down'}
                else:
                    emsg = ('Generic levels dimension {} not supported by '
                            'CubeCreator')
                    raise Exception(emsg.format(dim))
            else:
                raise

        if dim_spec["units"].startswith("days since "):
            dim_spec["units"] = set_time_units
        return dim_spec

    def _create_coord_from_spec(self, table, dim, set_time_units, frequency):
        dim_spec = self._get_coord_spec(table, dim, set_time_units)
        dim_spec['frequency'] = frequency

        if dim_spec["value"]:
            return self._construct_scalar_coord(dim_spec)

        return self._construct_array_coord(dim_spec)

    def _construct_array_coord(self, dim_spec):
        if dim_spec["units"].startswith("days since "):
            values = self._get_time_values(dim_spec)
            unit = Unit(dim_spec["units"], calendar='360_day')
        else:
            values = self._get_values(dim_spec)
            unit = Unit(dim_spec["units"])
        # Set up attributes dictionary
        coord_atts = {'stored_direction': dim_spec['stored_direction'],
                      'positive': dim_spec['positive']}
#                     'axis': dim_spec['axis'],}
        coord = iris.coords.DimCoord(values,
                                     standard_name=dim_spec["standard_name"],
                                     long_name=dim_spec["long_name"],
                                     var_name=dim_spec["out_name"],
                                     attributes=coord_atts,
                                     units=unit,
                                     )
        return coord

    def _get_values(self, dim_spec):
        valid_min = dim_spec['valid_min']
        if valid_min:
            valid_min = float(valid_min)
        else:
            valid_min = 0.0
        valid_max = dim_spec['valid_max']
        if valid_max:
            valid_max = float(valid_max)
        else:
            valid_max = 100.0
        decreasing = dim_spec['stored_direction'] == 'decreasing'
        endpoint = not dim_spec['standard_name'] == 'longitude'
        if decreasing:
            values = numpy.linspace(valid_max, valid_min, 20,
                                    endpoint=endpoint)
        else:
            values = numpy.linspace(valid_min, valid_max, 20,
                                    endpoint=endpoint)
        values = numpy.array(values)
        if dim_spec['requested']:
            requested = [float(val) for val in dim_spec['requested']]
            requested.sort(reverse=decreasing)
            for j in range(len(requested)):
                values[j] = requested[j]
            if decreasing:
                extra_values = numpy.linspace(len(requested), valid_min,
                                              20 - len(requested))
            else:
                extra_values = numpy.linspace(len(requested), valid_max,
                                              20 - len(requested))

            for j in range(len(requested), 20):
                values[j] = extra_values[j - len(requested)]

        return values

    def _get_time_values(self, dim_spec):
        frequency = dim_spec['frequency']
        if frequency == 'mon':
            delta = 30
        elif frequency == 'day':
            delta = 1
        elif frequency.ends_with('hr'):
            delta = float(frequency[:-2])/24
        else:
            raise Exception('Frequency {} not supported'.format(frequency))
        start = 0
        end = start + delta * 20
        return numpy.arange(start, end, step=delta)

    def _safely_join_dicts(self, d1, d2):
        """
        Returns a join of two dictionaries but raises exception
        if repeated keys found."
        :param d1: first dictionary
        :param d2: second dictionary
        :return: joint dictionary
        """
        d = d1

        for key, value in d2.items():
            if key in d:
                emsg = "Repeated key in dictionaries when merging: {}"
                raise Exception(emsg.format(key))
            d[key] = value

        return d

    def _parse_mip_table(self, table_name, get_var=None):
        "Reads MIP table and returns a dictionary."
        table_path = os.path.join(MIP_TABLE_DIR,
                                  'CMIP6_{}.json'.format(table_name))
        d = json.load(open(table_path))

        if get_var:
            header_dict = d["Header"]
            var_dict = d["variable_entry"][get_var]
            d = self._safely_join_dicts(header_dict, var_dict)

        return d

if __name__ == "__main__":
    unittest.main()
