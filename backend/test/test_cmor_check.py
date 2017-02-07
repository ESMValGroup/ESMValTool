"""
test_cmor_check.py
==================

Unit tests for the CMORCheck class.

"""

# Standard library imports
import os, sys
import unittest
import json

# Third-party imports
import numpy
import iris


# Temporary path fix
sys.path.extend([".", "..", "backend"])
# MIP_TABLE_DIR = os.environ.get("CMIP6_CMOR_TABLES", None)
MIP_TABLE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'cmip6-cmor-tables', 'Tables')

# Local imports
import cmor_check


def _safely_join_dicts(d1, d2):
    "Returns a join of two dictionaries but raises exception if repeated keys found."
    d = d1

    for key, value in d2.items():
        if key in d:
            raise Exception("Repeated key in dictionaries when merging: %s" % key)
        d[key] = value

    return d


def _parse_mip_table(table_name, get_var=None):
    "Reads MIP table and returns a dictionary."
    table_path = os.path.join(MIP_TABLE_DIR, table_name)
    d = json.load(open(table_path))

    if get_var:
        header_dict = d["Header"]
        var_dict = d["variable_entry"][get_var] 
        d = _safely_join_dicts(header_dict, var_dict)

    return d    


# Functions to create iris cubes to test
def _create_good_cube(get_var, set_time_units="days since 1850-1-1 00:00:00"):
    "Creates a cube based on a specification for ``get_var``."
    # Get a specification and build a cube from it. 
    spec = _parse_mip_table("CMIP6_Amon.json", get_var=get_var)
    coord_spec = _parse_mip_table("CMIP6_coordinate.json")

    coords = []
    scalar_coords = []
    index = 0
    for i, dim in enumerate(spec["dimensions"].split()):
        coord = create_coord_from_spec(coord_spec, dim, set_time_units)
        if isinstance(coord, iris.coords.DimCoord):
            coords.append((coord, index))
            index += 1
        else:
            scalar_coords.append(coord)

    if spec['valid_min']:
        valid_min = float(spec['valid_min'])
    else:
        valid_min = 0

    if spec['valid_max']:
        valid_max = float(spec['valid_max'])
    else:
        valid_max = valid_min + 100
         
    var_data = numpy.ones(len(coords) * [20], 'f') * (valid_min + (valid_max - valid_min)/2)
    cb = iris.cube.Cube(var_data,
                        standard_name=spec["standard_name"],
                        long_name=spec["long_name"],
                        var_name=spec["out_name"],
                        units=spec["units"],
                        attributes=None,
                        cell_methods=spec["cell_methods"])
    if spec['positive']:
        cb.attributes['positive'] = spec['positive']

    for coord, i in coords:
        cb.add_dim_coord(coord, i)

    for coord in scalar_coords:
        cb.add_aux_coord(coord)

    return cb


def _create_bad_cube(get_var, set_time_units="days since 1950-01-01 00:00:00"):

    coords = []
    scalar_coords = []
    coord_spec = _parse_mip_table("CMIP6_coordinate.json")
    spec = _parse_mip_table("CMIP6_Amon.json", get_var=get_var)

    index = 0
    for i, dim in enumerate(['time', 'height2m']):
        coord = create_coord_from_spec(coord_spec, dim, set_time_units)
        if isinstance(coord, iris.coords.DimCoord):
            coords.append((coord, index))
            index += 1
        else:
            scalar_coords.append(coord)

    coords.reverse()
    var_data = numpy.ones(len(coords) * [20], 'f')
    cb = iris.cube.Cube(var_data,
                        standard_name=spec["standard_name"],
                        long_name=spec["long_name"],
                        var_name=spec["out_name"],
                        units=spec["units"],
                        attributes=None,
                        cell_methods=spec["cell_methods"])

    for coord, i in coords:
        cb.add_dim_coord(coord, i)
    for coord in scalar_coords:
        cb.add_aux_coord(coord)

    return cb


def construct_scalar_variable(coord_spec):
    return iris.coords.AuxCoord(float(coord_spec['value']),
                                standard_name=coord_spec["standard_name"],
                                long_name=coord_spec["long_name"],
                                var_name=coord_spec["out_name"],
                                units=coord_spec["units"],
                                attributes=None)


def create_coord_from_spec(coord_spec, dim, set_time_units):
    dim_spec = coord_spec["axis_entry"][dim]
    # Check for time axis
    if dim_spec["units"].startswith("days since "):
        dim_spec["units"] = set_time_units

    if dim_spec["value"]:
        return construct_scalar_variable(dim_spec)

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

    if decreasing:
        values = numpy.linspace(valid_max, valid_min, 20)
    else:
        values = numpy.linspace(valid_min, valid_max, 20)
    values = numpy.array(values)
    if dim_spec['requested']:
        requested = [float(val) for val in dim_spec['requested']]
        requested.sort(reverse=decreasing)
        for j in range(len(requested)):
            values[j] = requested[j]
        if decreasing:
            extra_values = numpy.linspace(len(requested), valid_min, 20 - len(requested))
        else:
            extra_values = numpy.linspace(len(requested), valid_max, 20 - len(requested))

        for j in range(len(requested), 20):
            values[j] = extra_values[j - len(requested)]


    # Set up attributes dictionary
    coord_atts = {'stored_direction': dim_spec['stored_direction']}
    coord = iris.coords.DimCoord(values,
                                 standard_name=dim_spec["standard_name"],
                                 long_name=dim_spec["long_name"],
                                 var_name=dim_spec["out_name"],
                                 attributes=coord_atts,
                                 units=dim_spec["units"])
    return coord


class TestCMORCheckErrorReporting(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.table = 'Amon'
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var=self.varid)
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_bad_cube(self.varid)

    def test_failed_on_error_true(self):
        self.checker = cmor_check.CMORCheck(self.cube, self.table, fail_on_error=True)
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker._check_rank()
       
    def test_failed_on_error_false(self):
        self.checker = cmor_check.CMORCheck(self.cube, self.table, fail_on_error=False)
        self.checker._check_rank()

        self.checker._check_dim_names()
        assert len(self.checker._errors) > 1

    def test_report_error(self):
        checker = cmor_check.CMORCheck(self.cube, self.table,)
        self.assertFalse(checker.has_errors())
        checker.report_error('New error: {}', 'something failed')
        self.assertTrue(checker.has_errors())

    def test_report_raises(self):
        checker = cmor_check.CMORCheck(self.cube, self.table, fail_on_error=True)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.report_error('New error: {}', 'something failed')


class TestCMORCheckGoodCube(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.table = "Amon"
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var=self.varid) 
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_good_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube, self.table, fail_on_error=True)

    def test_read_tas_spec(self):
        assert self.tas_spec["units"] == "K"

    def test_check(self):
        self.checker.check()

    def test_postive_attributte(self):
        cube = _create_good_cube('tauu')
        checker = cmor_check.CMORCheck(cube, self.table)
        checker.check()


class TestCMORCheckBadCube(unittest.TestCase):

    def setUp(self):
        self.varid = "ta"
        self.table = "Amon"
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_bad_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube, self.table, fail_on_error=True)

    def test_check_rank_fails(self):
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker._check_rank()

    def test_check_dim_names_fails(self):
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker._check_dim_names()

    def test_check_coords_fails_reversed_direction(self):
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker._check_coords()

    def test_non_increasing(self):
        cube = _create_good_cube(self.varid)
        checker = cmor_check.CMORCheck(cube, self.table)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[-1], coord.points[0], len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_non_decreasing(self):
        cube = _create_good_cube(self.varid)
        checker = cmor_check.CMORCheck(cube, self.table)
        coord = cube.coord('air_pressure')
        values = numpy.linspace(coord.points[-1], coord.points[0], len(coord.points))
        self._update_coordinate_values(cube, coord, values, 2)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def _update_coordinate_values(self, cube, coord, values, position):
        cube.remove_coord(coord)
        new_coord = iris.coords.DimCoord(values,
                                         standard_name=coord.standard_name,
                                         long_name=coord.long_name,
                                         var_name=coord.var_name,
                                         units=coord.units)
        cube.add_dim_coord(new_coord, position)

    def test_not_valid_min(self):
        cube = _create_good_cube(self.varid)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0]-1, coord.points[-1], len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_not_valid_max(self):
        cube = _create_good_cube(self.varid)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0], coord.points[-1]+1, len(coord.points))
        self._update_coordinate_values(cube, coord, values, 1)
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_units(self):
        cube = _create_good_cube(self.varid)
        cube.coord('latitude').units = 'degrees_n'
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_time(self):
        cube = _create_good_cube(self.varid)
        checker = cmor_check.CMORCheck(cube, self.table)
        cube.coord('time').units = 'days'
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_standard_name(self):
        cube = _create_good_cube(self.varid)
        cube.coord('time').standard_name = 'region'
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_data_units(self):
        cube = _create_good_cube(self.varid)
        cube.units = 'hPa'
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_data_standard_name(self):
        cube = _create_good_cube(self.varid)
        cube.standard_name = 'wind_speed'
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_bad_data_var_name(self):
        cube = _create_good_cube(self.varid)
        cube.var_name = 'wind_speed'
        with self.assertRaises(KeyError):
            checker = cmor_check.CMORCheck(cube, self.table)
            checker.check()

    def test_bad_postive_attributte(self):
        cube = _create_good_cube('tauu')
        cube.attributes['positive'] = 'up'
        checker = cmor_check.CMORCheck(cube, self.table)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    # For the moment, we don't have a variable definition with these values to test

    # def test_data_not_valid_max(self):
    #     cube = _create_good_cube(self.varid)
    #     checker = cmor_check.CMORCheck(cube, self.table)
    #     cube.data[0] = 100000000000
    #     with self.assertRaises(cmor_check.CMORCheckError):
    #         checker.check()
    #
    # def test_data_not_valid_min(self):
    #     cube = _create_good_cube(self.varid)
    #     checker = cmor_check.CMORCheck(cube, self.table)
    #     cube.data[0] = -100000000000
    #     with self.assertRaises(cmor_check.CMORCheckError):
    #         checker.check()


if __name__ == "__main__":
    unittest.main()
