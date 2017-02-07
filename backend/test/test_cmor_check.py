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

if not MIP_TABLE_DIR:
   raise Exception("Cannot find location of MIP TABLES - please set CMIP6_CMOR_TABLES environment variable.")

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
def _create_good_cube(get_var, set_time_units="days since 1950-01-01 00:00:00"): 
    "Creates a cube based on a specification for ``get_var``."
    # Get a specification and build a cube from it. 
    spec = _parse_mip_table("CMIP6_Amon.json", get_var="tas")
    coord_spec = _parse_mip_table("CMIP6_coordinate.json")


    coords = []
    for i, dim in enumerate(spec["dimensions"].split()):
        dim_spec = coord_spec["axis_entry"][dim] 

        # Check for time axis
        if dim_spec["units"].startswith("days since "):
            dim_spec["units"] = set_time_units
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

        coord = iris.coords.DimCoord(numpy.linspace(valid_min, valid_max, 20),
                                     standard_name=dim_spec["standard_name"],
                                     long_name=dim_spec["long_name"],
                                     var_name=dim_spec["out_name"],
                                     units=dim_spec["units"])
        coords.append((coord, i)) 
         
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

    return cb


def _create_bad_cube(get_var, set_time_units="days since 1950-01-01 00:00:00"):

    coords = []
    coord_spec = _parse_mip_table("CMIP6_coordinate.json")

    for i, dim in enumerate(['time', 'longitude']):
        dim_spec = coord_spec["axis_entry"][dim]

        # Check for time axis
        if dim_spec["units"].startswith("days since "):
            dim_spec["units"] = set_time_units
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

        # Set up attributes dictionary
        coord_atts = {'stored_direction': dim_spec['stored_direction']}
        values = numpy.linspace(valid_min, valid_max, 20)

        # Reverse values so they should contradict the 'stored_direction' 
        values = values[::-1] 
        coord = iris.coords.DimCoord(values,
                                     standard_name=dim_spec["standard_name"],
                                     long_name=dim_spec["long_name"],
                                     var_name=dim_spec["out_name"],
                                     attributes=coord_atts,
                                     units=dim_spec["units"])
        coords.append((coord, i))

    coords.reverse()
    var_data = numpy.ones(len(coords) * [20], 'f')
    cb = iris.cube.Cube(var_data,
                        standard_name="geopotential_height",
                        long_name="Air Temperature",
                        var_name="tas",
                        units="K",
                        attributes=None,
                        cell_methods="")

    for coord, i in coords:
        cb.add_dim_coord(coord, i)

    return cb


class TestCMORCheckErrorReporting(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var=self.varid)
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_bad_cube(self.varid)

    def test_failed_on_error_true(self):
        self.checker = cmor_check.CMORCheck(self.cube, fail_on_error=True)
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker._check_rank()
       
    def test_failed_on_error_false(self):
        self.checker = cmor_check.CMORCheck(self.cube, fail_on_error=False)
        self.checker._check_rank()

        self.checker._check_dim_names()
        assert len(self.checker._errors) > 1

    def test_report_error(self):
        checker = cmor_check.CMORCheck(self.cube)
        self.assertFalse(checker.has_errors())
        checker.report_error('New error: {}', 'something failed')
        self.assertTrue(checker.has_errors())

    def test_report_raises(self):
        checker = cmor_check.CMORCheck(self.cube, fail_on_error=True)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.report_error('New error: {}', 'something failed')


class TestCMORCheckGoodCube(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var=self.varid) 
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_good_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube, fail_on_error=True)

    def test_read_tas_spec(self):
        assert self.tas_spec["units"] == "K"

    def test_check_rank(self):
        # Will raise exception if it fails
        self.checker._check_rank()

    def test_check_dim_names(self):
        # Will raise exception if it fails
        self.checker._check_dim_names()

    def test_check_coords(self):
        # Will raise exception if it fails
        self.checker._check_coords()


class TestCMORCheckBadCube(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_bad_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube, fail_on_error=True)

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
        checker = cmor_check.CMORCheck(cube)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[-1], coord.points[0], len(coord.points))
        self._update_latitude_values(cube, coord, values)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def _update_latitude_values(self, cube, coord, values):
        cube.remove_coord(coord)
        new_coord = iris.coords.DimCoord(values,
                                         standard_name=coord.standard_name,
                                         long_name=coord.long_name,
                                         var_name=coord.var_name,
                                         units=coord.units)
        cube.add_dim_coord(new_coord, 1)

    def test_not_valid_min(self):
        cube = _create_good_cube(self.varid)
        checker = cmor_check.CMORCheck(cube)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0]-1, coord.points[-1], len(coord.points))
        self._update_latitude_values(cube, coord, values)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()

    def test_not_valid_max(self):
        cube = _create_good_cube(self.varid)
        checker = cmor_check.CMORCheck(cube)
        coord = cube.coord('latitude')
        values = numpy.linspace(coord.points[0], coord.points[-1]+1, len(coord.points))
        self._update_latitude_values(cube, coord, values)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.check()


if __name__ == "__main__":
    unittest.main()
