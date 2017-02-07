"""
test_cmor_check.py
==================

Unit tests for the CMORCheck class.

"""

# Standard library imports
import os, re, sys
import unittest
import json

# Third-party imports
import numpy
import iris
from iris.cube import Cube


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

    axis_data = numpy.array(range(-50, 51, 5), 'f')

    coords = []
    for i, dim in enumerate(spec["dimensions"].split()):
        dim_spec = coord_spec["axis_entry"][dim] 

        # Check for time axis
        if dim_spec["units"].startswith("days since "):
            dim_spec["units"] = set_time_units
            
        coord = iris.coords.DimCoord(axis_data,
                                     standard_name=dim_spec["standard_name"],
                                     long_name=dim_spec["long_name"],
                                     var_name=dim_spec["out_name"],
                                     units=dim_spec["units"])
        coords.append((coord, i)) 
         
    var_data = numpy.ones(len(coords) * [len(axis_data)], 'f')
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

    axis_data = numpy.array(range(-50, 51, 5), 'f')

    coords = []
    coord_spec = _parse_mip_table("CMIP6_coordinate.json")
    for i, dim in enumerate(['time', 'longitude']):
        dim_spec = coord_spec["axis_entry"][dim]

        # Check for time axis
        if dim_spec["units"].startswith("days since "):
            dim_spec["units"] = set_time_units

        coord = iris.coords.DimCoord(axis_data,
                                     standard_name=dim_spec["standard_name"],
                                     long_name=dim_spec["long_name"],
                                     var_name=dim_spec["out_name"],
                                     units=dim_spec["units"])
        coords.append((coord, i))

    var_data = numpy.ones(len(coords) * [len(axis_data)], 'f')
    cb = iris.cube.Cube(var_data,
                        standard_name="air_temperature",
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
        self.cube = _create_good_cube(self.varid)

    def test_report_error(self):
        checker = cmor_check.CMORCheck(self.cube)
        self.assertFalse(checker.has_errors())
        checker.report_error('New error: {}', 'something failed')
        self.assertTrue(checker.has_errors())

    def test_report_raises(self):
        checker = cmor_check.CMORCheck(self.cube)
        with self.assertRaises(cmor_check.CMORCheckError):
            checker.report_error('New error: {}', 'something failed')


class TestCMORCheckGoodCube(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var=self.varid) 
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_good_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube)

    def test_read_tas_spec(self):
        assert self.tas_spec["units"] == "K"

    def test_check(self):
        # Will raise exception if it fails
        self.checker.check()


class TestCMORCheckBadCube(unittest.TestCase):

    def setUp(self):
        self.varid = "tas"
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube = _create_bad_cube(self.varid)
        self.checker = cmor_check.CMORCheck(self.cube)

    def test_check(self):
        # Will raise exception if it fails
        with self.assertRaises(cmor_check.CMORCheckError):
            self.checker.check()


if __name__ == "__main__":

    unittest.main()
