"""
test_cmor_check.py
==================

Unit tests for the CMORCheck class.

"""

# Standard library imports
import os, re, sys
import unittest
import json


# Temporary path fix
sys.path.extend([".", "..", "backend"])
MIP_TABLE_DIR = os.environ.get("CMIP6_CMOR_TABLES", None)

if not MIP_TABLE_DIR:
   raise Exception("Cannot find location of MIP TABLES - please set CMIP6_CMOR_TABLES environment variable.")

# Local imports
import cmor_check


# Functions to create example cubes to test
def _create_good_cube(): pass
def _create_bad_cube(): pass

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
    d = json.load(open(table_path).read())

    if get_var:
        header_dict = d["Header"]
        var_dict = d["variable_entry"][get_var] 
        d = _safely_join_dicts(header_dict, var_dict)

    return d    


class TestCMORCheck(unittest.TestCase):

    def setUp(self):
        self.tas_spec = _parse_mip_table("CMIP6_Amon.json", get_var="tas") 
        self.coords_dict = _parse_mip_table("CMIP6_coordinate.json")["axis_entry"]
        self.cube_good_1 = _create_good_cube()
        self.cube_bad_1 = _create_bad_cube()


    def test_read_tas_spec(self):
        assert tas_spec["units"] == "K"

    def test_check_rank(self):
        assert 

"""
_check_rank() - checks rank

_check_dim_names() - ...

_check_coords() - from JSON tables.

_check_time_coord() - from JSON tables.

_check_var_metadata() - from JSON tables.

_check_fill_value() - ...

_check_data_range() - using valid_min and valid_max.
"""

if __name__ == "__main__":

    unittest.main()
