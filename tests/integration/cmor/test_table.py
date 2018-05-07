"""Integration tests for the variable_info module"""

import os
import unittest

from esmvaltool.cmor.table import CMIP5Info, CMIP6Info


class TestCMIP6Info(unittest.TestCase):
    """Test for the CMIP6 info class"""

    @classmethod
    def setUpClass(cls):
        """
        Set up tests

        We read CMIP6Info once to keep tests times manageable
        """
        cls.variables_info = CMIP6Info()

    def test_custom_tables_location(self):
        """Test constructor with custom tables location"""
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, '..', '..', '..', 'esmvaltool',
                                        'cmor', 'tables', 'cmip6')
        cmor_tables_path = os.path.abspath(cmor_tables_path)
        CMIP6Info(cmor_tables_path)

    def test_get_variable_tas(self):
        """Get tas variable"""
        var = self.variables_info.get_variable('Amon', 'tas')
        self.assertEqual(var.short_name, 'tas')

    def test_get_variable_from_alias(self):
        """Get a variable from a known alias"""
        var = self.variables_info.get_variable('SImon', 'sic')
        self.assertEqual(var.short_name, 'siconc')

    def test_get_bad_variable(self):
        """Get none if a variable is not in the given table"""
        self.assertIsNone(self.variables_info.get_variable('Omon', 'tas'))


class TestCMIP5Info(unittest.TestCase):
    """Test for the CMIP5 info class"""

    @classmethod
    def setUpClass(cls):
        """
        Set up tests

        We read CMIP5Info once to keep testing times manageable
        """
        cls.variables_info = CMIP5Info()

    def test_custom_tables_location(self):
        """Test constructor with custom tables location"""
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, '..', '..', '..', 'esmvaltool',
                                        'cmor', 'tables', 'cmip5')
        cmor_tables_path = os.path.abspath(cmor_tables_path)
        CMIP5Info(cmor_tables_path)

    def test_get_variable_tas(self):
        """Get tas variable"""
        var = self.variables_info.get_variable('Amon', 'tas')
        self.assertEqual(var.short_name, 'tas')

    def test_get_bad_variable(self):
        """Get none if a variable is not in the given table"""
        self.assertIsNone(self.variables_info.get_variable('Omon', 'tas'))
