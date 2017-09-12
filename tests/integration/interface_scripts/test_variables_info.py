"""

Unit tests for the variable_info module.

"""

import os
import unittest

from esmvaltool.interface_scripts.variable_info import CMIP5Info, CMIP6Info


class TestCMIP6Info(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.variables_info = CMIP6Info()

    def test_constructor_optional_parameter(self):
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, '..', '..', '..',
                                        'esmvaltool', 'interface_scripts', 'cmip6-cmor-tables')
        cmor_tables_path = os.path.abspath(cmor_tables_path)
        CMIP6Info(cmor_tables_path)

    def test_get_variable_tas(self):
        var = self.variables_info.get_variable('Amon', 'tas')
        self.assertEqual(var.short_name, 'tas')

    def test_get_variable_with_changed_name(self):
        var = self.variables_info.get_variable('SImon', 'sic')
        self.assertEqual(var.short_name, 'siconc')

    def test_get_bad_variable(self):
        self.assertIsNone(self.variables_info.get_variable('Omon', 'tas'))


class TestCMIP5Info(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.variables_info = CMIP5Info()

    def test_constructor_optional_parameter(self):
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, '..', '..', '..',
                                        'esmvaltool', 'interface_scripts', 'cmip5-cmor-tables')
        cmor_tables_path = os.path.abspath(cmor_tables_path)
        CMIP5Info(cmor_tables_path)

    def test_get_variable_tas(self):
        var = self.variables_info.get_variable('Amon', 'tas')
        self.assertEqual(var.short_name, 'tas')

    def test_get_bad_variable(self):
        self.assertIsNone(self.variables_info.get_variable('Omon', 'tas'))
