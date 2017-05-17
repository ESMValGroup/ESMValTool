"""

Unit tests for the CMORCheck class.

"""

import unittest
import os

from backend.variable_info import VariableInfo, CoordinateInfo, VariablesInfo


class TestVariablesInfo(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.variables_info = VariablesInfo()

    def test_constructor_optional_parameter(self):
        cwd = os.path.dirname(os.path.realpath(__file__))
        cmor_tables_path = os.path.join(cwd, '..', '..', 'cmip6-cmor-tables')
        VariablesInfo(cmor_tables_path)

    def test_get_variable_tas(self):
        var = self.variables_info.get_variable('Amon', 'tas')
        self.assertEqual(var.short_name, 'tas')

    def test_get_variable_with_changed_name(self):
        var = self.variables_info.get_variable('SImon', 'sic')
        self.assertEqual(var.short_name, 'siconc')

    def test_get_bad_variable(self):
        self.assertIsNone(self.variables_info.get_variable('Omon', 'tas'))


class TestVariableInfo(unittest.TestCase):

    def setUp(self):
        self.value = 'value'

    def test_constructor(self):
        info = VariableInfo('var')
        self.assertEqual('var', info.short_name)

    def test_read_empty_dictionary(self):
        info = VariableInfo('var')
        info.read_json({})
        self.assertEquals('', info.standard_name)

    def test_read_standard_name(self):
        info = VariableInfo('var')
        info.read_json({'standard_name': self.value})
        self.assertEqual(info.standard_name, self.value)

    def test_read_long_name(self):
        info = VariableInfo('var')
        info.read_json({'long_name': self.value})
        self.assertEqual(info.long_name, self.value)

    def test_read_units(self):
        info = VariableInfo('var')
        info.read_json({'units': self.value})
        self.assertEqual(info.units, self.value)

    def test_read_valid_min(self):
        info = VariableInfo('var')
        info.read_json({'valid_min': self.value})
        self.assertEqual(info.valid_min, self.value)

    def test_read_valid_max(self):
        info = VariableInfo('var')
        info.read_json({'valid_max': self.value})
        self.assertEqual(info.valid_max, self.value)

    def test_read_positive(self):
        info = VariableInfo('var')
        info.read_json({'positive': self.value})
        self.assertEqual(info.positive, self.value)


class TestCoordinateInfo(unittest.TestCase):

    def setUp(self):
        self.value = 'value'

    def test_constructor(self):
        info = CoordinateInfo('var')
        self.assertEqual('var', info.name)

    def test_read_empty_dictionary(self):
        info = CoordinateInfo('var')
        info.read_json({})
        self.assertEquals('', info.standard_name)

    def test_read_standard_name(self):
        info = CoordinateInfo('var')
        info.read_json({'standard_name': self.value})
        self.assertEqual(info.standard_name, self.value)

    def test_read_var_name(self):
        info = CoordinateInfo('var')
        info.read_json({'var_name': self.value})
        self.assertEqual(info.var_name, self.value)

    def test_read_out_name(self):
        info = CoordinateInfo('var')
        info.read_json({'out_name': self.value})
        self.assertEqual(info.out_name, self.value)

    def test_read_units(self):
        info = CoordinateInfo('var')
        info.read_json({'units': self.value})
        self.assertEqual(info.units, self.value)

    def test_read_valid_min(self):
        info = CoordinateInfo('var')
        info.read_json({'valid_min': self.value})
        self.assertEqual(info.valid_min, self.value)

    def test_read_valid_max(self):
        info = CoordinateInfo('var')
        info.read_json({'valid_max': self.value})
        self.assertEqual(info.valid_max, self.value)

    def test_read_value(self):
        info = CoordinateInfo('var')
        info.read_json({'value': self.value})
        self.assertEqual(info.value, self.value)

    def test_read_requested(self):
        value = ['value1', 'value2']
        info = CoordinateInfo('var')
        info.read_json({'requested': value})
        self.assertEqual(info.requested, value)

