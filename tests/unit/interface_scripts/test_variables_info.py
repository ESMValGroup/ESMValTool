"""

Unit tests for the variable_info module.

"""

import unittest

from esmvaltool.interface_scripts.variable_info import (CoordinateInfo,
                                                        VariableInfo)


class TestVariableInfo(unittest.TestCase):
    def setUp(self):
        self.value = 'value'

    def test_constructor(self):
        info = VariableInfo('var')
        self.assertEqual('var', info.short_name)

    def test_read_empty_dictionary(self):
        info = VariableInfo('var')
        info.read_json({})
        self.assertEqual('', info.standard_name)

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
        self.assertEqual('', info.standard_name)

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
