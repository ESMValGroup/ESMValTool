"""Unit tests for the variable_info module."""

import unittest

from esmvaltool.cmor.table import CoordinateInfo, VariableInfo


class TestVariableInfo(unittest.TestCase):
    """Variable info tests"""

    def setUp(self):
        """Prepare for testing"""
        self.value = 'value'

    def test_constructor(self):
        """Test basic constructor"""
        info = VariableInfo('var')
        self.assertEqual('var', info.short_name)

    def test_read_empty_dictionary(self):
        """Test read empty dict"""
        info = VariableInfo('var')
        info.read_json({})
        self.assertEqual('', info.standard_name)

    def test_read_standard_name(self):
        """Test standard_name"""
        info = VariableInfo('var')
        info.read_json({'standard_name': self.value})
        self.assertEqual(info.standard_name, self.value)

    def test_read_long_name(self):
        """Test long_name"""
        info = VariableInfo('var')
        info.read_json({'long_name': self.value})
        self.assertEqual(info.long_name, self.value)

    def test_read_units(self):
        """Test units"""
        info = VariableInfo('var')
        info.read_json({'units': self.value})
        self.assertEqual(info.units, self.value)

    def test_read_valid_min(self):
        """Test valid_min"""
        info = VariableInfo('var')
        info.read_json({'valid_min': self.value})
        self.assertEqual(info.valid_min, self.value)

    def test_read_valid_max(self):
        """Test valid_max"""
        info = VariableInfo('var')
        info.read_json({'valid_max': self.value})
        self.assertEqual(info.valid_max, self.value)

    def test_read_positive(self):
        """Test positive"""
        info = VariableInfo('var')
        info.read_json({'positive': self.value})
        self.assertEqual(info.positive, self.value)


class TestCoordinateInfo(unittest.TestCase):
    """Tests for CoordinataInfo"""

    def setUp(self):
        """Prepare for testing"""
        self.value = 'value'

    def test_constructor(self):
        """Test constructor"""
        info = CoordinateInfo('var')
        self.assertEqual('var', info.name)

    def test_read_empty_dictionary(self):
        """Test empty dict"""
        info = CoordinateInfo('var')
        info.read_json({})
        self.assertEqual('', info.standard_name)

    def test_read_standard_name(self):
        """Test standard_name"""
        info = CoordinateInfo('var')
        info.read_json({'standard_name': self.value})
        self.assertEqual(info.standard_name, self.value)

    def test_read_var_name(self):
        """Test var_name"""
        info = CoordinateInfo('var')
        info.read_json({'var_name': self.value})
        self.assertEqual(info.var_name, self.value)

    def test_read_out_name(self):
        """Test out_name"""
        info = CoordinateInfo('var')
        info.read_json({'out_name': self.value})
        self.assertEqual(info.out_name, self.value)

    def test_read_units(self):
        """Test units"""
        info = CoordinateInfo('var')
        info.read_json({'units': self.value})
        self.assertEqual(info.units, self.value)

    def test_read_valid_min(self):
        """Test valid_min"""
        info = CoordinateInfo('var')
        info.read_json({'valid_min': self.value})
        self.assertEqual(info.valid_min, self.value)

    def test_read_valid_max(self):
        """Test valid_max"""
        info = CoordinateInfo('var')
        info.read_json({'valid_max': self.value})
        self.assertEqual(info.valid_max, self.value)

    def test_read_value(self):
        """Test value"""
        info = CoordinateInfo('var')
        info.read_json({'value': self.value})
        self.assertEqual(info.value, self.value)

    def test_read_requested(self):
        """Test requested"""
        value = ['value1', 'value2']
        info = CoordinateInfo('var')
        info.read_json({'requested': value})
        self.assertEqual(info.requested, value)
