"""Unit tests for the variable_info module."""

import unittest

import mock

from esmvaltool.cmor.fix import Fix, fix_data, fix_file, fix_metadata


class TestFixFile(unittest.TestCase):
    """Fix file tests"""

    def setUp(self):
        """Prepare for testing"""
        self.filename = 'filename'
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_file.return_value = 'new_filename'

    def test_fix(self):
        """Check that the returned fix is applied"""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                return_value=[self.mock_fix]):
            file_returned = fix_file('filename', 'short_name', 'project',
                                     'model', 'output_dir')
            self.assertNotEqual(file_returned, self.filename)
            self.assertEqual(file_returned, 'new_filename')

    def test_nofix(self):
        """Check that the same file is returned if no fix is available"""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes', return_value=[]):
            file_returned = fix_file('filename', 'short_name', 'project',
                                     'model', 'output_dir')
            self.assertEqual(file_returned, self.filename)


class TestGetCube(unittest.TestCase):
    """Test get cube by var_name method"""

    def setUp(self):
        """Prepare for testing"""
        self.cube_1 = mock.Mock()
        self.cube_1.var_name = 'cube1'
        self.cube_2 = mock.Mock()
        self.cube_2.var_name = 'cube2'
        self.cubes = [self.cube_1, self.cube_2]
        self.fix = Fix()

    def test_get_first_cube(self):
        """Test selecting first cube"""
        self.assertIs(self.cube_1,
                      self.fix.get_cube_from_list(self.cubes, "cube1"))

    def test_get_second_cube(self):
        """Test selecting second cube."""
        self.assertIs(self.cube_2,
                      self.fix.get_cube_from_list(self.cubes, "cube2"))

    def test_get_default_raises(self):
        """Check that the default raises (Fix is not a cube)."""
        with self.assertRaises(Exception):
            self.fix.get_cube_from_list(self.cubes)

    def test_get_default(self):
        """Check that the default raises (Fix is a cube)."""
        self.cube_1.var_name = 'Fix'
        self.assertIs(self.cube_1, self.fix.get_cube_from_list(self.cubes))


class TestFixMetadata(unittest.TestCase):
    """Fix metadata tests."""

    def setUp(self):
        """Prepare for testing."""
        self.cube = mock.Mock()
        self.cube.attributes = {'source_file': 'source_file'}
        self.fixed_cube = mock.Mock()
        self.fixed_cube.attributes = {'source_file': 'source_file'}
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_metadata.return_value = [self.fixed_cube]

    def test_fix(self):
        """Check that the returned fix is applied."""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                return_value=[self.mock_fix]):
            cube_returned = fix_metadata([self.cube], 'short_name', 'project',
                                         'model')[0]
            self.assertTrue(cube_returned is not self.cube)
            self.assertTrue(cube_returned is self.fixed_cube)

    def test_nofix(self):
        """Check that the same cube is returned if no fix is available."""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes', return_value=[]):
            cube_returned = fix_metadata([self.cube], 'short_name', 'project',
                                         'model')[0]
            self.assertTrue(cube_returned is self.cube)
            self.assertTrue(cube_returned is not self.fixed_cube)

    def test_cmor_checker_called(self):
        """Check that the cmor check is done."""
        checker = mock.Mock()
        checker.return_value = mock.Mock()
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes', return_value=[]):
            with mock.patch(
                    'esmvaltool.cmor.fix._get_cmor_checker',
                    return_value=checker) as get_mock:
                fix_metadata([self.cube], 'short_name', 'project', 'model',
                             'cmor_table', 'mip', 'frequency')
                get_mock.assert_called_once_with(
                    automatic_fixes=True,
                    fail_on_error=False,
                    frequency='frequency',
                    mip='mip',
                    short_name='short_name',
                    table='cmor_table')
                checker.assert_called_once_with(self.cube)
                checker.return_value.check_metadata.assert_called_once_with()


class TestFixData(unittest.TestCase):
    """Fix data tests."""

    def setUp(self):
        """Prepare for testing."""
        self.cube = mock.Mock()
        self.fixed_cube = mock.Mock()
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_data.return_value = self.fixed_cube

    def test_fix(self):
        """Check that the returned fix is applied."""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                return_value=[self.mock_fix]):
            cube_returned = fix_data(self.cube, 'short_name', 'project',
                                     'model')
            self.assertTrue(cube_returned is not self.cube)
            self.assertTrue(cube_returned is self.fixed_cube)

    def test_nofix(self):
        """Check that the same cube is returned if no fix is available."""
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes', return_value=[]):
            cube_returned = fix_data(self.cube, 'short_name', 'project',
                                     'model')
            self.assertTrue(cube_returned is self.cube)
            self.assertTrue(cube_returned is not self.fixed_cube)

    def test_cmor_checker_called(self):
        """Check that the cmor check is done"""
        checker = mock.Mock()
        checker.return_value = mock.Mock()
        with mock.patch(
                'esmvaltool.cmor._fixes.fix.Fix.get_fixes', return_value=[]):
            with mock.patch(
                    'esmvaltool.cmor.fix._get_cmor_checker',
                    return_value=checker) as get_mock:
                fix_data(self.cube, 'short_name', 'project', 'model',
                         'cmor_table', 'mip', 'frequency')
                get_mock.assert_called_once_with(
                    automatic_fixes=True,
                    fail_on_error=False,
                    frequency='frequency',
                    mip='mip',
                    short_name='short_name',
                    table='cmor_table')
                checker.assert_called_once_with(self.cube)
                checker.return_value.check_data.assert_called_once_with()
