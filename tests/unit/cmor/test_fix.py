"""Unit tests for the variable_info module."""

import unittest

from esmvaltool.cmor.fix import fix_file, fix_data, fix_metadata
import mock


class TestFixFile(unittest.TestCase):
    """Fix file tests"""

    def setUp(self):
        """Prepare for testing"""
        self.filename = 'filename'
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_file.return_value = 'new_filename'

    def test_fix(self):
        """Check that the returned fix is applied"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[self.mock_fix]):
            file_returned = fix_file('filename', 'short_name', 'project',
                                     'model', 'output_dir')
            self.assertNotEqual(file_returned, self.filename)
            self.assertEqual(file_returned, 'new_filename')

    def test_nofix(self):
        """Check that the same file is returned if no fix is available"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[]):
            file_returned = fix_file('filename', 'short_name', 'project',
                                     'model', 'output_dir')
            self.assertEqual(file_returned, self.filename)


class TestFixMetadata(unittest.TestCase):
    """Fix metadata tests"""

    def setUp(self):
        """Prepare for testing"""
        self.cube = mock.Mock()
        self.fixed_cube = mock.Mock()
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_metadata.return_value = self.fixed_cube

    def test_fix(self):
        """Check that the returned fix is applied"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[self.mock_fix]):
            cube_returned = fix_metadata(self.cube, 'short_name', 'project',
                                         'model')
            self.assertTrue(cube_returned is not self.cube)
            self.assertTrue(cube_returned is self.fixed_cube)

    def test_nofix(self):
        """Check that the same cube is returned if no fix is available"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[]):
            cube_returned = fix_metadata(self.cube, 'short_name', 'project',
                                         'model')
            self.assertTrue(cube_returned is self.cube)
            self.assertTrue(cube_returned is not self.fixed_cube)

    def test_cmor_checker_called(self):
        """Check that the cmor check is done"""
        checker = mock.Mock()
        checker.return_value = mock.Mock()
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[]):
            with mock.patch('esmvaltool.cmor.fix._get_cmor_checker',
                            return_value=checker) as get_mock:
                fix_metadata(self.cube, 'short_name', 'project',
                             'model', 'cmor_table', 'mip')
                get_mock.assert_called_once_with(automatic_fixes=True,
                                                 fail_on_error=False,
                                                 mip='mip',
                                                 short_name='short_name',
                                                 table='cmor_table')
                checker.assert_called_once_with(self.cube)
                checker.return_value.check_metadata.assert_called_once_with()


class TestFixData(unittest.TestCase):
    """Fix data tests"""

    def setUp(self):
        """Prepare for testing"""
        self.cube = mock.Mock()
        self.fixed_cube = mock.Mock()
        self.mock_fix = mock.Mock()
        self.mock_fix.fix_data.return_value = self.fixed_cube

    def test_fix(self):
        """Check that the returned fix is applied"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[self.mock_fix]):
            cube_returned = fix_data(self.cube, 'short_name', 'project',
                                     'model')
            self.assertTrue(cube_returned is not self.cube)
            self.assertTrue(cube_returned is self.fixed_cube)

    def test_nofix(self):
        """Check that the same cube is returned if no fix is available"""
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[]):
            cube_returned = fix_data(self.cube, 'short_name', 'project',
                                     'model')
            self.assertTrue(cube_returned is self.cube)
            self.assertTrue(cube_returned is not self.fixed_cube)

    def test_cmor_checker_called(self):
        """Check that the cmor check is done"""
        checker = mock.Mock()
        checker.return_value = mock.Mock()
        with mock.patch('esmvaltool.cmor._fixes.fix.Fix.get_fixes',
                        return_value=[]):
            with mock.patch('esmvaltool.cmor.fix._get_cmor_checker',
                            return_value=checker) as get_mock:
                fix_data(self.cube, 'short_name', 'project',
                         'model', 'cmor_table', 'mip')
                get_mock.assert_called_once_with(automatic_fixes=True,
                                                 fail_on_error=False,
                                                 mip='mip',
                                                 short_name='short_name',
                                                 table='cmor_table')
                checker.assert_called_once_with(self.cube)
                checker.return_value.check_data.assert_called_once_with()
