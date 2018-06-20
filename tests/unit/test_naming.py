"""Checks to ensure that files follow the naming convention"""

import os
import unittest


class TestNaming(unittest.TestCase):
    """Test naming of files and folders"""

    def setUp(self):
        folder = os.path.join(__file__, '..', '..', '..')
        self.esmvaltool_folder = os.path.abspath(folder)

    def test_windows_reserved_names(self):
        """
        Check that no file or folder uses a Windows reserved name

        Files can not differ from a reserved name by the extension only
        """
        reserved_names = {'CON', 'PRN', 'AUX', 'NUL', 'COM1', 'COM2', 'COM3',
                          'COM4', 'COM5', 'COM6', 'COM7', 'COM8', 'COM9',
                          'LPT1', 'LPT2', 'LPT3', 'LPT4', 'LPT5', 'LPT6',
                          'LPT7', 'LPT8', 'LPT9'}

        for dirpath, dirnames, filenames in os.walk(self.esmvaltool_folder):
            error_msg = 'Reserved name found at {}'.format(dirpath)
            self.assertTrue(reserved_names.isdisjoint(dirnames), error_msg)
            self.assertTrue(reserved_names.isdisjoint(filenames), error_msg)
            without_extensions = (os.path.splitext(filename)[0]
                                  for filename in filenames)
            self.assertTrue(reserved_names.isdisjoint(without_extensions),
                            error_msg)

    def test_avoid_case_related_collisions(self):
        """
        Check that there are no names differing only in the capitalization

        This includes folders differing from files
        """
        for dirpath, dirnames, filenames in os.walk(self.esmvaltool_folder):
            self.assertEqual(len(filenames) + len(dirnames),
                             len({name.lower() for name
                                  in filenames + dirnames}),
                             'Colliding names found at {0}'.format(dirpath))
