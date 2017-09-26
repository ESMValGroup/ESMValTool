"""
Test script for "MyVar"
"""

import shutil
import tempfile
import unittest

from .esmvaltool_testlib import ESMValToolTest


class TestDiagnostic(unittest.TestCase):
    def setUp(self):
        """ Run before the test"""
        self.output_directory = tempfile.mkdtemp()

    def tearDown(self):
        """ Clean up after the test"""
        shutil.rmtree(self.output_directory)

    def test_diagnostic(self):
        """
        this script runs the actual test

        You can set the following options below:
        execute : bool
            True: executes the namelist you specified to generate actual output
                  data which is the compared against reference data.
            False : If you only want to test the testing environment then it
                    could be useful to set this parameter to False
        graphics : bool
            check similarity of graphics
        checksum_files : [list, None, 'all']
            specifies the files for which the MD5 checksum should be calculated
            'all' : use all files found in specified reference directory
            None : do nothing
            [list] : list with existing filename to be checked
        files : [list, None, 'all']
            specifies the files for which their existance checked.
            'all' : use all files found in specified reference directory
            None : do nothing
            [list] : list with existing filename to be checked
        """
        T = ESMValToolTest(
            namelist='namelist_MyVar.yml',
            output_directory=self.output_directory)
        T.run_tests(
            execute=True, graphics=None, checksum_files='all', files='all')
        self.assertTrue(T.sucess)


if __name__ == "__main__":
    unittest.main()
