"""
Test script for "Diurnal cycle of convection (first harmonic for precipitation)"
"""

import os
import unittest

from esmvaltool_testlib import ESMValToolTest

# set ESMValTool root directory
esmval_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')


class MyDiagnosticTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self):

        # 1) specify here the full path of the namelist to be tested
        #    (relative to ESMValTool root)
        xmlfile = "namelist_DiurnalCycle_harmonic.xml"
        nml = os.path.join(os.environ["ESMValTool_easytest_nmldirs"], xmlfile)
        # <<<<<<<<<

        # 2) define here the location of the reference directory
        #    note that it is expeced that the directory has the same name
        #    as the namelist
        refdir = os.path.join(os.environ["ESMValTool_easytest_refdirs"], xmlfile[:-4], "plots")
        # <<<<<<<<<<<<<<<
        super(MyDiagnosticTest, self).__init__(nml=nml,
                                               refdirectory=refdir,
                                               esmval_dir=esmval_dir)


class TestDiagnostic(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

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
        T = MyDiagnosticTest()
        T.run_nml()
        T.run_tests(execute=False, graphics=None, checksum_files='all', files='all')
        self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()
