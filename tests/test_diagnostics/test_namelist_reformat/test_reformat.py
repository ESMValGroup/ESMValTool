# -*- coding: utf-8 -*-

# This file is part of ESMValTool

"""
Test script for "python dummy"
"""

import os
import unittest
import sys
sys.path.append('../..')

from esmvaltool_testlib import ESMValToolTest


class ReformatTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "namelist_reformat.xml"  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(ReformatTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        r.update({'tas' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        return r


class TestDiagnostic(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        #specify list with expected output files as tuples with directory and filename
        reffiles=[('work','tsline/tsline_tas_nomask_noanom_nodetr_historical_1980-2013.nc'),('plot','tsline/tsline_tas_nomask_noanom_nodetr_historical_1980-2013.ps')]

        T = ReformatTest(files=reffiles)  # provide subdirectory within plot_dir
        T.run_nml()
        T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()


