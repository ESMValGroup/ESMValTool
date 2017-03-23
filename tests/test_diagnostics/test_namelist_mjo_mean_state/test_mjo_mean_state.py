# -*- coding: utf-8 -*-

# This file is part of ESMValTool

"""
Test script for "python dummy"
"""

import os
import unittest
import sys
sys.path.append('../..')

from esmvaltool_testlib import ESMValToolTest, ESMValTestDiagnostic


class MJOMeanTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "namelist_mjo_mean_state.xml"  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(MJOMeanTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        # TODO: document variable field
        r.update({'pr' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'ua' : {'method' : 'constant', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        return r





class TestDiagnostic(ESMValTestDiagnostic):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        #specify list with expected output files as tuples with directory and filename
        reffiles=[('plot','test_image.png'),('plot','python_test_out.tab')]

        T = MJOMeanTest(files=reffiles)
        T.run_nml()
        #T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        #self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()


