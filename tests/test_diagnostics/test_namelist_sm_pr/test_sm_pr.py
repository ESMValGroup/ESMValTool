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


class SMPRTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "namelist_sm_pr.xml"  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(SMPRTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        r.update({'pr' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'mrsos' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        #~ r.update({'evspsbl' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        return r


class TestDiagnostic(ESMValTestDiagnostic):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        #specify list with expected output files as tuples with directory and filename
        reffiles=[('plot','ro_coefficient-rel-pr_biases.png'),('plot','ro-et_coefficient_biases.png')]
        #~ models = ['MPIESM', 'MPIESM-LR'] # models as specified in testing namelist
        #~ for m in models:
            #~ reffiles.append(('plot',m + '_bias_ET_catchments.txt'))
            #~ reffiles.append(('plot',m + '_bias_precip_catchments.txt'))
            #~ reffiles.append(('plot',m + '_bias_runoff_catchments.txt'))
            #~ reffiles.append(('plot',m + '_ET_catchments.txt'))
            #~ reffiles.append(('plot',m + '_precip_catchments.txt'))
            #~ reffiles.append(('plot',m + '_runoff_catchments.txt'))
            #~ reffiles.append(('plot',m + '_sep-bias-plot-ET.png'))
            #~ reffiles.append(('plot',m + '_sep-bias-plot-precip.png'))
            #~ reffiles.append(('plot',m + '_sep-bias-plot-runoff.png'))

        T = SMPRTest(files=reffiles)  # provide subdirectory within plot_dir
        T.run_nml()
        #~ T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        #~ self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()


