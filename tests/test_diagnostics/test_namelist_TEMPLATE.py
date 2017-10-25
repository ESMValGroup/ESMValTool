# -*- coding: utf-8 -*-

# This file is part of ESMValTool

"""
Test script 
"""

import os
import unittest
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + '..' + os.sep + '..')
from esmvaltool_testlib import ESMValToolTest, ESMValTestDiagnostic


class DiagnosticTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xml_single = 'namelist_myfirsttest.xml'
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + xml_single  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(DiagnosticTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        r.update({'pr' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'mrsos' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        return r


class TestDiagnostic(ESMValTestDiagnostic):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        # specify list with expected output files as tuples with directory and filename
        reffiles=[('plot','ro_coefficient-rel-pr_biases.png'),('plot','ro-et_coefficient_biases.png')]
        # ALTERNATIVE 1: example if output files based on model names need to be processed
        # models = ['MPIESM', 'MPIESM-LR'] # models as specified in testing namelist
        # for m in models:
            # reffiles.append(('plot',m + '_bias_ET_catchments.txt'))
            # reffiles.append(('plot',m + '_bias_precip_catchments.txt'))
            # reffiles.append(('plot',m + '_bias_runoff_catchments.txt'))
            # reffiles.append(('plot',m + '_ET_catchments.txt'))
            # reffiles.append(('plot',m + '_precip_catchments.txt'))
            # reffiles.append(('plot',m + '_runoff_catchments.txt'))
            # reffiles.append(('plot',m + '_sep-bias-plot-ET.png'))
            # reffiles.append(('plot',m + '_sep-bias-plot-precip.png'))
            # reffiles.append(('plot',m + '_sep-bias-plot-runoff.png'))
        # ALTERNATIVE 2: read reference files directly from a separate ASCII file
        # reffiles = self.read_reffiles('myreffiles.txt')

        T = DiagnosticTest(files=reffiles)  
        T.run_nml()
        T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()


