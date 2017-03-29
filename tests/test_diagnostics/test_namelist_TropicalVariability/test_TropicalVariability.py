# -*- coding: utf-8 -*-

# This file is part of ESMValTool

"""
Test script for "python dummy"
"""

import os
import unittest
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep + '..' + os.sep + '..')
from esmvaltool_testlib import ESMValToolTest, ESMValTestDiagnostic


class TropicalVariabilityTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "namelist_TropicalVariability.xml"  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(TropicalVariabilityTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        r.update({'ts' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'pr' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'ua' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'va' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})

        return r




class TestDiagnostic(ESMValTestDiagnostic):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        #specify list with expected output files as tuples with directory and filename
        reffiles=[]


        # generate list of expected reference files
        models = ['MPI-ESM-LR', 'HadGEM2-ES', 'IPSL-CM5A-MR'] # models as specified in testing namelist
        start_year = 1998
        stop_year = 2004

        # subdir TropicalVariability
        subdir = 'TropicalVariability'
        fname = subdir + os.sep + 'TropicalVariability_tspr-mmday_T2Ms_Atlantic-seasonal-mean.png'
        reffiles.append(('plot', fname))
        fname = subdir + os.sep + 'TropicalVariability_tspr-mmday_T2Ms_Indian-seasonal-mean.png'
        reffiles.append(('plot', fname))
        fname = subdir + os.sep + 'TropicalVariability_tspr-mmday_T2Ms_Pacific-seasonal-mean.png'
        reffiles.append(('plot', fname))

        for m in models:
            fname = subdir + os.sep + 'TropicalVariability_tspr-mmday_T2Ms_scatter_CMIP5_' + m + '_Amon_amip_r1i1p1_' + str(start_year) + '-' + str(stop_year) + '.png'
            reffiles.append(('plot', fname))

        # subdir TropicalVariability_EQ
        subdir = 'TropicalVariability_EQ'
        fname = subdir + os.sep + 'TropicalVariability_EQ_ua-1000_T2Ms_Atlantic-equatorial-mean.png'
        reffiles.append(('plot', fname))
        fname = subdir + os.sep + 'TropicalVariability_EQ_ua-1000_T2Ms_Indian-equatorial-mean.png'
        reffiles.append(('plot', fname))
        fname = subdir + os.sep + 'TropicalVariability_EQ_ua-1000_T2Ms_Pacific-equatorial-mean.png'
        reffiles.append(('plot', fname))


        # subdir TropicalVariability_wind
        subdir = 'TropicalVariability_wind'
        for m in models:
            for r in ['Atlantic', 'Indian', 'Pacific']:
                fname = subdir + os.sep + 'TropicalVariability_wind_ua-925va-925_T2Ms_' + r + '_CMIP5_' + m + '_Amon_amip_r1i1p1_' + str(start_year) + '-' + str(stop_year) + '.png'
                reffiles.append(('plot', fname))


        T = TropicalVariabilityTest(files=reffiles)  # provide subdirectory within plot_dir
        T.run_nml()
        T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        self.assertTrue(T.sucess)

if __name__ == "__main__":
    unittest.main()


