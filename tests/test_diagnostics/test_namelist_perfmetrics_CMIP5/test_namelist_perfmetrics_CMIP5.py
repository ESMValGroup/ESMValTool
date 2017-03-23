# -*- coding: utf-8 -*-

# This file is part of ESMValTool

"""
Test script for "Performance metrics for essential climate parameters (CMIP5)"
"""

import os
import unittest
import sys
sys.path.append('../..')

from esmvaltool_testlib import ESMValToolTest, ESMValTestDiagnostic


class PerfmetricTest(ESMValToolTest):
    """
    Define class for Test metric here
    """
    def __init__(self, **kwargs):
        xmlfile = os.path.dirname(os.path.realpath(__file__)) + os.sep + "namelist_perfmetrics_CMIP5.xml"  # the realpath ensures that the file is also foudn wehn nosetests runs from an upper directory
        super(PerfmetricTest, self).__init__(nml=xmlfile, **kwargs)

    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        """
        r = {}
        rpath = self._default_input_dir
        r.update({'ta' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'tas' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'ts' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'pr' : {'method' : 'uniform', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'ua' : {'method' : 'uniform', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'va' : {'method' : 'uniform', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'sm' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'sic' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'zg' : {'method' : 'uniform', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'rlut' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'LW_CRE' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'SW_CRE' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'hus' : {'method' : 'uniform', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        r.update({'clt' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'rsut' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'od550aer' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'od870aer' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'abs550aer' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'od550lt1aer' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})
        r.update({'toz' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 2}})

        return r



class TestPerfmetric(ESMValTestDiagnostic):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_general_output(self):
        #specify list with expected output files as tuples with directory and filename
        reffiles=[('plot','test_image.png'),('plot','python_test_out.tab')]

        T = PerfmetricTest(files=reffiles)
        T.run_nml()
        #T.run_tests(execute=False, graphics=None, checksum_files=None, files='all', check_size_gt_zero=True)
        #self.assertTrue(T.sucess)


if __name__ == "__main__":
    unittest.main()
