# -*- coding: utf-8 -*-

# This file is part of ESMValTool



import os
import unittest

from esmvaltool_testlib import ESMValToolTest
from xml.etree import ElementTree as et



class Test(unittest.TestCase):

    def setUp(self):
        # implement here everything you would like to see happen BEFORE a test is executed

        # to allow that test find the ESMValTool modules, we add here pathes to the system path
        #~ esmval_path = os.path.dirname(os.path.realpath(__file__)) + os.sep + '..' + os.sep
        #~ sys.path.append(esmval_path)
        pass

    def tearDown(self):
        # implement here everything you would like to see happen AFTER a test was executed
        pass

    @unittest.expectedFailure
    def test_init(self):
        nml = 'namelist_dummy_python.xml'
        E = ESMValToolTest(nml='./test_diagnostics/test_python_dummy/' + nml, files=[('plot','abcd.txt')])
        self.assertTrue(os.path.exists('.' + os.sep + 'data' + os.sep + nml + os.sep + 'output'))

    @unittest.expectedFailure
    def test_modify_nml(self):
        nml = 'namelist_dummy_python.xml'
        E = ESMValToolTest(nml='./test_diagnostics/test_python_dummy/' + nml, files=[('plot','abcd.txt')])

        # check names and output file existence
        self.assertEqual(E.nml, E._default_dir + os.sep + os.path.basename(E._nml)[:-4] + '_new.xml')
        self.assertTrue(os.path.exists(E.nml))

        # check that attributes were changed properly
        tree = et.parse(E.nml)
        self.assertEqual(tree.find('GLOBAL/wrk_dir').text, E.wrk_dir)
        self.assertEqual(tree.find('GLOBAL/plot_dir').text, E.plot_dir)
        self.assertEqual(tree.find('GLOBAL/climo_dir').text, E.climo_dir)


if __name__ == "__main__":
    unittest.main()


