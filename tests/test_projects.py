# -*- coding: utf-8 -*-

# This file is part of ESMValTool


"""
Tests are implemented using *assert* statements
"""

import sys
import os
import glob

import unittest

class Test(unittest.TestCase):

    def setUp(self):
        # implement here everything you would like to see happen BEFORE a test is executed

        # to allow that test find the ESMValTool modules, we add here pathes to the system path
        esmval_path = os.path.dirname(os.path.realpath(__file__)) + os.sep + '..' + os.sep
        sys.path.append(esmval_path)

    def tearDown(self):
        # implement here everything you would like to see happen AFTER a test was executed
        pass

    def test_project_init(self):
        from interface_scripts.projects import Project
        # most simple test of init
        P = Project()

    def test_cmip5_init(self):
        from interface_scripts.projects import CMIP5
        P = CMIP5()
        self.assertEqual(P.basename, 'CMIP5')


if __name__ == "__main__":
    unittest.main()


