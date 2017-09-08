# -*- coding: utf-8 -*-

# This file is part of ESMValTool


"""
Tests are implemented using *assert* statements
"""

import unittest


class Test(unittest.TestCase):

    def setUp(self):
        # implement here everything you would like to see happen BEFORE a test is executed
        pass

    def tearDown(self):
        # implement here everything you would like to see happen AFTER a test was executed
        pass

    def test_model_get_line(self):
        from esmvaltool.interface_scripts.model import Model
        M = Model('modelname', 'diagname')
        self.assertEqual(M.get_model_line(), 'modelname')

if __name__ == "__main__":
    unittest.main()


