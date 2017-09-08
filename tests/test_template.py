# -*- coding: utf-8 -*-

# This file is part of ESMValTool


"""
Tests are implemented using *assert* statements
"""

import sys
import os

import unittest


class Test(unittest.TestCase):

    def setUp(self):
        # implement here everything you would like to see happen BEFORE a test is executed
        pass

    def tearDown(self):
        # implement here everything you would like to see happen AFTER a test was executed
        pass

    def test_basic(self):
        assert True

if __name__ == "__main__":
    unittest.main()


