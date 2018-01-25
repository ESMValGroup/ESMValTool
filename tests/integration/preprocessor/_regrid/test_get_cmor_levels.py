"""
Integration tests for the :func:
`esmvaltool.preprocessor.regrid.get_cmor_levels`
function.

"""

from __future__ import absolute_import, division, print_function

import unittest

from esmvaltool.preprocessor import _regrid


class TestGetCmorLevels(unittest.TestCase):

    def test_cmip6_alt40(self):
        self.assertListEqual(_regrid.get_cmor_levels('CMIP6_alt40'),
                             [240.0, 720.0, 1200.0, 1680.0, 2160.0, 2640.0,
                              3120.0, 3600.0, 4080.0, 4560.0, 5040.0, 5520.0,
                              6000.0, 6480.0, 6960.0, 7440.0, 7920.0, 8400.0,
                              8880.0, 9360.0, 9840.0, 10320.0, 10800.0,
                              11280.0, 11760.0, 12240.0, 12720.0, 13200.0,
                              13680.0, 14160.0, 14640.0, 15120.0, 15600.0,
                              16080.0, 16560.0, 17040.0, 17520.0, 18000.0,
                              18480.0, 18960.0])

    def test_cmip5_alt40(self):
        self.assertListEqual(_regrid.get_cmor_levels('CMIP5_Amon_plevs'),
                             [240.0, 720.0, 1200.0, 1680.0, 2160.0, 2640.0,
                              3120.0, 3600.0, 4080.0, 4560.0, 5040.0, 5520.0,
                              6000.0, 6480.0, 6960.0, 7440.0, 7920.0, 8400.0,
                              8880.0, 9360.0, 9840.0, 10320.0, 10800.0,
                              11280.0, 11760.0, 12240.0, 12720.0, 13200.0,
                              13680.0, 14160.0, 14640.0, 15120.0, 15600.0,
                              16080.0, 16560.0, 17040.0, 17520.0, 18000.0,
                              18480.0, 18960.0])
