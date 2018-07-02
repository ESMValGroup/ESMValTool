"""Unit tests for :func:`esmvaltool.preprocessor._time_area.season_slice`."""

from __future__ import absolute_import, division, print_function

import iris
from iris.cube import Cube
import iris.coords
import iris.coord_categorisation
from cf_units import Unit
import numpy as np

import tests
from esmvaltool.preprocessor._time_area import extract_month


class Test(tests.Test):
    def setUp(self):
        self.cube = Cube(np.arange(1, 25), var_name='co2', units='J')
        self.cube.add_dim_coord(
            iris.coords.DimCoord(
                np.arange(15., 720., 30.),
                standard_name='time',
                units=Unit('days since 1950-01-01', calendar='gregorian')),
            0)
        iris.coord_categorisation.add_month_number(self.cube, 'time')

    def test_get_january(self):
        sliced = extract_month(self.cube, 1)
        print(sliced)
        self.assertTrue((np.array([1, 1]) ==
                         sliced.coord('month_number').points).all())
