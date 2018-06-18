import unittest
import tempfile

import numpy as np
from cf_units import Unit
import iris
from iris.cube import Cube
import netCDF4
import os

from esmvaltool.cmor._fixes.CMIP5.MPI_ESM_LR import pctisccp, so


class TestPctisccp2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='pctisccp', units='J')
        self.fix = pctisccp()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))


class TestSo(unittest.TestCase):
    def setUp(self):
        file_handler, self.filename = tempfile.mkstemp('.nc')
        os.close(file_handler)
        handler = netCDF4.Dataset(self.filename, 'w')
        handler.createDimension('time', 2)
        var = handler.createVariable('so', float, ('time',))
        var.units = 'psu'
        var[:] = np.array((1.0, 1.0))
        handler.close()
        self.fix = so()

    def tearDown(self):
        os.remove(self.filename)

    def test_fix_data(self):
        cube = iris.load_cube(self.filename)
        cube = self.fix.fix_metadata(cube)
        self.assertEqual(cube.data[0], 1)
        self.assertEqual(cube.data[1], 1)
        self.assertEqual(cube.units, Unit('1'))
        self.assertNotIn('invalid_units', cube.attributes)
