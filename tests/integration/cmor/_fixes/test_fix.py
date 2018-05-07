import unittest
import tempfile
import os
import shutil

from iris.cube import Cube

from esmvaltool.cmor.fix import Fix


class TestFix(unittest.TestCase):
    def setUp(self):
        """Set up temp folder"""
        self.temp_folder = tempfile.mkdtemp()

    def tearDown(self):
        """Remove temp folder"""
        shutil.rmtree(self.temp_folder)

    def test_get_fix(self):
        from esmvaltool.cmor._fixes.CMIP5.CanESM2 import fgco2
        self.assertListEqual(
            Fix.get_fixes('CMIP5', 'CanESM2', 'fgco2'), [fgco2()])

    def test_get_fixes_with_replace(self):
        from esmvaltool.cmor._fixes.CMIP5.BNU_ESM import ch4
        self.assertListEqual(Fix.get_fixes('CMIP5', 'BNU-ESM', 'ch4'), [ch4()])

    def test_get_fixes_with_generic(self):
        from esmvaltool.cmor._fixes.CMIP5.CESM1_BGC import (
            allvars, co2)
        self.assertListEqual(
            Fix.get_fixes('CMIP5', 'CESM1-BGC', 'co2'), [allvars(),
                                                         co2()])

    def test_get_fix_no_project(self):
        self.assertListEqual(
            Fix.get_fixes('BAD_PROJECT', 'BNU-ESM', 'ch4'), [])

    def test_get_fix_no_model(self):
        self.assertListEqual(Fix.get_fixes('CMIP5', 'BAD_MODEL', 'ch4'), [])

    def test_get_fix_no_var(self):
        self.assertListEqual(Fix.get_fixes('CMIP5', 'BNU-ESM', 'BAD_VAR'), [])

    def test_fix_metadata(self):
        cube = Cube([0])
        reference = Cube([0])

        self.assertEqual(Fix().fix_metadata(cube), reference)

    def test_fix_data(self):
        cube = Cube([0])
        reference = Cube([0])

        self.assertEqual(Fix().fix_data(cube), reference)

    def test_fix_file(self):
        filepath = 'sample_filepath'
        self.assertEqual(Fix().fix_file(filepath, 'preproc'), filepath)

    def test_fixed_filenam(self):
        filepath = os.path.join(self.temp_folder, 'file.nc')
        output_dir = os.path.join(self.temp_folder, 'fixed')
        os.makedirs(output_dir)
        fixed_filepath = Fix().get_fixed_filepath(output_dir, filepath)
        self.assertTrue(fixed_filepath,
                        os.path.join(output_dir, 'file.nc'))
