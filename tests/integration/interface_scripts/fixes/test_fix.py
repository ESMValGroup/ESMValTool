import unittest

from iris.cube import Cube

from esmvaltool.interface_scripts.fixes.fix import Fix


class TestFix(unittest.TestCase):
    def test_get_fix(self):
        from esmvaltool.interface_scripts.fixes.CMIP5.CanESM2 import fgco2
        self.assertListEqual(
            Fix.get_fixes('CMIP5', 'CanESM2', 'fgco2'), [fgco2()])

    def test_get_fixes_with_replace(self):
        from esmvaltool.interface_scripts.fixes.CMIP5.BNU_ESM import ch4
        self.assertListEqual(Fix.get_fixes('CMIP5', 'BNU-ESM', 'ch4'), [ch4()])

    def test_get_fixes_with_generic(self):
        from esmvaltool.interface_scripts.fixes.CMIP5.CESM1_BGC import allvars, co2
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
        self.assertEqual(Fix().fix_file(filepath), filepath)
