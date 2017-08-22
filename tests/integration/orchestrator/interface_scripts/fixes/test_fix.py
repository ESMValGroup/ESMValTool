import unittest
from orchestrator.interface_scripts.fixes.fix import Fix


class TestFix(unittest.TestCase):

    def test_get_fix(self):
        from orchestrator.interface_scripts.fixes.CMIP5.BNU_ESM import ch4
        self.assertIs(Fix.get_fix('CMIP5', 'BNU-ESM', 'ch4'), ch4)

    def test_get_fix_no_project(self):
        self.assertIsNone(Fix.get_fix('BAD_PROJECT', 'BNU-ESM', 'ch4'))

    def test_get_fix_no_model(self):
        self.assertIsNone(Fix.get_fix('CMIP5', 'BAD_MODEL', 'ch4'))

    def test_get_fix_no_var(self):
        self.assertIsNone(Fix.get_fix('CMIP5', 'BNU-ESM', 'BAD_VAR'))


