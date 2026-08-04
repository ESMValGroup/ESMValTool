"""Checks to ensure that files follow the naming convention"""

import os
import unittest

IGNORE = {
    ".git",
    ".github",
    ".eggs",
    ".pixi",
    "ESMValTool.egg-info",
    "__pycache__",
    "test-reports",
}


class TestNaming(unittest.TestCase):
    """Test naming of files and folders"""

    def setUp(self):
        """Prepare tests"""
        folder = os.path.join(__file__, "..", "..", "..")
        self.esmvaltool_folder = os.path.abspath(folder)

    def test_no_namelist(self):
        """
        Check that there are no namelist references in file and folder names

        This will help us to avoid bad merges with stale branches
        """
        exclude_paths = ["esmvaltool/diag_scripts/cvdp/cvdp"]

        for dirpath, dirnames, filenames in os.walk(self.esmvaltool_folder):
            # we need to modify in-place dirnames so that we don't walk
            # over the contents of the dirs that need be ignored
            dirnames[:] = [dirn for dirn in dirnames if dirn not in IGNORE]
            print(dirnames)
            if any([item in dirpath for item in exclude_paths]):
                continue
            self.assertFalse(
                any(
                    "namelist" in name.lower() for name in filenames + dirnames
                ),
                f'Namelist reference found at {dirpath}. Please use "recipe" instead',
            )
