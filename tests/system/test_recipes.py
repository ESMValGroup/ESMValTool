"""Test script to compare the output of ESMValTool against previous runs."""

import pytest

from .esmvaltool_testlib import RECIPES, ESMValToolTest


@pytest.mark.parametrize("recipe", RECIPES)
def test_recipe(tmpdir, recipe):  # noqa
    """Create a test for each recipe in RECIPES and run those."""
    test = ESMValToolTest(
        recipe=recipe,
        output_directory=str(tmpdir),
        ignore=['tmp/*/*', '*log*.txt', '*.log'],
        checksum_exclude=['pdf', 'ps', 'png', 'eps', 'epsi', 'nc'])

    test.run(graphics=None,
             files='all',
             check_size_gt_zero=True,
             checksum_files='all',
             check_file_content=['nc'])

    assert test.success
