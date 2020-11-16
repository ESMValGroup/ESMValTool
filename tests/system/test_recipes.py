"""Test script to compare the output of ESMValTool against previous runs."""

import shutil
import tempfile

import pytest

from .esmvaltool_testlib import RECIPES, ESMValToolTest


@pytest.fixture
def output_directory():
    """Create a directory for storing ESMValTool output."""
    tmp = tempfile.mkdtemp()
    yield tmp
    shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.parametrize("recipe", RECIPES)
def test_recipe(output_directory, recipe):  # noqa
    """Create a test for each recipe in RECIPES and run those."""
    test = ESMValToolTest(
        recipe=recipe,
        output_directory=output_directory,
        ignore=['tmp/*/*', '*log*.txt', '*.log'],
        checksum_exclude=['pdf', 'ps', 'png', 'eps', 'epsi', 'nc'])

    test.run(
        graphics=None,
        files='all',
        check_size_gt_zero=True,
        checksum_files='all',
        check_file_content=['nc'])

    assert test.sucess
