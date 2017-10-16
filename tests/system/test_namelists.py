"""Test script to compare the output of ESMValTool against previous runs."""

import shutil
import tempfile

import pytest

from .esmvaltool_testlib import NAMELISTS, ESMValToolTest


@pytest.fixture
def output_directory():
    """Create a directory for storing ESMValTool output."""
    tmp = tempfile.mkdtemp()
    yield tmp
    shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.parametrize("namelist", NAMELISTS)
def test_namelist(output_directory, namelist):  # noqa
    """Create a test for each namelist in NAMELISTS and run those."""
    test = ESMValToolTest(
        namelist=namelist,
        output_directory=output_directory,
        ignore=['work/interface_data/*/*/*', '*log*.txt', '*.log'],
        checksum_exclude=['pdf', 'ps', 'png', 'eps', 'epsi', 'nc'])

    test.run(
        graphics=None,
        files='all',
        check_size_gt_zero=True,
        checksum_files='all',
        check_file_content=['nc'])

    assert test.sucess
