"""Provides tests for the cvdp diagnostic."""

import os

import pytest

from esmvaltool.diag_scripts.cvdp.cvdp_wrapper import create_link


# def test_setup_driver():
#     assert False
#
@pytest.fixture(scope='session')
def test_create_links(tmpdir_factory):
    """Test create_link function."""
    cfg = dict()
    link_dir = tmpdir_factory.mkdir("link")
    cfg['run_dir'] = os.path.join(link_dir.dirname)

    testfile = tmpdir_factory.mkdir("sub").join("file_2009-2010.nc")
    testfile.write("Test")
    filepath = os.path.join(testfile.dirname, testfile.basename)

    link = create_link(cfg, filepath)
    if not os.path.islink(link):
        raise AssertionError()


#
# def test_setup_namelist():
#     assert False
#
# def test_log_functions():
#     assert False
#
# def test_cvdp_available():
#     assert False
#
# def test_nco_available():
#     assert False
