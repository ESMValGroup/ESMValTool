"""
Provides tests for the cvdp diagnostic

"""
import os
import pytest

from esmvaltool.diag_scripts.cvdp.cvdp_wrapper import *

#def test_setup_driver():
#    assert False
#
@pytest.fixture(scope='session')
def test_create_links(tmpdir_factory):

    cfg = dict()
    l = tmpdir_factory.mkdir("link")
    cfg['run_dir'] = os.path.join(l.dirname)

    p = tmpdir_factory.mkdir("sub").join("file_2009-2010.nc")
    p.write("Test")
    filepath =  os.path.join(p.dirname, p.basename)

    link = create_link(cfg, filepath)
    assert os.path.islink(link)
#
#def test_setup_namelist():
#    assert False
#
#def test_log_functions():
#    assert False
#
#def test_cvdp_available():
#    assert False
#
#def test_nco_available():
#    assert False

