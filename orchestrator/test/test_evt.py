import pytest
import yaml
import os
from parser.evt import Namelist

TESTDATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'testdata',
)


def test_init_namelist():
    s = """
--- !Namelist
GLOBAL: True
PREPROCESS: True
MODELS: True
DIAGNOSTICS: True
"""
    n = yaml.load(s)
    assert n.GLOBAL is True
    assert n.PREPROCESS is True
    assert n.MODELS is True
    assert n.DIAGNOSTICS is True


def test_load_namelist():
    s = file('{}/namelist_test002.yaml'.format(TESTDATA_DIR), 'r')
    n = yaml.load(s)
    assert isinstance(n, Namelist)
    assert n.GLOBAL["write_plots"] == True
    assert n.PREPROCESS["pptag1"] == True
    assert n.MODELS == ["models1", "models2", "models3"]
    assert n.DIAGNOSTICS == ["diag1", "diag2", "diag3"]


def test_load_namelist_full():
    s = file('{}/namelist_test005.yml'.format(TESTDATA_DIR), 'r')
    n = yaml.load(s)
    assert isinstance(n, Namelist)
    assert n.GLOBAL["write_plots"] == True
    assert n.GLOBAL["write_netcdf"] == True
    assert n.GLOBAL["verbosity"] == 1
    assert n.GLOBAL["exit_on_warning"] == False
    assert n.GLOBAL["output_file_type"] == "ps"
    assert 'select_level' in n.PREPROCESS.keys()
    assert isinstance(n.MODELS, list)
    assert isinstance(n.DIAGNOSTICS, dict)
    for k, v in n.DIAGNOSTICS.items():
        assert 'description' in v
        assert 'scripts' in v
