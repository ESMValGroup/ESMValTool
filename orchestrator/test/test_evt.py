import pytest
import yaml
from evt import Namelist

def test_init_namelist():
    s = """
--- !Namelist
GLOBAL: True
PREPROCESS: True
MODELS: True
DIAGNOSTICS: True
"""
    ans = yaml.load(s)
    assert ans.GLOBAL is True
    assert ans.PREPROCESS is True
    assert ans.MODELS is True
    assert ans.DIAGNOSTICS is True

def test_load_namelist():
    s = file('test/testdata/namelist_test002.yaml', 'r')
    ans = yaml.load(s)
    assert isinstance(ans,Namelist) 
    assert ans.GLOBAL["write_plots"] == True
    assert ans.PREPROCESS["pptag1"] == True
    assert ans.MODELS == ["models1", "models2", "models3"]
    assert ans.DIAGNOSTICS == ["diag1", "diag2", "diag3"]

