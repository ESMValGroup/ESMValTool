import pytest
from evt import Namelist

def test_Namelist_init():
    n = Namelist()
    assert n.__name__ == "Namelist"

