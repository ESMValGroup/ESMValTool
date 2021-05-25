"""Tests for the module :mod:`esmvaltool.diag_scripts.shared._diag`."""
import pytest
from packaging import version

from esmvaltool import __version__ as current_version
from esmvaltool.diag_scripts.shared import Datasets, Variable, Variables


@pytest.mark.parametrize('cls', [Datasets, Variable, Variables])
def test_removal_deprecated_classes(cls):
    """Test if deprecated classes should be removed."""
    if version.parse(current_version) >= version.parse('2.4'):
        assert False, (f"{cls} has been deprecated in version 2.2 and should "
                       f"be removed in version 2.4 or higher (current "
                       f"version: {current_version})")
