import logging
import os
from distutils.version import LooseVersion

import iris

from ._version import __version__

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def get_script_root():
    """Return the location of the ESMValTool installation."""
    return os.path.abspath(os.path.dirname(__file__))


def use_legacy_iris():
    """Return True if legacy iris is used."""
    return LooseVersion(iris.__version__) < LooseVersion("2.0.0")
