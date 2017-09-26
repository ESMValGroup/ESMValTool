import logging
import os

from main import __version__

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def get_script_root():
    """ Return the location of the ESMValTool installation."""
    return os.path.dirname(__file__)
