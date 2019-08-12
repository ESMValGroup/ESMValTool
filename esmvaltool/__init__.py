"""ESMValTool diagnostics package."""
import os

__version__ = '2.0.0b0'


def get_script_root():
    """Return the location of the ESMValTool installation."""
    return os.path.abspath(os.path.dirname(__file__))
