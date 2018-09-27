"""
Look at this module for guidance how to write your own.

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvaltool.preprocessor._regrid import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean

Pipe output through logger;
"""


def i_does_diags(my_files_dict):
    """
    Example of personal diagnostic function.

    Arguments:
        run - dictionary of data files

    Returns:
        whatever you want

    """
    metrics = my_files_dict.keys()

    return metrics
