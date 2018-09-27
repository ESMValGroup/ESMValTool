"""Module for snow metrics."""

import os
import numpy as np
import iris

from esmvaltool.preprocessor._regrid import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean


def i_does_diags(my_files_dict):
    """
    Example of personal diagnostic function.

    Arguments:
        run - dictionary of data files

    Returns:
        whatever you want

    """
    metrics = my_files_dict.keys()
    print(metrics)

    return metrics
