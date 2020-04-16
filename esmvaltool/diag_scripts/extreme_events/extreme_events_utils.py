#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for the calculations of the ETCCDI Climate Change Indices."""
import logging
import os
from pprint import pformat
import numpy as np
import iris

#from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
#                                            select_metadata, sorted_metadata)
#from esmvaltool.diag_scripts.shared._base import (
#    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
#from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

def event_count_time(cube, threshold, lower=True, aggregate="full"):
    """Compute the number of events."""
    
    return

def __binary_translation(cube, threshold, lower=True):
    """Compute boolean for exceeding threshold of data in cube."""
    if lower:
        res_cube = cube < threshold
    else:
        res_cube = cube > threshold
    return res_cube