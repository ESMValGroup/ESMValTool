#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculations of the ETCCDI Climate Change Indices. (http://etccdi.pacificclimate.org/list_27_indices.shtml)"""
import logging
import os
from pprint import pformat
import numpy as np
import iris
from extreme_events_utils import event_count_time
from cf_units import Unit

#from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
#                                            select_metadata, sorted_metadata)
#from esmvaltool.diag_scripts.shared._base import (
#    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
#from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def fdETCCDI_yr(cube):
    """FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 degC."""
    logger.info("Computing yearly number of frost days.")
    
    # set aggregation level
    agg = 'year'
    
    # add year auxiliary coordinate
    iris.coord_categorisation.add_year(cube, 'time', name=agg)
    
    # get unit
    c_unit = cube.units
    logger.info("The cube's unit is {}.".format(c_unit))
    
    # set threshold depending on unit
    if c_unit == "K":
        threshold = 273.15
        
    elif c_unit == "deg_c":
        threshold = 0
        
    elif c_unit == "deg_f":
        threshold = 32
        
    else:
        logger.error("The unit {} can not be interpreted as a temperature unit!".format(c_unit))
        return
    
    logger.info("Threshold is {} {}.".format(threshold, c_unit))
        
    # compute index 
    res_cube = event_count_time(cube, threshold, aggregate=agg)
    
    # adjust variable name and unit
    res_cube.rename('frost_days')
    res_cube.units = Unit('days per year')
    
    return res_cube