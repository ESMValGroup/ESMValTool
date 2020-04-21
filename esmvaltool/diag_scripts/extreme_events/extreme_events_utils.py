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

# aggregation wrapper for iris.Cube.collapsed
agg_wrapper = {
        'full': (lambda cube: 
            iris.cube.Cube.collapsed(cube,
                                     'time',
                                     iris.analysis.SUM)),
        'year': (lambda cube: 
            iris.cube.Cube.aggregated_by(cube,
                                         'year',
                                         iris.analysis.SUM)),
        }

    
def __event_count_time__(cube, threshold, logic='lt', aggregate='full'):
    """Compute the number of events."""
    
    # compute the binarized version of the cube
    bin_cube = __binary_translation__(cube, threshold, logic)
    
    # aggregate within the aggregate horizon
    res_cube = agg_wrapper[aggregate](bin_cube)
    
    return res_cube


def __binary_translation__(cube, threshold, logic='lt'):
    """Compute boolean for exceeding threshold of data in cube."""
    
    #test cube against threshold and write into res_cube
    if logic == 'lt':
        thresh_data = cube.core_data() < threshold
    elif logic == 'gt':
        thresh_data = cube.core_data() > threshold
    elif logic == 'le':
        thresh_data = cube.core_data() <= threshold   
    elif logic == 'ge':
        thresh_data = cube.core_data() >= threshold  
    else:
        logger.error('The logic {} is not available!'.format(logic))
        
    #copy cube
    res_cube = cube.copy(data = thresh_data)
    
    return res_cube


def numdaysyear_base(cube, threshold=273.15, logic='l'):
    """Compute number of days per year for specific logic and threshold"""
    
    # set aggregation level
    agg = 'year'
    
    # add year auxiliary coordinate
    if agg not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name=agg)
        
    # calculate event count
    res_cube = __event_count_time__(cube, threshold, logic='lt', aggregate=agg)
    
    return res_cube
    