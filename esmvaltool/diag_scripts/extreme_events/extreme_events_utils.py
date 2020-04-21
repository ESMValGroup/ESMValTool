#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for the calculations of the ETCCDI Climate Change Indices."""
import logging
import os
from pprint import pformat
import numpy as np
import iris
from cf_units import Unit

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
    logger.info("assessing the logic '{}'".format(logic))
    
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


def __numdaysyear_base__(cube, threshold=273.15, logic='lt'):
    """Compute number of days per year for specific logic and threshold"""
    
    # set aggregation level
    agg = 'year'
    
    # add year auxiliary coordinate
    if agg not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name=agg)
        
    # calculate event count
    res_cube = __event_count_time__(cube, threshold, logic=logic, aggregate=agg)
    
    return res_cube


def numdaysyear_wrapper(cubes, specs):
    """Wrapper function for several number of days per year indices"""
    
    # test if required variable is in cubes
    fdcube = None
    for var, cube in cubes.items():
        if var in specs['required']:
            fdcube = cube
        
    if fdcube is None:
        logger.error('Cannot calculate number of {} for any of the following variables: {}'.format(
                specs['name'], cubes.keys()))
        return
    
    logger.info('Computing yearly number of {}.'.format(specs['name']))
    
    # get cube unit
    c_unit = fdcube.units
    logger.info("The cube's unit is {}.".format(c_unit))
    
    # get threshold
    threshold = specs['threshold']['value']
    
    # convert depending on unit of cube
    if not c_unit == specs['threshold']['unit']:
        threshold = Unit(specs['threshold']['unit']).convert(threshold, c_unit)
    
    logger.info('Threshold is {} {}.'.format(threshold, c_unit))
        
    # compute index 
    res_cube = __numdaysyear_base__(fdcube,
                                    threshold,
                                    logic=specs['threshold']['logic'])
    
    # adjust variable name and unit
    res_cube.rename(specs['cf_name'])
    res_cube.units = Unit('days per year')
    
    return res_cube
    