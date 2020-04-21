#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculations of the ETCCDI Climate Change Indices. (http://etccdi.pacificclimate.org/list_27_indices.shtml)"""
import logging
import os
import sys
from pprint import pformat
import numpy as np
import iris
from extreme_events_utils import numdaysyear_base
from cf_units import Unit
import yaml

#from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
#                                            select_metadata, sorted_metadata)
#from esmvaltool.diag_scripts.shared._base import (
#    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
#from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

index_definition = yaml.load("""
annual_number_of_frost_days:
    name: frost days
    required:
        - tasmin
    threshold:
        value: 273.15
        unit: K
        logic: lt
    cf_name: number_of_days_with_air_temperature_below_freezing_point
annual_number_of_summer_days:
    name: summer days
    required:
        - tasmax
    threshold:
        value: 298.15
        unit: K
        logic: gt
    cf_name: number_of_days_with_air_temperature_above_25_degree_Celsius
annual_number_of_icing_days:
    name: icing days
    required:
        - tasmax
    threshold:
        value: 273.15
        unit: K
        logic: lt
    cf_name: number_of_days_where_air_temperature_remains_below_freezing_point
annual_number_of_tropical_nights:
    name: tropical nights
    required:
        - tasmin
    threshold:
        value: 293.15
        unit: K
        logic: gt
    cf_name: number_of_days_where_air_temperature_remains_above_20_degre_Celsius
annual_number_of_days_where_cumulative_precipitation_is_above_10_mm:
    name: R10mm
    required:
        - pr
    threshold:
        value: 10
        unit: mm day-1
        logic: ge
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_10_mm
annual_number_of_days_where_cumulative_precipitation_is_above_20_mm:
    name: R20mm
    required:
        - pr
    threshold:
        value: 20
        unit: mm day-1iris.coord_categorisation.add_year(cube, 'time', name=agg)
        logic: ge
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_20_mm
""")
print("INDEX_DEFINITION:")
print(yaml.dump(index_definition))

index_method = {
        "annual_number_of_frost_days": "fdETCCDI_yr",
        "annual_number_of_summer_days": "suETCCDI_yr",
        "annual_number_of_icing_days": "idETCCDI_yr",
        "annual_number_of_tropical_nights": "trETCCDI_yr",
        "annual_number_of_days_where_cumulative_precipitation_is_above_10_mm": 
            "r10mmETCCDI_yr",
        "annual_number_of_days_where_cumulative_precipitation_is_above_20_mm": 
            "r20mmETCCDI_yr",
        }

method_index = {}
for index, method in index_method.items():
    method_index[method] = index

def fdETCCDI_yr(cubes):
    """FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 degC."""
    
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    logger.info(specs)
    
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
        threshold = Unit.conform(threshold, specs['threshold']['unit'], c_unit)
    
    logger.info('Threshold is {} {}.'.format(threshold, c_unit))
        
    # compute index 
    res_cube = numdaysyear_base(fdcube,
                                threshold,
                                logic=specs['threshold']['logic'])
    
    # adjust variable name and unit
    res_cube.rename(specs['cf_name'])
    res_cube.units = Unit('days per year')
    
    return res_cube


def suETCCDI_yr(cube):
    """SU, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degC."""
    logger.info('Computing yearly number of summer days.')
    
    # set aggregation level
    agg = 'year'
    
    # add year auxiliary coordinate
    if agg not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name=agg)
    
    # get unit
    c_unit = cube.units
    logger.info("The cube's unit is {}.".format(c_unit))
    
    # set threshold depending on unit
    if c_unit == "K":
        threshold = 273.15 + 25
        
    elif c_unit == "deg_c":
        threshold = 25
        
    elif c_unit == "deg_f":
        threshold = 77
        
    else:
        logger.error('The unit {} can not be interpreted as a temperature unit!'.format(c_unit))
        return
    
    logger.info('Threshold is {} {}.'.format(threshold, c_unit))
        
    # compute index 
    res_cube = event_count_time(cube, threshold, logic='gt', aggregate=agg)
    
    # adjust variable name and unit
    res_cube.rename('summer_days')
    res_cube.units = Unit('days per year')
    
    return res_cube


def idETCCDI_yr(cube):
    """ID, Number of icing days: Annual count of days when TX (daily maximum temperature) < 0 degC."""
    logger.info('Computing yearly number of icing days.')
    
    # set aggregation level
    agg = 'year'
    
    # add year auxiliary coordinate
    if agg not in [cc.long_name for cc in cube.coords()]:
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
        logger.error('The unit {} can not be interpreted as a temperature unit!'.format(c_unit))
        returniris.coord_categorisation.add_year(cube, 'time', name=agg)
    
    logger.info('Threshold is {} {}.'.format(threshold, c_unit))
        
    # compute index 
    res_cube = event_count_time(cube, threshold, logic='lt', aggregate=agg)
    
    # adjust variable name and unit
    res_cube.rename('icing_days')
    res_cube.units = Unit('days per year')
    
    return res_cube