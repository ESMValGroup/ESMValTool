#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculations of the ETCCDI Climate Change Indices. (http://etccdi.pacificclimate.org/list_27_indices.shtml)"""
import logging
import os
import sys
#from pprint import pformat
import numpy as np
import iris
from extreme_events_utils import sum_perc_ex_wd, numdaysyear_wrapper, var_perc_ex, select_value, __nonzero_mod__, _check_required_variables, gsl_check_specs, gsl_aggregator, merge_SH_NH_cubes#, __nonzero_mod__
from cf_units import Unit
import yaml
#import dask.dataframe as dd
import dask.array as da
#import datetime
#from calendar import monthrange

import esmvalcore.preprocessor as prep

#from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
#                                            select_metadata, sorted_metadata)
#from esmvaltool.diag_scripts.shared._base import (
#    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
#from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

index_definition = yaml.load("""
annual_number_of_frost_days:
    name: frost days
    period: annual
    required:
        - tasmin
    threshold:
        value: 273.15
        unit: K
        logic: less
    cf_name: number_of_days_with_air_temperature_below_freezing_point
annual_number_of_summer_days:
    name: summer days
    period: annual
    required:
        - tasmax
    threshold:
        value: 298.15
        unit: K
        logic: greater
    cf_name: number_of_days_with_air_temperature_above_25_degree_Celsius
annual_number_of_icing_days:
    name: icing days
    period: annual
    required:
        - tasmax
    threshold:
        value: 273.15
        unit: K
        logic: less
    cf_name: number_of_days_where_air_temperature_remains_below_freezing_point
annual_number_of_tropical_nights:
    name: tropical nights
    period: annual
    required:
        - tasmin
    threshold:
        value: 293.15
        unit: K
        logic: greater
    cf_name: number_of_days_where_air_temperature_remains_above_20_degre_Celsius
annual_number_of_days_where_cumulative_precipitation_is_above_10_mm:
    name: R10mm
    period: annual
    required:
        - pr
    threshold:
        value: 10
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_10_mm
annual_number_of_days_where_cumulative_precipitation_is_above_20_mm:
    name: R20mm
    period: annual
    required:
        - pr
    threshold:
        value: 20
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_20_mm
annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm:
    name: R{}mm
    period: annual
    required:
        - pr
    threshold:
        unit: mm day-1
        logic: greater_equal
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_{}_mm
monthly_maximum_value_of_daily_maximum_temperature:
    name: TXx
    period: monthly
    required:
        - tasmax
    logic: max
    cf_name: monthly_maximum_value_of_daily_maximum_temperature
monthly_maximum_value_of_daily_minimum_temperature:
    name: TNx
    period: monthly
    required:
        - tasmin
    logic: max
    cf_name: monthly_maximum_value_of_daily_minimum_temperature
monthly_minimum_value_of_daily_maximum_temperature:
    name: TXn
    period: monthly
    required:
        - tasmax
    logic: min
    cf_name: monthly_minimum_value_of_daily_maximum_temperature
monthly_minimum_value_of_daily_minimum_temperature:
    name: TNn
    period: monthly
    required:
        - tasmin
    logic: min
    cf_name: monthly_minimum_value_of_daily_minimum_temperature
monthly_maximum_1day_precipitation:
    name: Rx1day
    period: monthly
    required:
        - pr
    spell:
        value: 1
        unit: day
    logic: max
    cf_name: monthly_maximum_1day_precipitation
monthly_maximum_5day_precipitation:
    name: Rx5day
    period: monthly
    required:
        - pr
    spell:
        value: 5
        unit: day
    logic: max
    cf_name: monthly_maximum_5day_precipitation
annual_total_precipitation_in_wet_days:
    name: PRCPTOT
    period: annual
    required:
        - pr
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days
annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%:
    name: R95pTOT
    period: annual
    required:
        - pr
    threshold:
        value: 95
        unit: percent
        logic: greater
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%
annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%:
    name: R99pTOT
    period: annual
    required:
        - pr
    threshold:
        value: 99
        unit: percent
        logic: greater
    logic: sum
    cf_name: annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%
monthly_number_of_days_where_daily_minimum_temperature_below_10%:
    name: TN10p
    period: monthly
    required:
        - tasmin
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: monthly_number_of_days_where_daily_minimum_temperature_below_10%
monthly_number_of_days_where_daily_minimum_temperature_above_90%:
    name: TN90p
    period: monthly
    required:
        - tasmin
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: monthly_number_of_days_where_daily_minimum_temperature_above_90%
monthly_number_of_days_where_daily_maximum_temperature_below_10%:
    name: TX10p
    period: monthly
    required:
        - tasmax
    threshold:
        value: 10
        unit: percent
        logic: less
    cf_name: monthly_number_of_days_where_daily_maximum_temperature_below_10%
monthly_number_of_days_where_daily_maximum_temperature_above_90%:
    name: TX90p
    period: monthly
    required:
        - tasmax
    threshold:
        value: 90
        unit: percent
        logic: greater
    cf_name: monthly_number_of_days_where_daily_maximum_temperature_above_90%
daily_temperature_range:
    name: DTR
    period: daily
    required:
        - tasmin
        - tasmax
    logic: diff
    cf_name: daily_temperature_range
annual_growing_season_length:
    name: GSL
    required:
        - tas
    start:
        threshold:
            value: 278.15
            unit: K
            logic: greater
        spell:
            value: 6
            unit: days
            logic: equal
        time:
            delay: 0
            len: 6
            unit: month
    end:
        threshold:
            value: 278.15
            unit: K
            logic: less
        spell:
            value: 6
            unit: days
            logic: equal
        time:
            delay: 6
            len: 6
            unit: month
    spatial_subsets:
        NH:
            latitude: [0, 90]
            longitude: [0, 360]
        SH:
            latitude: [0, -90]
            longitude: [0, 360]
    cf_name: annual_growing_season_length
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
        "annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm":
            "rnnmmETCCDI_yr",
        "monthly_maximum_value_of_daily_maximum_temperature":
            "txxETCCDI_m",
        "monthly_maximum_value_of_daily_minimum_temperature":
            "tnxETCCDI_m",
        "monthly_minimum_value_of_daily_maximum_temperature":
            "txnETCCDI_m",
        "monthly_minimum_value_of_daily_minimum_temperature":
            "tnnETCCDI_m",
        "monthly_maximum_1day_precipitation":
            "rx1dayETCCDI_m",
        "monthly_maximum_5day_precipitation":
            "rx5dayETCCDI_m",
        "annual_total_precipitation_in_wet_days":
            "prcptot",
        "daily_temperature_range":
            "dtr",
        "annual_growing_season_length": "gslETCCDI_yr",
        "monthly_number_of_days_where_daily_minimum_temperature_below_10%":
            "tn10pETCCDI_m",
        "monthly_number_of_days_where_daily_maximum_temperature_above_90%":
            "tx90pETCCDI_m",
        "monthly_number_of_days_where_daily_minimum_temperature_above_90%":
            "tn90pETCCDI_m",
        "monthly_number_of_days_where_daily_maximum_temperature_below_10%":
            "tx10pETCCDI_m",
        "annual_total_precipitation_in_wet_days_where_daily_precipitation_above_99%":
            "r99ptotETCCDI_yr",
        "annual_total_precipitation_in_wet_days_where_daily_precipitation_above_95%":
            "r95ptotETCCDI_yr",
        }

method_index = {}
for index, method in index_method.items():
    method_index[method] = index

def fdETCCDI_yr(cubes, **kwargs):
    """FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 degC."""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube


def suETCCDI_yr(cubes, **kwargs):
    """SU, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degC."""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube


def idETCCDI_yr(cubes, **kwargs):
    """ID, Number of icing days: Annual count of days when TX (daily maximum temperature) < 0 degC."""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube


def trETCCDI_yr(cubes, **kwargs):
    """TR, Number of tropical nights: Annual count of days when TN (daily minimum temperature) > 20 degC"""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube


def r20mmETCCDI_yr(cubes, **kwargs):
    """R20mm, Annual count of days when PRCP≥ 20mm"""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube


def r10mmETCCDI_yr(cubes, **kwargs):
    """R10mm, Annual count of days when PRCP≥ 10mm"""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube

def rnnmmETCCDI_yr(cubes, **kwargs):
    """Rnnmm, Annual count of days when PRCP≥ nnmm"""

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    threshold_v = kwargs['cfg'].pop('cumprec_threshold_nn', None)
    specs['threshold']['value'] = kwargs['cfg']['cumprec_threshold_nn'] = threshold_v
    specs['name'] = specs['name'].format(threshold_v)
    specs['cf_name'] = specs['cf_name'].format(threshold_v)

    # actual calculation
    res_cube = numdaysyear_wrapper(cubes, specs)

    return res_cube

def txxETCCDI_m(alias_cubes, **kwargs):
    """TXx, monthly maximum value of daily maximum temperature."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def tnxETCCDI_m(alias_cubes, **kwargs):
    """TNx, monthly maximum value of daily minimum temperature."""
    # TODO: Here is a lot of repetition that should be cleaned.
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def txnETCCDI_m(alias_cubes, **kwargs):
    """TNx, monthly minimum value of daily maximum temperature."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def tnnETCCDI_m(alias_cubes, **kwargs):
    """TNx, monthly minimum value of daily minimum temperature."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def rx1dayETCCDI_m(alias_cubes, **kwargs):
    """Calulates the Rx1day climate index: Monthly_maximum_1day_precipitation."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def rx5dayETCCDI_m(alias_cubes, **kwargs):
    """Calulates the Rx5day climate index: Monthly_maximum_5day_precipitation."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    alias_cubes = {key: cube.rolling_window('time', iris.analysis.SUM, specs['spell']['value'])
            for key, cube in alias_cubes.items()}
    res_cube = select_value(alias_cubes, specs)
    return res_cube

def prcptot(alias_cubes, **kwargs):
    """Calculates the PRCPTOT climate index: Annual total precipitation in wet days."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    return select_value(alias_cubes, specs)

def dtr(alias_cubes, **kwargs):
    """Calculates the DTR climate index: Daily temperature range."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    _check_required_variables(specs['required'], [item.var_name for _,item in alias_cubes.items()])
    result_cube = alias_cubes['tasmax'] - alias_cubes['tasmin']
    result_cube.rename(specs['cf_name'])
    return result_cube

def tn10pETCCDI_m(alias_cubes, **kwargs):
    """TN10p, Percentage of days when TN < 10th percentile"""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    result_cube = var_perc_ex(alias_cubes, specs, kwargs['cfg'])
    return result_cube

def tx10pETCCDI_m(alias_cubes, **kwargs):
    """TX10p, Percentage of days when TX < 10th percentile"""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    result_cube = var_perc_ex(alias_cubes, specs, kwargs['cfg'])
    return result_cube

def tn90pETCCDI_m(alias_cubes, **kwargs):
    """TN90p, Percentage of days when TN > 90th percentile"""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    result_cube = var_perc_ex(alias_cubes, specs, kwargs['cfg'])
    return result_cube

def tx90pETCCDI_m(alias_cubes, **kwargs):
    """TX90p, Percentage of days when TX > 90th percentile"""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    result_cube = var_perc_ex(alias_cubes, specs, kwargs['cfg'])
    return result_cube

def r95ptotETCCDI_yr(alias_cubes, **kwargs):
    """Annual total PRCP when RR > 95th percentile."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    result_cube = sum_perc_ex_wd(alias_cubes, specs, kwargs['cfg'])
    
    return result_cube

def r99ptotETCCDI_yr(alias_cubes, **kwargs):
    """Annual total PRCP when RR > 99th percentile."""
    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]

    result_cube = sum_perc_ex_wd(alias_cubes, specs, kwargs['cfg'])
    
    return result_cube

def gslETCCDI_yr(cubes, **kwargs):
    """GSL, Growing season length: Annual (1st Jan to 31st Dec in Northern Hemisphere (NH), 1st July to 30th June in Southern Hemisphere (SH)) count between first span of at least 6 days with daily mean temperature TG>5 degC and first span after July 1st (Jan 1st in SH) of 6 days with TG<5 degC. """

    logger.info('Loading ETCCDI specifications...')
    specs = index_definition[method_index[sys._getframe().f_code.co_name]]
    
    specs = gsl_check_specs(cubes, specs)
    
    # hemispheric split  
    res_cubes = []
    for hemisphere in ['NH', 'SH']:
        logger.info('Computing hemisphere: {}'.format(hemisphere))
        cube = cubes[specs['required'][0]]
        # constrain hemispheres
        regional_constraints = specs['spatial_subsets'][hemisphere]
        cube = cube.extract(iris.Constraint(
                latitude = lambda cell: 
                    np.min(regional_constraints['latitude'])
                    <= cell <=
                    np.max(regional_constraints['latitude']), 
                longitude = lambda cell:
                    np.min(regional_constraints['longitude'])
                    <= cell <=
                    np.max(regional_constraints['longitude']),
                ))
        
        thresh_specs = {'start': specs['start'], 'end': specs['end']}
        
        # add year auxiliary coordinate
        agg = 'year'
        
        if hemisphere == 'NH':
            # possibly using time information for start and end
            if agg not in [cc.long_name for cc in cube.coords()]:
                iris.coord_categorisation.add_year(cube, 'time', name=agg)
        else:
            if agg not in [cc.long_name for cc in cube.coords()]:
                iris.coord_categorisation.add_season_year(cube, 'time', name=agg, seasons=('jasondj','fmamj'))
        
        # actual calculation
        res_cube = gsl_aggregator(cube, thresh_specs)
        
        res_cubes.append(res_cube)
        
    # combine hemispheres
    res_cube = merge_SH_NH_cubes(res_cubes)
    
    # adjust cube information
    res_cube.rename(specs['cf_name'])
    res_cube.units = Unit('days per year')

    return res_cube
