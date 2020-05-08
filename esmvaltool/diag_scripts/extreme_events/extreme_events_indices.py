#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculations of the ETCCDI Climate Change Indices. (http://etccdi.pacificclimate.org/list_27_indices.shtml)"""
import logging
import os
import sys
from pprint import pformat
import numpy as np
import iris
from extreme_events_utils import numdaysyear_wrapper, select_value
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
    period: annual
    required:
        - tasmin
    threshold:
        value: 273.15
        unit: K
        logic: lt
    cf_name: number_of_days_with_air_temperature_below_freezing_point
annual_number_of_summer_days:
    name: summer days
    period: annual
    required:
        - tasmax
    threshold:
        value: 298.15
        unit: K
        logic: gt
    cf_name: number_of_days_with_air_temperature_above_25_degree_Celsius
annual_number_of_icing_days:
    name: icing days
    period: annual
    required:
        - tasmax
    threshold:
        value: 273.15
        unit: K
        logic: lt
    cf_name: number_of_days_where_air_temperature_remains_below_freezing_point
annual_number_of_tropical_nights:
    name: tropical nights
    period: annual
    required:
        - tasmin
    threshold:
        value: 293.15
        unit: K
        logic: gt
    cf_name: number_of_days_where_air_temperature_remains_above_20_degre_Celsius
annual_number_of_days_where_cumulative_precipitation_is_above_10_mm:
    name: R10mm
    period: annual
    required:
        - pr
    threshold:
        value: 10
        unit: mm day-1
        logic: ge
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_10_mm
annual_number_of_days_where_cumulative_precipitation_is_above_20_mm:
    name: R20mm
    period: annual
    required:
        - pr
    threshold:
        value: 20
        unit: mm day-1
        logic: ge
    cf_name: annual_number_of_days_where_cumulative_precipitation_is_above_20_mm
annual_number_of_days_where_cumulative_precipitation_is_above_nn_mm:
    name: R{}mm
    period: annual
    required:
        - pr
    threshold:
        unit: mm day-1
        logic: ge
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
            "prcptot"
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

