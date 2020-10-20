#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for the calculations of the ETCCDI Climate Change Indices."""
import logging
import os
import sys
from pprint import pformat
import numpy as np
import iris
from cf_units import Unit
import esmvalcore.preprocessor
import dask.array as da
import dask.dataframe as dd
import datetime
import collections
from scipy.stats.mstats import mquantiles


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


def __event_count_time__(cube, threshold, logic='less', aggregate='full'):
    """Compute the number of events."""

    # compute the binarized version of the cube
    bin_cube = __boolean_translation__(cube, threshold, logic)

    # aggregate within the aggregate horizon
    res_cube = agg_wrapper[aggregate](bin_cube)

    return res_cube


def __boolean_translation__(cube, threshold, logic='less'):
    """Compute boolean for exceeding threshold of data in cube."""
    logger.info("assessing the logic '{}'".format(logic))

    #test cube against threshold and write into res_cube
    thresh_data = getattr(da, logic)(cube.core_data(),threshold)
    thresh_data = thresh_data.astype(int)

    #copy cube
    res_cube = cube.copy(data = thresh_data)

    return res_cube


def __cumsum_of_boolean__(cube, threshold, logic='less'):
    """Compute cumsum for exceeding threshold of data in cube."""
    # binarize cube
    bin_cube = __boolean_translation__(cube, threshold, logic=logic)

    # calculate cumsum with setback after first 0
    cs_data = __cumsum_with_setback_multiD__(bin_cube.core_data(), 0)
    cs_cube = bin_cube.copy(data=cs_data)

    return cs_cube


def __cumsum_with_setback__(array):
    """Compute cumsum with setback after first 0"""
    cs_not_array = da.logical_not(array).cumsum()
    dataframe = dd.from_dask_array(da.stack([array, cs_not_array],
                                            axis=1),
                                   columns=['data', 'groups'])
    
    cs_df = dataframe.groupby('groups').cumsum()
    return cs_df


def __cumsum_without_setback__(array):
    """Compute cumsum"""
    return array.cumsum(axis=0)


def __cumsum_with_setback_multiD__(arrayMD, axis):
    """Compute cumsum with setback after first 0 along an axis"""
    return da.apply_along_axis(__cumsum_with_setback__,
                               axis, arrayMD,
                               shape=(arrayMD.shape[axis],),
                               dtype=int)


def __cumsum_without_setback_multiD__(arrayMD, axis):
    """Compute cumsum along an axis"""
    return da.apply_along_axis(__cumsum_without_setback__,
                               axis, arrayMD,
                               shape=(arrayMD.shape[axis],),
                               dtype=int)

def threshold_span(cube, threshold_specs, span_specs):
    """Gets the binary for any exceeding of thresholds within the respective span"""

    # cumulative summation of the booleans depending of threshold
    cum_bin = __cumsum_of_boolean__(cube, threshold=threshold_specs['value'], logic=threshold_specs['logic'])

    # boolean translation of the cummation depending on span
    span_bin = __boolean_translation__(cum_bin, threshold=span_specs['value'], logic=span_specs['logic'])

    return span_bin

def max_span_yr(cube, threshold_specs, agg):
    """Gets the maximum value for a cube per year exceeding the threshold"""

    threshold_specs['value'] = __adjust_threshold__(threshold_specs, cube.units)

    if agg != 'year':
        raise Exception(f"Period {agg} not implemented.")
    
    cum_bin_per_agg = []
    # cumulative summation of the booleans depending of threshold for each agg increment
    for years in np.unique(cube.coord('year').points):
        loc_cube = cube.extract(iris.Constraint(year = lambda cell: cell == years))
        cum_bin_per_agg.append(__cumsum_of_boolean__(loc_cube, threshold=threshold_specs['value'], logic=threshold_specs['logic']))
        
    cum_bin = iris.cube.CubeList(cum_bin_per_agg).concatenate_cube()

    # calculation of maximum span    
    max_span = cum_bin.aggregated_by(agg, iris.analysis.MAX)

    return max_span


def __first_appearence_data__(data, fill_value):
    """Get the first appearence of True and assign the respective time"""

    time_rec = __argmax_wrapper__(data, fill_value)

    time_rec = da.ma.masked_where(time_rec==fill_value, time_rec)

    return time_rec


def __argmax_wrapper__(array, fill_value):
    """wrapper with check on false data"""
    ret = da.argmax(array)
    if not array[ret]:
        ret=fill_value
    
    return ret


def __threshold_over_span__(array, threshold, span, fill_value):
    """find threhsold over span"""
    array_thresh = getattr(da, threshold['logic'])(array, threshold['value'])
    array_sum = __cumsum_with_setback__(array_thresh)
    array_span = getattr(da, span['logic'])(array_sum, span['value'])
    array_span = array_span * 1
    
    return array_span


def __adjust_masked_values__(d1, d2, year_length):
    """adjust for nonending or nonstarting growing seasons"""
    
    # no start scenario
    if da.ma.getmaskarray(d1) and not da.ma.getmaskarray(d2):
        d1 = d2.copy()+1
    # no end scenario
    if da.ma.getmaskarray(d2) and not da.ma.getmaskarray(d1):
        d2 = da.array(year_length)
        
    return d1, d2


def __gsl_per_year_per_ts__(array, split_pnt, specs, fill_value):
    """calculate the growing season per pixel for a annual timeseries"""
    
    start_span = __threshold_over_span__(array[:split_pnt], specs['start']['threshold'], specs['start']['spell'], fill_value).to_dask_array(lengths=True).squeeze()
    end_span = __threshold_over_span__(array[split_pnt:], specs['end']['threshold'], specs['end']['spell'], fill_value).to_dask_array(lengths=True).squeeze()

    # first appearence with shift to first day of first appearence
    first_appearence = __first_appearence_data__(start_span, fill_value)+2-specs['start']['spell']['value']
    
    end_span_cumusm = __cumsum_without_setback__(end_span)
    
    # first appearence with shift to first day of first appearence
    last_appearence = __first_appearence_data__(end_span_cumusm, fill_value)+2+split_pnt#-span_end['value']
    
    first_appearence, last_appearence = __adjust_masked_values__(first_appearence, last_appearence, len(array))
    
    res = last_appearence - first_appearence + 1
    
    #TODO: decide on tradeoff between memory and processing time
#    return res 
    return res.compute()
    

def __gsl_applied_MD__(array_MD, axis, split_pnt, specs, fill_value):
    """calculate the growing season for a multipixel annual timeseries"""
    
    annual_gsl = da.apply_along_axis(__gsl_per_year_per_ts__, axis, array_MD, split_pnt, specs, fill_value, shape = (1,), dtype = int)
    
    return da.reshape(annual_gsl, annual_gsl.shape[0:2])


def merge_SH_NH_cubes(loCubes, fill_value = -9999):
    """merges a list of 2 cubes (northern and southern hemisphere) with overlaping timeseries"""
    which_longer = np.argmax([lC.shape[0] for lC in loCubes])
    
    # search the longer time series
    longer = loCubes.pop(which_longer)
    shorter = loCubes[0]
    
    # adjust time coordinate and attach missing year segments
    longer_time = longer.coord('year')
    shorter_time = shorter.coord('year')
    
    longer_time_pts = longer_time.points
    shorter_time_pts = shorter_time.points
    
    add_segments = longer_time[[i for i,ltp in enumerate(longer_time_pts) if not ltp in shorter_time_pts]]
    
    additional_segments = [shorter]
    for i in add_segments:
        loc_slice = shorter[0,:,:].copy() * 0 + fill_value
        loc_slice.metadata = shorter.metadata
        loc_slice.replace_coord(i)
        loc_slice = iris.util.new_axis(loc_slice, 'year')
        additional_segments.append(loc_slice)
        
    shorter = iris.cube.CubeList(additional_segments).concatenate_cube()
    
    # adjust latitude coordinate and delete duplicate latitude segments
    longer_lat = longer.coord('latitude')
    shorter_lat = shorter.coord('latitude')
    
    longer_lat_pts = longer_lat.points
    shorter_lat_pts = shorter_lat.points
    
    del_segments_ind = [i for i,llp in enumerate(longer_lat_pts) if llp in shorter_lat_pts]
    
    
    if len(del_segments_ind) > 0:
        shorter = shorter.extract(iris.Constraint(latitude = lambda cell: cell not in del_segments_ind))
    
    return iris.cube.CubeList([shorter, longer]).concatenate_cube()


def gsl_aggregator(cube, specs):
    """ calculate the growing season for a cube per year"""
    GSL = iris.analysis.Aggregator('gsl_aggregator',
                                    __gsl_applied_MD__,
                                    lazy_func=__gsl_applied_MD__
                                    )
    
    if 'mon' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_month_number(cube, 'time', name='mon')
    
    res_cubes = []
    
    for years in np.unique(cube.coord('year').points):
        loc_cube = cube.extract(iris.Constraint(year = lambda cell: cell == years))
        months = loc_cube.coord('mon').points
        if np.all([ m in months for m in np.arange(1,13)]):
            mid_mon = months[0] + specs['end']['time']['delay'] - 1
            mid_mon_index = [i for i, e in enumerate(months) if e == mid_mon][-1]
            res_cube = loc_cube.collapsed('year', GSL, split_pnt=mid_mon_index, specs=specs, fill_value=-9999)
            res_cube.remove_coord('time')
            res_cube.remove_coord('mon')
            res_cubes.append(res_cube)
    
    return iris.cube.CubeList(res_cubes).merge_cube()
    

def __numdaysyear_base__(cube, threshold=0.0, logic='less'):
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

    logger.info('Computing annual number of {}.'.format(specs['name']))

    # get cube unit
    c_unit = fdcube.units
    logger.info("The cube's unit is {}.".format(c_unit))

    # get threshold
    threshold = __adjust_threshold__(specs['threshold'], c_unit)
    
    logger.info('Threshold is {} {}.'.format(threshold, c_unit))

    # compute index
    res_cube = __numdaysyear_base__(fdcube,
                                    threshold,
                                    logic=specs['threshold']['logic'])

    # adjust variable name and unit
    res_cube.rename(specs['cf_name'])
    res_cube.units = Unit('days per year')

    return res_cube


def __adjust_threshold__(specs_threshold, unit):
    """adjusts a threshold to the requested unit (if possible) and returns value only"""
    threshold = specs_threshold['value']
    
    # convert depending on unit
    if not unit == specs_threshold['unit']:
        threshold = Unit(specs_threshold['unit']).convert(threshold, unit)
        
    return threshold


def convert_ETCCDI_units(cube):
    """Converting cubes' units if necessary"""
    if cube.standard_name == 'precipitation_flux' and cube.units == 'kg m-2 s-1':
        cube.data = cube.data * 24. * 3600.
        cube.units = Unit('mm day-1')
    if cube.units == 'K':
        cube.convert_units('celsius')
    return cube


def select_value(alias_cubes, specs):
    """Select value per period according to given logical operator."""
    logger.info(f"Computing the {specs['cf_name'].replace('_',' ')}.")
    _check_required_variables(specs['required'], [item.var_name for _,item in alias_cubes.items()])
    # Here we assume that only one variable is required.
    cube = [item for _,item in alias_cubes.items() if item.var_name in specs['required']].pop()
    if 'period' not in specs.keys():
        raise Exception(f"Period needs to be specified.")
    statistic_function = getattr(esmvalcore.preprocessor, f"{specs['period']}_statistics", None)
    if statistic_function:
        result_cube = statistic_function(cube, specs['logic'])
    else:
        raise Exception(f"Period {specs['period']} not implemented.")
    result_cube.rename(specs['cf_name'])
    return result_cube

def _check_required_variables(required, available):
    missing = [item for item in required if item not in available]
    if len(missing):
        raise Exception(f"Missing required variable {' and '.join(missing)}.")

def regional_constraint_from_specs(spec):
    ###DOES NOT WORK###
    """Converting regional specs into iris constraints"""
    constraints = {}
    for st_set, st_def in spec['spatial_subsets'].items():
        logger.info(st_set)
        constraints[st_set] = iris.Constraint(
                latitude = lambda cell: np.min(st_def['latitude']) <= cell <= np.max(st_def['latitude']), 
                longitude = lambda cell: np.min(st_def['longitude']) <= cell <= np.max(st_def['longitude']),
                )
        logger.info(np.min(st_def['latitude']))
        logger.info(np.max(st_def['latitude']))
    logger.info(constraints)
    return constraints


def __nonzero_mod__(x, mod):
    """Compute modulo without 0 (e.g. for months)"""
    res = x%mod
    if res == 0:
        return mod
    else:
        return res


def gsl_check_specs(cubes, specs):
    """Check the units in specs for gsl conformity"""

    if len(specs['required'])!=1:
        logger.error('Searching too many cubes (should be one):'.format(
                specs['required']))
        raise Exception(f'Wrong data.')

    if specs['required'][0] not in cubes.keys():
        logger.error('Required cube not available. Looking for {}.'.format(
                specs['required'][0]))
        raise Exception(f'Wrong data.')
    else:
        c_unit = cubes[specs['required'][0]].units

    if not cubes[specs['required'][0]].attributes['frequency'] == Unit('day'):
        logger.error('Requires cube frequency in days, got {}.'.format(
                cubes[specs['required'][0]].attributes['frequency']))
        raise Exception(f'Wrong data.')
        

    for season in ['start', 'end']:
        specs[season]['threshold']['value'] = __adjust_threshold__(
                specs['start']['threshold'], c_unit)
        specs[season]['threshold']['unit'] =  c_unit
        if not specs[season]['spell']['unit'] == Unit('day'):
            logger.error('Requires span in days, got {}.'.format(
                    specs[season]['spell']))
            raise Exception(f'Wrong specifications.')
        if not specs[season]['time']['unit'] == Unit('month'):
            logger.error('Requires time in months, got {}.'.format(
                    specs[season]['time']))
            raise Exception(f'Wrong specifications.')

    return specs


def lazy_percentiles(cube, percentiles, dims="time"):
    """Calculate percentiles lazy"""
    perc_dict = collections.OrderedDict()
    
    pre_slice = cube.collapsed(dims,iris.analysis.MEAN)
    
    dask_data = cube.core_data()
    dask_data[da.ma.getmaskarray(dask_data)] = np.nan
    
    if dims == "time":
        dask_perc = da.apply_along_axis(mquantiles,
                                        0,
                                        dask_data,
                                        percentiles,
                                        alphap=1/3,
                                        betap=1/3)
        
    elif np.all(np.sort(dims) == np.sort(["latitude","longitude"])):
        dask_perc = da.apply_along_axis(mquantiles,
                                        1,
                                        dask_data.reshape(dask_data.shape[0],
                                                          -1),
                                        percentiles,
                                        alphap=1/3,
                                        betap=1/3)
    else:
        logger.error("other dimensions not implemented yet")
    
    dask_perc = da.ma.masked_array(dask_perc, da.isnan(dask_perc))
    
    for act_p in np.sort(percentiles):
        
        sub_cube = pre_slice.copy()
        
        if dims == "time":
            sub_data = (dask_perc[
                    [act_p == p for p in percentiles],:,:]).squeeze()
        elif np.all(np.sort(dims) == np.sort(["latitude","longitude"])):
            sub_data = (dask_perc[
                    :,[act_p == p for p in percentiles]]).squeeze()
    
        sub_cube.data = sub_data
        
        perc_dict.update({act_p:sub_cube})
        
    return perc_dict

def var_perc_ex(alias_cubes, specs, cfg):
    """Calculate the number of exceeding values above/below percentual threshold"""
    _check_required_variables(specs['required'], [item.var_name for _,item in alias_cubes.items()])
    if len(specs['required'])!=1:
        logger.error('Searching too many cubes (should be one):'.format(
                specs['required']))
        raise Exception(f'Wrong data.')
    if specs['threshold']['unit'] != 'percent':
        logger.error('Threshold has wrong unit. Expected "percent":'.format(
                specs['threshold']['unit']))
        raise Exception(f'Wrong unit.')
        
    cube = alias_cubes[specs['required'][0]]
    if 'doy' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='doy')
    if 'year' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name='year')
        
    base_range = cfg.copy().pop('base_range', None)
    analysis_range = cfg.copy().pop('analysis_range', None)
    
    if not base_range is None:
        base_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(base_range) <= cell.point.year <= np.max(base_range)))
    else: 
        base_cube = cube.copy()
        
    if not analysis_range is None:
        analysis_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(analysis_range) <= cell.point.year <= np.max(analysis_range)))
    else: 
        analysis_cube = cube.copy()
        
    base_cube_years = list(set(base_cube.coord('year').points))
    analysis_cube_years = list(set(analysis_cube.coord('year').points))
        
    overlap = [None]
    if (base_range is None) or \
       (analysis_range is None):
           overlap = list(set(analysis_cube_years).intersection(base_cube_years))
    elif (min(base_range)<min(analysis_range)<max(base_range)) or \
         (min(base_range)<max(analysis_range)<max(base_range)):
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled besides noticed overlaps.")
    else:
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled.")
    
    percentile_cubes = {}
    for y in analysis_cube_years:
        if y in overlap:
            all_percentiles_doy = []
            reduced_base_cube_years = base_cube_years.copy()
            reduced_base_cube_years.remove(y)
            for by in reduced_base_cube_years:
                loc_base_cube = __resample_cube__(base_cube, y, by)
                percentiles_doy = __calculate_base_percentiles__(loc_base_cube, specs['threshold']['value'])
                all_percentiles_doy.append(percentiles_doy)
            percentile_cubes.update({str(y): all_percentiles_doy})                
        else:
            percentiles_doy = __calculate_base_percentiles__(base_cube, specs['threshold']['value'])
            percentile_cubes.update({str(y): [percentiles_doy]})
            
    analysis_percentiles = []

    for cs in analysis_cube.slices(['latitude', 'longitude']):
        loc_cs_thresholds = []
        for loc_percentiles in percentile_cubes[str(cs.coord('year').points[0])]:
            thresh_data = getattr(da, specs['threshold']['logic'])(
                                          cs.core_data(),
                                          loc_percentiles.extract(iris.Constraint(
                                                  doy = lambda cell:
                                                      cell == cs.coord(
                                                              'doy').points[0])).core_data())
            cs_threshold = cs.copy(data=thresh_data*1.)
            cs_threshold.remove_coord('doy')
            cs_threshold = iris.util.new_axis(cs_threshold, 'time')
            cs_threshold.remove_coord('year')
            loc_cs_thresholds.append(cs_threshold)
            
        loc_cs_threshold = loc_cs_thresholds[0]
        if len(loc_cs_thresholds) > 1:
            for cube in loc_cs_thresholds[1:]:
                loc_cs_threshold += cube
            loc_cs_threshold = loc_cs_threshold/len(loc_cs_thresholds)
            
        analysis_percentiles.append(loc_cs_threshold)
        for ap_cube in analysis_percentiles:
            ap_cube.metadata = cube.metadata
            ap_cube.rename(specs['cf_name'])
            ap_cube.units = Unit('days')
    
    percentiles_threshold = iris.cube.CubeList(analysis_percentiles).concatenate_cube()
    statistic_function = getattr(esmvalcore.preprocessor, f"{specs['period']}_statistics", None)
    if statistic_function:
        result_cube = statistic_function(percentiles_threshold, 'sum')
    else:
        raise Exception(f"Period {specs['period']} not implemented.")

    # adjust cube information
    result_cube.rename(specs['cf_name'])
    result_cube.units = Unit('days')

    return result_cube

def __resample_cube__(base_cube, year, newyear):
    """exchange year by new year from base cube"""
    new_cube_l = base_cube.extract(iris.Constraint(year = lambda cell:cell < year))
    new_cube_h = base_cube.extract(iris.Constraint(year = lambda cell:cell > year))
    old_year = base_cube.extract(iris.Constraint(year = lambda cell:cell == year))
    extra = base_cube.extract(iris.Constraint(year = lambda cell:cell == newyear))
#    new_time = [datetime.datetime(year, 1, 1) + datetime.timedelta(int(doy) - 1)  for doy in extra.coord('doy').points]
    iris.util.demote_dim_coord_to_aux_coord(extra, 'time')
    time_delta = old_year.coord('time').points[0] - extra.coord('time').points[0]
    new_time = extra.coord('time')+time_delta
    extra.remove_coord('time')
    extra.add_dim_coord(new_time,0)
    extra.remove_coord('year')
    iris.coord_categorisation.add_year(extra, 'time', name='year')
    extra = extra.extract(iris.Constraint(year = lambda cell:cell == year))
    cube_list = [new_cube_l, extra, new_cube_h]
    new_cube = iris.cube.CubeList(list(filter(None, cube_list))).concatenate_cube()
    return new_cube

def __calculate_base_percentiles__(base_cube, percentilelist):
    """Calculate percentiles from base cube"""
    percentiles = []
    for doy in np.arange(1,367):
        doy_set = [__nonzero_mod__(doy-2, 365), __nonzero_mod__(doy+2, 365)]
        if doy_set[0] < doy_set[1]:
            constraint_5day = iris.Constraint(doy = lambda cell:
                doy_set[0] <= cell <= doy_set[1])
        else:
            constraint_5day = iris.Constraint(doy = lambda cell:
                doy_set[0] <= cell <= 366 or 0 <= cell <= doy_set[1])
        loc_cube = lazy_percentiles(
                base_cube.extract(constraint_5day),
                [percentilelist], dims='time')
        loc_cube = loc_cube[percentilelist]
        loc_cube.remove_coord('doy')
        loc_cube.remove_coord('height')
        new_coord = iris.coords.AuxCoord(doy, long_name='doy', units='1')
        loc_cube.add_aux_coord(new_coord)
        loc_cube = iris.util.new_axis(loc_cube, 'doy')
        loc_cube.remove_coord('time')
        percentiles.append(loc_cube)
        
    percentiles_doy = iris.cube.CubeList(percentiles).concatenate_cube()
    return percentiles_doy

def sum_perc_ex_wd(alias_cubes, specs, cfg):
    """Calculate the sum of exceeding values above/below percentual threshold [precipitation only, wet days]"""
    _check_required_variables(specs['required'], [item.var_name for _,item in alias_cubes.items()])
    if len(specs['required'])!=1:
        logger.error('Searching too many cubes (should be one):'.format(
                specs['required']))
        raise Exception(f'Wrong data.')
    if specs['threshold']['unit'] != 'percent':
        logger.error('Threshold has wrong unit. Expected "percent":'.format(
                specs['threshold']['unit']))
        raise Exception(f'Wrong unit.')
        
    cube = alias_cubes[specs['required'][0]]
    if 'doy' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='doy')
    if 'year' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name='year')
        
    base_range = cfg.copy().pop('base_range', None)
    analysis_range = cfg.copy().pop('analysis_range', None)
    
    if not base_range is None:
        base_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(base_range) <= cell.point.year <= np.max(base_range)))
    else: 
        base_cube = cube.copy()
        
    if not analysis_range is None:
        analysis_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(analysis_range) <= cell.point.year <= np.max(analysis_range)))
    else: 
        analysis_cube = cube.copy()
        
    base_cube_years = list(set(base_cube.coord('year').points))
    analysis_cube_years = list(set(analysis_cube.coord('year').points))
    base_cube.data = da.ma.masked_less(base_cube.core_data(), 1)
    
    overlap = [None]
    if (base_range is None) or \
       (analysis_range is None):
           overlap = list(set(analysis_cube_years).intersection(base_cube_years))
    elif (min(base_range)<min(analysis_range)<max(base_range)) or \
         (min(base_range)<max(analysis_range)<max(base_range)):
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled besides noticed overlaps.")
    else:
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled.")
    
    percentile_cubes = {}
    for y in analysis_cube_years:
        if y in overlap:
            all_percentiles_doy = []
            reduced_base_cube_years = base_cube_years.copy()
            reduced_base_cube_years.remove(y)
            for by in reduced_base_cube_years:
                loc_base_cube = __resample_cube__(base_cube, y, by)
                percentiles_doy = __calculate_base_percentiles__(loc_base_cube, specs['threshold']['value'])
                all_percentiles_doy.append(percentiles_doy)
            percentile_cubes.update({str(y): all_percentiles_doy})                
        else:
            percentiles_doy = __calculate_base_percentiles__(base_cube, specs['threshold']['value'])
            percentile_cubes.update({str(y): [percentiles_doy]})
            
    
    
    analysis_percentiles = []
    for cs in analysis_cube.slices(['latitude', 'longitude']):
        loc_cs_thresholds = []
        for loc_percentiles in percentile_cubes[str(cs.coord('year').points[0])]:
            thresh_data = getattr(da, specs['threshold']['logic'])(
                                      cs.core_data(),
                                      loc_percentiles[specs['threshold']['value']].core_data())
            cs_threshold = cs.copy(data=da.where(da.logical_or(da.logical_not(da.ma.getdata(thresh_data)), da.ma.getmaskarray(thresh_data)),0,cs.core_data()))
            cs_threshold = iris.util.new_axis(cs_threshold, 'time')
            cs_threshold.remove_coord('year')
            loc_cs_thresholds.append(cs_threshold)
        
        loc_cs_threshold = loc_cs_thresholds[0]
        if len(loc_cs_thresholds) > 1:
            for cube in loc_cs_thresholds[1:]:
                loc_cs_threshold += cube
            loc_cs_threshold = loc_cs_threshold/len(loc_cs_thresholds)
            
        analysis_percentiles.append(loc_cs_threshold)
        for ap_cube in analysis_percentiles:
            ap_cube.metadata = cube.metadata
            ap_cube.rename(specs['cf_name'])
            ap_cube.units = Unit('mm per year')
    
    percentiles_threshold = iris.cube.CubeList(analysis_percentiles).concatenate_cube()
    statistic_function = getattr(esmvalcore.preprocessor, f"{specs['period']}_statistics", None)
    if statistic_function:
        result_cube = statistic_function(percentiles_threshold, specs['logic'])
    
    # adjust cube information
    result_cube.rename(specs['cf_name'])
    result_cube.units = Unit('mm per year')
    
    return result_cube

def spell_perc_ex_thresh(alias_cubes, specs, cfg):
    """Calculate the number of days if 6+ days have exceeding values above/below percentual threshold"""

    _check_required_variables(specs['required'], [item.var_name for _,item in alias_cubes.items()])
    if len(specs['required'])!=1:
        logger.error('Searching too many cubes (should be one):'.format(
                specs['required']))
        raise Exception(f'Wrong data.')
    if specs['threshold']['unit'] != 'percent':
        logger.error('Threshold has wrong unit. Expected "percent":'.format(
                specs['threshold']['unit']))
        raise Exception(f'Wrong unit.')
        
    cube = alias_cubes[specs['required'][0]]
    if 'doy' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='doy')
    if 'year' not in [cc.long_name for cc in cube.coords()]:
        iris.coord_categorisation.add_year(cube, 'time', name='year')
        
    base_range = cfg.copy().pop('base_range', None)
    analysis_range = cfg.copy().pop('analysis_range', None)
    
    if not base_range is None:
        base_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(base_range) <= cell.point.year <= np.max(base_range)))
    else: 
        base_cube = cube.copy()
        
    if not analysis_range is None:
        analysis_cube = cube.extract(iris.Constraint(time = lambda cell:
            np.min(analysis_range) <= cell.point.year <= np.max(analysis_range)))
    else: 
        analysis_cube = cube.copy()
        
    base_cube_years = list(set(base_cube.coord('year').points))
    analysis_cube_years = list(set(analysis_cube.coord('year').points))
        
    overlap = [None]
    if (base_range is None) or \
       (analysis_range is None):
           overlap = list(set(analysis_cube_years).intersection(base_cube_years))
    elif (min(base_range)<min(analysis_range)<max(base_range)) or \
         (min(base_range)<max(analysis_range)<max(base_range)):
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled besides noticed overlaps.")
    else:
           logger.warning("Both, base temporal range and analysis temporal" +
                       " range are specified. \n" +
                       "Bootstrapping is thus canceled.")
    
    percentile_cubes = {}
    for y in analysis_cube_years:
        if y in overlap:
            all_percentiles_doy = []
            reduced_base_cube_years = base_cube_years.copy()
            reduced_base_cube_years.remove(y)
            for by in reduced_base_cube_years:
                loc_base_cube = __resample_cube__(base_cube, y, by)
                percentiles_doy = __calculate_base_percentiles__(loc_base_cube, specs['threshold']['value'])
                all_percentiles_doy.append(percentiles_doy)
            percentile_cubes.update({str(y): all_percentiles_doy})                
        else:
            percentiles_doy = __calculate_base_percentiles__(base_cube, specs['threshold']['value'])
            percentile_cubes.update({str(y): [percentiles_doy]})
            
    analysis_percentiles = []

    for cs in analysis_cube.slices(['latitude', 'longitude']):
        loc_cs_thresholds = []
        for loc_percentiles in percentile_cubes[str(cs.coord('year').points[0])]:
            thresh_data = getattr(da, specs['threshold']['logic'])(
                                          cs.core_data(),
                                          loc_percentiles.extract(iris.Constraint(
                                                  doy = lambda cell:
                                                      cell == cs.coord(
                                                              'doy').points[0])).core_data())
            cs_threshold = cs.copy(data=thresh_data*1.)
            cs_threshold.remove_coord('doy')
            cs_threshold = iris.util.new_axis(cs_threshold, 'time')
            cs_threshold.remove_coord('year')
            loc_cs_thresholds.append(cs_threshold)
            
        loc_cs_threshold = loc_cs_thresholds[0]
        if len(loc_cs_thresholds) > 1:
            for cube in loc_cs_thresholds[1:]:
                loc_cs_threshold += cube
            loc_cs_threshold = loc_cs_threshold/len(loc_cs_thresholds)
            
        analysis_percentiles.append(loc_cs_threshold)
        for ap_cube in analysis_percentiles:
            ap_cube.metadata = cube.metadata
            ap_cube.rename(specs['cf_name'])
            ap_cube.units = Unit('days per year')
    
    percentiles_threshold = iris.cube.CubeList(analysis_percentiles).concatenate_cube()
    
    perc_thresh_cube = percentiles_threshold.copy(data=__cumsum_with_setback_multiD__(percentiles_threshold.core_data(), 0))
    
    if 'year' not in [cc.long_name for cc in perc_thresh_cube.coords()]:
        iris.coord_categorisation.add_year(perc_thresh_cube, 'time', name='year')
    
    spells_logic = __boolean_translation__(perc_thresh_cube, specs['spell']['value'], logic=specs['spell']['logic'])
    spells_equal = __boolean_translation__(perc_thresh_cube, specs['spell']['value'], logic='equal')
    
    spells_inc = spells_logic + 5 * spells_equal
    
    statistic_function = getattr(esmvalcore.preprocessor, f"{specs['period']}_statistics", None)
    if statistic_function:
        result_cube = statistic_function(spells_inc, 'sum')
    else:
        raise Exception(f"Period {specs['period']} not implemented.")
        
    # adjust cube information
    result_cube.rename(specs['cf_name'])
    result_cube.units = Unit('days per year')
    
    return result_cube

def add_filename(cube, fun):
    """ adds a filename to the cube attributes based on the cube attributes """
    # add attributes
    
#    try:
#        dty_start = cube.coords("year")[0].points[0]
#        dty_end = cube.coords("year")[0].points[-1]
#    except:
#        dts = cube.coords("time")[0]
#        dty_start = dts.units.num2date(dts.points[0]).year
#        dty_end = dts.units.num2date(dts.points[-1]).year

    cube.attributes.update(
        {"FI_index": fun,
#         "FI_model": cube.attributes['model_id'],
#         "FI_experiment": cube.attributes['experiment_id'],
#         "FI_ensemble": cube.attributes['parent_experiment_rip'],
#         "FI_temporal_coverage": "{}-{}".format(dty_start, dty_end),
                })
#    cube.attributes.update({"FI_filename":
#        "{}_{}_{}_{}_{}".format(cube.attributes["FI_index"],
#         cube.attributes["FI_model"],
#         cube.attributes["FI_experiment"],
#         cube.attributes["FI_ensemble"],
#         cube.attributes["FI_temporal_coverage"],)})
    