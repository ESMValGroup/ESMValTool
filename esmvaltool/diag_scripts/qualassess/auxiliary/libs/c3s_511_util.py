#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:17:48 2018

"""

import imp
import os
import sys
import numpy as np
import dask
from scipy import stats
import cf_units
#import matplotlib.pyplot as plt
#import time
import iris
import collections
#from memory_profiler import profile

# sys.path.insert(0,
#                os.path.abspath(os.path.join(os.path.join(
#                        os.path.dirname(os.path.abspath(
#                                __file__)), os.pardir),
#                                os.pardir)))

#from .TempStab.TempStab import TempStab as TS
import logging

logger = logging.getLogger(os.path.basename(__file__))


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout


def __getInfoFromFile__(filename):
    """
    routine to read cfg file
    """
    f = open(filename)
    __config__ = imp.load_source('cfg', '', f)
    f.close()
    return __config__


def __remove_all_aux_coords__(cube):
    """
    remove all auxiliary coordinates from cube
    """
    for dim in cube.coords():
        if isinstance(dim, iris.coords.AuxCoord):
             cube.remove_coord(dim)
            
    return

def __minmeanmax__(array):
    """
    calculate minimum, maximum and average of array
    """
    return (np.nanmin(array), np.nanmean(array), np.nanmax(array))


def temporal_trend(cube, pthres=1.01):
    """
    calculate temporal trend of the data over time
    the slope of the temporal trend has unit [dataunit/day]
    Parameters
    ----------
    return_object : bool
        specifies if a C{Data} object shall be returned [True]
        or if a numpy array shall be returned [False]
    pthres : float
        specifies significance threshold; all values above this threshold will be masked
    Returns
    -------
    The following variables are returned:
    correlation, slope, intercept, p-value
    """
    # time difference in days and decades
    dt = np.asarray(cube.coord('time').points - cube.coord('time').points[0]
                    ).astype('float') / 365.2425 / 10
    # ensure that a masked array is used
    x = np.ma.array(dt, mask=dt != dt)

    R, S, I, P, C = __corr_single__(cube, x, pthres=pthres)

    R.long_name = cube.long_name + ' (correlation)'
    S.long_name = cube.long_name + ' ($\partial x / \partial t$)'
    I.long_name = cube.long_name + ' (offset)'
    P.long_name = cube.long_name + ' (p-value)'

    if S.units in [None, 'no_unit', '1', 'unknown']:
        S.units = cf_units.Unit('0.1 year-1')
    else:
        S.units = cf_units.Unit(str(S.units) + ' 0.1 year-1')
    R.units = '-'
    I.units = cube.units
    P.units = 1
    return R, S, I, P


def __corr_single__(cube, x, pthres=1.01):
    """
    The routine correlates a data vector with all data of the
    current object. It returns several correlation measures
    Parameters
    ----------
    cube :
    x : ndarray
        the data vector correlations should be calculated with,
        numpy array [time]
    pthres : float
        significance threshold. All values below this
        threshold will be returned as valid
    """

    if cube.ndim != 3:
        raise ValueError('Invalid geometry!')

    nt, ny, nx = cube.shape

    if nt != len(x):
        raise ValueError('Inconsistent geometries')

    # check if 'x' is a masked array where the mask has the same
    # size as the data. If this is not the case, then set it
    if isinstance(x, np.ma.core.MaskedArray):
        if isinstance(x.mask, np.ndarray):
            if x.mask.shape != x.data.shape:
                raise ValueError('Invalid mask geometry!')
        else:
            raise ValueError(
                'The mask of the dataset needs to be an array!')
    else:
        raise ValueError('Expect masked array as input in corr_single')

    # get data with at least one valid value
    dat, msk = __get_valid_data__(cube, mode='thres', thres=3)
    xx, n = dat.shape
    logger.info(('   Number of grid points: ' + str(n)))

    R = np.ones((ny, nx)) * np.nan  # output matrix for correlation
    P = np.ones((ny, nx)) * np.nan  # output matrix for p-value
    S = np.ones((ny, nx)) * np.nan  # output matrix for slope
    I = np.ones((ny, nx)) * np.nan  # output matrix for intercept
    CO = np.ones((ny, nx)) * np.nan  # output matrix for covariance

    R.shape = (-1)
    S.shape = (-1)
    P.shape = (-1)
    I.shape = (-1)
    CO.shape = (-1)

    logger.info('   Calculating correlation ...')
    res = [stats.mstats.linregress(x, dat[:, i]) for i in range(n)]

    res = np.asarray(res)
    slope = res[:, 0]
    intercept = res[:, 1]
    r_value = res[:, 2]
    p_value = res[:, 3]
#    std_err = res[:, 4]

    R[msk] = r_value
    P[msk] = p_value
    I[msk] = intercept
    S[msk] = slope
    R.shape = (ny, nx)
    P.shape = (ny, nx)
    I.shape = (ny, nx)
    S.shape = (ny, nx)

    #--- prepare output data objects
    Rout = cube[0, :, :].copy()  # copy object to get coordinates
    Rout.long_name = 'correlation'
    msk = (P > pthres) | (np.isnan(R))
    #msk = np.zeros_like(R).astype('bool')
    Rout.data = np.ma.array(R, mask=msk).copy()
    Rout.units = None

    Sout = cube[0, :, :].copy()  # copy object to get coordinates
    Sout.long_name = 'slope'
    Sout.data = np.ma.array(S, mask=msk).copy()
    Sout.units = cube.units

    Iout = cube[0, :, :].copy()  # copy object to get coordinates
    Iout.long_name = 'intercept'
    Iout.data = np.ma.array(I, mask=msk).copy()
    Iout.units = cube.units

    Pout = cube[0, :, :].copy()  # copy object to get coordinates
    Pout.long_name = 'p-value'
    Pout.data = np.ma.array(P, mask=msk).copy()
    Pout.units = 1

    Cout = cube[0, :, :].copy()  # copy object to get coordinates
    Cout.long_name = 'covariance'
    # currently not supported: covariance!
    Cout.data = np.ma.array(np.ones(P.shape) * np.nan, mask=msk).copy()
    Cout.units = None

    return Rout, Sout, Iout, Pout, Cout


def __get_valid_data__(cube, mode='all', thres=-99):
    """
    this routine calculates from the masked array
    only the valid data and returns it together with its
    coordinate as vector
    valid means that ALL timestamps need to be valid!
    Parameters
    ----------
    return_mask : bool
        specifies if the mask applied to the original data should
        be returned as well
    mode : str
        analysis mode ['all','one','thres']
        'all': all timestamps need to be valid
        'one': at least a single dataset needs to be valid
        'thres' : number of valid timesteps needs to be abovt a threshold
    thres : int
        threshold for minimum number of valid values (needed when mode=='thres')
    """

    if mode == 'thres':
        assert thres > 0, 'Threshold needs to be > 0!'

    if cube.ndim == 3:

        n = len(cube.coord("time").points)

        # vectorize the data

        data = cube.data.reshape(n, -1)

        if not np.ma.is_masked(data):
            data = np.ma.masked_array(data, np.full_like(data, False))
        # set pixels with NaN to invalid
        data.mask[np.isnan(data.data)] = True

        # extract only valid (not masked data)
        if mode == 'all':
            # identify all ngrid cells where all timesteps are valid
            msk = np.sum(~data.mask, axis=0) == n
        elif mode == 'one':
            # identify ONE grid cell where all timesteps are valid
            msk = np.sum(~data.mask, axis=0) > 0
        elif mode == 'thres':
            msk = np.sum(~data.mask, axis=0) > thres
        else:
            raise ValueError('Invalid option in get_valid_data() %s' %
                             mode)

        data = data[:, msk]

    elif cube.ndim == 2:
        data = cube.data.reshape(-1)
        msk = ~data.mask
        data = data[msk]
    else:
        raise ValueError('Unsupported dimension!')

    return data, msk


#def __loc_TSA_fun__(array, **kwargs):
#
#    breakpoint_method = kwargs.get('breakpoint_method', 'CUMSUMADJ')
#    max_num_period = kwargs.get('max_num_periods', 3)
#    periods_method = kwargs.get('periods_method', 'autocorr')
#    temporal_resolution = kwargs.get('temporal_resolution', 1.)
#    minimum_available_data_points = kwargs.get('min_avail_pts', 1)
#
#    if minimum_available_data_points < 2:
#        assert False, "No trend calculation possible for " + \
#            "less than 2 data points"
#
#    RES = None
#    done = -2
#
#    timearray = kwargs.get('dates', None)
#
#    if timearray is not None:
#        if array.mask.sum() < len(array.mask) - (minimum_available_data_points - 1):
#            try:
#                with HiddenPrints():
#                    #                    print array
#                    TSA = TS(timearray, array,
#                             breakpoint_method=breakpoint_method,
#                             detrend=True,
#                             deseason=True,
#                             max_num_periods=max_num_period,
#                             periods_method=periods_method,
#                             temporal_resolution=temporal_resolution)
#                    RES = TSA.analysis(homogenize=True)
#                    done = 2
#            except BaseException:
#                try:
#                    with HiddenPrints():
#                        TSA = TS(timearray, array,
#                                 breakpoint_method=breakpoint_method,
#                                 detrend=True,
#                                 deseason=False,
#                                 max_num_periods=max_num_period,
#                                 periods_method=periods_method,
#                                 temporal_resolution=temporal_resolution)
#                        RES = TSA.analysis(homogenize=True)
#                        done = 1
#                except BaseException:
#                    try:
#                        with HiddenPrints():
#                            TSA = TS(timearray, array,
#                                     breakpoint_method=breakpoint_method,
#                                     detrend=False,
#                                     deseason=False,
#                                     max_num_periods=max_num_period,
#                                     periods_method=periods_method,
#                                     temporal_resolution=temporal_resolution)
#                            RES = TSA.analysis(homogenize=True)
#                            done = 0
#                    except BaseException:
#                        done = -9
#        else:
#            done = -1
#    else:
#        "Error in timearray."
#
#    if RES is not None:
#        slope_diff = RES["homogenized_trend"]["slope"]
#        fin_res = np.atleast_1d(np.append(np.array([slope_diff,
#                                                    len(RES["breakpoints"]),
#                                                    done]), TSA.homogenized))
#
#    else:
#        fin_res = np.atleast_1d(np.append(np.array([np.nan, np.nan, done]),
#                                          np.ones(array.shape) * np.nan))
#
#    return fin_res


#def __TS_of_cube__(cube, **kwargs):
#
#    breakpoint_method = kwargs.get('breakpoint_method', 'CUMSUMADJ')
#    max_num_period = kwargs.get('max_num_periods', 3)
#    periods_method = kwargs.get('periods_method', 'autocorr')
#    temporal_resolution = kwargs.get('temporal_resolution', 1.)
#    minimum_available_data_points = kwargs.get('min_avail_pts', 1)
#
#    min_trend = cube[0, :, :].copy()
#    num_bp = cube[0, :, :].copy()
#    version = cube[0, :, :].copy()
#    homogenized = cube.copy()
#
#    timearray = kwargs.get('dates', None)
#
#    if timearray is None:
#        timearray = cube.coord("time").points
#
#    res = np.apply_along_axis(__loc_TSA_fun__, 0, cube.data,
#                              dates=timearray,
#                              breakpoint_method=breakpoint_method,
#                              max_num_period=max_num_period,
#                              periods_method=periods_method,
#                              min_avail_pts=minimum_available_data_points,
#                              temporal_resoution=temporal_resolution)
#
#    mask = np.isnan(res[0, :, :]) + min_trend.data.mask
#    min_trend.data = np.ma.array(data=res[0, :, :] * 365.2425 * 10, mask=mask)
#    if min_trend.units in [None, 'no_unit', '1', 'unknown']:
#        min_trend.units = cf_units.Unit('0.1 year-1')
#    else:
#        min_trend.units += cf_units.Unit(str(min_trend.units) + ' 0.1 year-1')
#
#    num_bp.data = np.ma.array(data=res[1, :, :], mask=mask)
#    num_bp.units = cf_units.Unit("1")
#
#    homogenized.data = np.ma.array(
#        data=res[3:, :, :], mask=np.isnan(res[3:, :, :]))
#
#    version.data = np.ma.array(data=res[2, :, :], mask=mask)
#    version.units = cf_units.Unit("1")
#
#    return({"slope": min_trend, "number_breakpts": num_bp,
#            "version": version, "homogenized": homogenized})


def weighted_STD_DEV(cube, dim, weights=None):

    if weights is None:
        return (cube.collapsed(dim, iris.analysis.RMS)**2 -
                cube.collapsed(dim, iris.analysis.MEAN)**2)**0.5
    else:
        return (
            cube.collapsed(
                dim,
                iris.analysis.RMS,
                weights=weights)**2 -
            cube.collapsed(
                dim,
                iris.analysis.MEAN,
                weights=weights)**2)**0.5
                    

def dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating and merging keys. The ``merge_dct`` is 
    merged into ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            if not k in dct:
                dct[k] = merge_dct[k]
            else:
                if isinstance(dct[k],list):
                    if isinstance(merge_dct[k],list):
                        dct[k] = dct[k] + merge_dct[k]
                    else:
                        dct[k] = dct[k] + [merge_dct[k]]
                else:
                    if isinstance(merge_dct[k],list):
                        dct[k] = [dct[k]] + merge_dct[k]
                    else:
                        dct[k] = [dct[k], merge_dct[k]]


def dask_weighted_mean_wrapper(cube, spatial_weights, dims=None):
    
    if not isinstance(dims,list):
        dims = [dims]
    
    # weights only necessary if latitude is aggregated
    if "latitude" in dims:
        latlon_list = []
        
        SPW = spatial_weights.compute()
        
        latc = cube.coord("latitude").points
        cpos = (len(latc) - 1.) / 2.
        latc = (latc[int(np.ceil(cpos))] + latc[int(np.floor(cpos))])/2
        latb = cube.coord("latitude").bounds
        latb = [latb[0][0],latb[-1][-1]]
        
        for latlon in cube.slices(["latitude","longitude"]):
            
            for rcoord in ["day_of_month", "day_of_year", "month_number", "year"]:
                if rcoord in [coord.name() for coord in latlon.coords()]:
                    latlon.remove_coord(rcoord)
                
            latlon_master = next(latlon.slices_over("latitude"))
                
            latlon_master.coord("latitude").points = latc
            latlon_master.coord("latitude").bounds = latb
            
            data = latlon_master.core_data() 
            
            
            data = dask.array.sum(latlon.core_data() * SPW /
                                  dask.array.sum((latlon.core_data() *
                                                  0. + 1.) * SPW, axis=0),
                                                  axis=0)
            
            latlon_master.data = data
            
            latlon_list.append(latlon_master)
                
            latlon = None
        
        dims.remove("latitude")
        
        cube_list = iris.cube.CubeList(latlon_list)
        
        new_cube = cube_list.merge_cube()
        
        del latlon_list[:]
        del latlon_list
        del cube_list
        del SPW
        
        if len(dims):
            new_cube = new_cube.collapsed(dims, iris.analysis.MEAN)
            
        dims.append("latitude")
            
    else:
        new_cube = cube.collapsed(dims, iris.analysis.MEAN)
        
    return new_cube


def dask_weighted_stddev_wrapper(cube, spatial_weights, dims=None):    
    
    weighted_mean = dask_weighted_mean_wrapper(cube,
                                               spatial_weights,
                                               dims=dims)
    
    for rcoord in ["day_of_month", "day_of_year", "month_number", "year"]:
        if (rcoord in [coord.name() for coord in weighted_mean.coords()] or
            rcoord in [coord.name() for coord in cube.coords()]):
            try:
                cube.remove_coord(rcoord)
            except:
                cube.remove_aux_factory(rcoord)
    
    if "air_pressure" in [coord.name() for coord in cube.coords()]:
        if np.all((dask.array.flip(
                weighted_mean.coord("air_pressure").core_points(),0) == 
            cube.coord("air_pressure").core_points()).compute()):
            latlontim_list = []
            for latlontim in cube.slices(["time","latitude","longitude"]):
                latlontim_list.append(latlontim)
               
            cube_list = iris.cube.CubeList(latlontim_list)
            cube = cube_list.merge_cube()

    sqdifference = (cube - weighted_mean) ** 2
    
    variance = dask_weighted_mean_wrapper(sqdifference,
                                          spatial_weights,
                                          dims=dims)
    
    std_dev = variance ** 0.5

    return std_dev


def lazy_climatology(cube, t_coord):
    
    t_coord_dict = collections.OrderedDict()
    
    for act_t in np.sort(np.unique(cube.coord(t_coord).points)):
        
        sub_cube = cube.extract(
                iris.Constraint(
                        coord_values={t_coord:lambda point: point == act_t}))
#        sub_cube = cube.extract(iris.Constraint(month_number = lambda point: point == act_t))

        sub_mean = sub_cube.collapsed("time",iris.analysis.MEAN)
            
        t_coord_dict.update({act_t:sub_mean})
    
    return t_coord_dict


def minmax_cubelist(cubelist,perc,symmetric=False):

    try: 
        vminmax = dask.array.concatenate([__dask_perc_masked__(plotcube, perc)
                   for plotcube in cubelist])
    except:
        vminmax = dask.array.concatenate([dask.array.percentile(
                plotcube.core_data().ravel(),perc)
                   for plotcube in cubelist])
    
    if symmetric:
        vminmax = dask.array.concatenate([vminmax, -vminmax])
        
    vmin = dask.array.atleast_1d(vminmax.min())
    vmax = dask.array.atleast_1d(vminmax.max())
    vminmax = dask.array.concatenate([vmin,vmax])

    return vminmax.compute()


def __dask_perc_masked__(cube,perc):
    ravelcube = cube.core_data().ravel()
    percentile = dask.array.percentile(
            ravelcube[~dask.array.ma.getmaskarray(ravelcube)],
            perc)
    return percentile


def lazy_percentiles(cube, percentiles, dims="time"):
    
    perc_dict = collections.OrderedDict()
    
    pre_slice = cube.collapsed(dims,iris.analysis.MEAN)
    
    dask_data = cube.core_data()
    dask_data[dask.array.ma.getmaskarray(dask_data)] = np.nan
    
    if dims == "time":
        dask_perc = dask.array.apply_along_axis(np.nanpercentile,
                                                0,dask_data,percentiles)
    elif np.all(np.sort(dims) == np.sort(["latitude","longitude"])):
        dask_perc = dask.array.apply_along_axis(np.nanpercentile,
                                                1,dask_data.reshape(
                                                        dask_data.shape[0],-1),
                                                        percentiles)  
    else:
        logger.error("other dimensions not implemented yet")
    
    dask_perc = dask.array.ma.masked_array(dask_perc,
                                              dask.array.isnan(dask_perc))
    
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
