"""
*********************************************************************
APPLICATE/TRR Ocean Diagnostics
*********************************************************************
"""
import logging
import os
import joblib
from collections import OrderedDict
import iris

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

import inspect
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib as mpl
mpl.use('agg')  #noqa
import matplotlib.pylab as plt
import math
from matplotlib import cm
from netCDF4 import num2date
#import seawater as sw
from collections import OrderedDict
from cdo import Cdo
import cmocean.cm as cmo
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
import pandas as pd
import pyresample
from scipy.interpolate import interp1d
import ESMF
import pyproj

from esmvaltool.diag_scripts.ocean3d.utils import genfilename
from esmvaltool.diag_scripts.ocean3d.regions import hofm_regions, transect_points


def hofm_data(model_filenames, mmodel, cmor_var, areacello_fx,
              max_level, region, diagworkdir):

    logger.info("Extract  {} data for {}, region {}".format(cmor_var, mmodel, region))
    
    datafile = Dataset(model_filenames[mmodel])
    print(model_filenames[mmodel])
    datafile_area = Dataset(areacello_fx[mmodel])
    lon = datafile.variables['lon'][:]
    lat = datafile.variables['lat'][:]
    lev = datafile.variables['lev'][:]
    time = num2date(datafile.variables['time'][:],
                    datafile.variables['time'].units)
    if lon.ndim == 2:
        lon2d, lat2d = lon, lat
    elif lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)

    areacello = datafile_area.variables['areacello'][:]

    lev_limit = lev[lev <= max_level].shape[0] + 1

    indexesi, indexesj = hofm_regions(region, lon2d, lat2d)

    # Fix for climatology
    # ESMValTool reduces the dimentions if one of the
    # dimentions is "empty"
    if datafile.variables[cmor_var].ndim < 4:
        series_lenght = 1
    else:
        series_lenght = datafile.variables[cmor_var].shape[0]
    
    oce_hofm = np.zeros((lev[0:lev_limit].shape[0], series_lenght))
    for mon in range(series_lenght):
        # print(mon)
        for ind, depth in enumerate(lev[0:lev_limit]):
            # fix for climatology
            if datafile.variables[cmor_var].ndim < 4:
                level_pp = datafile.variables[cmor_var][ind, :, :] 
            else:
                level_pp = datafile.variables[cmor_var][mon, ind, :, :] 

            ## This is fix fo make models with 0 as missing values work,
            ## should be fixed in fixes that do not work for now in the new backend
            if not isinstance(level_pp, np.ma.MaskedArray):
                # print(mmodel)
                level_pp = np.ma.masked_equal(level_pp, 0)
            data_mask = level_pp[indexesi,indexesj].mask
            area_masked = np.ma.masked_where(data_mask, areacello[indexesi,indexesj])
            result = (area_masked*level_pp[indexesi,indexesj]).sum()/area_masked.sum()
            oce_hofm[ind, mon] = result

    ofilename = genfilename(diagworkdir, cmor_var,
                            mmodel, region, 'hofm')
    ofilename_levels = genfilename(diagworkdir, cmor_var,
                                   mmodel, region, 'levels')
    ofilename_time = genfilename(diagworkdir, cmor_var,
                                 mmodel, region, 'time')
    print(ofilename)
    np.save(ofilename, oce_hofm)
    np.save(ofilename_levels, lev[0:lev_limit])
    np.save(ofilename_time, time)
    datafile.close()

def transect_data(mmodel,  cmor_var, 
              max_level, region, diagworkdir, mult = 2, observations='PHC'):
    logger.info("Extract  {} data for {}, region {}".format(cmor_var, mmodel, region))
    ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')
    print(ifilename)
    datafile = Dataset(ifilename)
    grid = ESMF.Grid(filename=ifilename,
                 filetype=ESMF.FileFormat.GRIDSPEC
               )
    # lon = datafile_t.variables['lon'][:]
    # lat = datafile_t.variables['lat'][:]
    lev = datafile.variables['lev'][:]
    
    # if lon.ndim == 2:
    #     lon2d, lat2d = lon, lat
    # elif lon.ndim == 1:
    #     lon2d, lat2d = np.meshgrid(lon, lat)

    lev_limit = lev[lev <= max_level].shape[0] + 1
    # indexesi, indexesj = hofm_regions(region, lon2d, lat2d)
    lon_s4new, lat_s4new = transect_points(region, mult = mult)

    coord_sys=ESMF.CoordSys.SPH_DEG
    domask=True

    locstream = ESMF.LocStream(lon_s4new.shape[0], name="Atlantic Inflow Section", coord_sys=coord_sys)
    # appoint the section locations
    locstream["ESMF:Lon"] = lon_s4new
    locstream["ESMF:Lat"] = lat_s4new
    if domask:
        locstream["ESMF:Mask"] = np.array(np.ones(lon_s4new.shape[0]), dtype=np.int32)
    

    secfield = np.zeros((lon_s4new.shape[0],datafile.variables[cmor_var].shape[1]))
    if mmodel != observations:
        ndepths = datafile.variables[cmor_var].shape[1]
    else:
        ndepths = datafile.variables[cmor_var].shape[0]

    for kind in range(0, ndepths):
        print('indice = ', kind)
        # Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
    #     srcgrid = ESMF.Grid(filename=gridfile,
    #                      filetype=ESMF.FileFormat.GRIDSPEC)
        sourcefield = ESMF.Field(grid, staggerloc=ESMF.StaggerLoc.CENTER, name = 'MPI',)

        if mmodel != observations:
            thetao = datafile.variables[cmor_var][0,kind,:,:]
        else:
            thetao = datafile.variables[cmor_var][0,kind,:,:]
        print(mmodel)
        print(thetao.shape)
        if isinstance(thetao, np.ma.core.MaskedArray):
            sourcefield.data[...] = thetao.filled(0).T
        else:
            sourcefield.data[...] = thetao.T
        
        # create a field on the locstream
        dstfield = ESMF.Field(locstream, name='dstfield')
        dstfield.data[:] = 0.0

        # create an object to regrid data from the source to the destination field
        dst_mask_values=None
        if domask:
                dst_mask_values=np.array([0])

        regrid = ESMF.Regrid(sourcefield, dstfield,
                            regrid_method=ESMF.RegridMethod.NEAREST_STOD,
                            #regrid_method=ESMF.RegridMethod.BILINEAR,
                            unmapped_action=ESMF.UnmappedAction.IGNORE,
                            dst_mask_values=dst_mask_values
                            )

        # do the regridding from source to destination field
        dstfield = regrid(sourcefield, dstfield)
        secfield[:,kind] = dstfield.data

    g = pyproj.Geod(ellps='WGS84')
    (az12, az21, dist) = g.inv(lon_s4new[0:-1], lat_s4new[0:-1], lon_s4new[1:], lat_s4new[1:])
    dist = dist.cumsum()/1000
    dist = np.insert(dist, 0, 0)

    ofilename = genfilename(diagworkdir, cmor_var,
                            mmodel, region, 'transect')
    ofilename_depth = genfilename(diagworkdir, 'depth',
                            mmodel, region, 'transect')
    ofilename_dist = genfilename(diagworkdir, 'distance',
                            mmodel, region, 'transect')

    np.save(ofilename, secfield)
    np.save(ofilename_depth, lev)
    np.save(ofilename_dist, dist)

    datafile.close()

def tsplot_data(mmodel, 
              max_level, region, diagworkdir, observations = 'PHC'):
    logger.info("Extract  TS data for {}, region {}".format(mmodel, region))
    ifilename_t = genfilename(diagworkdir, 'thetao',
                             mmodel,  data_type='timmean', extension='.nc')
    ifilename_s = genfilename(diagworkdir, 'so',
                             mmodel,  data_type='timmean', extension='.nc')

    datafile_t = Dataset(ifilename_t)
    datafile_s = Dataset(ifilename_s)

    lon = datafile_t.variables['lon'][:]
    lat = datafile_t.variables['lat'][:]
    lev = datafile_t.variables['lev'][:]
    
    if lon.ndim == 2:
        lon2d, lat2d = lon, lat
    elif lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)

    lev_limit = lev[lev <= max_level].shape[0] + 1
    indexesi, indexesj = hofm_regions(region, lon2d, lat2d)

    temp = np.array([])
    salt = np.array([])
    depth_model = np.array([])
    for ind, depth in enumerate(lev[0:lev_limit]):
        if mmodel != observations:
            level_pp = datafile_t.variables['thetao'][0, ind, :, :]
            level_pp_s = datafile_s.variables['so'][0, ind, :, :]
        else:
            level_pp = datafile_t.variables['thetao'][0, ind, :, :]
            level_pp_s = datafile_s.variables['so'][0, ind, :, :]
        ## This is fix fo make models with 0 as missing values work,
        ## should be fixed in fixes that do not work for now in the new backend
        if not isinstance(level_pp, np.ma.MaskedArray):
            level_pp = np.ma.masked_equal(level_pp, 0)
            level_pp_s = np.ma.masked_equal(level_pp_s, 0)
        temp = np.hstack((temp, level_pp[indexesi,indexesj].compressed()))
        salt = np.hstack((salt, level_pp_s[indexesi,indexesj].compressed()))
        depth_temp = np.zeros_like(level_pp[indexesi,indexesj].compressed())
        depth_temp[:] = depth
        depth_model = np.hstack((depth_model, depth_temp))

    ofilename_t = genfilename(diagworkdir, 'thetao',
                            mmodel, region, 'tsplot')
    ofilename_s = genfilename(diagworkdir, 'so',
                            mmodel, region, 'tsplot')
    ofilename_depth = genfilename(diagworkdir, 'depth',
                            mmodel, region, 'tsplot')
    np.save(ofilename_t, temp)
    np.save(ofilename_s, salt)
    np.save(ofilename_depth, depth_model)

    datafile_t.close()
    datafile_s.close()

def meanminmax_data(model_filenames, mmodel, cmor_var, areacello_fx,
              min_level, max_level, region, diagworkdir):
    logger.info("Extract  {} data for {}, region {}".format(cmor_var, mmodel, region))

    datafile = Dataset(model_filenames[mmodel])
    datafile_area = Dataset(areacello_fx[mmodel])
    lon = datafile.variables['lon'][:]
    lat = datafile.variables['lat'][:]
    lev = datafile.variables['lev'][:]
    time = num2date(datafile.variables['time'][:],
                    datafile.variables['time'].units)
    if lon.ndim == 2:
        lon2d, lat2d = lon, lat
    elif lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)

    areacello = datafile_area.variables['areacello'][:]

    levels_range = lev[(lev <= max_level) & (lev >= min_level)]

    indexesi, indexesj = hofm_regions(region, lon2d, lat2d)

    series_lenght = datafile.variables[cmor_var].shape[0]


    mean_property_mon = []
    max_property_mon = []
    min_property_mon = []
    for mon in range(series_lenght):
            # print(mon)
            mean_property = []
            max_property = []
            min_property = []
            for ind, depth in enumerate(levels_range):
                level_pp = datafile.variables[cmor_var][mon, ind, :, :]
                ## This is fix fo make models with 0 as missing values work,
                ## should be fixed in fixes that do not work for now in the new backend
                if not isinstance(level_pp, np.ma.MaskedArray):
                    level_pp = np.ma.masked_equal(level_pp, 0)
                data_mask = level_pp[indexesi,indexesj].mask
                area_masked = np.ma.masked_where(data_mask, areacello[indexesi,indexesj])
                result_mean = (area_masked*level_pp[indexesi,indexesj]).sum()/area_masked.sum()
                mean_property.append(result_mean)
                result_max = np.nanmax(level_pp[indexesi,indexesj])
                max_property.append(result_max)
                result_min = np.nanmin(level_pp[indexesi,indexesj])
                min_property.append(result_min)
            mean_property_mon.append(np.nanmean(np.array(mean_property)))
            max_property_mon.append(np.nanmax(np.array(max_property)))
            min_property_mon.append(np.nanmin(np.array(min_property)))

    meanminmax = np.vstack((mean_property_mon, min_property_mon, max_property_mon ))
    
    ofilename = genfilename(diagworkdir, cmor_var,
                            mmodel, region, 'meanminmax')
    
    # ofilename_levels = genfilename(diagworkdir, cmor_var,
    #                                mmodel, region, 'meanminmax_levels')
    ofilename_time = genfilename(diagworkdir, cmor_var,
                                 mmodel, region, 'meanminmax_time')
     
    np.save(ofilename, meanminmax)
    # np.save(ofilename_levels, lev[0:lev_limit])
    np.save(ofilename_time, time)
    datafile.close()

def tsplot_data_clim(max_level, region, diagworkdir):
    logger.info("Extract  TS data from PHC climatology")
    phc = Dataset('/mnt/lustre01/work/ab0995/a270088/ESMV/obsdata/Tier2/PHC/phc3.0_annual.nc')

    temp_phc = phc.variables['temp'][:]
    salt_phc = phc.variables['salt'][:]
    lon_phc = phc.variables['lon'][:]
    lat_phc = phc.variables['lat'][:]
    depth_phc = phc.variables['depth'][:]
    lonc, latc = np.meshgrid(lon_phc, lat_phc)

    depth3d_phc = np.zeros(salt_phc.shape)
    for i in range(depth_phc.shape[0]):
        depth3d_phc[i, :, :] = depth_phc[i]

    ptemp_phc = sw.ptmp(salt_phc, temp_phc, depth3d_phc)
    temp_phc = ptemp_phc
    
    indexesi, indexesj = hofm_clim_regions(region, lonc, latc)
    
    lev_limit = depth_phc[depth_phc <= max_level].shape[0] + 1

    temp_clim = np.array([])
    salt_clim = np.array([])
    depth_clim = np.array([])
    for ind, depth in enumerate(depth_phc[0:lev_limit]):
        level_pp = np.ma.masked_invalid(temp_phc[ind, :, :])
        level_pp_s = np.ma.masked_invalid(salt_phc[ind, :, :])
        temp_clim = np.hstack((temp_clim, level_pp[indexesi,indexesj].compressed()))
        salt_clim = np.hstack((salt_clim, level_pp_s[indexesi,indexesj].compressed()))
        depth_temp = np.zeros_like(level_pp[indexesi,indexesj].compressed())
        depth_temp[:] = depth
        depth_clim = np.hstack((depth_clim, depth_temp))

    ofilename_t = genfilename(diagworkdir, 'thetao',
                            'PHC3.0', region, 'tsplot')
    ofilename_s = genfilename(diagworkdir, 'so',
                            'PHC3.0', region, 'tsplot')
    ofilename_depth = genfilename(diagworkdir, 'depth',
                            'PHC3.0', region, 'tsplot')
    np.save(ofilename_t, temp_clim)
    np.save(ofilename_s, salt_clim)
    np.save(ofilename_depth, depth_clim)

    phc.close()

    
def aw_core(model_filenames, diagworkdir, region, cmor_var):
    logger.info("Calculate AW core statistics")
    aw_core_parameters = {}

    for i, mmodel in enumerate(model_filenames):
        aw_core_parameters[mmodel] = {}
        logger.info("Plot profile {} data for {}, region {}".format(cmor_var,
                                                            mmodel,
                                                            region))
        ifilename = genfilename(diagworkdir, cmor_var, 
                                mmodel, region, 'hofm', '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var,
                                    mmodel, region, 'levels', '.npy')

        hofdata = np.load(ifilename)
        lev = np.load(ifilename_levels)
        
        profile = (hofdata)[:, :].mean(axis=1)
        maxvalue = np.max(profile[(lev >= 200) & (lev <= 1000)])
        maxvalue_index = np.where(profile==maxvalue)[0][0]
        maxvalue_depth = lev[maxvalue_index]

        aw_core_parameters[mmodel]['maxvalue'] = maxvalue
        aw_core_parameters[mmodel]['maxvalue_index'] = maxvalue_index
        aw_core_parameters[mmodel]['maxvalue_depth'] = maxvalue_depth
    
    # Hardcoded values for PHC
    # aw_core_parameters['PHC'] = {}
    # aw_core_parameters['PHC']['maxvalue'] = 0.991 + 273.15
    # aw_core_parameters['PHC']['maxvalue_index'] = 11
    # aw_core_parameters['PHC']['maxvalue_depth'] = 300.0

    return aw_core_parameters