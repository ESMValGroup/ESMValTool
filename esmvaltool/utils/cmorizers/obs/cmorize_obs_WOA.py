#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###############################################################################
## REFORMAT SCRIPT FOR WOA Oxygen climatology
###############################################################################
##
## Tier
##    2 (freely available data set other than obs4MIPs and ana4MIPs)
##
## Source
##    More Information:
##      https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
##    Reference:
##
## Last access
##    16/10/2018
##
## Download and processing instructions
##    Download:
##      https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
##    Processing:
##        this script (reformat_obs_woa_o2.py)
##
## Caveats
##    ...
##
## ############################################################################
"""

import errno
import os
import numpy as np
from cf_units import Unit
import iris


def makeLonSafe(lon):
    """
    Makes sure that the value is between 0. and 360 (up to 360,)
    """
    while True:
        if 0. <= lon < 360.:
                    return lon
        if lon < 0:
                    lon+=360.
        if lon >= 360.:
                    lon -= 360.

def makeLonSafeArr(lon):
    """
    Makes sure that the entire array is between -180 and 180.
    """
    out_lon = np.zeros_like(lon)
    if lon.ndim == 3:
      for (l,ll,lll,) , lo in np.ndenumerate(lon):
        out_lon[l,ll,lll] = makeLonSafe(lo)
      return out_lon
    if lon.ndim == 2:
      for l,lon1 in enumerate(lon):
        for ll,lon2 in enumerate(lon1):
          out_lon[l,ll] = makeLonSafe(lon2)
      return out_lon
    if lon.ndim == 1:
      for l,lon1 in enumerate(lon):
        out_lon[l] = makeLonSafe(lon1)
      return out_lon
    assert False


def fix_longitude(cube):
    """Longitude needs to be 0 < 360 and monotonic"""
    lon = cube.coord('longitude')
    data = cube.data
    #mask = cube.mask
    newlon = []
    lonpoints = makeLonSafeArr(lon.points)

    lonbounds = makeLonSafeArr(lon.bounds)
    #lon.points = lonpoints
    #lon.bounds = lonbounds

    newlon = np.arange(0., 360., 1.) + 0.5
    newlonbounds = np.zeros_like(lonbounds)
    newdata = np.zeros_like(data)
    #newmask = np.zeros_like(mask)

    londict = {}
    for l,lo in enumerate(newlon):
        a = int(np.argwhere(lonpoints == lo)[0][0])
        londict[l] = a

    print('bounds:', lonbounds.shape, newlonbounds.shape)
    for l, a in londict.items(): # l : new location, a: old location
        print('fixing:', l, '-->', a)
        newlonbounds[l] = lonbounds[a]
        newdata[...,l] = data[...,a]
        #newmask[...,l] = mask[...,a]

    print(newlonbounds)
    newlonbounds[-1,-1] = 360.
    cube.data = newdata
    #cube.mask = newmask
    lon.points = newlon
    lon.bounds = newlonbounds
    cube.attributes['geospatial_lon_min'] = 0.
    cube.attributes['geospatial_lon_max'] = 360.


def fix_time(cube, time_res):
    """Fix the time units and values to something sensible."""
    time = cube.coord('time')
    print(time.units)
    newunits = Unit('days since 2000-01-01 00:00:00', calendar='gregorian')
    time.units = newunits
    if time_res == 'Annual':
        time.points = np.array([182.5, ], dtype=np.float32)
        time.bounds = np.array([[ 0., 365., ]], dtype=np.float32)
    else:
        assert 0


def remove_coord_system(cube):
    """
    Remove the coordinate system from WOA data.

    Needed, because models don't have a coord system, and model and data
    need to have the same coord system.
    """
    cube.coord('latitude').coord_system = None
    cube.coord('longitude').coord_system = None


def apply_standard_transformations(cube, time_res):
    #cube.data.set_fill_value(1.e20)
    depth = cube.coord('depth')
    depth.standard_name = 'depth'
    depth.long_name = 'ocean depth coordinate'
    depth.var_name = 'lev'
    lat = cube.coord('latitude')
    lat.var_name = 'lat'
    lat.long_name = 'latitude'
    lon = cube.coord('longitude')
    lon.var_name = 'lon'
    lon.long_name = 'longitude'
    fix_time(cube, time_res)
    fix_longitude(cube)
    remove_coord_system(cube)
    print(cube)


def save_variable(cube, var, outdir):
        #o2_OBS_WOA_clim_TO3Y_o2_20000101-20001231.nc
    filename = var+'_WOA_L3_clim_TO3Y_'+var+'_20000101-20001231.nc'
    outpath = os.path.join(outdir, filename)
    print('save_variable:', outpath)
    iris.save(cube, outpath, ) #fill_value=1.e20)


def extract_thetao(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    t_an = 'Objectively analyzed mean fields for sea_water_temperature at standard depth levels.'
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == t_an:
                print('found t_an')
                break
    print(cube)
    std = 'sea_water_potential_temperature'
    long_name = 'Sea Water Potential Temperature'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'thetao'
    cube.units = Unit('celsius')
    cube.convert_units(Unit('kelvin'))
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'thetao', outdir)

#
def extract_so(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    s_an = 'Objectively analyzed mean fields for salinity at standard depth levels.'
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == s_an:
                print('found s_an')
                break
    print(cube)
    std = 'sea_water_salinity'
    long_name = 'Sea Water Salinity'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'so'
    cube.units = Unit('Unknown')
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'so', outdir)


def extract_o2(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    o_mn = 'Average of all unflagged interpolated values at each standard depth level for volume_fraction_of_oxygen_in_sea_water in each grid-square which contain at least one measurement.'
    o_an = 'Objectively analyzed mean fields for volume_fraction_of_oxygen_in_sea_water at standard depth levels.'
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == o_an:
                print('found o_an')
                break
    print(cube)
    std = 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
    long_name = 'Dissolved Oxygen Concentration'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'o2'
    cube.data = cube.data * 44.661 / 1000. # Convert from ml/l to mol/m^3
    cube.units = Unit('mol m-3')
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'o2', outdir)


def extract_no3(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    n_an = 'Objectively analyzed mean fields for moles_concentration_of_nitrate_in_sea_water at standard depth levels.'
    field = n_an
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == field:
                print('found field in netcdf', field)
                break
    print(cube)
    std = 'mole_concentration_of_nitrate_in_sea_water'
    long_name = 'Dissolved Nitrate Concentration'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'no3'
    cube.data = cube.data / 1000. # Convert from ml/l to mol/m^3
    cube.units = Unit('mol m-3')
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'no3', outdir)


def extract_po4(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    p_an = 'Objectively analyzed mean fields for moles_concentration_of_phosphate_in_sea_water at standard depth levels.'
    field = p_an
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == field:
                print('found field in netcdf', field)
                break
    print(cube)
    std = 'mole_concentration_of_phosphate_in_sea_water'
    long_name = 'Dissolved Phosphate Concentration'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'po4'
    cube.data = cube.data / 1000. # Convert from ml/l to mol/m^3
    cube.units = Unit('mol m-3')
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'po4', outdir)


def extract_si(inpath, outdir, time_res):
    print(inpath, outdir)
    cubes = iris.load_raw(inpath, ) #
    i_an = 'Objectively analyzed mean fields for moles_concentration_of_silicate_in_sea_water at standard depth levels.'
    field = i_an
    for cube in cubes:
        print(cube.long_name)
        if cube.long_name == field:
                print('found field in netcdf', field)
                break
    print(cube)
    std = 'mole_concentration_of_silicate_in_sea_water'
    long_name = 'Dissolved Silicate Concentration'
    cube.standard_name = std
    cube.long_name = long_name
    cube.var_name = 'si'
    cube.data = cube.data / 1000. # Convert from ml/l to mol/m^3
    cube.units = Unit('mol m-3')
    apply_standard_transformations(cube, time_res)
    save_variable(cube, 'si', outdir)


def mkdirs(outdir):
        try:
            os.makedirs(outdir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def main(project_info):
    # data can be downloaded from: https://data.nodc.noaa.gov/woa/WOA13/DATAv2/
    fields_to_make = [ 'thetao', 'so']#'no3', 'po4', 'si'] # 'o2',
    time_res = 'Annual'
    #
    if 'thetao' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Temperature_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_decav81B0_t00_01.nc')
        mkdirs(outdir)
        extract_thetao(inpath, outdir, time_res)

    #
    if 'so' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Salinity_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_decav81B0_s00_01.nc')
        mkdirs(outdir)
        extract_so(inpath, outdir, time_res)

    if 'o2' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Oxygen_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_all_o00_01.nc')
        mkdirs(outdir)
        extract_o2(inpath, outdir, time_res)
    #
    if 'no3' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Nitrate_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_all_n00_01.nc')
        mkdirs(outdir)
        extract_no3(inpath, outdir, time_res)

    #
    if 'po4' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Phosphate_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_all_p00_01.nc')
        mkdirs(outdir)
        extract_po4(inpath, outdir, time_res)

    #
    if 'si' in fields_to_make:
        basedir = '/data/sthenno1/scratch/ledm/Observations/WOA/Silicate_All_Periods'
        outdir = os.path.join(basedir, 'obsdata', 'Tier2', 'WOA')

        inpath = os.path.join(basedir, 'woa13_all_i00_01.nc')
        mkdirs(outdir)
        extract_si(inpath, outdir, time_res)


if __name__ == '__main__':
    main(None)
