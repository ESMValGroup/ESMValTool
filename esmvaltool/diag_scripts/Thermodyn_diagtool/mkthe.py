#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:06:30 2018

@author: Valerio2

!
!     This script computes equivalent potential temperatures and 
!     temperatures representative of the sensible and latent heat
!     exchanges in the lower layers of the troposphere. Some
!     estimate of the boundary layer height and lifting condensation
!     level is also provided.
!
!     Authors: Frank Lunkeit and Valerio Lembo (University of Hamburg)
!

"""

from netCDF4 import Dataset
import numpy as np
import os
from cdo import *
cdo = Cdo()

#Temporary definitions
#modelname = 'MPI-ESM-LR'


alv     = 2.5008e6    #Latent heat of vaporization
p0      = 100000.     #reference pressure
rv      = 461.51      #Gas constant for water vapour
tmelt   = 273.15      #freezing temp.
akap    = 0.286       #Kappa (Poisson constant R/Cp)
gascon  = 287.0       #Gas constant
cpv     = 719.         #specific heat capacity for water vapour
ra1     = 610.78      #Parameter for Magnus-Teten-Formula
ra2     = 17.2693882  #for saturation vapor pressure
ra4     = 35.86       #over liquid water
hs      = 300.        #stable boundary layer height (m)
hu      = 1000.       #unstable boundary layer height (m)
ricrs   = 0.39        #Critical Richardson number for stable layer
ricru   = 0.28        #Critical Richardson number for unstable layer
bp      = 0.03        #buoyancy parameter (m*s^-2*K^-1)
g       = 9.81        #Gravity acceleration (m*s^-2)

#folder   = '/Users/Valerio2/LEC_python/data/'
#ts_file  = folder+'CMIP5_Amon_historical_{}_r1i1p1_T2Ms_ts_1990-1995.nc'.format(modelname) #fort.10
#hus_file = folder+'CMIP5_Amon_historical_{}_r1i1p1_T3M_hus_1990-1995.nc'.format(modelname) #fort.11
#tas_file = folder+'{}_tas.nc'.format(modelname) #fort.12
#ps_file  = folder+'CMIP5_Amon_historical_MPI-ESM-LR_r1i1p1_T2Ms_ps_1990-1995.nc' #fort.13
#V_file   = folder+'MPI-ESM-LR_V.nc' #fort.14
#hfss_file= folder+'CMIP5_Amon_historical_MPI-ESM-LR_r1i1p1_T2Ms_hfss_1990-1995.nc' #fort.15
#te_file  = folder+'MPI-ESM-LR_te.nc' #fort.16


def mkthe(wdir,ts_file,hus_file,tas_file,ps_file,uas_file,vas_file,hfss_file,te_file,modelname):
    
    ts_miss_file=wdir+'/ts.nc'
    removeif(ts_miss_file)
    cdo.setctomiss('0',input=ts_file,output = ts_miss_file)
    hus_miss_file=wdir+'/hus.nc'
    removeif(hus_miss_file)
    cdo.setctomiss('0',input=hus_file,output = hus_miss_file)
    tas_miss_file=wdir+'/tas.nc'
    removeif(tas_miss_file)
    cdo.setctomiss('0',input='-monmean {}'.format(tas_file),output = tas_miss_file)
    ps_miss_file=wdir+'/ps.nc'
    removeif(ps_miss_file)
    cdo.setctomiss('0',input=ps_file,output = ps_miss_file)
    V_miss_file=wdir+'/V.nc'
    removeif(V_miss_file)
    V_file=wdir+'/{}_V.nc'.format(modelname)
    removeif(V_file)
    cdo.sqrt(input='-add -sqr {} -sqr {}'.format(uas_file,vas_file),options='-b F32',output = V_file)
    cdo.setctomiss('0',input=V_file,output = V_miss_file)
    hfss_miss_file=wdir+'/hfss.nc'
    removeif(hfss_miss_file)
    cdo.setctomiss('0',input=hfss_file,output = hfss_miss_file)
    te_miss_file=wdir+'/te.nc'
    removeif(te_miss_file)
    cdo.setctomiss('0',input=te_file,output = te_miss_file)
    
    dataset0 = Dataset(ts_miss_file)
    ts    = dataset0.variables['ts'][:, :, :]
    lats  = dataset0.variables['lat'][:]
    nlat  = len(lats)
    lons  = dataset0.variables['lon'][:]
    nlon  = len(lons)
    time  = dataset0.variables['time'][:]
    ntime = len(time)
    
    dataset = Dataset(hus_miss_file)
    hus     = dataset.variables['hus'][:, :, :, :]
    lev     = dataset.variables['plev'][:]
    nlev    = len(lev)
    dataset = Dataset(tas_miss_file)
    tas     = dataset.variables['tas'][:, :, :]
    dataset = Dataset(ps_miss_file)
    ps      = dataset.variables['ps'][:, :, :]
    dataset = Dataset(V_miss_file)
    V       = dataset.variables['uas'][:, :, :]
    dataset = Dataset(hfss_miss_file)
    hfss    = dataset.variables['hfss'][:, :, :]
    dataset = Dataset(te_miss_file)
    te    = dataset.variables['rlut'][:, :, :]
    
    #print(lev)
    #print(np.shape(hus))
    huss = hus[:,0,:,:]
    huss = np.where(lev[0] >= ps, huss, 0.)  
    #ps_m = np.where(lev[0] <= ps, ps, 0.)  
    #print(lev[0])
    for l in range(nlev):
        y=np.where(lev[l] <= ps)  
        aux = hus[:,l,:,:]
        #print(lev[l])
        #aux = np.where((huss == 0.), aux, 0.)
        aux = np.where((ps >= lev[l]), aux, 0.)       
        huss = huss + aux
    
    ricr = ricru
    h    = hu
    ricr = np.where(hfss>=0.75,ricr,ricrs)
    h    = np.where(hfss>=0.75,h   ,hs)
    
#  !
#  !get tlcl from Magnus formula (as in PlaSim)
#  !
    e      = huss*ps/(huss+gascon/rv)          #!Water vapour pressure
    td_inv = (1/tmelt)-(rv/alv)*np.log(e/ra1)  #!Dewpoint temperature
    td     = 1/td_inv
    hlcl   = 125.*(ts-td)                     #!Empirical formula for LCL height
                     
#  !
#  !Negative heights are replaced by the height of the stable
#  !boundary layer (lower constraint to the height of the cloud layer)
#  !
    hlcl   = np.where(hlcl >= 0.,hlcl,h)
    
    cp=gascon/akap
    ztlcl=ts-(g/cp)*hlcl
    
#  !
#  !Compute the pseudo-adiabatic lapse rate to obtain the height of cloud top knowing emission temperature
#  !
    gw=(g/cp)*(1+((alv*huss)/(gascon*ztlcl))/(1+((alv**2*huss*0.622)/(cp*gascon*ztlcl**2))))
#   hlcl=-(ztlcl-ts)/gw
#   hlcl   = np.where(hlcl >= 0.,hlcl,h)
    htop=-(te-ztlcl)/gw+hlcl
#  !
#  !Compute equivalent potential temperature (optional output)
#  !
   # thes=ths*np.exp((alv*huss)/(cp*ztlcl))
#  !
#  !Use potential temperature and critical Richardson number to compute
#  !temperature and height of the boundary layer top
#  !
    ths=ts*(p0/ps)**akap
    thz=ths+0.03*ricr*(V)**2/h
    pz=ps*np.exp((-g*h)/(gascon*ts))  # Barometric equation 
    tz=thz*(p0/pz)**(-akap)
    
    nc_attrs, nc_dims, nc_vars = ncdump(dataset0,'ts',True)
    
    nc_f = wdir+'/tlcl.nc'.format(modelname)
    removeif(nc_f)
    w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean LCL temperature from {} model. Calculated by Thermodynamics model diagnostics \
                            in ESMValTool. Author Valerio Lembo, Meteorologisches Institut, Universität Hamburg.".format(modelname)
    w_nc_fid.createDimension('time', None)
    w_nc_dim = w_nc_fid.createVariable('time', dataset0.variables['time'].dtype,\
                                       ('time',))
    for ncattr in dataset0.variables['time'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['time'].getncattr(ncattr))
    # Assign the dimension data to the new NetCDF file.
    w_nc_fid.variables['time'][:] = time
    w_nc_fid.createDimension('lat', len(lats))
    w_nc_dim = w_nc_fid.createVariable('lat',dataset0.variables['lat'].dtype,\
                                       ('lat',))
    for ncattr in dataset0.variables['lat'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lat'].getncattr(ncattr))
    w_nc_fid.variables['lat'][:] = lats
    w_nc_fid.createDimension('lon', len(lons))
    w_nc_dim = w_nc_fid.createVariable('lon', dataset0.variables['lon'].dtype,\
                                       ('lon',))
    for ncattr in dataset0.variables['lon'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lon'].getncattr(ncattr))
    w_nc_fid.variables['lon'][:] = lons
    w_nc_var = w_nc_fid.createVariable('tlcl', 'f8', ('time','lat','lon'))
    w_nc_var.setncatts({'long_name': u"LCL Temperature",'units': u"K", 'level_desc': u"surface",\
                    'var_desc': u"LCL temperature from LCL height (Magnus formulas \
                    and dry adiabatic lapse ratio",'statistic': 'monthly mean'})
    w_nc_fid.variables['tlcl'][:] = ztlcl
    w_nc_fid.close()  # close the new file
    
    nc_f = wdir+'/tabl.nc'.format(modelname)
    removeif(nc_f)
    w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean temperature at BL top for {} model. Calculated by Thermodynamics model diagnostics \
                            in ESMValTool. Author Valerio Lembo, Meteorologisches Institut, Universität Hamburg.".format(modelname)
    w_nc_fid.createDimension('time', None)
    w_nc_dim = w_nc_fid.createVariable('time', dataset0.variables['time'].dtype,\
                                       ('time',))
    for ncattr in dataset0.variables['time'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['time'].getncattr(ncattr))
    # Assign the dimension data to the new NetCDF file.
    w_nc_fid.variables['time'][:] = time
    w_nc_fid.createDimension('lat', len(lats))
    w_nc_dim = w_nc_fid.createVariable('lat', dataset0.variables['lat'].dtype,\
                                       ('lat',))
    for ncattr in dataset0.variables['lat'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lat'].getncattr(ncattr))
    w_nc_fid.variables['lat'][:] = lats
    w_nc_fid.createDimension('lon', len(lons))
    w_nc_dim = w_nc_fid.createVariable('lon', dataset0.variables['lon'].dtype,\
                                       ('lon',))
    for ncattr in dataset0.variables['lon'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lon'].getncattr(ncattr))
    w_nc_fid.variables['lon'][:] = lons
    w_nc_var = w_nc_fid.createVariable('tabl', 'f8', ('time','lat','lon'))
    w_nc_var.setncatts({'long_name': u"Temperature at BL top",'units': u"K", 'level_desc': u"surface",\
                    'var_desc': u"Temperature at the Boundary Layer top, from boundary layer thickness \
                    and barometric equation",'statistic': u'monthly mean'})
    w_nc_fid.variables['tabl'][:] = tz
    w_nc_fid.close()  # close the new file
    
    nc_f = wdir+'/htop.nc'.format(modelname)
    removeif(nc_f)
    w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean height of the BL top for {} model. Calculated by Thermodynamics model diagnostics \
                            in ESMValTool. Author Valerio Lembo, Meteorologisches Institut, Universität Hamburg.".format(modelname)
    w_nc_fid.createDimension('time', None)
    w_nc_dim = w_nc_fid.createVariable('time', dataset0.variables['time'].dtype,\
                                       ('time',))
    for ncattr in dataset0.variables['time'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['time'].getncattr(ncattr))
    # Assign the dimension data to the new NetCDF file.
    w_nc_fid.variables['time'][:] = time
    w_nc_fid.createDimension('lat', len(lats))
    w_nc_dim = w_nc_fid.createVariable('lat', dataset0.variables['lat'].dtype,\
                                       ('lat',))
    for ncattr in dataset0.variables['lat'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lat'].getncattr(ncattr))
    w_nc_fid.variables['lat'][:] = lats
    w_nc_fid.createDimension('lon', len(lons))
    w_nc_dim = w_nc_fid.createVariable('lon', dataset0.variables['lon'].dtype,\
                                       ('lon',))
    for ncattr in dataset0.variables['lon'].ncattrs():
        w_nc_dim.setncattr(ncattr, dataset0.variables['lon'].getncattr(ncattr))
    w_nc_fid.variables['lon'][:] = lons
    w_nc_var = w_nc_fid.createVariable('htop', 'f8', ('time','lat','lon'))
    w_nc_var.setncatts({'long_name': u"Height at BL top",'units': u"m", 'level_desc': u"surface",\
                    'var_desc': u"Height at the Boundary Layer top, from boundary layer thickness \
                    and barometric equation",'statistic': u'monthly mean'})
    w_nc_fid.variables['htop'][:] = htop
    w_nc_fid.close()  # close the new file
    
    #return ztlcl,tz,htop

def ncdump(nc_fid,key,verb):
    """
    Prints the NetCDF file attributes for a given key

    Parameters
    ----------
    key : unicode
        a valid netCDF4.Dataset.variables key
    """

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    
    return nc_attrs, nc_dims, nc_vars

    
def removeif(filename):
    try:
        os.remove(filename)
    except OSError:
	pass