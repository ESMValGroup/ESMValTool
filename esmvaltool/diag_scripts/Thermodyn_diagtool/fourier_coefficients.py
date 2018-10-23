#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 15:40:08 2018

@author: Valerio Lembo
"""

import numpy as np
from netCDF4 import Dataset

tainput = 'inputen.nc'
tasinput = 'tas.nc'
gpres=np.array([16,32,48,64,96,128,256,384,512,1024,2048,4096])
fcres=np.array([5,10,15,21,31,43,85,127,171,341,683,1365])

class Fourier_coeff():
    
    from fourier_coefficients import *
    
    def fourier_coeff(self, outfile, tainput):
        
        fourcoeff = Fourier_coeff()
        
        fileo    = outfile
        
        dataset = Dataset(tainput)
        lon   = dataset.variables['lon'][:]
        lat   = dataset.variables['lat'][:]
        lev   = dataset.variables['plev'][:]
        time  = dataset.variables['time'][:]
        nlon  = len(lon)
        nlat  = len(lat)
        nlev  = len(lev)
        ntime = len(time)
        i = np.min(np.where(2*nlat <= gpres))
        trunc=fcres[i]+1
        wave2 = np.linspace(0,trunc-1,trunc)
    
        #print(f.field.keys())
        ta=dataset.variables['ta'][:,:,:,:]
        ua=dataset.variables['ua'][:,:,:,:] 
        va=dataset.variables['va'][:,:,:,:]
        wap=dataset.variables['wap'][:,:,:,:]
        #dataset = Dataset(tasinput)
        #tas=dataset.variables['tas'][:,:,:]
    
        #print(np.shape(ta))
        tafft_p = np.fft.fft(ta,axis=3)[:,:,:,:trunc/2]/(nlon)
        tafft_p = np.transpose(tafft_p, (0, 1, 3, 2))
        uafft_p = np.fft.fft(ua,axis=3)[:,:,:,:trunc/2]/(nlon)
        uafft_p = np.transpose(uafft_p, (0, 1, 3, 2))
        vafft_p = np.fft.fft(va,axis=3)[:,:,:,:trunc/2]/(nlon)
        vafft_p = np.transpose(vafft_p, (0, 1, 3, 2))
        wapfft_p = np.fft.fft(wap,axis=3)[:,:,:,:trunc/2]/(nlon)
        wapfft_p = np.transpose(wapfft_p, (0, 1, 3, 2))
        
        tafft = np.zeros([ntime,nlev,trunc,nlat])
        uafft = np.zeros([ntime,nlev,trunc,nlat])
        vafft = np.zeros([ntime,nlev,trunc,nlat])
        wapfft = np.zeros([ntime,nlev,trunc,nlat])
        tafft[:,:,0::2,:]=np.real(tafft_p)
        tafft[:,:,1::2,:]=np.imag(tafft_p)
        uafft[:,:,0::2,:]=np.real(uafft_p)
        uafft[:,:,1::2,:]=np.imag(uafft_p)
        vafft[:,:,0::2,:]=np.real(vafft_p)
        vafft[:,:,1::2,:]=np.imag(vafft_p)
        wapfft[:,:,0::2,:]=np.real(wapfft_p)
        wapfft[:,:,1::2,:]=np.imag(wapfft_p)
            
        fourcoeff.pr_output(tafft,uafft,vafft,wapfft,tainput,fileo,wave2, 
                            'ta','ua','va','wap',verb=True)
        
        
    def pr_output(self, var1, var2, var3, var4, nc_f, fileo, wave2, name1, name2, name3, name4, verb=True):
        '''
        NAME
            NetCDF with Python
        PROGRAMMER(S)
            Chris Slocum
        REVISION HISTORY
            20140320 -- Initial version created and posted online
            20140722 -- Added basic error handling to ncdump
                        Thanks to K.-Michael Aye for highlighting the issue
        '''
        
        from fourier_coefficients import *
        
        fourcoeff = Fourier_coeff()
    
        nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                                             # and create an instance of the ncCDF4 class
        nc_attrs, nc_dims, nc_vars = fourcoeff.ncdump(nc_fid,'ta',verb)
        
        # Extract data from NetCDF file
        time = nc_fid.variables['time'][:]  # extract the coordinate
        plev = nc_fid.variables['plev'][:]  # extract the coordinate
        lats = nc_fid.variables['lat'][:]  # extract the coordinate
        #lon = nc_fid.variables['lon'][:]  # extract the coordinate
    
        #print(nc_fid.variables['plev'])
            
        # Writing NetCDF files
        var_nc_fid = Dataset(fileo, 'w', format='NETCDF4')
        var_nc_fid.description = "Fourier coefficients"
    #    var2_nc_fid = Dataset(fileo, name2, format='NETCDF4')
    #    var2_nc_fid.description = "Eastward wind"
    #    var3_nc_fid = Dataset(fileo, name3, format='NETCDF4')
    #    var3_nc_fid.description = "Northward wind"
    #    var4_nc_fid = Dataset(fileo, name4, format='NETCDF4')
    #    var4_nc_fid.description = "Lagrangian tendency of air pressure"
        
        # Using our previous dimension info, we can create the new dimensions.
        
        var_nc_fid.createDimension('time', len(time))
        var_nc_dim = var_nc_fid.createVariable('time', nc_fid.variables['time'].dtype,\
                                           ('time',))
        for ncattr in nc_fid.variables['time'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['time'].getncattr(ncattr))
        var_nc_fid.variables['time'][:] = time
        
        var_nc_fid.createDimension('plev', len(plev))
        var_nc_dim = var_nc_fid.createVariable('plev', nc_fid.variables['plev'].dtype,\
                                           ('plev',))
        for ncattr in nc_fid.variables['plev'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['plev'].getncattr(ncattr))
        var_nc_fid.variables['plev'][:] = plev
        
        var_nc_fid.createDimension('wave', len(wave2))
        var_nc_dim = var_nc_fid.createVariable('wave', nc_fid.variables['plev'].dtype,\
                                           ('wave',))
        var_nc_fid.variables['wave'][:] = wave2
    
    #    var_nc_fid.createDimension('lon', len(lon))
    #    var_nc_dim = var_nc_fid.createVariable('lon', nc_fid.variables['lon'].dtype,\
    #                                       ('lon',))
    #    for ncattr in nc_fid.variables['lon'].ncattrs():
    #        var_nc_dim.setncattr(ncattr, nc_fid.variables['lon'].getncattr(ncattr))
    #    var_nc_fid.variables['lon'][:] = lon
            
        var_nc_fid.createDimension('lat', len(lats))
        var_nc_dim = var_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                           ('lat',))
        for ncattr in nc_fid.variables['lat'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
        var_nc_fid.variables['lat'][:] = lats
    
        nc_fid.close()
        
        var1_nc_var = var_nc_fid.createVariable(name1, 'f8', ('time','plev','wave','lat'))
        fourcoeff.varatts(var1_nc_var,name1)
        #print(np.shape(var1))
        #print(np.shape(var1_nc_var))
        var_nc_fid.variables[name1][:,:,:,:] = var1
        var2_nc_var = var_nc_fid.createVariable(name2, 'f8', ('time','plev','wave','lat'))
        fourcoeff.varatts(var2_nc_var,name2)
        var_nc_fid.variables[name2][:,:,:,:] = var2
        var3_nc_var = var_nc_fid.createVariable(name3, 'f8', ('time','plev','wave','lat'))
        fourcoeff.varatts(var3_nc_var,name3)
        var_nc_fid.variables[name3][:,:,:,:] = var3
        var4_nc_var = var_nc_fid.createVariable(name4, 'f8', ('time','plev','wave','lat'))
        fourcoeff.varatts(var4_nc_var,name4)
        var_nc_fid.variables[name4][:,:,:,:] = var4
        
        var_nc_fid.close()  # close the new file
        
    
    
    def ncdump(self, nc_fid,key,verb):
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
    
    def varatts(self, w_nc_var,varname):
        
        if varname == 'ta':
            w_nc_var.setncatts({'long_name': u"Air temperature",'units': u"K", 'level_desc': 'pressure levels'})
        elif varname == 'ua':
            w_nc_var.setncatts({'long_name': u"Eastward wind",'units': u"m s-1", 'level_desc': 'pressure levels'})
        elif varname == 'va':
            w_nc_var.setncatts({'long_name': u"Northward wind",'units': u"m s-1", 'level_desc':'pressure levels'})
        elif varname == 'wap':
            w_nc_var.setncatts({'long_name': u"Lagrangian tendency of air pressure",'units': u"Pa s-1", 'level_desc': 'pressure levels'})
            
    def print_ncattr(self, key):
            """
            Prints the NetCDF file attributes for a given key
    
            Parameters
            ----------
            key : unicode
                a valid netCDF4.Dataset.variables key
            """
            try:
                print "\t\ttype:", repr(nc_fid.variables[key].dtype)
                for ncattr in nc_fid.variables[key].ncattrs():
                    print '\t\t%s:' % ncattr,\
                          repr(nc_fid.variables[key].getncattr(ncattr))
            except KeyError:
                print "\t\tWARNING: %s does not contain variable attributes" % key