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
g0 = 9.81               # Gravity acceleration
gam = 0.0065            # Standard atmosphere lapse rate
gascon  = 287.0         # Gas constant
rho0 = 1.2              # Mean air density
p0 = 10000              # Reference tropospheric pressure

class Fourier_coeff():
    
    from fourier_coefficients import *
    
    def fourier_coeff(self, tadiagfile, outfile, tainput, tasinput):
        
        fourcoeff = Fourier_coeff()
        
        fileo    = outfile
        fileta   = tadiagfile
        fileps   = '/work/um0005/u234097/ESMV/esmvaltool_output/ps.nc'
        
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
        dataset = Dataset(tasinput)
        tas=dataset.variables['tas'][:,:,:]
        tas=tas[:,::-1,:]
        
        ta1_fx = np.array(ta)
        deltat = np.zeros([ntime,nlev,nlat,nlon])
        ps = np.full([ntime,nlat,nlon],p0)
        for i in np.arange(nlev-1,0,-1):
            #h1 = np.ma.masked_where(ta1_fx[:,i-1,:,:]==0, ta1_fx[:,i-1,:,:])
            h1 = np.ma.masked_where(ta1_fx[:,i,:,:]!=0, ta1_fx[:,i,:,:])
            #mask=h1.mask
            if np.any(h1.mask>0):
                #deltat[:,i-1,:,:] = np.where(ta1_fx[:,i-1,:,:]!=0,deltat[:,i-1,:,:], \
                #      (ta1_fx[:,i,:,:] - tas))
                deltat[:,i-1,:,:] = np.where(ta1_fx[:,i-1,:,:]!=0,deltat[:,i-1,:,:], \
                      (ta1_fx[:,i,:,:] - tas))
                deltat[:,i-1,:,:] = (1*np.array(h1.mask))*np.array(deltat[:,i-1,:,:])
                dp = -((p0*g0/(gam*gascon))*deltat[:,i-1,:,:]/tas)
                ps = np.where(ta1_fx[:,i-1,:,:]!=0,ps,lev[i-1] + dp)
                for k in np.arange(0,nlev-i-1,1):
                    h3 = np.ma.masked_where(ta1_fx[:,i+k,:,:]!=0, ta1_fx[:,i+k,:,:])
                    if np.any(h3.mask>0):
                        deltat[:,i-1,:,:] = np.where(ta1_fx[:,i+k,:,:]!=0, \
                              deltat[:,i-1,:,:], (ta1_fx[:,i+k+1,:,:] - tas))
                        dp = -((p0*g0/(gam*gascon))*deltat[:,i-1,:,:]/tas)
                        ps = np.where(ta1_fx[:,i+k,:,:]!=0,ps,lev[i+k] + dp)
                    else:
                        pass
                #dp = -((p0*g0/(gam*gascon))*deltat[:,i-1,:,:]/tas)
                #ps = np.where(deltat[:,i-1,:,:]==0,ps,lev[i-1] + dp)
            else:
                pass
        ta2_fx = np.array(ta)
        mask = np.zeros([nlev,ntime,nlat,nlon])
        dat = np.zeros([nlev,ntime,nlat,nlon])
        tafr_bar = np.zeros([nlev,ntime,nlat,nlon])
        deltap = np.zeros([ntime,nlev,nlat,nlon])
        #ta[:,0,:,:] = tas 
        for i in np.arange(nlev):
            deltap[:,i,:,:] = ps - lev[i]
            h2 = np.ma.masked_where(ta2_fx[:,i,:,:]==0, ta2_fx[:,i,:,:])
            mask[i,:,:,:] = np.array(h2.mask)       
            tafr_bar[i,:,:,:] = 1*np.array(mask[i,:,:,:])*(tas - \
                    gam*gascon/(g0*ps)*deltap[:,i,:,:]*tas)
            dat[i,:,:,:] = ta2_fx[:,i,:,:] * (1-1*np.array(mask[i,:,:,:]))
            ta[:,i,:,:] = dat[i,:,:,:]+tafr_bar[i,:,:,:]
        fourcoeff.pr_output_diag(ta, tainput, fileta, 'ta', verb=True)
        fourcoeff.pr_output_ps(ps, tainput, fileps, 'ps', verb=True)
        
        #print(np.shape(ta))
        tafft_p = np.fft.fft(ta,axis=3)[:,:,:,:trunc/2]/(nlon)
        #tafft_p = np.transpose(tafft_p, (0, 1, 3, 2))
        uafft_p = np.fft.fft(ua,axis=3)[:,:,:,:trunc/2]/(nlon)
        #uafft_p = np.transpose(uafft_p, (0, 1, 3, 2))
        vafft_p = np.fft.fft(va,axis=3)[:,:,:,:trunc/2]/(nlon)
        #vafft_p = np.transpose(vafft_p, (0, 1, 3, 2))
        wapfft_p = np.fft.fft(wap,axis=3)[:,:,:,:trunc/2]/(nlon)
        #wapfft_p = np.transpose(wapfft_p, (0, 1, 3, 2))
        
        tafft = np.zeros([ntime,nlev,nlat,trunc])
        uafft = np.zeros([ntime,nlev,nlat,trunc])
        vafft = np.zeros([ntime,nlev,nlat,trunc])
        wapfft = np.zeros([ntime,nlev,nlat,trunc])
        tafft[:,:,:,0::2]=np.real(tafft_p)
        tafft[:,:,:,1::2]=np.imag(tafft_p)
        uafft[:,:,:,0::2]=np.real(uafft_p)
        uafft[:,:,:,1::2]=np.imag(uafft_p)
        vafft[:,:,:,0::2]=np.real(vafft_p)
        vafft[:,:,:,1::2]=np.imag(vafft_p)
        wapfft[:,:,:,0::2]=np.real(wapfft_p)
        wapfft[:,:,:,1::2]=np.imag(wapfft_p)
            
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
                    
        # Writing NetCDF files
        var_nc_fid = Dataset(fileo, 'w', format='NETCDF4')
        var_nc_fid.description = "Fourier coefficients"
        
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
            
        var_nc_fid.createDimension('lat', len(lats))
        var_nc_dim = var_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                           ('lat',))
        for ncattr in nc_fid.variables['lat'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
        var_nc_fid.variables['lat'][:] = lats
    
        nc_fid.close()
        
        var1_nc_var = var_nc_fid.createVariable(name1, 'f8', ('time','plev','lat','wave'))
        fourcoeff.varatts(var1_nc_var,name1)
        var_nc_fid.variables[name1][:,:,:,:] = var1
        var2_nc_var = var_nc_fid.createVariable(name2, 'f8', ('time','plev','lat','wave'))
        fourcoeff.varatts(var2_nc_var,name2)
        var_nc_fid.variables[name2][:,:,:,:] = var2
        var3_nc_var = var_nc_fid.createVariable(name3, 'f8', ('time','plev','lat','wave'))
        fourcoeff.varatts(var3_nc_var,name3)
        var_nc_fid.variables[name3][:,:,:,:] = var3
        var4_nc_var = var_nc_fid.createVariable(name4, 'f8', ('time','plev','lat','wave'))
        fourcoeff.varatts(var4_nc_var,name4)
        var_nc_fid.variables[name4][:,:,:,:] = var4
        
        var_nc_fid.close()  # close the new file
        
    
    def pr_output_diag(self, var1, nc_f, fileo, name1, verb=True):
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
        lons = nc_fid.variables['lon'][:]  # extract the coordinate
            
        # Writing NetCDF files
        var_nc_fid = Dataset(fileo, 'w', format='NETCDF4')
        var_nc_fid.description = "Fourier coefficients"
        
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
        
        var_nc_fid.createDimension('lon', len(lons))
        var_nc_dim = var_nc_fid.createVariable('lon', nc_fid.variables['plev'].dtype,\
                                           ('lon',))
        var_nc_fid.variables['lon'][:] = lons
            
        var_nc_fid.createDimension('lat', len(lats))
        var_nc_dim = var_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                           ('lat',))
        for ncattr in nc_fid.variables['lat'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
        var_nc_fid.variables['lat'][:] = lats
    
        nc_fid.close()
        
        var1_nc_var = var_nc_fid.createVariable(name1, 'f8', ('time','plev','lat','lon'))
        fourcoeff.varatts(var1_nc_var,name1)
        var_nc_fid.variables[name1][:,:,:,:] = var1
        
        var_nc_fid.close()  # close the new file
        
    
    def pr_output_ps(self, var1, nc_f, fileo, name1, verb=True):
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
        lats = nc_fid.variables['lat'][:]  # extract the coordinate
        lons = nc_fid.variables['lon'][:]  # extract the coordinate
            
        # Writing NetCDF files
        var_nc_fid = Dataset(fileo, 'w', format='NETCDF4')
        var_nc_fid.description = "Fourier coefficients"
        
        # Using our previous dimension info, we can create the new dimensions.        
        var_nc_fid.createDimension('time', len(time))
        var_nc_dim = var_nc_fid.createVariable('time', nc_fid.variables['time'].dtype,\
                                           ('time',))
        for ncattr in nc_fid.variables['time'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['time'].getncattr(ncattr))
        var_nc_fid.variables['time'][:] = time
        
        var_nc_fid.createDimension('lon', len(lons))
        var_nc_dim = var_nc_fid.createVariable('lon', nc_fid.variables['plev'].dtype,\
                                           ('lon',))
        var_nc_fid.variables['lon'][:] = lons
            
        var_nc_fid.createDimension('lat', len(lats))
        var_nc_dim = var_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                           ('lat',))
        for ncattr in nc_fid.variables['lat'].ncattrs():
            var_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
        var_nc_fid.variables['lat'][:] = lats
    
        nc_fid.close()
        
        var1_nc_var = var_nc_fid.createVariable(name1, 'f8', ('time','lat','lon'))
        fourcoeff.varatts(var1_nc_var,name1)
        var_nc_fid.variables[name1][:,:,:] = var1
        
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
