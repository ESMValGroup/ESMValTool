# -*- coding: utf-8 -*-

import os
from shutil import move
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
from scipy import interpolate
import math
import logging
#from cartopy import config
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from PyAstronomy import pyaC
from cdo import *
#from auxiliary import verbosity

###############################################################################

class Plot_script():
    
    from plot_script import *
        
    def latwgt(self,lat,tr):
        
        #Routine to weigh fields with cosines of latitude
        pi    = math.pi
        conv  = 2*pi/360
        
        dlat = np.zeros(len(lat))
        for i in range(len(lat)-1):
            dlat[i]=abs(lat[i+1]-lat[i])
        dlat[len(lat)-1] = dlat[len(lat)-2]
        
        latr  = conv*lat
        dlatr = conv*dlat
        
        tr2=np.zeros((np.shape(tr)[0],np.shape(tr)[1]))
        for j in range(len(lat)):
            tr2[:,j] = tr[:,j]*np.cos(latr[j])*dlatr[j]/2
        return tr2
    
###############################################################################
    
    
    def hemean(self,hem,lat,inp):
        
        #Routine for hemispheric means:
        #hem=1 indicates SH, hem=else indicates NH
        
        plotsmod = Plot_script()
        
        j_end = np.shape(inp)[1]
        zmn = plotsmod.latwgt(lat,inp)
        hmean = []
        if hem == 1:
            if j_end%2 == 0:
                hmean = np.nansum(zmn[:,j_end/2+1:j_end],axis=1)
            else:
                hmean = (np.nansum(zmn[:,(j_end+3)/2:j_end],axis=1)+0.5*zmn[:,(j_end+3)/2-1])
        else:
            if j_end%2 == 0:
                hmean = np.nansum(zmn[:,1:j_end/2],axis=1)
            else:
                hmean = (np.nansum(zmn[:,1:(j_end-1)/2],axis=1)+0.5*zmn[:,(j_end-1)/2+1])
            
        return hmean
    
    ###############################################################################
    
    
    def transport(self,zmean,gmean,lat):
        
        #Routine for the computation of meridional transports (of energy,
        # as well as other quantities)
        plotsmod = Plot_script()
        
        pi  = math.pi
        
        dlat             = np.zeros(len(lat))
        for i in range(len(lat)-1):
            dlat[i] = abs(lat[i+1]-lat[i])
        dlat[len(lat)-1] = dlat[len(lat)-2]
        
        zmn_ub           = np.zeros((np.shape(zmean)[0],np.shape(zmean)[1]))
        for i in range(len(gmean)):
            for j in range(np.shape(zmean)[1]):
                zmn_ub[i,j] = zmean[i,j]-gmean[i]
        
        zmn_ub[np.isnan(zmn_ub)]    = 0
        cumb                        = np.zeros((np.shape(zmean)[0],np.shape(zmean)[1]))
        transp                      = np.zeros((np.shape(zmean)[0],np.shape(zmean)[1]))
        
        for j in range(len(lat)-1):
            cumb[:,j] = -2*np.nansum(plotsmod.latwgt(lat[j:len(lat)],zmn_ub[:,j:len(lat)]),axis=1)
        
        R      = 6.371*10**6
        transp = 2*pi*cumb*R*R
    
    
        return [zmn_ub,transp]
    
    ###############################################################################
    
    def transp_max(self,lat,transp,lim):
        
        deriv  = np.gradient(transp)
        #print(deriv)
        xc, xi = pyaC.zerocross1d(lat, deriv, getIndices=True)
        #print(xc)
        yi     = np.zeros(2)
        xc_cut = np.zeros(2)
        j=0
        for i in range(len(xc)):
                if abs(xc[i])<=lim:
                    xc_cut[j] = xc[i]
                    yi[j]     = interpolate.interp1d(lat,transp,kind='cubic')(xc[i])
                    #print(xc_cut[j],yi[j])
                    j=j+1
                    if j==2:
                        break
                else:
                    pass
        #print(xc_cut,yi)
        
        return [xc_cut,yi]
    
    ###############################################################################
    
    def postsealand(self,path,file,filemask,model,name,transp_name,transpy,time,lat,lon,opt):
        
        plotsmod = Plot_script()
        
        dataset = netcdf_dataset(filemask)
        lsmsk = dataset.variables['sftlf'][:, :]
        lsmsk            = lsmsk/100
        lsm_ext=np.zeros((1,len(lat),len(lon)))
        lsm_ext[0,:,:] = lsmsk
        if opt in {'sea','oc','ocean'}:
            zmnlsm   = np.nanmean(1-lsm_ext,axis=2)
            mytitle  = '{} meridional transports - Oceans'.format(transp_name)
        elif opt in {'land','continents'}:
            zmnlsm   = np.nanmean(lsm_ext,axis=2)
            mytitle  = '{} meridional transports - Land'.format(transp_name)
        else:
            logger.debug('No meaningful option')
        glob       = np.nansum(plotsmod.latwgt(lat,zmnlsm),axis=1)
        dlon       = lon[1]-lon[0]
        NHlsm      = plotsmod.hemean(0,lat,zmnlsm)*360/dlon
        SHlsm      = plotsmod.hemean(1,lat,zmnlsm)*360/dlon
        dataset = netcdf_dataset(file)
        var   = dataset.variables[name][:, :, :]
        dataset = netcdf_dataset(filemask)
        var_r = np.reshape(var,(np.shape(var)[0]/12,12,len(lat),len(lon)))
        var   = np.nanmean(var_r,axis=1)
        #Compute means over land and oceans
        zmean=np.zeros((len(time),len(lat)))
        for i in range(len(lat)):
            zmean[:,i]   = np.nanmean(var[:,i,:],axis=1)*zmnlsm[0,i]
        zmean[np.isnan(zmean)] = 0
        zmean_w  = plotsmod.latwgt(lat,zmean)
        gmean    = np.nansum(zmean_w,axis=1)*glob
        nhgmean  = plotsmod.hemean(0,lat,zmean_w)*NHlsm
        shgmean  = plotsmod.hemean(1,lat,zmean_w)*SHlsm
        timeser = np.column_stack((gmean,shgmean,nhgmean))
        #Compute metrics for the transports
        transp = []
        transp  = plotsmod.transport(zmean,gmean,lat)
        #zmnba   = transp[0]
        transpp = transp[1]
        transp_mn= np.nanmean(transpp,axis=0)
        yr_ext = []
        lat_max = list()
        tr_max  = list()
        lim=45
        for i in range(len(time)):
            yr_ext       = plotsmod.transp_max(lat,transpp[i,:],lim)
            lat_max.append(yr_ext[0])
            tr_max.append(yr_ext[1])
        fig = plt.figure()
        ax  = plt.subplot(111)
        plt.plot(lat, transp_mn)
        plt.title(mytitle)
        plt.xlabel('Latitude [deg]')
        plt.ylabel('[W]')
        plt.tight_layout()
        plt.xlim(-90, 90)
        plt.ylim(transpy)
        plt.grid()
        plt.savefig(path+'/{}_{}_{}_transp.png'.format(model,name,opt))
        #plt.show(fig)
        plt.close(fig)
    
    #Import variables over land and oceans separately
    
    ###############################################################################
    
    def balances(self,workdir,plotpath,filena,name,model_name,lsm):
    
        plotsmod = Plot_script()
        cdo      = Cdo()
        
        nsub       = len(filena)
        sep        = '.nc'
        path       = plotpath
        model      = model_name
    
        timesery = np.zeros([nsub,2])
        if nsub == 3:
            ext_name=['TOA Energy Budget','Atmospheric Energy Budget','Surface Energy Budget']
            timesery[0,:]    = (-2,2)
            rangect      = [-100,100]
            transpty     = (-6E15, 6E15)
            timesery[1,:]   =(-1,1)
            timesery[2,:]   =(-3,3)
        elif nsub == 2:
            ext_name=['Water mass budget','Latent heat budget']
            timesery[0,:]   =(-3E-6,3E-6)
            rangecw       = [-1E-4,1E-4]
            transpwy     = (-2E9, 2E9)
            timesery[1,:]   =(-20,20)
            rangecl       = [-150,150]
            transply     = (-6E15, 6E15)
        else:           
            logger.debug('No option recognized')
            logger.debug('I will now quit')
            quit()
        #Import files
        filena[0]     = filena[0].split(sep, 1)[0]
        filename   = filena[0]+'.nc'
        dataset = netcdf_dataset(filename)
        lats  = dataset.variables['lat'][:]
        lons  = dataset.variables['lon'][:]
        time  = dataset.variables['time'][:]
        nlats = len(lats)
        nlons = len(lons)
        ntime = len(time)
        yr    = len(time)/12
        timey = np.linspace(0,yr-1,num=yr)
        var=np.zeros([nsub,ntime,nlats,nlons])
        for i in np.arange(nsub):
            filena[i]     = filena[i].split(sep, 1)[0]
            filename   = filena[i]+'.nc'
            #print(filename)
            dataset = netcdf_dataset(filename)
            var[i,:,:,:]   = dataset.variables[name[i]][:, :, :]
        #Compute annual mean values
        var_r = np.reshape(var,(nsub,np.shape(var)[1]/12,12,nlats,nlons))
        vary   = np.nanmean(var_r,axis=2)
        
        #Compute the zonal mean
        zmean = np.nanmean(vary, axis=3)  
        #Compute the climatological mean map
        tmean = np.nanmean(vary, axis=1)
        
        #Compute global and hemispheric means as function of years   
        transp_mean=np.zeros([nsub,nlats])
        lat_maxm=np.zeros([nsub,2,len(timey)])
        tr_maxm=np.zeros([nsub,2,len(timey)])
        lim=[55,55,30]
        for i in np.arange(nsub):
            zmean_w = plotsmod.latwgt(lats,zmean[i,:,:])
            gmean   = np.nansum(zmean_w, axis=1)
            shmean  = plotsmod.hemean(0,lats,zmean[i,:,:])
            nhmean  = plotsmod.hemean(1,lats,zmean[i,:,:])
            timeser = np.column_stack((gmean,shmean,nhmean))
            #Compute transports
            transp  = plotsmod.transport(zmean[i,:,:],gmean,lats) 
            transpp = transp[1]
            transp_mean[i,:] = np.nanmean(transpp,axis=0)
            yr_ext = []
            lat_max = list()
            tr_max  = list()
            for t in range(len(timey)):
                yr_ext       = plotsmod.transp_max(lats,transpp[t,:],lim[i])
                lat_max.append(yr_ext[0])
                tr_max.append(yr_ext[1])       
            for t in range(len(timey)):
                lat_maxm[i,:,t]=lat_max[t]
                tr_maxm[i,:,t] =tr_max[t]
        
            tgmean  = np.nanmean(gmean)
            
            fig = plt.figure()
            ax  = plt.subplot(111)
            ax.plot(timey,timeser[:,0],'k',label='Global')
            ax.plot(timey,timeser[:,1],'r',label='SH')
            ax.plot(timey,timeser[:,2],'b',label='NH')
            plt.title('Annual mean {}'.format(ext_name[i]))
            plt.xlabel('Years')
            plt.ylabel('[W/m2]')
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),  shadow=True, ncol=3)
            plt.tight_layout()
            plt.ylim(timesery[i,:])
            plt.grid()
            plt.savefig(path+'/{}_{}_timeser.png'.format(model,name[i]))
            plt.close(fig)
        
        cm='bwr'
        if nsub == 3:
            fig = plt.figure(figsize=(12, 22))
            ax = plt.subplot(311,projection=ccrs.PlateCarree())
            ax.coastlines()
            plt.contourf(lons, lats, tmean[0,:,:],60,transform=ccrs.PlateCarree())
            plt.pcolor(lons, lats, tmean[0,:,:], vmin=rangect[0], vmax=rangect[1],cmap=cm,antialiaseds='True')
            plt.colorbar()
            plt.title('Climatological Mean {}'.format(ext_name[0]))
            #plt.tight_layout()
            plt.grid()
            ax = plt.subplot(312,projection=ccrs.PlateCarree())
            ax.coastlines()
            plt.contourf(lons, lats, tmean[1,:,:],60,transform=ccrs.PlateCarree())
            plt.pcolor(lons, lats, tmean[1,:,:], vmin=rangect[0], vmax=rangect[1],cmap=cm,antialiaseds='True')
            plt.colorbar()
            plt.title('Climatological Mean {}'.format(ext_name[1]))
            #plt.tight_layout()
            plt.grid()
            ax = plt.subplot(313,projection=ccrs.PlateCarree())
            ax.coastlines()
            plt.contourf(lons, lats, tmean[2,:,:],60,transform=ccrs.PlateCarree())
            plt.pcolor(lons, lats, tmean[2,:,:], vmin=rangect[0], vmax=rangect[1],cmap=cm,antialiaseds='True')
            plt.colorbar()
            plt.title('Climatological Mean {}'.format(ext_name[2]))
            #plt.tight_layout()
            plt.grid()
            plt.savefig(path+'/{}_energy_climap.png'.format(model))
            #plt.show()
            plt.close(fig)
        elif nsub == 2:
            fig = plt.figure()
            ax = plt.subplot(111,projection=ccrs.PlateCarree())
            ax.coastlines()
            plt.contourf(lons, lats, tmean[0,:,:],60,transform=ccrs.PlateCarree())
            plt.pcolor(lons, lats, tmean[0,:,:], vmin=rangecw[0], vmax=rangecw[1],cmap=cm,antialiaseds='True')
            plt.colorbar()
            #plt.colorbar(ax,fraction=0.046, pad=0.04)
            plt.title('Climatological Mean {}'.format(ext_name[0]))
            #plt.tight_layout()
            plt.grid()
            plt.savefig(path+'/{}_{}_climap.png'.format(model,name[0]))
            plt.close(fig)
            fig = plt.figure()
            #fig.set_size_inches(12, 22)
            #ax = plt.subplot(111)
            ax = plt.subplot(111,projection=ccrs.PlateCarree())
            ax.coastlines()
            plt.contourf(lons, lats, tmean[1,:,:],60,transform=ccrs.PlateCarree())
            plt.pcolor(lons, lats, tmean[1,:,:], vmin=rangecl[0], vmax=rangecl[1],cmap=cm,antialiaseds='True')
            #plt.colorbar(ax,fraction=0.046, pad=0.04)
            plt.colorbar()
            plt.title('Climatological Mean {}'.format(ext_name[1]))
            #plt.tight_layout()
            plt.grid()
            plt.savefig(path+'/{}_{}_climap.png'.format(model,name[1]))
            plt.close(fig)
        if nsub == 3:
            fig = plt.figure()
            ax  = plt.subplot(111)
            for i in np.arange(nsub):
                filename=filena[i]+'.nc'
                if name[i]=='toab':
                    nameout='total'
                elif name[i]=='atmb':
                    nameout='atmos'
                elif name[i]=='surb':
                    nameout='ocean'
                nc_f=workdir+'/{}_transp_mean_{}.nc'.format(nameout,model_name)
                plotsmod.removeif(nc_f)
                plotsmod.pr_output(transp_mean[i,:], filename, nc_f, nameout, verb=True)
                name_model = '{}_{}'.format(nameout, model_name)
                lat_model  = 'lat_{}'.format(model_name)
                cdo.chname('{},{}'.format(nameout,name_model), input=nc_f, 
                           output='aux.nc')
                move('aux.nc', nc_f)
                cdo.chname('lat,{}'.format(lat_model), input=nc_f, 
                           output='aux.nc')
                move('aux.nc', nc_f)
                plt.plot(lats, transp_mean[i,:])
            plt.title('Meridional heat transports')
            plt.xlabel('Latitude [deg]',fontsize=10)
            plt.ylabel('[W]',fontsize=10)
            plt.tight_layout()
            plt.ylim(transpty)
            plt.xlim(-90, 90)
            plt.grid()
            plt.savefig(path+'/{}_transp.png'.format(model))
            #plt.show(fig)
            plt.close(fig)
        elif nsub == 2:
            fig = plt.figure()
            ax  = plt.subplot(111)
            nc_f=workdir+'/{}_transp_mean_{}.nc'.format('wmass',model)
            plotsmod.removeif(nc_f)
            plotsmod.pr_output(transp_mean[0,:],filename, nc_f, 'wmass', verb=True)
            plt.plot(lats, transp_mean[0,:])
            plt.title('Water mass transports',fontsize=10)
            plt.xlabel('Latitude [deg]',fontsize=10)
            plt.ylabel('[W]')
            plt.tight_layout()
            plt.ylim(transpwy)
            plt.xlim(-90, 90)
            plt.grid()
            plt.savefig(path+'/{}_wmass_transp.png'.format(model))
            plt.close(fig)
            fig = plt.figure()
            ax  = plt.subplot(111)
            nc_f=workdir+'/{}_transp_mean_{}.nc'.format('latent',model)
            plotsmod.removeif(nc_f)
            plotsmod.pr_output(transp_mean[1,:],filename, nc_f, 'latent', verb=True)
            plt.plot(lats, transp_mean[1,:])
            plt.title('Latent heat transports',fontsize=10)
            plt.xlabel('Latitude [deg]',fontsize=10)
            plt.ylabel('[W]')
            plt.tight_layout()
            plt.ylim(transply)
            plt.xlim(-90, 90)
            plt.grid()
            plt.savefig(path+'/{}_latent_transp.png'.format(model))
            #plt.show(fig)
            plt.close(fig)
        
        colors = (0,0,0)
        if nsub == 3:           
            fig = plt.figure()
            fig.set_size_inches(12, 22)
            fig.set_size_inches(12, 22)
            ax  = plt.subplot(321)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[0,0,:], tr_maxm[0,0,:], c=colors, alpha=1)
            plt.title('Scatter plot of total peak magnitude against peak position - SH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_sh)
            #plt.ylim(scattery_sh)
            #plt.tight_layout()
            plt.grid()
            
            ax  = plt.subplot(322)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[0,1,:], tr_maxm[0,1,:], c=colors, alpha=1)
            plt.title('Scatter plot of total peak magnitude against peak position - NH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_nh)
            #plt.ylim(scattery_nh)
            #plt.tight_layout()
            plt.grid()
            
            ax  = plt.subplot(323)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[1,0,:], tr_maxm[1,0,:], c=colors, alpha=1)
            plt.title('Scatter plot of atmospheric peak magnitude against peak position - SH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_sh)
            #plt.ylim(scattery_sh)
            #plt.tight_layout()
            plt.grid()
            
            ax  = plt.subplot(324)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[1,1,:], tr_maxm[1,1,:], c=colors, alpha=1)
            plt.title('Scatter plot of atmospheric peak magnitude against peak position - NH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_nh)
            #plt.ylim(scattery_nh)
            #plt.tight_layout()
            plt.grid()
            
            ax  = plt.subplot(325)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[2,0,:], tr_maxm[2,0,:], c=colors, alpha=1)
            plt.title('Scatter plot of oceanic peak magnitude against peak position - SH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_sh)
            #plt.ylim(scattery_sh)
            #plt.tight_layout()
            plt.grid()
           
            ax  = plt.subplot(326)
            ax.set_figsize=(50,50)
            plt.scatter(lat_maxm[2,1,:], tr_maxm[2,1,:], c=colors, alpha=1)
            plt.title('Scatter plot of oceanic peak magnitude against peak position - NH',fontsize=10,y=1.04)
            plt.xlabel('Peak position [degrees of latitude]',fontsize=10)
            plt.ylabel('Peak magnitude [W]',fontsize=10)
            #plt.xlim(scatterx_nh)
            #plt.ylim(scattery_nh)
            #plt.tight_layout()
            plt.grid()
            #plt.show(fig)
            plt.savefig(path+'/{}_scatpeak.png'.format(model))
            
            plt.close(fig)
        
        #Import land-sea mask and perform operations on land and oceans
    #    if lsm in {'y','yes'}:
    #        #filemask = "/Users/Valerio2/ESMValTool-private/modeldata/sftlf_fx_MPI-ESM-LR_historical_r0i0p0.nc"
    #        opt      = 'oc'
    #        file     = "{}_ocean.nc".format(filena)
    #        postsealand(path,file,filemask,model,name,transp_name,transpy,timey,lats,lons,opt)
    #        opt      = 'land'
    #        file     = "{}_land.nc".format(filena)
    #        postsealand(path,file,filemask,model,name,transp_name,transpy,timey,lats,lons,opt)
    #    else:
    #        pass
        
        #map2d(np.nanmean(varla,axis=0),lons,lats)
    
    def entropy(self,plotpath,filename,name,ext_name,model_name):
    
        path       = plotpath
        model      = model_name
        
        if ext_name == 'Vertical entropy production':
            rangec      = [-0.01,0.1]
            cm         = 'YlOrBr'
        elif ext_name == 'Horizontal entropy production':
            rangec       = [-0.5,0.5]
            cm           = 'bwr'
        elif ext_name == 'Sensible Heat entropy production':
            rangec       = [-0.01,0.01]
            cm           = 'YlOrBr'
        elif ext_name == 'Evaporation entropy production':
            rangec       = [0,1]
            cm           = 'YlOrBr'
        elif ext_name == 'Rainfall precipitation entropy production':
            rangec       = [0,1]
            cm           = 'YlOrBr'
        elif ext_name == 'Snowfall precipitation entropy production':
            rangec       = [0,0.25]
            cm           = 'YlOrBr'
        elif ext_name == 'Phase changes ice -> rain entropy production':
            rangec       = [0,0.05]
            cm           = 'YlOrBr'
        elif ext_name == 'Phase changes vapor -> snow entropy production':
            rangec       = [0,0.001]
            cm           = 'YlOrBr'
        elif ext_name == 'Potential energy entropy production':
            rangec       = [0,0.1]
            cm           = 'YlOrBr'
        else:           
            logger.debug('No name recognized')
            logger.debug('I will now quit')
            quit()
            #Import files
            #filename="/Users/Valerio2/ESMValTool-private/TRR181_valerio/MPI-ESM-LR_toab_ymm.nc"
        dataset = netcdf_dataset(filename)
        var   = dataset.variables[name][:, :, :]
        lats  = dataset.variables['lat'][:]
        lons  = dataset.variables['lon'][:]
        time  = dataset.variables['time'][:]
        
        #Compute the climatological mean map
        tmean = np.nanmean(var, axis=0)
        fig = plt.figure()
        #fig.set_size_inches(12, 22)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        plt.contourf(lons, lats, tmean,60,transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean, vmin=rangec[0], vmax=rangec[1],cmap=cm,antialiaseds='True')    
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name))
        plt.tight_layout()
        plt.grid()
        plt.savefig(path+'/{}_{}_climap.png'.format(model,name))
        #plt.show()
        plt.close(fig)
    
    def plot_ellipse(self,semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e3,
                     ax=None,plot_kwargs=None,fill=False,fill_kwargs=None,
                     data_out=False,cov=None,mass_level=0.68):
        '''
            An easy to use function for plotting ellipses in Python 2.7!
        '''
    
        plotsmod = Plot_script()
        
        # Get Ellipse Properties from cov matrix
        if cov is not None:
            eig_vec,eig_val,u = np.linalg.svd(cov)
            # Make sure 0th eigenvector has positive x-coordinate
            if eig_vec[0][0] < 0:
                eig_vec[0] *= -1
            semimaj = np.sqrt(eig_val[0])
            semimin = np.sqrt(eig_val[1])
            if mass_level is None:
                multiplier = np.sqrt(2.279)
            else:
                distances = np.linspace(0,20,20001)
                chi2_cdf = plotsmod.chi2.cdf(distances,df=2)
                multiplier = np.sqrt(distances[np.where(np.abs(chi2_cdf-mass_level)==np.abs(chi2_cdf-mass_level).min())[0][0]])
            semimaj *= multiplier
            semimin *= multiplier
            phi = np.arccos(np.dot(eig_vec[0],np.array([1,0])))
            if eig_vec[0][1] < 0 and phi > 0:
                phi *= -1
    
        # Generate data for ellipse structure
        theta = np.linspace(0,2*np.pi,theta_num)
        r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        data = np.array([x,y])
        S = np.array([[semimaj,0],[0,semimin]])
        R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
        T = np.dot(R,S)
        data = np.dot(T,data)
        data[0] += x_cent
        data[1] += y_cent
    
        # Output data?
        if data_out == True:
            return data
    
        # Plot!
        return_fig = False
        if ax is None:
            return_fig = True
            fig,ax = plt.subplots()
    
        if plot_kwargs is None:
            ax.plot(data[0],data[1],color='b',linestyle='-')
        else:
            ax.plot(data[0],data[1],**plot_kwargs)
    
        plot_kwargs = {'color':'black'}
        if fill == True:
            ax.fill(data[0],data[1],**fill_kwargs)
    
        if return_fig == True:
            return fig
        
    def pr_output(self,varout, filep, nc_f, nameout, verb=True):
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
    
        plotsmod = Plot_script()
        
        nc_fid = netcdf_dataset(filep, 'r')  # Dataset is the class behavior to open the file
                                             # and create an instance of the ncCDF4 class
        nc_attrs, nc_dims, nc_vars = plotsmod.ncdump(nc_fid,nameout,verb)
        # Extract data from NetCDF file
        lats = nc_fid.variables['lat'][:]  # extract the coordinate
        
        # Writing NetCDF files
        w_nc_fid = netcdf_dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "Total, atmospheric and oceanic annual mean meridional heat transports"
        
        # Using our previous dimension info, we can create the new dimensions.
        
        w_nc_fid.createDimension('lat', len(lats))
        w_nc_dim = w_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                           ('lat',))
        for ncattr in nc_fid.variables['lat'].ncattrs():
            w_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
        w_nc_fid.variables['lat'][:] = lats
        w_nc_var = w_nc_fid.createVariable(nameout, 'f8', ('lat'))
        plotsmod.varatts(w_nc_var,nameout)
        w_nc_fid.variables[nameout][:] = varout
        w_nc_fid.close()  # close the new file
       
        nc_fid.close()
        
    
    
    def ncdump(self,nc_fid,key,verb):
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
    
    def varatts(self,w_nc_var,varname):
        
        if varname == 'total':
            w_nc_var.setncatts({'long_name': u"Total meridional heat transport",'units': u"W", 'level_desc': 'TOA'})
        elif varname == 'atmos':
            w_nc_var.setncatts({'long_name': u"Atmospheric meridional heat transport",'units': u"W", 'level_desc': 'Vertically integrated'})
        elif varname == 'ocean':
            w_nc_var.setncatts({'long_name': u"Oceanic meridional heat transport",'units': u"W", 'level_desc':'sfc'})
        elif varname == 'wmass':
            w_nc_var.setncatts({'long_name': u"Meridional water mass transport",'units': u"W", 'level_desc': 'sfc'})
        elif varname == 'latent':
             w_nc_var.setncatts({'long_name': u"Meridional latent heat transport",'units': u"W", 'level_desc': 'sfc'})
             
    def removeif(self,filename):
        try:
            os.remove(filename)
        except OSError:
            pass