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
#import seaborn as sns
import palettable
#plt.style.context('seaborn-talk')
#sns.set_context("paper")
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from esmvaltool.diag_scripts.ocean3d.utils import genfilename, dens_back
from esmvaltool.diag_scripts.ocean3d.interpolation import interpolate_pyresample, interpolate_esmf, closest_depth
from esmvaltool.diag_scripts.ocean3d.getdata import  transect_points

def hofm_plot(model_filenames, cmor_var,
              max_level, region, diagworkdir, diagplotdir,
              levels, ncols=3, cmap=cm.Spectral_r, observations='PHC'):
    # I am not sure why it happens, but 
    # the del below delete instance from the original dict
    # copy works.
    model_filenames = model_filenames.copy()
    if observations:
        del model_filenames[observations]
    ncols = 3
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(8*ncols,2*nrows*ncols))

    for mmodel in model_filenames:
        logger.info("Plot  {} data for {}, region {}".format(cmor_var,
                                                      mmodel,
                                                      region))

        ifilename = genfilename(diagworkdir, cmor_var, 
                                mmodel, region, 'hofm', '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var,
                                       mmodel, region, 'levels', '.npy')
        ifilename_time = genfilename(diagworkdir, cmor_var,
                                     mmodel, region,'time', '.npy')
        print(ifilename)

        hofdata = np.load(ifilename)
        lev = np.load(ifilename_levels)
        time = np.load(ifilename_time)
        if cmor_var == 'thetao':
            hofdata = hofdata-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        lev_limit = lev[lev <= max_level].shape[0]+1

        series_lenght = time.shape[0]

        months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        plt.subplot(nrows, ncols, nplot)

        plt.contourf(months, depth, hofdata, cmap=cmap,
                     levels=levels,
                     extend='both')
        plt.yticks(size=15)
        plt.ylabel('m', size=15, rotation='horizontal')
        plt.ylim(max_level, 0)

        cb = plt.colorbar(pad=0.01)
        cb.set_label(cb_label, rotation='vertical', size=15)
        # cb.set_ticks(size=15)
        cb.ax.tick_params(labelsize=15)
        ygap = int((np.round(series_lenght/12.)/5)*12)
        year_indexes = range(series_lenght)[::ygap]
        year_value = []
        for index in year_indexes:
            year_value.append(time[index].year)
        plt.xticks(year_indexes, year_value, size=15)

        plt.title(mmodel, size=20)
        nplot=nplot+1

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var,
                             region= region, data_type='hofm')
    plt.savefig(pltoutname, dpi=100)


def hofm_plot2(model_filenames, cmor_var,
              max_level, region, diagworkdir, diagplotdir,
              levels, ncols=3, cmap=cm.Spectral_r, observations='PHC'):
    model_filenames = model_filenames.copy()
    if observations:
        del model_filenames[observations]
    ncols = 3
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(8*ncols,2*nrows*ncols))

    for mmodel in model_filenames:
        ifilename = genfilename(diagworkdir, cmor_var, 
                                mmodel, region, 'hofm', '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var,
                                       mmodel, region, 'levels', '.npy')
        ifilename_time = genfilename(diagworkdir, cmor_var,
                                     mmodel, region,'time', '.npy')
        logger.info("Plot  {} data for {}, region {}".format(cmor_var,
                                                      mmodel,
                                                      region))
        hofdata = np.load(ifilename)
        
        lev = np.load(ifilename_levels)
        time = np.load(ifilename_time)
        if cmor_var == 'thetao':
            hofdata = hofdata-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var =='so':
            cb_label = 'psu'

        hofdata_mean = hofdata[:, :].mean(axis=1)
        hofdata_anom = (hofdata[:, :].transpose()-hofdata_mean).transpose()

        lev_limit = lev[lev <= max_level].shape[0]+1

        series_lenght = time.shape[0]

        months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        plt.subplot(nrows, ncols, nplot)

        plt.contourf(months, depth, hofdata_anom, cmap=cmap,
                     levels=levels,
                     extend='both')
        plt.yticks(size=15)
        plt.ylabel('m', size=15, rotation='horizontal')
        plt.ylim(max_level, 0)

        cb = plt.colorbar(pad =0.01)
        cb.set_label(cb_label, rotation='horizontal', size=15)
        # cb.set_ticks(size=15)
        cb.ax.tick_params(labelsize=15)
        ygap = int((np.round(series_lenght/12.)/5)*12)
        year_indexes = range(series_lenght)[::ygap]
        year_value = []
        for index in year_indexes:
            year_value.append(time[index].year)
        plt.xticks(year_indexes, year_value, size=15)

        plt.title(mmodel, size=20)
        nplot=nplot+1

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var,
                             region= region, data_type='hofm2')
    # pltoutname = os.path.join(diagplotdir,
    #                              'arctic_ocean_{}_{}_hofm2'.format(cmor_var,
    #                                                               region))
    plt.savefig(pltoutname, dpi=100)

def tsplot_plot(model_filenames, max_level, region, diagworkdir, diagplotdir,
                ncols=3, cmap=cm.Set1, observations = 'PHC'):
    #ncols = 3
    # phc_dict = OrderedDict([('PHC3.0','./')])
    # phc_dict.update(model_filenames)
    # model_filenames = phc_dict

    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(8*ncols,2*nrows*ncols))

    for mmodel in model_filenames:
        logger.info("Plot  tsplot data for {}, region {}".format(mmodel,
                                                      region))

        ifilename_t = genfilename(diagworkdir, 'thetao',
                                mmodel, region, 'tsplot',  '.npy')
        ifilename_s = genfilename(diagworkdir, 'so',
                                mmodel, region, 'tsplot',  '.npy')
        ifilename_depth = genfilename(diagworkdir, 'depth',
                                mmodel, region, 'tsplot',  '.npy')

        temp = np.load(ifilename_t)
        salt = np.load(ifilename_s)
        depth = np.load(ifilename_depth)

        #lev_limit = lev[lev <= max_level].shape[0]+1

        #series_lenght = time.shape[0]

        #months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        plt.subplot(nrows, ncols, nplot)

        si2, ti2, dens = dens_back(33, 36.,-2, 6)
        if mmodel == 'PHC3.0':
            temp = temp
        else:
            temp = temp-273.15

        cs = plt.contour(si2, ti2, dens, colors='k', levels = np.linspace(dens.min(),dens.max(),15), alpha=0.3)
        plt.scatter(salt[::], temp[::], c=depth, s=3.0,  cmap=cmap, edgecolors='none', vmax=max_level)
        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.1f')
        plt.xlim(33, 36.)
        plt.ylim(-2.1, 6)
        plt.xlabel('Salinity', size=20)
        plt.ylabel('Temperature, $^{\circ}$C', size=20)
        plt.xticks(size=15)
        plt.yticks(size=15)
        cb = plt.colorbar(pad=0.03)
        cb.ax.get_yaxis().labelpad = 15
        cb.set_label('depth, m', rotation=270, size=20)
        cb.ax.tick_params(labelsize=15)
#        plt.tight_layout()

        plt.title(mmodel, size=20)
        nplot=nplot+1

    # ifilename_t = genfilename(diagworkdir, 'thetao',
    #                         'PHC3.0', region, 'tsplot',  '.npy')
    # ifilename_s = genfilename(diagworkdir, 'so',
    #                         'PHC3.0', region, 'tsplot',  '.npy')
    # ifilename_depth = genfilename(diagworkdir, 'depth',
    #                         'PHC3.0', region, 'tsplot',  '.npy')
    

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, 'tsplot',
                             region= region, data_type='tsplot')
    plt.savefig(pltoutname, dpi=100)

def plot_profile(model_filenames, cmor_var,
                 max_level, region, diagworkdir, diagplotdir,
                 cmap=cm.Set2, dpi=100, observations='PHC'):
    print(model_filenames)
    level_clim = Dataset(model_filenames[observations]).variables['lev'][:]
    # plt.style.use('seaborn-paper')
    plt.figure(figsize=(5, 6))
    ax = plt.subplot(111)

    color = iter(cmap(np.linspace(0, 1, len(model_filenames))))
    lev_limit_clim = level_clim[level_clim <= max_level].shape[0]+1
    mean_profile = np.zeros((level_clim[:lev_limit_clim].shape[0],
                             len(model_filenames)-1))
    mean_profile_counter = 0
    for i, mmodel in enumerate(model_filenames):
        logger.info("Plot profile {} data for {}, region {}".format(cmor_var,
                                                             mmodel,
                                                             region))
        ifilename = genfilename(diagworkdir, cmor_var, 
                                mmodel, region, 'hofm', '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var,
                                       mmodel, region, 'levels', '.npy')

        hofdata = np.load(ifilename)
        lev = np.load(ifilename_levels)
        

        if cmor_var == 'thetao':
            hofdata = hofdata-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        #datafile = Dataset(model_filenames[mmodel])

        #lev = datafile.variables['lev'][:]
        lev_limit = lev[lev <= max_level].shape[0]+1
        profile = (hofdata)[:, :].mean(axis=1)
        
        if mmodel != observations:
            c = next(color)
        else:
            c = 'k'
        plt.plot(profile,
                 lev[0:lev_limit],
                 label=mmodel,
                 c=c)

        profile_interpolated = np.interp(level_clim[:lev_limit_clim],
                                         lev[0:lev_limit],
                                         profile)
        if mmodel != observations:

            print('include {} in to the mean'.format(mmodel))
            mean_profile[:, mean_profile_counter] = profile_interpolated
            mean_profile_counter += 1

    mean_profile_mean = np.nanmean(mean_profile, axis=1)

    # plt.plot(clim_prof[:lev_limit_clim],
    #          level_clim[:lev_limit_clim],
    #          label='PHC',
    #          linestyle='--',
    #          color='r',
    #          lw=3)
    # ofilename = genfilename(diagworkdir, cmor_var, 
    #                             'PHC-CLIMATOLOGY', region, 'profile', '.npy')
    # ofilename_levels = genfilename(diagworkdir, cmor_var, 
    #                             'PHC-CLIMATOLOGY', region, 'profile-levels', '.npy')
    # ofilename = os.path.join(diagworkdir,
    #                          'arctic_ocean_{}_{}_{}_profile.npy'.format(cmor_var,
    #                                                                  'PHC-CLIMATOLOGY',
    #                                                                  region))
    # ofilename_levels = os.path.join(diagworkdir,
    #                              'arctic_ocean_{}_{}_{}_profile-levels.npy'.format(cmor_var,
    #                                                               'PHC-CLIMATOLOGY',
    #                                                               region))
    # np.save(ofilename, clim_prof[:lev_limit_clim])
    # np.save(ofilename_levels, level_clim[:lev_limit_clim])

    plt.plot(mean_profile_mean,
             level_clim[:lev_limit_clim],
             label='MODEL-MEAN',
             linestyle='--',
             color='k',
             lw=3)
    # ofilename = genfilename(diagworkdir, cmor_var, 
    #                             'MODEL-MEAN', region, 'profile', '.npy')
    # ofilename_levels = genfilename(diagworkdir, cmor_var, 
    #                             'MODEL-MEAN', region, 'profile-levels', '.npy')
    # ofilename = os.path.join(diagworkdir,
    #                          'arctic_ocean_{}_{}_{}_profile.npy'.format(cmor_var,
    #                                                                  'MODEL-MEAN',
    #                                                                  region))
    # ofilename_levels = os.path.join(diagworkdir,
    #                              'arctic_ocean_{}_{}_{}_profile_levels.npy'.format(cmor_var,
    #                                                                         'MODEL-MEAN',
    #                                                                          region))
    # np.save(ofilename, mean_profile_mean)
    # np.save(ofilename_levels, level_clim[:lev_limit_clim])


    plt.xticks(size=12)
    plt.yticks(size=12)

    plt.xlabel(cb_label, size=12, rotation='horizontal')
    plt.ylabel('m', size=12, rotation='horizontal')

    plt.ylim(0, max_level)

    # plt.legend(loc=0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    
    plt.gca().invert_yaxis()
    # plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var,
                             'MULTIMODEL', region, 'profile')
    # pltoutname = os.path.join(diagplotdir,
    #                           'arctic_ocean_{}_{}_profile'.format(cmor_var,
    #                                                               region))
    plt.savefig(pltoutname, dpi=dpi, bbox_inches='tight')

def plot2d_original_grid(model_filenames, cmor_var, 
                         depth, region, diagworkdir, diagplotdir, dpi=100):

    ncols = 4
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    # nplot = 1
    figsize=(5*ncols,1.5*nrows*ncols)
    fig, ax = plt.subplots(nrows,ncols, 
                           figsize=figsize,
#                           subplot_kw=dict(projection=proj),
                           constrained_layout=True)
    ax = ax.flatten()
    
    for indx, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_original_grid {} for {}".format(cmor_var, mmodel))
        # ifilename = os.path.join(diagworkdir,
        #                     'arctic_ocean_{}_{}_timmean.nc'.format(cmor_var,
        #                                                           mmodel))
        ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')
        print(ifilename)
        mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l', ax=ax[indx])
        datafile = Dataset(ifilename)
        lon = datafile.variables['lon'][:]
        lat = datafile.variables['lat'][:]
        lev = datafile.variables['lev'][:]
        depth_target, level_target = closest_depth(lev, depth)
        if lon.ndim == 2:
            lon2d, lat2d = lon, lat
        elif lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        #areacello = datafile.variables['areacello'][:]

        if datafile.variables[cmor_var].ndim < 4:
            data = datafile.variables[cmor_var][level_target,:,:]
        else:
            data = datafile.variables[cmor_var][0,level_target,:,:]

        if cmor_var == 'thetao':
            data = data-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'
        
        # plt.subplot(nrows, ncols, nplot)
        
        

        xx, yy = mm(lon2d, lat2d)
        mm.drawmapboundary(fill_color='0.9')
        mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
        mm.drawcoastlines(linewidth=0.1)
        image = ax[indx].contourf(xx,yy, data,
                    levels=np.round(np.linspace(-2,4,20),1),
                    cmap = cmo.balance, 
                    extend='both',
                  )
        ax[indx].set_title("{}, {} m".format(mmodel, np.round(lev[level_target],1)), size=18)
        ax[indx].set_rasterization_zorder(-1)
        

        # plt.title("{}, {} m".format(mmodel, np.round(lev[level],1)), size=18)
        # nplot=nplot+1
    for delind in range(indx+1, len(ax)):
        fig.delaxes(ax[delind])
    
    cb = fig.colorbar(image, orientation='horizontal', ax=ax, pad=0.01, shrink = 0.9)
    cb.set_label(cb_label, rotation='horizontal', size=12)
    cb.ax.tick_params(labelsize=12)

    pltoutname = genfilename(diagplotdir, cmor_var, "MULTIMODEL", data_type='plot2d_{}_depth'.format(str(depth)))
    plt.savefig(pltoutname, dpi=dpi)

def plot_aw_core_stat(aw_core_parameters, diagplotdir):
    logger.info("Plot AW core statistics")
    df = pd.DataFrame(aw_core_parameters).transpose()
    df['maxvalue'] = df.maxvalue - 273.15
    plt.figure()
    df.maxvalue.plot(kind='barh')
    plt.xlabel('$^{\circ}$C')
    pltoutname = genfilename(diagplotdir, variable='aw-core-temp', region='EB', data_type='awiCoreTemp')
    # pltoutname = os.path.join(diagplotdir,
    #                           'arctic_ocean_aw-core-temp_EB_awiCoreTemp.png')
    plt.tight_layout()
    plt.savefig(pltoutname, dpi=100)
    
    plt.figure()
    df.maxvalue_depth.plot(kind='barh')
    plt.xlabel('m')
    pltoutname = genfilename(diagplotdir, variable='aw-core-depth', region='EB', data_type='awiCoreTemp')
    # pltoutname = os.path.join(diagplotdir,
    #                           'arctic_ocean_aw-core-depth_EB_awiCoreTemp.png')
    plt.tight_layout()
    plt.savefig(pltoutname, dpi=100)



def plot2d_original_grid_AWdepth(model_filenames, cmor_var, 
                         aw_core_parameters, region, diagworkdir, diagplotdir, dpi=100, observations='PHC'):
    ncols = 4
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(5*ncols,1.5*nrows*ncols))
    mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    
    for mmodel in model_filenames:
        print(mmodel)
        logger.info("Plot plot2d_original_grid_AWdepth {} for {}".format(cmor_var, mmodel))
        ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')
        # ifilename = os.path.join(diagworkdir,
        #                     'arctic_ocean_{}_{}_timmean.nc'.format(cmor_var,
        #                                                           mmodel))
        datafile = Dataset(ifilename)
        lon = datafile.variables['lon'][:]
        lat = datafile.variables['lat'][:]
        lev = datafile.variables['lev'][:]
        if lon.ndim == 2:
            lon2d, lat2d = lon, lat
        elif lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        #areacello = datafile.variables['areacello'][:]
        maxvalue_index = aw_core_parameters[mmodel]['maxvalue_index']
        if mmodel != observations:
            data = datafile.variables[cmor_var][0, maxvalue_index, :, :]
        else:
            data = datafile.variables[cmor_var][maxvalue_index, :, :]
            

        if cmor_var == 'thetao':
            data = data-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'
        
        plt.subplot(nrows, ncols, nplot)

        xx, yy = mm(lon2d, lat2d)
        mm.drawmapboundary(fill_color='0.9')
        mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
        mm.drawcoastlines(linewidth=0.1)
        mm.contourf(xx,yy, data,
                    levels=np.round(np.linspace(-2, 2.3, 41), 1),
                    cmap = cmo.balance, 
                    extend='both',
                  )
        cb = plt.colorbar(orientation='horizontal', pad=0.03, shrink=0.95)
        cb.set_label(cb_label, rotation='horizontal', size=12)
        cb.ax.tick_params(labelsize=12)

        plt.title("{}, {} m".format(mmodel, np.round(lev[maxvalue_index], 1)), size=18)
        nplot = nplot+1

    # phc = Dataset('/mnt/lustre01/work/ab0995/a270088/ESMV/obsdata/Tier2/PHC/phc3.0_annual.nc')
    # temp_phc = phc.variables['temp'][:]
    # salt_phc = phc.variables['salt'][:]
    # lon_phc = phc.variables['lon'][:]
    # lat_phc = phc.variables['lat'][:]
    # depth_phc = phc.variables['depth'][:]

    # depth3d_phc = np.zeros(salt_phc.shape)
    # for i in range(depth_phc.shape[0]):
    #     depth3d_phc[i, :, :] = depth_phc[i]

    # ptemp_phc = sw.ptmp(salt_phc, temp_phc, depth3d_phc)
    # temp_phc = ptemp_phc

    # lonc, latc = np.meshgrid(lon_phc, lat_phc)

    # maxvalue_index = aw_core_parameters['PHC3']['maxvalue_index']

    # data = temp_phc[maxvalue_index, :, :]

    # if cmor_var == 'thetao':
    #     data = data
    #     cb_label = '$^{\circ}$C'
    # elif cmor_var == 'so':
    #     cb_label = 'psu'
    
    # plt.subplot(nrows, ncols, nplot)

    # xx, yy = mm(lonc, latc)
    # mm.drawmapboundary(fill_color='0.9')
    # mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
    # mm.drawcoastlines(linewidth=0.1)
    # mm.contourf(xx,yy, data,
    #             levels=np.round(np.linspace(-2, 2.3, 41),1),
    #             cmap = cmo.balance, 
    #             extend='both',
    #           )
    # cb = plt.colorbar(orientation='horizontal', pad=0.03, shrink=0.95)
    # cb.set_label(cb_label, rotation='horizontal', size=12)
    # cb.ax.tick_params(labelsize=12)

    # plt.title("{}, {} m".format('PHC3', np.round(depth_phc[maxvalue_index],1)), size=18)
    nplot=nplot+1

    plt.tight_layout()
    pltoutname = os.path.join(diagplotdir,
                              'arctic_ocean_{}_{}_plot2dAWdepth'.format(cmor_var,
                                                                  region))
    plt.savefig(pltoutname, dpi=dpi)


#####################################
######################################

def plot2d_speed(model_filenames, 
                 depth, vmax, diagworkdir, diagplotdir, dpi=100):
    ncols = 4
    nplots = len(model_filenames)+1
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(5*ncols,1.5*nrows*ncols))
    mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    for mmodel in model_filenames:
        logger.info("Plot plot2d_speed for {}".format(mmodel))
        ifilename_u = genfilename(diagworkdir, 'uo',
                             mmodel,  data_type='timmean', extension='.nc')
        ifilename_v = genfilename(diagworkdir, 'vo',
                             mmodel,  data_type='timmean', extension='.nc')
        datafile_u = Dataset(ifilename_u)
        datafile_v = Dataset(ifilename_v)
        lon = datafile_u.variables['lon'][:]
        lat = datafile_u.variables['lat'][:]
        lev = datafile_u.variables['lev'][:]
        if lon.ndim == 2:
            lon2d, lat2d = lon, lat
        elif lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        #areacello = datafile.variables['areacello'][:]
        iz=abs(abs(lev)-abs(depth)).argmin()
        target_depth = lev[iz]
        #maxvalue_index = aw_core_parameters[mmodel]['maxvalue_index']
        data_u = datafile_u.variables['uo'][0, iz, :, :]
        data_v = datafile_v.variables['vo'][0, iz, :, :]
        speed = np.hypot(data_u, data_v)
        # if cmor_var == 'thetao':
        #     data = data-273.15
        cb_label = 'm/s'
        # elif cmor_var == 'so':
        #     cb_label = 'psu'
        
        plt.subplot(nrows, ncols, nplot)
        
        

        xx, yy = mm(lon2d, lat2d)
        mm.drawmapboundary(fill_color='0.9')
        mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
        mm.drawcoastlines(linewidth=0.1)
        mm.pcolormesh(xx,yy, speed,
              vmin =0,
              vmax = vmax,
            #levels=np.round(np.linspace(0,0.05,20),4),
            cmap = cmo.speed, 
            #extend='both',
          )
        # mm.contourf(xx,yy, speed,
        #             levels=np.round(np.linspace(-2, 2.3, 41), 1),
        #             cmap = cmo.balance, 
        #             extend='both',
        #           )
        cb = plt.colorbar(orientation='horizontal', pad=0.03, shrink=0.95, format='%1.3g')
        cb.set_label(cb_label, rotation='horizontal', size=12)
        cb.ax.tick_params(labelsize=10)
        # for l in cb.ax.yaxis.get_ticklabels():
        #     #l.set_weight("bold")
        #     l.set_fontsize(10)
        plt.title("{}, {} m".format(mmodel, np.round(target_depth, 1)), size=18)
        nplot = nplot+1

        plt.tight_layout()
    # pltoutname = os.path.join(diagplotdir,
    #                           'arctic_ocean_{}_{}_plot2dAWdepth'.format(cmor_var,
    #                                                               region))
    pltoutname = genfilename(diagplotdir, 'speed', "MULTIMODEL", data_type='plot2d_speed_{}_level'.format(str(int(depth))))
    plt.savefig(pltoutname, dpi=dpi)

#######################################
#######################################

def plot2d_bias(model_filenames, cmor_var, 
                         depth, region, diagworkdir, diagplotdir, contours=[-3, 3, 21], dpi=100, observations = 'PHC'):

    ncols = 4
    nplots = len(model_filenames)+1
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(5*ncols,1.5*nrows*ncols))

    ifilename_obs =  genfilename(diagworkdir, cmor_var,
                             observations,  data_type='timmean', extension='.nc')

    mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    model_filenames = model_filenames.copy()
    del model_filenames[observations]
    
    for mmodel in model_filenames:
        logger.info("Plot plot2d_bias {} for {}".format(cmor_var, mmodel))

        ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')

        # lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_pyresample(ifilename_obs, ifilename, depth, cmor_var)
        lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_esmf(ifilename_obs, ifilename, depth, cmor_var)

        if cmor_var == 'thetao':
            data_onlev_obs_cyc = data_onlev_obs_cyc-273.15
            interpolated  = interpolated-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        plt.subplot(nrows, ncols, nplot)
        
        xx, yy = mm(lonc, latc)
        mm.drawmapboundary(fill_color='0.9')
        mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
        mm.drawcoastlines(linewidth=0.1)
        print(np.linspace(contours[0], contours[1], contours[2]))
        mm.contourf(xx,yy, interpolated - data_onlev_obs_cyc,
                    levels=np.round(np.linspace(contours[0], contours[1], contours[2]),2),
                    cmap = cmo.balance, 
                    extend='both',
                  )
        cb = plt.colorbar(orientation='horizontal', pad=0.03, shrink=0.95)
        cb.set_label(cb_label, rotation='horizontal', size=12)
        cb.ax.tick_params(labelsize=12)

        plt.title("{}, {} m".format(mmodel, np.round(target_depth,1)), size=18)
        nplot=nplot+1

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var, "MULTIMODEL", data_type='plot2d_bias_{}_level'.format(str(int(target_depth))))
    print(pltoutname)

    plt.savefig(pltoutname, dpi=dpi)

def plot2d_bias2(model_filenames, cmor_var, 
                         depth, region, diagworkdir, diagplotdir, contours=[-3, 3, 21], dpi=100, observations = 'PHC'):

    ncols = 4
    nplots = len(model_filenames)+1
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    proj = ccrs.Robinson()
    fig, ax = plt.subplots(nrows,ncols, 
                       figsize=(5*ncols,1.5*nrows*ncols),
                       subplot_kw=dict(projection=proj), constrained_layout=True)
    ax = ax.flatten()
    # plt.figure(figsize=(5*ncols,1.5*nrows*ncols))

    ifilename_obs =  genfilename(diagworkdir, cmor_var,
                             observations,  data_type='timmean', extension='.nc')

    # mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    model_filenames = model_filenames.copy()
    del model_filenames[observations]
    
    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_bias {} for {}".format(cmor_var, mmodel))

        ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')

        # lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_pyresample(ifilename_obs, ifilename, depth, cmor_var)
        lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_pyresample(ifilename_obs, ifilename, depth, cmor_var)
        joblib.dump(lonc, 'lonc')
        joblib.dump(latc, 'latc')
        joblib.dump(data_onlev_obs_cyc, 'data_onlev_obs_cyc')
        joblib.dump(interpolated, 'interpolated')

        if cmor_var == 'thetao':
            data_onlev_obs_cyc = data_onlev_obs_cyc-273.15
            interpolated  = interpolated-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        ax[ind].set_extent([-180,180,-80,90], crs=ccrs.PlateCarree())
        ax[ind].coastlines(resolution = '50m',lw=0.5)
        # plt.subplot(nrows, ncols, nplot)
        
        # xx, yy = mm(lonc, latc)
        # mm.drawmapboundary(fill_color='0.9')
        # mm.fillcontinents(color='#f0f0f0',lake_color='#f0f0f0')
        # mm.drawcoastlines(linewidth=0.1)
        print(np.linspace(contours[0], contours[1], contours[2]))
        image = ax[ind].pcolormesh(lonc,\
                 latc,\
                 interpolated - data_onlev_obs_cyc,
                 vmin = contours[0],
                 vmax = contours[1],
                #  levels=np.round(np.linspace(contours[0], contours[1], contours[2]),2),
                 transform=ccrs.PlateCarree(), 
                 cmap=cmo.balance , 
                #  extend='both'
                 )
        ax[ind].set_title("{}, {} m".format(mmodel, np.round(target_depth,1)), size=18)
        ax[ind].add_feature(cfeature.GSHHSFeature(levels=[1],scale='coarse', facecolor='lightgray'))

        # mm.contourf(xx,yy, interpolated - data_onlev_obs_cyc,
        #             levels=np.round(np.linspace(contours[0], contours[1], contours[2]),2),
        #             cmap = cmo.balance, 
        #             extend='both',
        #           )
        # cb = fig.colorbar(image, orientation='horizontal', ax=ax, pad=0.02, shrink = 0.9)
        # cb.set_label(label, size=20)
        # cb.ax.tick_params(labelsize=15)
        # fig.set_constrained_layout_pads(w_pad=4./72., h_pad=4./72.,
        # hspace=0.04, wspace=0.)

        # cb = plt.colorbar(orientation='horizontal', pad=0.03, shrink=0.95)
        # cb.set_label(cb_label, rotation='horizontal', size=12)
        # cb.ax.tick_params(labelsize=12)

        # plt.title("{}, {} m".format(mmodel, np.round(target_depth,1)), size=18)
        # nplot=nplot+1

    # plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var, "MULTIMODEL", data_type='plot2d_bias_{}_level'.format(str(int(target_depth))))
    print(pltoutname)

    plt.savefig(pltoutname, dpi=dpi)

def transect_plot(model_filenames, cmor_var,
              max_level, region, diagworkdir, diagplotdir,
              levels, ncols=3, cmap=cm.Spectral_r):

    ncols = 3
    nplots = len(model_filenames)+1
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    
    plt.figure(figsize=(8*ncols,2*nrows*ncols))
    plt.subplot(nrows, ncols, nplot)
    lon_s4new, lat_s4new = transect_points(region, mult = 2)
    g = pyproj.Geod(ellps='WGS84')
    (az12, az21, dist) = g.inv(lon_s4new[0:-1], lat_s4new[0:-1], lon_s4new[1:], lat_s4new[1:])
    dist = dist.cumsum()/1000
    dist = np.insert(dist, 0, 0)

    m = Basemap(width=8000000,height=8000000,
            resolution='l',projection='stere',
            lat_ts=40,lat_0=90,lon_0=0.)
    m.drawcoastlines()
    
    xpt,ypt = m(lon_s4new,lat_s4new)
    m.scatter(xpt,ypt,c=dist, s=10, cmap=cm.Spectral, edgecolors='none' )
    cb = plt.colorbar()
    cb.set_label('Along-track distance, km', rotation='vertical', size=15)
    nplot=nplot+1
    for mmodel in model_filenames:
        logger.info("Plot  {} data for {}, region {}".format(cmor_var,
                                                      mmodel,
                                                      region))
        ifilename = genfilename(diagworkdir, cmor_var,
                            mmodel, region, 'transect', '.npy')
        ifilename_depth = genfilename(diagworkdir, 'depth',
                                mmodel, region, 'transect', '.npy')
        ifilename_dist = genfilename(diagworkdir, 'distance',
                                mmodel, region, 'transect', '.npy')
        print(ifilename)

        data = np.load(ifilename)
        data = np.ma.masked_equal(data.T, 0)
        lev = np.load(ifilename_depth)
        dist = np.load(ifilename_dist)

        if cmor_var == 'thetao':
            data = data-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        lev_limit = lev[lev <= max_level].shape[0]+1

        plt.subplot(nrows, ncols, nplot)

        plt.contourf(dist, lev[:lev_limit], data[:lev_limit,:], levels=levels, extend='both', cmap=cmo.thermal)
        plt.gca().invert_yaxis()
        cb = plt.colorbar(pad=0.01)
        cb.set_label(cb_label, rotation='horizontal', size=15)
        plt.yticks(size=15)
        plt.ylabel('Depth, m', size=15, rotation='vertical')
        plt.ylim(max_level, 0)
        plt.xticks(size=15)
        plt.xlabel('Along-track distance, km', size=15, rotation='horizontal')
        cb.ax.tick_params(labelsize=15)

        plt.title(mmodel, size=20)
        nplot=nplot+1

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var,
                             region= region, data_type='transect')
    plt.savefig(pltoutname, dpi=100)

def plot_hist(model_filenames, cmor_var,min_level, max_level, region, 
                 bins, diagworkdir, diagplotdir, ncols = 4, dpi=100):
    nplots = len(model_filenames)+1
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(5*ncols,1.5*nrows*ncols))
    # mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    for mmodel in model_filenames:

        logger.info("Plot hist {} for {}".format(cmor_var, mmodel))
        ifilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')
                                                     
        datafile = Dataset(ifilename)
        lon = datafile.variables['lon'][:]
        lat = datafile.variables['lat'][:]
        lev = datafile.variables['lev'][:]
        if lon.ndim == 2:
            lon2d, lat2d = lon, lat
        elif lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        #areacello = datafile.variables['areacello'][:]
        
        indexesi, indexesj = hofm_regions(region, lon2d, lat2d)
        iz_min=abs(abs(lev)-abs(min_level)).argmin()
        iz_max=abs(abs(lev)-abs(max_level)).argmin()

        data = datafile.variables[cmor_var][0, iz_min:iz_max, indexesi, indexesj]
        if not isinstance(data, np.ma.MaskedArray):
                data = np.ma.masked_equal(data, 0)

        if cmor_var == 'thetao':
            data = data-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'
        
        plt.subplot(nrows, ncols, nplot)

        
        ax = sns.distplot(data.compressed(), bins = bins)
        plt.ylim(0,0.5)
        plt.xlim(bins[0]-1, bins[-1]+1)
        
        plt.title("{}".format(mmodel), size=18)
        plt.xticks(size=15)
        plt.yticks(size=15)
        plt.xlabel(cb_label, size = 15)
        plt.ylabel('pdf', size = 15)
        nplot = nplot+1

        plt.tight_layout()

    pltoutname = genfilename(diagplotdir, cmor_var,
                             region= region, data_type='hist')
    plt.savefig(pltoutname, dpi=dpi)


def meanminmax_plot(model_filenames, cmor_var,
              min_tollerance, max_tollerance, min_optimum, max_optimum, region, diagworkdir, diagplotdir,
              ylimits, ncols=3):
    ncols = 3
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots/ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(8*ncols,2*nrows*ncols))

    for mmodel in model_filenames:
        logger.info("Plot  {} data for {}, region {}".format(cmor_var,
                                                      mmodel,
                                                      region))

        ifilename = genfilename(diagworkdir, cmor_var,
                            mmodel, region, 'meanminmax', '.npy')
        # ifilename_levels = genfilename(diagworkdir, cmor_var,
        #                                mmodel, region, 'levels', '.npy')
        ifilename_time = genfilename(diagworkdir, cmor_var,
                                 mmodel, region, 'meanminmax_time', '.npy')
        print(ifilename)

        data = np.load(ifilename)
        # lev = np.load(ifilename_levels)
        time = np.load(ifilename_time)
        if cmor_var == 'thetao':
            data = data-273.15
            cb_label = '$^{\circ}$C'
        elif cmor_var == 'so':
            cb_label = 'psu'

        # lev_limit = lev[lev <= max_level].shape[0]+1

        # series_lenght = time.shape[0]

        # months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        plt.subplot(nrows, ncols, nplot)
        plt.plot(time, data[0,:])
        plt.plot(time, data[1,:])
        plt.plot(time, data[2,:])

        y1 = np.zeros(time.shape[0])
        y2 = np.zeros(time.shape[0])
        y3 = np.zeros(time.shape[0])
        y4 = np.zeros(time.shape[0])
        y1[:] = min_optimum
        y2[:] = max_optimum
        y3[:] = min_tollerance
        y4[:] = max_tollerance

        plt.fill_between(time, y1, y2, alpha = 0.2, color = 'b')
        plt.fill_between(time, y3, y1, alpha = 0.2, color='r')
        plt.fill_between(time, y2, y4, alpha = 0.2, color='r')


        # plt.contourf(months, depth, hofdata, cmap=cmap,
        #              levels=levels,
        #              extend='both')
        # plt.yticks(size=15)
        # plt.ylabel('m', size=15, rotation='horizontal')
        plt.ylim(ylimits)

        # cb = plt.colorbar(pad=0.01)
        # cb.set_label(cb_label, rotation='horizontal', size=15)
        # # cb.set_ticks(size=15)
        # cb.ax.tick_params(labelsize=15)
        # ygap = int((np.round(series_lenght/12.)/5)*12)
        # year_indexes = range(series_lenght)[::ygap]
        # year_value = []
        # for index in year_indexes:
        #     year_value.append(time[index].year)
        # plt.xticks(year_indexes, year_value, size=15)

        plt.title(mmodel, size=20)
        nplot=nplot+1

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir, cmor_var,
                             region= region, data_type='meanminmax')
    plt.savefig(pltoutname, dpi=100)