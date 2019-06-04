"""
*********************************************************************
APPLICATE/TRR Ocean Diagnostics
*********************************************************************
"""
import logging
import os
# import joblib
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
from esmvaltool.diag_scripts.ocean3d.getdata import transect_points, load_meta


def create_plot(model_filenames, ncols=3, projection=None):
    '''Creates matplotlib figure and set of axis that corespond to
    the number of models that should be plotted.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    ncols: int
        Number of columns in the plot. The number of rows
        will be calculated authomatically.
    projection: cartopy projection istance

    Returns
    -------
    fig: matplotlib figure
    ax: list with matplotlib axis, flattened
    '''
    # Calcualte number of plots
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots / float(ncols))
    ncols, nrows = int(ncols), int(nrows)

    # different projections will have to have different
    # coefficints for creating good looking figsize,
    # now the numbers are working well with NorthPolarStereo
    # and PlateCaree
    if projection:
        figsize = (5 * ncols, 1.55 * nrows * ncols)
        fig, ax = plt.subplots(nrows,
                               ncols,
                               figsize=figsize,
                               subplot_kw=dict(projection=projection),
                               constrained_layout=True)
    # this workd well for usual plots
    else:
        figsize = (8 * ncols, 2 * nrows * ncols)
        fig, ax = plt.subplots(nrows, ncols, figsize=figsize)
    if isinstance(ax, np.ndarray):
        ax = ax.flatten()
    else:
        ax = [ax]
    return fig, ax

def label_and_conversion(cmor_var, data):
    ''' Converts data if needed and returns
    formatted version of the colorbar label.

    Parameters
    ----------
    cmor_var: str
        name of the cmor variable
    data: numpy array
        array with the data

    Returns
    -------
    cb_label: str
        formatted units for the cmor_var
    data: numpy array
        data, converted if needed.
    '''
    if cmor_var == 'thetao':
        data = data - 273.15
        cb_label = '$^{\circ}$C'
    elif cmor_var == 'so':
        cb_label = 'psu'
    return cb_label, data

def hofm_plot(model_filenames,
              cmor_var,
              max_level,
              region,
              diagworkdir,
              diagplotdir,
              levels,
              ncols=3,
              cmap=cm.Spectral_r,
              observations='PHC'):
    '''Plot Hovmoeller diagram from data at diagworkdir

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the CMOR variable
    max_level: float
        maximum depth level the Hovmoeller diagrams should go to.
    region: str
        name of the region predefined in `hofm_regions` function.
    diagworkdir: str
        path to work directory.
    diagplotdir: str
        path to plotting directory.
    levels: list
        levels for the contour plot.
    ncols: str
        number of columns in the resulting plot
        (raws will be calculated from total number of plots)
    cmap: matplotlib.cmap object
        color map
    observations: str
        name of the dataset with observations

    Returns
    -------
    None
    '''
    model_filenames = model_filenames.copy()
    if observations:
        del model_filenames[observations]

    fig, ax = create_plot(model_filenames)

    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot  %s data for %s, region %s", cmor_var, mmodel,
                    region)

        ifilename = genfilename(diagworkdir, cmor_var, mmodel, region, 'hofm',
                                '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var, mmodel, region,
                                       'levels', '.npy')
        ifilename_time = genfilename(diagworkdir, cmor_var, mmodel, region,
                                     'time', '.npy')

        hofdata = np.load(ifilename, allow_pickle=True)
        lev = np.load(ifilename_levels, allow_pickle=True)
        time = np.load(ifilename_time, allow_pickle=True)

        cb_label, hofdata = label_and_conversion(cmor_var, hofdata)

        lev_limit = lev[lev <= max_level].shape[0] + 1
        series_lenght = time.shape[0]
        months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        image = ax[ind].contourf(months,
                                 depth,
                                 hofdata,
                                 cmap=cmap,
                                 levels=levels,
                                 extend='both')

        ygap = int((np.round(series_lenght / 12.) / 5) * 12)
        year_indexes = list(range(series_lenght)[::ygap])
        year_value = []
        for index in year_indexes:
            year_value.append(time[index].year)

        ax[ind].set_xticks(year_indexes)
        ax[ind].set_xticklabels(year_value, size=15)
        ax[ind].set_title(mmodel, size=20)
        ax[ind].set_ylabel('m', size=15, rotation='horizontal')
        ax[ind].invert_yaxis()
        ax[ind].tick_params(axis='y', labelsize=15)

        cb = fig.colorbar(image, ax=ax[ind], pad=0.01)
        cb.set_label(cb_label, rotation='vertical', size=15)
        cb.ax.tick_params(labelsize=15)

    # delete unused axis
    for delind in range(ind + 1, len(ax)):
        fig.delaxes(ax[delind])

    plt.tight_layout()
    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             region=region,
                             data_type='hofm')
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
        hofdata = np.load(ifilename, allow_pickle=True)

        lev = np.load(ifilename_levels, allow_pickle=True)
        time = np.load(ifilename_time, allow_pickle=True)
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

        temp = np.load(ifilename_t, allow_pickle=True)
        salt = np.load(ifilename_s, allow_pickle=True)
        depth = np.load(ifilename_depth, allow_pickle=True)

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


def plot_profile(model_filenames,
                 cmor_var,
                 max_level,
                 region,
                 diagworkdir,
                 diagplotdir,
                 cmap=cm.Set2,
                 dpi=100,
                 observations='PHC'):
    '''Plot profiles from previously calculated data for
    Hovmoeller diagrams.
    '''
    level_clim = Dataset(model_filenames[observations]).variables['lev'][:]
    plt.figure(figsize=(5, 6))
    ax = plt.subplot(111)

    color = iter(cmap(np.linspace(0, 1, len(model_filenames))))
    lev_limit_clim = level_clim[level_clim <= max_level].shape[0] + 1

    mean_profile = np.zeros(
        (level_clim[:lev_limit_clim].shape[0], len(model_filenames) - 1))
    mean_profile_counter = 0

    for i, mmodel in enumerate(model_filenames):
        logger.info("Plot profile %s data for %s, region %s",
                    cmor_var, mmodel, region)
        # construct input filenames
        ifilename = genfilename(diagworkdir, cmor_var, mmodel, region, 'hofm',
                                '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var, mmodel, region,
                                       'levels', '.npy')
        # load data
        hofdata = np.load(ifilename, allow_pickle=True)
        lev = np.load(ifilename_levels, allow_pickle=True)

        # convert data if needed and set labeles
        cb_label, hofdata = label_and_conversion(cmor_var, hofdata)

        # set index for maximum level (max_level+1)
        lev_limit = lev[lev <= max_level].shape[0] + 1

        # calculate mean profile
        profile = (hofdata)[:, :].mean(axis=1)

        if mmodel != observations:
            c = next(color)
        else:
            c = 'k'

        plt.plot(profile, lev[0:lev_limit], label=mmodel, c=c)

        # interpolate to standard levels and add to mean profile
        profile_interpolated = np.interp(level_clim[:lev_limit_clim],
                                         lev[0:lev_limit], profile)
        if mmodel != observations:

            print('include {} in to the mean'.format(mmodel))
            mean_profile[:, mean_profile_counter] = profile_interpolated
            mean_profile_counter += 1

    mean_profile_mean = np.nanmean(mean_profile, axis=1)

    plt.plot(mean_profile_mean,
             level_clim[:lev_limit_clim],
             label='MODEL-MEAN',
             linestyle='--',
             color='k',
             lw=3)

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
    pltoutname = genfilename(diagplotdir, cmor_var, 'MULTIMODEL', region,
                             'profile')
    plt.savefig(pltoutname, dpi=dpi, bbox_inches='tight')


def plot2d_original_grid(model_filenames,
                         cmor_var,
                         depth,
                         levels,
                         diagworkdir,
                         diagplotdir,
                         dpi=100,
                         explicit_depths=None,
                         projection=ccrs.NorthPolarStereo(),
                         bbox=[-180, 180, 60, 90]
                        ):
    ''' Plot 2d maps on original grid using cartopy.

    Parameters
    ----------
    model_filenames:OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the variable
    depth: int
        we will plot the data on the model level that is closest to the `depth`.
        Ignored if explicit_depths is provided.
    levels: tuple
        values to be used for vmin and vmax in the form of (vmin, vmax)
    diagworkdir: str
        path to the working directory
    diagplotdir: str
        path to the plot directory
    dpi: int
        the dpi values to save the figure
    explicit_depths: dict
        Output of the `aw_core` function.
        It's a dictionary where for each model there is a maximum temperature,
        depth level in the model, index of the depth level in the model.
        If provided the `depth` parameter is excluded.
    projection: instance of cartopy projection (ccrs)
    bbox: list
        bounding box. It will be the input for cartopy `set_extent`.

    Retuns
    ------
    None
    '''

    fig, ax = create_plot(model_filenames,
                          ncols=4,
                          projection=projection)

    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_original_grid {} for {}".format(
            cmor_var, mmodel))

        ifilename = genfilename(diagworkdir,
                                cmor_var,
                                mmodel,
                                data_type='timmean',
                                extension='.nc')

        metadata = load_meta(datapath=model_filenames[mmodel], fxpath=None)
        datafile, lon2d, lat2d, lev, time, areacello = metadata

        if not explicit_depths:
            depth_target, level_target = closest_depth(lev, depth)
        else:
            level_target = explicit_depths[mmodel]['maxvalue_index']
            depth_target = lev[level_target]

        if datafile.variables[cmor_var].ndim < 4:
            data = datafile.variables[cmor_var][level_target, :, :]
        else:
            data = datafile.variables[cmor_var][0, level_target, :, :]

        cb_label, data = label_and_conversion(cmor_var, data)

        left, right, down, up = bbox

        ax[ind].set_extent([left, right, down, up], crs=ccrs.PlateCarree())
        # Only pcolormesh is working for now with cartopy,
        # contourf is failing to plot curvilinear meshes,
        # let along the unstructures ones.
        image = ax[ind].pcolormesh(
            lon2d,
            lat2d,
            data,
            vmin=levels[0],
            vmax=levels[-1],
            transform=ccrs.PlateCarree(),
            cmap=cmo.balance,
        )

        ax[ind].add_feature(
            cfeature.GSHHSFeature(levels=[1],
                                  scale="low",
                                  facecolor="lightgray"))
        ax[ind].set_title("{}, {} m".format(mmodel,
                                            np.round(lev[level_target], 1)),
                          size=18)
        ax[ind].set_rasterization_zorder(-1)

    # delete unused axis
    for delind in range(ind+1, len(ax)):
        fig.delaxes(ax[delind])

    # set common colorbar
    cb = fig.colorbar(image,
                      orientation='horizontal',
                      ax=ax,
                      pad=0.01,
                      shrink=0.9)
    cb.set_label(cb_label, rotation='horizontal', size=18)
    cb.ax.tick_params(labelsize=18)

    if not explicit_depths:
        plot_type = 'plot2d_{}_depth'.format(str(depth))
    else:
        plot_type = "plot2d_different_levels"

    # save the figure
    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             "MULTIMODEL",
                             data_type=plot_type)
    plt.savefig(pltoutname, dpi=dpi)


def plot2d_bias(model_filenames,
                cmor_var,
                depth,
                diagworkdir,
                diagplotdir,
                levels,
                dpi=100,
                observations='PHC',
                projection=ccrs.NorthPolarStereo(),
                bbox=[-180, 180, 60, 90]
                ):
    '''Plot 2d maps of the bias relative to climatology.

    Parameters
    ----------
    model_filenames:OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the variable
    depth: int
    diagworkdir: str
        path to the working directory
    diagplotdir: str
        path to the plot directory
        we will plot the data on the model level that is closest to the `depth`.
    levels: tuple
        values to be used for vmin and vmax in the form of (vmin, vmax)
    dpi: int
        the dpi values to save the figure
    observations: str
        name of the observations
    projection: instance of cartopy projection (ccrs)
    bbox: list
        bounding box. It will be the input for cartopy `set_extent`.

    Retuns
    ------
    None
    '''

    fig, ax = create_plot(model_filenames,
                          ncols=4,
                          projection=projection)

    ifilename_obs = genfilename(diagworkdir,
                                cmor_var,
                                observations,
                                data_type='timmean',
                                extension='.nc')

    # mm = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
    model_filenames = model_filenames.copy()
    del model_filenames[observations]

    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_bias {} for {}".format(cmor_var, mmodel))

        ifilename = genfilename(diagworkdir,
                                cmor_var,
                                mmodel,
                                data_type='timmean',
                                extension='.nc')

        # lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_pyresample(ifilename_obs, ifilename, depth, cmor_var)
        lonc, latc, target_depth, data_onlev_obs_cyc, interpolated = interpolate_esmf(
            ifilename_obs, ifilename, depth, cmor_var)

        metadata = load_meta(datapath=model_filenames[mmodel], fxpath=None)
        datafile, lon2d, lat2d, lev, time, areacello = metadata

        cb_label, data_onlev_obs_cyc = label_and_conversion(
            cmor_var, data_onlev_obs_cyc)
        cb_label, interpolated = label_and_conversion(cmor_var, interpolated)

        left, right, down, up = bbox

        ax[ind].set_extent([left, right, down, up], crs=ccrs.PlateCarree())
        # Only pcolormesh is working for now with cartopy,
        # contourf is failing to plot curvilinear meshes,
        # let along the unstructures ones.
        image = ax[ind].contourf(
            lonc,
            latc,
            interpolated - data_onlev_obs_cyc,
            levels=levels,
            extend='both',
            # vmin=contours[0],
            # vmax=contours[-1],
            transform=ccrs.PlateCarree(),
            cmap=cmo.balance,
        )

        ax[ind].add_feature(
            cfeature.GSHHSFeature(levels=[1],
                                  scale="low",
                                  facecolor="lightgray"))

        ax[ind].set_title("{}, {} m".format(mmodel, int(target_depth)),
                          size=18)
        ax[ind].set_rasterization_zorder(-1)

    for delind in range(ind + 1, len(ax)):
        fig.delaxes(ax[delind])

    # set common colorbar
    cb = fig.colorbar(image,
                      orientation='horizontal',
                      ax=ax,
                      pad=0.01,
                      shrink=0.9)
    cb.set_label(cb_label, rotation='horizontal', size=18)
    cb.ax.tick_params(labelsize=18)

    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             "MULTIMODEL",
                             data_type='plot2d_bias_{}_level'.format(
                                 str(int(target_depth))))
    print(pltoutname)

    plt.savefig(pltoutname, dpi=dpi)

def plot_aw_core_stat(aw_core_parameters, diagplotdir):
    ''' Generate 2 plots: depth of the AW core and
    temperature of the AW core. Use pandas plot functionality.

    Parameters
    ----------
    aw_core_parameters: dict
        Output of the `aw_core` function.
        It's a dictionary where for each model there is maximum temperature,
        depth level in the model, index of the depth level in the model.
    diagplotdir: str
        plot folder

    Returns
    -------
    None
    '''

    logger.info("Plot AW core statistics")
    # Convert dictionary to pandas Dataframe
    df = pd.DataFrame(aw_core_parameters).transpose()
    df['maxvalue'] = df.maxvalue - 273.15
    plt.figure()
    df.maxvalue.plot(kind='barh')
    plt.xlabel('$^{\circ}$C')
    pltoutname = genfilename(diagplotdir,
                             variable='aw-core-temp',
                             region='EB',
                             data_type='awiCoreTemp')

    plt.tight_layout()
    plt.savefig(pltoutname, dpi=100)

    plt.figure()
    df.maxvalue_depth.plot(kind='barh')
    plt.xlabel('m')
    pltoutname = genfilename(diagplotdir,
                             variable='aw-core-depth',
                             region='EB',
                             data_type='awiCoreTemp')

    plt.tight_layout()
    plt.savefig(pltoutname, dpi=100)


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

        data = np.load(ifilename, allow_pickle=True)
        data = np.ma.masked_equal(data.T, 0)
        lev = np.load(ifilename_depth, allow_pickle=True)
        dist = np.load(ifilename_dist, allow_pickle=True)

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
