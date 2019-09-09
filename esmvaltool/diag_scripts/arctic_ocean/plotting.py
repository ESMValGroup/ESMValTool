# -*- coding: utf-8 -*-
import logging
import math
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from netCDF4 import Dataset

from esmvaltool.diag_scripts.arctic_ocean.getdata import (load_meta,
                                                          transect_points)
from esmvaltool.diag_scripts.arctic_ocean.interpolation import (
    closest_depth, interpolate_esmf)
from esmvaltool.diag_scripts.arctic_ocean.utils import (dens_back, genfilename,
                                                        point_distance)

logger = logging.getLogger(os.path.basename(__file__))

mpl.use('agg')


def create_plot(model_filenames, ncols=3, projection=None, nplots_increment=0):
    """Create base figure for multipanel plot.

    Creates matplotlib figure and set of axis that corespond to
    the number of models that should be plotted.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    ncols: int
        Number of columns in the plot. The number of rows
        will be calculated authomatically.
    projection: cartopy projection istance
    nplots_increment: int
        allows to increase or decrease number of plots.

    Returns
    -------
    fig: matplotlib figure
    ax: list with matplotlib axis, flattened
    """
    # Calcualte number of plots on the figure
    nplots = len(model_filenames) + nplots_increment
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
        figsize = (10 * ncols, 2.5 * nrows * ncols)
        fig, ax = plt.subplots(nrows, ncols, figsize=figsize)
    # if you have more than one axis, flatten the array
    # this way it is easier to handle it.
    if isinstance(ax, np.ndarray):
        ax = ax.flatten()
    # if only one axis is created wrap it in a list.
    else:
        ax = [ax]
    return fig, ax


def label_and_conversion(cmor_var, data):
    """Convert data if needed.

    And returns formatted version of the colorbar label.

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
    """
    if cmor_var == 'thetao':
        # Check if we use K (CMIP5)
        # or degC (CMIP6)
        if data.min() < 100:
            data = data
        else:
            data = data - 273.15
        cb_label = r'$^{\circ}$C'
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
    """Plot Hovmoeller diagram from data at diagworkdir.

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
    """
    # remove "model" that contain observations,
    # since there will be no monthly data
    model_filenames = model_filenames.copy()
    if observations:
        del model_filenames[observations]
    # creata basis for the muli panel figure
    fig, ax = create_plot(model_filenames, ncols=ncols)

    # plot data on the figure, axis by axis
    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot  %s data for %s, region %s", cmor_var, mmodel,
                    region)
        # generate input filenames that
        # the data prepared by `hofm_data` function
        ifilename = genfilename(diagworkdir, cmor_var, mmodel, region, 'hofm',
                                '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var, mmodel, region,
                                       'levels', '.npy')
        ifilename_time = genfilename(diagworkdir, cmor_var, mmodel, region,
                                     'time', '.npy')
        # load the data
        hofdata = np.load(ifilename, allow_pickle=True)
        lev = np.load(ifilename_levels, allow_pickle=True)
        time = np.load(ifilename_time, allow_pickle=True)

        # convert data if needed and get labeles for colorbars
        cb_label, hofdata = label_and_conversion(cmor_var, hofdata)

        # find index of the model level closes to the `max_level`
        # and add 1, to make a plot look better
        lev_limit = lev[lev <= max_level].shape[0] + 1

        # get the length of the time series
        series_lenght = time.shape[0]

        # create 2d arrays with coordinates of time and depths
        months, depth = np.meshgrid(range(series_lenght), lev[0:lev_limit])

        # plot an image for the model on the axis (ax[ind])
        image = ax[ind].contourf(months,
                                 depth,
                                 hofdata,
                                 cmap=cmap,
                                 levels=levels,
                                 extend='both')
        # Generate tick marks with years that looks ok
        ygap = int((np.round(series_lenght / 12.) / 5) * 12)
        year_indexes = list(range(series_lenght)[::ygap])
        year_value = []
        for index in year_indexes:
            year_value.append(time[index].year)

        # set properties of the axis
        ax[ind].set_xticks(year_indexes)
        ax[ind].set_xticklabels(year_value, size=15)
        ax[ind].set_title(mmodel, size=20)
        ax[ind].set_ylabel('m', size=15, rotation='horizontal')
        ax[ind].invert_yaxis()
        ax[ind].tick_params(axis='y', labelsize=15)

        # Add a colorbar
        cb = fig.colorbar(image, ax=ax[ind], pad=0.01)
        cb.set_label(cb_label, rotation='vertical', size=15)
        cb.ax.tick_params(labelsize=15)

    # delete unused axis
    for delind in range(ind + 1, len(ax)):
        fig.delaxes(ax[delind])
    # tighten the layout
    plt.tight_layout()
    # generate the path to the output file
    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             region=region,
                             data_type='hofm')
    # save figure
    plt.savefig(pltoutname, dpi=100)


def tsplot_plot(model_filenames,
                max_level,
                region,
                diagworkdir,
                diagplotdir,
                ncols=3,
                cmap=cm.Set1):
    """Plot a TS diagram.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    max_level: float
        maximum depth level the TS plot shoud go.
    region: str
        name of the region predefined in `hofm_regions` function.
    diagworkdir: str
        path to work directory.
    diagplotdir: str
        path to plotting directory.
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
    """
    # Setup a figure
    nplots = len(model_filenames)
    ncols = float(ncols)
    nrows = math.ceil(nplots / ncols)
    ncols = int(ncols)
    nrows = int(nrows)
    nplot = 1
    plt.figure(figsize=(8 * ncols, 2 * nrows * ncols))

    # loop over models
    for mmodel in model_filenames:
        logger.info("Plot  tsplot data for %s, region %s", mmodel, region)
        # load mean data created by `tsplot_data`
        ifilename_t = genfilename(diagworkdir, 'thetao', mmodel, region,
                                  'tsplot', '.npy')
        ifilename_s = genfilename(diagworkdir, 'so', mmodel, region, 'tsplot',
                                  '.npy')
        ifilename_depth = genfilename(diagworkdir, 'depth', mmodel, region,
                                      'tsplot', '.npy')

        temp = np.load(ifilename_t, allow_pickle=True)
        salt = np.load(ifilename_s, allow_pickle=True)
        depth = np.load(ifilename_depth, allow_pickle=True)
        # Still old fashioned way to setup a plot, works best for now.
        plt.subplot(nrows, ncols, nplot)
        # calculate background with density isolines
        si2, ti2, dens = dens_back(33, 36., -2, 6)

        # convert form Kelvin if needed
        if temp.min() > 100:
            temp = temp - 273.15

        # plot the background
        cs = plt.contour(si2,
                         ti2,
                         dens,
                         colors='k',
                         levels=np.linspace(dens.min(), dens.max(), 15),
                         alpha=0.3)
        # plot the scatter plot
        plt.scatter(salt[::],
                    temp[::],
                    c=depth,
                    s=3.0,
                    cmap=cmap,
                    edgecolors='none',
                    vmax=max_level)
        # adjust the plot
        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.1f')
        plt.xlim(33, 36.)
        plt.ylim(-2.1, 6)
        plt.xlabel('Salinity', size=20)
        plt.ylabel(r'Temperature, $^{\circ}$C', size=20)
        plt.xticks(size=15)
        plt.yticks(size=15)
        # setup the colorbar
        cb = plt.colorbar(pad=0.03)
        cb.ax.get_yaxis().labelpad = 15
        cb.set_label('depth, m', rotation=270, size=20)
        cb.ax.tick_params(labelsize=15)

        plt.title(mmodel, size=20)
        nplot = nplot + 1

    plt.tight_layout()
    # save the plot
    pltoutname = genfilename(diagplotdir,
                             'tsplot',
                             region=region,
                             data_type='tsplot')
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
    """Plot profiles.

    From previously calculated data for Hovmoeller diagrams.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the CMOR variable
    region: str
        name of the region predefined in `hofm_regions` function.
    diagworkdir: str
        path to work directory.
    diagplotdir: str
        path to plotting directory.
    cmap: matplotlib.cmap object
        color map
    dpi: int
        dpi fro the output figure
    observations: str
        name of the dataset with observations

    Returns
    -------
    None
    """
    level_clim = Dataset(model_filenames[observations]).variables['lev'][:]
    plt.figure(figsize=(5, 6))
    ax = plt.subplot(111)

    color = iter(cmap(np.linspace(0, 1, len(model_filenames))))
    lev_limit_clim = level_clim[level_clim <= max_level].shape[0] + 1

    mean_profile = np.zeros(
        (level_clim[:lev_limit_clim].shape[0], len(model_filenames) - 1))
    mean_profile_counter = 0

    for mmodel in model_filenames:
        logger.info("Plot profile %s data for %s, region %s", cmor_var, mmodel,
                    region)
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

    # Here we are ploting the mean profile separately
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

    # we shift the legend and plot it
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
                         cmap,
                         dpi=100,
                         explicit_depths=None,
                         projection=ccrs.NorthPolarStereo(),
                         bbox=[-180, 180, 60, 90],
                         ncols=4):
    """Plot 2d maps on original grid using cartopy.

    Parameters
    ----------
    model_filenames:OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the variable
    depth: int
        we will plot the data on the model
        level that is closest to the `depth`.
        Ignored if explicit_depths is provided.
    levels: tuple
        values to be used for vmin and vmax in the form of (vmin, vmax)
    diagworkdir: str
        path to the working directory
    diagplotdir: str
        path to the plot directory
    cmap:  matplotlib colormap
        colormap
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
    ncols: int
        number of columns.

    Retuns
    ------
    None
    """
    fig, ax = create_plot(model_filenames, ncols=ncols, projection=projection)

    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_original_grid %s for %s", cmor_var, mmodel)

        ifilename = genfilename(diagworkdir,
                                cmor_var,
                                mmodel,
                                data_type='timmean',
                                extension='.nc')

        metadata = load_meta(ifilename, fxpath=None)
        datafile = metadata['datafile']
        lon2d = metadata['lon2d']
        lat2d = metadata['lat2d']
        lev = metadata['lev']

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
            cmap=cmap,
        )

        ax[ind].add_feature(
            cfeature.GSHHSFeature(levels=[1],
                                  scale="low",
                                  facecolor="lightgray"))
        ax[ind].set_title("{}, {} m".format(mmodel, np.round(depth_target, 1)),
                          size=18)
        ax[ind].set_rasterization_zorder(-1)

    # delete unused axis
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
                cmap,
                levels,
                dpi=100,
                observations='PHC',
                projection=ccrs.NorthPolarStereo(),
                bbox=(-180, 180, 60, 90),
                ncols=4):
    """Plot 2d maps of the bias relative to climatology.

    Parameters
    ----------
    model_filenames:OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the variable
    depth: int
        we will plot the data on the model level
        that is closest to the `depth`.
    diagworkdir: str
        path to the working directory
    diagplotdir: str
        path to the plot directory
    levels: tuple
        values to be used for vmin and vmax in the form of (vmin, vmax)
    dpi: int
        the dpi values to save the figure
    observations: str
        name of the observations
    projection: instance of cartopy projection (ccrs)
    bbox: list
        bounding box. It will be the input for cartopy `set_extent`.
    ncols: int
        number of columns.

    Retuns
    ------
    None
    """
    # setupa a base figure
    fig, ax = create_plot(model_filenames, ncols=ncols, projection=projection)
    # get the filename of observations
    ifilename_obs = genfilename(diagworkdir,
                                cmor_var,
                                observations,
                                data_type='timmean',
                                extension='.nc')
    # get the metadata for observations (we just need a size)
    metadata = load_meta(datapath=model_filenames[observations], fxpath=None)
    lon2d = metadata['lon2d']

    # Create an empty array to store the mean.
    # One point larger along long to acount for cyclic point
    model_mean = np.zeros((lon2d.shape[0], lon2d.shape[1] + 1))
    print("MODEL MEAN SHAPE {}".format(model_mean.shape))

    # delete observations from the model list
    model_filenames = model_filenames.copy()
    del model_filenames[observations]
    # loop over models
    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot plot2d_bias %s for %s", cmor_var, mmodel)
        # get the filename with the mean generated by the `timemean`
        ifilename = genfilename(diagworkdir,
                                cmor_var,
                                mmodel,
                                data_type='timmean',
                                extension='.nc')
        # do the interpolation to the observation grid
        # the output is
        lonc, latc, target_depth, data_obs, interpolated = interpolate_esmf(
            ifilename_obs, ifilename, depth, cmor_var)
        # get the label and convert data if needed
        cb_label, data_obs = label_and_conversion(cmor_var, data_obs)
        cb_label, interpolated = label_and_conversion(cmor_var, interpolated)
        # add to the mean model
        model_mean = model_mean + interpolated
        # set the map extent
        left, right, down, up = bbox
        ax[ind].set_extent([left, right, down, up], crs=ccrs.PlateCarree())
        # Only pcolormesh is working for now with cartopy,
        # contourf is failing to plot curvilinear meshes,
        # let along the unstructures ones.
        image = ax[ind].contourf(
            lonc,
            latc,
            interpolated - data_obs,
            levels=levels,
            extend='both',
            # vmin=contours[0],
            # vmax=contours[-1],
            transform=ccrs.PlateCarree(),
            cmap=cmap,
        )
        # fill continents
        ax[ind].add_feature(
            cfeature.GSHHSFeature(levels=[1],
                                  scale="low",
                                  facecolor="lightgray"))

        ax[ind].set_title("{}, {} m".format(mmodel, int(target_depth)),
                          size=18)
        ax[ind].set_rasterization_zorder(-1)
    # calculate the model mean and plot it
    model_mean = model_mean / len(model_filenames)
    ax[ind + 1].set_extent([left, right, down, up], crs=ccrs.PlateCarree())
    image = ax[ind + 1].contourf(
        lonc,
        latc,
        model_mean - data_obs,
        levels=levels,
        extend='both',
        # vmin=contours[0],
        # vmax=contours[-1],
        transform=ccrs.PlateCarree(),
        cmap=cmo.balance,
    )

    ax[ind + 1].add_feature(
        cfeature.GSHHSFeature(levels=[1], scale="low", facecolor="lightgray"))

    ax[ind + 1].set_title("Model mean bias, {} m".format(int(target_depth)),
                          size=18)
    ax[ind + 1].set_rasterization_zorder(-1)
    # delete the axis that are not needed
    for delind in range(ind + 2, len(ax)):
        fig.delaxes(ax[delind])
    # set common colorbar
    cb = fig.colorbar(image,
                      orientation='horizontal',
                      ax=ax,
                      pad=0.01,
                      shrink=0.9)
    cb.set_label(cb_label, rotation='horizontal', size=18)
    cb.ax.tick_params(labelsize=18)
    # save the picture
    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             "MULTIMODEL",
                             data_type='plot2d_bias_{}_level'.format(
                                 str(int(target_depth))))
    print(pltoutname)

    plt.savefig(pltoutname, dpi=dpi)


def plot_aw_core_stat(aw_core_parameters, diagplotdir):
    """Generate plotsa t AW core depth.

    Depth of the AW core and temperature of the AW core.
    Use pandas plot functionality.

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
    """
    logger.info("Plot AW core statistics")
    # Convert dictionary to pandas Dataframe
    df = pd.DataFrame(aw_core_parameters).transpose()

    plt.figure()
    df.maxvalue.plot(kind='barh')
    plt.xlabel(r'$^{\circ}$C')
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


def transect_map(region,
                 diagplotdir,
                 projection=ccrs.NorthPolarStereo(),
                 bbox=(-180, 180, 60, 90),
                 mult=2):
    """Plot the map with points of the transect overlayed.

    Parameters
    ----------
    region: str
        name of the region predefined in `transect_points` function.
    diagplotdir: str
        path to the plot directory
    projection: instance of cartopy ccrs
        cartopy progection
    bbox: list
        four values - [left, right, bottom, top]
    mult: int
        miltiplicator for the number of points.
        E.g. mult=2 increase the number of points 2 times.

    Returns
    -------
    None

    """
    logger.info("Create transect map for region %s", region)
    lon_s4new, lat_s4new = transect_points(region, mult=mult)
    dist = point_distance(lon_s4new, lat_s4new)
    fig, ax = plt.subplots(1,
                           1,
                           subplot_kw=dict(projection=projection),
                           constrained_layout=True)

    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    image = ax.scatter(lon_s4new,
                       lat_s4new,
                       s=10,
                       c=dist,
                       transform=ccrs.PlateCarree(),
                       cmap=cm.Spectral,
                       edgecolors='none')
    ax.coastlines(resolution="50m")

    cb = fig.colorbar(image, ax=ax)
    cb.set_label('Along-track distance, km', rotation='vertical', size=15)
    pltoutname = genfilename(diagplotdir,
                             'allvars',
                             region=region,
                             data_type='transect_map')
    plt.savefig(pltoutname, dpi=100)


def transect_plot(model_filenames,
                  cmor_var,
                  max_level,
                  region,
                  diagworkdir,
                  diagplotdir,
                  levels,
                  ncols=3,
                  cmap=cmo.thermal):
    """Plot transects.

    Parameters
    ----------
    model_filenames:OrderedDict
        OrderedDict with model names as keys and input files as values.
    cmor_var: str
        name of the variable
    region: str
        name of the region predefined in `transect_points` function.
    diagworkdir: str
        path to the working directory
    diagplotdir: str
        path to the plot directory
    levels: tuple
       contours - (minimum, maximum, number of levels)
    ncols: int
        number of columns in the plot
    cmap: matplotlib colormap instance

    Returns
    -------
    None

    """
    fig, ax = create_plot(model_filenames, ncols=ncols)

    # get transect positions and calculate distances between points
    lon_s4new, lat_s4new = transect_points(region, mult=2)
    dist = point_distance(lon_s4new, lat_s4new)

    # loop over models
    for ind, mmodel in enumerate(model_filenames):
        logger.info("Plot  %s data for %s, region %s", cmor_var, mmodel,
                    region)
        # construct file names and get the data
        ifilename = genfilename(diagworkdir, cmor_var, mmodel, region,
                                'transect', '.npy')
        ifilename_depth = genfilename(diagworkdir, 'depth', mmodel, region,
                                      'transect', '.npy')
        ifilename_dist = genfilename(diagworkdir, 'distance', mmodel, region,
                                     'transect', '.npy')

        data = np.load(ifilename, allow_pickle=True)
        data = np.ma.masked_equal(data.T, 0)
        lev = np.load(ifilename_depth, allow_pickle=True)
        dist = np.load(ifilename_dist, allow_pickle=True)
        # get labeles and convert the data
        cb_label, data = label_and_conversion(cmor_var, data)
        # index of the maximum depth
        lev_limit = lev[lev <= max_level].shape[0] + 1

        image = ax[ind].contourf(dist,
                                 lev[:lev_limit],
                                 data[:lev_limit, :],
                                 levels=levels,
                                 extend='both',
                                 cmap=cmap)
        # plot settings
        ax[ind].set_ylabel('Depth, m', size=15, rotation='vertical')
        ax[ind].set_xlabel('Along-track distance, km',
                           size=15,
                           rotation='horizontal')
        ax[ind].set_title(mmodel, size=20)
        ax[ind].set_ylim(max_level, 0)
        # ax[ind].invert_yaxis()
        ax[ind].tick_params(axis='both', labelsize=15)
        # color bar settings
        cb = fig.colorbar(image, ax=ax[ind], pad=0.01)
        cb.set_label(cb_label, rotation='vertical', size=15)
        cb.ax.tick_params(labelsize=15)
    # fig.set_constrained_layout_pads(w_pad=2./30., h_pad=2./30.,
    #     hspace=10, wspace=10)
    for delind in range(ind + 1, len(ax)):
        fig.delaxes(ax[delind])

    pltoutname = genfilename(diagplotdir,
                             cmor_var,
                             region=region,
                             data_type='transect')
    plt.savefig(pltoutname, dpi=100)
