"""FUNCTIONS FOR PLOTS.

Plotting module for Thermodyn_diagtool.

The module provides plots for a single model of:
- climatological mean maps of TOA, atmospheric and surface energy budgets;
- annual mean time series of TOA, atmospheric and surface energy budgets anom.;
- climatological mean maps of latent energy and water mass budgets;
- annual mean time series of latent energy and water mass budget anom.;
- meridional section of meridional enthalpy transports;
- meridional section of meridional water mass transports;
- scatter plots of atmospheric vs. oceani peak magnitudes in the two hem.;
- climatological mean maps of every component of the entropy budget.

@author: valerio.lembo@uni-hamburg.de, Valerio Lembo, Hamburg University, 2018.
"""
import math
import os
from shutil import move

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from cdo import Cdo
from matplotlib import rcParams
from netCDF4 import Dataset
from scipy import interpolate, stats

from esmvaltool.diag_scripts.shared import ProvenanceLogger
from esmvaltool.diag_scripts.thermodyn_diagtool import (fourier_coefficients,
                                                        provenance_meta)


def balances(cfg, wdir, plotpath, filena, name, model):
    """Plot everything related to energy and water mass budgets.

    This method provides climatological annal mean maps of TOA, atmospheric
    and surface energy budgets, time series of annual mean anomalies in the
    two hemispheres and meridional sections of meridional enthalpy
    transports. Scatter plots of oceanic vs. atmospheric meridional
    enthalpy transports are also provided.

    Arguments:
    - wdir: the working directory;
    - plotpath: the path where the plot has to be saved;
    - filena: the files containing input fields;
    - name: the name of the variable associated with the input field;
    - model: the name of the model to be analysed;
    """
    cdo = Cdo()
    provlog = ProvenanceLogger(cfg)
    nsub = len(filena)
    pdir = plotpath
    plotentname = pdir + '/{}_heat_transp.png'.format(model)
    plotwmbname = pdir + '/{}_wmb_transp.png'.format(model)
    plotlatname = pdir + '/{}_latent_transp.png'.format(model)

    # timesery = np.zeros([nsub, 2])
    dims, ndims, tmean, zmean, timeser = global_averages(nsub, filena, name)
    transp_mean = np.zeros([nsub, ndims[1]])
    lat_maxm = np.zeros([nsub, 2, len(dims[3])])
    tr_maxm = np.zeros([nsub, 2, len(dims[3])])
    lim = [55, 55, 25]
    for i_f in np.arange(nsub):
        transp = transport(zmean[i_f, :, :], timeser[i_f, :, 0], dims[1])
        transp_mean[i_f, :], list_peak = transports_preproc(
            dims[1], ndims[3], lim[i_f], transp)
        lat_maxm[i_f, :, :] = list_peak[0]
        tr_maxm[i_f, :, :] = list_peak[1]
    if nsub == 3:
        ext_name = [
            'TOA Energy Budget', 'Atmospheric Energy Budget',
            'Surface Energy Budget'
        ]
        transpty = (-6E15, 6E15)
        coords = [dims[0], dims[1]]
        plot_climap_eb(model, pdir, coords, tmean, ext_name)
        fig = plt.figure()
        strings = ['Meridional heat transports', 'Latitude [deg]', '[W]']
        lats = dims[1]
        for i in np.arange(nsub):
            filename = filena[i] + '.nc'
            if name[i] == 'toab':
                nameout = 'total'
            elif name[i] == 'atmb':
                nameout = 'atmos'
            elif name[i] == 'surb':
                nameout = 'ocean'
            nc_f = wdir + '/{}_transp_mean_{}.nc'.format(nameout, model)
            removeif(nc_f)
            lat_model = 'lat_{}'.format(model)
            pr_output(transp_mean[i, :], filename, nc_f, nameout, lat_model)
            name_model = '{}_{}'.format(nameout, model)
            cdo.chname('{},{}'.format(nameout, name_model),
                       input=nc_f,
                       output='aux.nc')
            move('aux.nc', nc_f)
            cdo.chname('lat,{}'.format(lat_model), input=nc_f, output='aux.nc')
            move('aux.nc', nc_f)
            attr = ['{} meridional enthalpy transports'.format(nameout), model]
            provrec = provenance_meta.get_prov_transp(attr, filename,
                                                      plotentname)
            provlog.log(nc_f, provrec)
            plot_1m_transp(lats, transp_mean[i, :], transpty, strings)
        plt.grid()
        plt.savefig(plotentname)
        plt.close(fig)
        plot_1m_scatter(model, pdir, lat_maxm, tr_maxm)
    elif nsub == 2:
        ext_name = ['Water mass budget', 'Latent heat budget']
        transpwy = (-2E9, 2E9)
        transply = (-6E15, 6E15)
        coords = [dims[0], dims[1]]
        plot_climap_wm(model, pdir, coords, tmean, ext_name, name)
        nc_f = wdir + '/{}_transp_mean_{}.nc'.format('wmb', model)
        removeif(nc_f)
        filena[0] = filena[0].split('.nc', 1)[0]
        filename = filena[0] + '.nc'
        pr_output(transp_mean[0, :], filename, nc_f, 'wmb', 'lat')
        attr = ['water mass transport', model]
        provrec = provenance_meta.get_prov_transp(attr, filename, plotwmbname)
        provlog.log(nc_f, provrec)
        nc_f = wdir + '/{}_transp_mean_{}.nc'.format('latent', model)
        removeif(nc_f)
        filena[1] = filena[1].split('.nc', 1)[0]
        filename = filena[1] + '.nc'
        pr_output(transp_mean[1, :], filename, nc_f, 'latent', 'lat')
        attr = ['latent energy transport', model]
        provrec = provenance_meta.get_prov_transp(attr, filename, plotlatname)
        provlog.log(nc_f, provrec)
        strings = ['Water mass transports', 'Latitude [deg]', '[kg*s-1]']
        fig = plt.figure()
        plot_1m_transp(dims[1], transp_mean[0, :], transpwy, strings)
        plt.grid()
        plt.savefig(plotwmbname)
        plt.close(fig)
        strings = ['Latent heat transports', 'Latitude [deg]', '[W]']
        fig = plt.figure()
        plot_1m_transp(dims[1], transp_mean[1, :], transply, strings)
        plt.grid()
        plt.savefig(plotlatname)
        plt.close(fig)
    for i_f in np.arange(nsub):
        fig = plt.figure()
        axi = plt.subplot(111)
        axi.plot(dims[3], timeser[i_f, :, 0], 'k', label='Global')
        axi.plot(dims[3], timeser[i_f, :, 1], 'r', label='SH')
        axi.plot(dims[3], timeser[i_f, :, 2], 'b', label='NH')
        plt.title('Annual mean {}'.format(ext_name[i_f]))
        plt.xlabel('Years')
        plt.ylabel('[W/m2]')
        axi.legend(loc='upper center',
                   bbox_to_anchor=(0.5, -0.07),
                   shadow=True,
                   ncol=3)
        plt.tight_layout()
        plt.grid()
        plt.savefig(pdir + '/{}_{}_timeser.png'.format(model, name[i_f]))
        plt.close(fig)


def entropy(plotpath, filename, name, ext_name, model):
    """Plot everything rleated to annual mean maps of mat. entr. prod.

    Arguments:
    - plotpath: the path where the plot has to be saved;
    - filename: the file containing input fields;
    - name: the name of the variable associated with the input field;
    - ext_name: the long name of the input field
    - model: the name of the model to be analysed;
    """
    pdir = plotpath
    if ext_name == 'Vertical entropy production':
        rangec = [-0.01, 0.1]
        c_m = 'YlOrBr'
    elif ext_name == 'Horizontal entropy production':
        rangec = [-0.5, 0.5]
        c_m = 'bwr'
    elif ext_name == 'Sensible Heat entropy production':
        rangec = [-0.01, 0.01]
        c_m = 'YlOrBr'
    elif ext_name == 'Evaporation entropy production':
        rangec = [0, 1]
        c_m = 'YlOrBr'
    elif ext_name == 'Rainfall precipitation entropy production':
        rangec = [0, 1]
        c_m = 'YlOrBr'
    elif ext_name == 'Snowfall precipitation entropy production':
        rangec = [0, 0.25]
        c_m = 'YlOrBr'
    elif ext_name == 'Snow melting entropy production':
        rangec = [0, 0.05]
        c_m = 'YlOrBr'
    elif ext_name == 'Potential energy entropy production':
        rangec = [0, 0.1]
        c_m = 'YlOrBr'
    else:
        quit()
    with Dataset(filename) as dataset:
        var = dataset.variables[name][:, :, :]
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
    tmean = np.nanmean(var, axis=0)
    fig = plt.figure()
    axi = plt.axes(projection=ccrs.PlateCarree())
    coords = [lons, lats]
    title = 'Climatological Mean {}'.format(ext_name)
    plot_climap(axi, coords, tmean, title, rangec, c_m)
    plt.savefig(pdir + '/{}_{}_climap.png'.format(model, name))
    plt.close(fig)


def global_averages(nsub, filena, name):
    """Compute zonal mean, global mean, time mean averages.

    Arguments:
    - nsub: the number of variables for which averages must be computed;
    - filena: the name of the file containing the variable (without extension);
    - name: the names of the variables;
    """
    sep = '.nc'
    filena[0] = filena[0].split(sep, 1)[0]
    filename = filena[0] + sep
    with Dataset(filename) as dataset:
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        time = dataset.variables['time'][:]
    nlats = len(lats)
    nlons = len(lons)
    ntime = len(time)
    yr_0 = int(len(time) / 12)
    timey = np.linspace(0, yr_0 - 1, num=yr_0)
    dims = [lons, lats, time, timey]
    ndims = [nlons, nlats, ntime, yr_0]
    var = np.zeros([nsub, ntime, nlats, nlons])
    for i in np.arange(nsub):
        filena[i] = filena[i].split(sep, 1)[0]
        filename = filena[i] + '.nc'
        with Dataset(filename) as dataset:
            dataset = Dataset(filename)
            var[i, :, :, :] = dataset.variables[name[i]][:, :, :]
    var_r = np.reshape(var,
                       (nsub, int(np.shape(var)[1] / 12), 12, nlats, nlons))
    vary = np.nanmean(var_r, axis=2)
    zmean = np.nanmean(vary, axis=3)
    tmean = np.nanmean(vary, axis=1)
    timeser = np.zeros([nsub, yr_0, 3])
    for i_f in np.arange(nsub):
        zmean_w = latwgt(lats, zmean[i_f, :, :])
        gmean = np.nansum(zmean_w, axis=1)
        shmean = hemean(0, lats, zmean[i_f, :, :])
        nhmean = hemean(1, lats, zmean[i_f, :, :])
        timeser[i_f, :, :] = np.column_stack((gmean, shmean, nhmean))
    return dims, ndims, tmean, zmean, timeser


def hemean(hem, lat, inp):
    """Compute hemispheric averages.

    Arguments:
    - hem: a parameter for the choice of the hemisphere (1 stands for SH);
    - lat: latitude (in degrees);
    - inp: input field;
    """
    j_end = np.shape(inp)[1]
    zmn = latwgt(lat, inp)
    hmean = []
    if hem == 1:
        if j_end % 2 == 0:
            hmean = 2 * np.nansum(zmn[:, int(j_end / 2):j_end], axis=1)
        else:
            hmean = 2 * np.nansum(zmn[:, int((j_end + 1) / 2):j_end], axis=1)
    else:
        if j_end % 2 == 0:
            hmean = 2 * np.nansum(zmn[:, 1:int(j_end / 2)], axis=1)
        else:
            hmean = 2 * np.nansum(zmn[:, 1:int((j_end - 1) / 2)], axis=1)
    return hmean


def init_plotentr(model, pdir, flist):
    """Define options for plotting maps of entropy production components.

    Arguments:
    - model: the name of the model;
    - path: the path to the plots directory;
    - flist: a list of files containing the components of the entropy
      production with the direct method;
    """
    entropy(pdir, flist[0], 'ssens', 'Sensible Heat entropy production', model)
    entropy(pdir, flist[1], 'sevap', 'Evaporation entropy production', model)
    entropy(pdir, flist[2], 'srain',
            'Rainfall precipitation entropy production', model)
    entropy(pdir, flist[3], 'ssnow',
            'Snowfall precipitation entropy production', model)
    entropy(pdir, flist[4], 'smelt', 'Snow melting entropy production', model)
    entropy(pdir, flist[5], 'spotp', 'Potential energy entropy production',
            model)


def latwgt(lat, t_r):
    """Compute weighted average over latitudes.

    Arguments:
    - lat: latitude (in degrees);
    - tr: the field to be averaged (time,lat);
    """
    p_i = math.pi
    conv = 2 * p_i / 360
    dlat = np.zeros(len(lat))
    for i in range(len(lat) - 1):
        dlat[i] = abs(lat[i + 1] - lat[i])
    dlat[len(lat) - 1] = dlat[len(lat) - 2]
    latr = conv * lat
    dlatr = conv * dlat
    tr2 = np.zeros((np.shape(t_r)[0], np.shape(t_r)[1]))
    for j in range(len(lat)):
        tr2[:, j] = t_r[:, j] * np.cos(latr[j]) * dlatr[j] / 2
    return tr2


def plot_climap_eb(model, pdir, coords, tmean, ext_name):
    """Plot climatological mean maps of TOA, atmospheric, oceanic energy budg.

    Arguments:
    - model: the name of the model;
    - pdir: a plots directory;
    - coords: the lon and lat coordinates;
    - tmean: the climatological mean (3,lat,lon) maps of the three budgets;
    - ext_name: the extended name of the budget, to be used for the title;
    """
    rangect = [-100, 100]
    fig = plt.figure(figsize=(12, 22))
    axi = plt.subplot(311, projection=ccrs.PlateCarree())
    title = 'Climatological Mean {}'.format(ext_name[0])
    plot_climap(axi, coords, tmean[0, :, :], title, rangect, 'bwr')
    axi = plt.subplot(312, projection=ccrs.PlateCarree())
    title = 'Climatological Mean {}'.format(ext_name[1])
    plot_climap(axi, coords, tmean[1, :, :], title, rangect, 'bwr')
    axi = plt.subplot(313, projection=ccrs.PlateCarree())
    title = 'Climatological Mean {}'.format(ext_name[2])
    plot_climap(axi, coords, tmean[2, :, :], title, rangect, 'bwr')
    plt.savefig(pdir + '/{}_energy_climap.png'.format(model))
    plt.close(fig)


def plot_climap_wm(model, pdir, coords, tmean, ext_name, name):
    """Plot climatological mean maps of water mass and latent energy budgets.

    Arguments:
    - model: the name of the model;
    - pdir: a plots directory;
    - coords: the lon and lat coordinates;
    - tmean: the climatological mean (3,lat,lon) maps of the three budgets;
    - ext_name: the extended name of the budget, to be used for the title;
    - name: the variable name, used for the file name of the figure;
    """
    rangecw = [-1E-4, 1E-4]
    rangecl = [-150, 150]
    fig = plt.figure()
    axi = plt.subplot(111, projection=ccrs.PlateCarree())
    title = 'Climatological Mean {}'.format(ext_name[0])
    plot_climap(axi, coords, tmean[0, :, :], title, rangecw, 'bwr')
    plt.savefig(pdir + '/{}_{}_climap.png'.format(model, name[0]))
    plt.close(fig)
    fig = plt.figure()
    axi = plt.subplot(111, projection=ccrs.PlateCarree())
    title = 'Climatological Mean {}'.format(ext_name[1])
    plot_climap(axi, coords, tmean[1, :, :], title, rangecl, 'bwr')
    plt.savefig(pdir + '/{}_{}_climap.png'.format(model, name[1]))
    plt.close(fig)


def plot_climap(axi, coords, fld, title, rrange, c_m):
    """Plot very colourful maps.

    Arguments:
    - axi: an axis identifier;
    - coords: the lon and lat coordinates;
    - fld: the field to be plotted;
    - title: the title to appear on the figure;
    - rrange: the range for the color bar;
    - c_m: a color map identifier;
    """
    axi.coastlines()
    lons = np.linspace(0, 360, len(coords[0])) - (coords[0][1] - coords[0][0])
    plt.contourf(lons, coords[1], fld, 60, transform=ccrs.PlateCarree())
    plt.pcolor(lons,
               coords[1],
               fld,
               vmin=rrange[0],
               vmax=rrange[1],
               cmap=c_m,
               antialiaseds='True')
    plt.colorbar()
    plt.title(title)
    plt.grid()


def plot_ellipse(semimaj, semimin, phi, x_cent, y_cent, a_x):
    """Plot ellipses in Python in a simple way.

    This method plots ellipses with matplotlib.

    Arguments:
    - semimaj: the length of the major axis;
    - semimin: the length of the minor axis;
    - phi: the tilting of the semimaj axis;
    - (x_cent, y_cent): the coordinates of the ellipse centre;
    - theta_num: the number of points to sample along ellipse from 0-2pi;
    - ax: an object containing the axis properties;
    - plot_kwargs: matplotlib.plot keyword arguments;
    - fill: a flag to fill the inside of the ellipse;
    - fill_kwargs: keyword arguments for matplotlib.fill;
    - data_out: a flag to return the ellipse samples without plotting;
    - cov: a 2x2 covariance matrix;
    - mass_level: a number defining the fractional probability enclosed, if
    cov is given;
    """
    theta = np.linspace(0, 2 * np.pi, 100)
    r_r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
    x_x = r_r * np.cos(theta)
    y_x = r_r * np.sin(theta)
    data = np.array([x_x, y_x])
    s_ax = np.array([[semimaj, 0], [0, semimin]])
    r_angle = np.array([[np.cos(phi), -np.sin(phi)],
                        [np.sin(phi), np.cos(phi)]])
    t_t = np.dot(r_angle, s_ax)
    data = np.dot(t_t, data)
    data[0] += x_cent
    data[1] += y_cent
    a_x.plot(data[0], data[1], color='b', linestyle='-')


def plot_1m_scatter(model, pdir, lat_maxm, tr_maxm):
    """Plot the scatter plots of atmospheric vs. oceanic peaks and locations.

    The function produces scatter plots for the atmospheric vs. oceanic peak
    magnitudes in the NH (a) and SH (b), atmospheric vs. ocean peak locations
    in the NH (c) and SH (d).

    Arguments:
    - model: the name of the model;
    - pdir: a plots directory;
    - lat_maxm: the positions of the peaks;
    - tr_maxm: the magnitudes of the peaks;
    """
    fig = plt.figure()
    fig.set_size_inches(12, 12)
    axi = plt.subplot(221)
    axi.set_figsize = (50, 50)
    plt.scatter(tr_maxm[1, 0, :], tr_maxm[2, 0, :], c=(0, 0, 0), alpha=1)
    plt.title('(a) Atm. vs ocean magnitude - SH', fontsize=13, y=1.02)
    plt.xlabel('Atmos. trans. [W]', fontsize=11)
    plt.ylabel('Oceanic trans. [W]', fontsize=11)
    plt.grid()
    axi = plt.subplot(222)
    axi.set_figsize = (50, 50)
    plt.scatter(tr_maxm[1, 1, :], tr_maxm[2, 1, :], c=(0, 0, 0), alpha=1)
    plt.title('(b) Atm. vs ocean magnitude - NH', fontsize=13, y=1.02)
    plt.xlabel('Atmos. trans. [W]', fontsize=11)
    plt.ylabel('Oceanic trans. [W]', fontsize=11)
    plt.grid()
    axi = plt.subplot(223)
    axi.set_figsize = (50, 50)
    plt.scatter(lat_maxm[1, 0, :], lat_maxm[2, 0, :], c=(0, 0, 0), alpha=1)
    plt.title('(c) Atm. vs ocean location - SH', fontsize=13, y=1.02)
    plt.xlabel('Atmos. trans. position [degrees of latitude]', fontsize=11)
    plt.ylabel('Oceanic trans. position [degrees of latitude]', fontsize=11)
    plt.grid()
    axi = plt.subplot(224)
    axi.set_figsize = (50, 50)
    plt.scatter(lat_maxm[1, 1, :], lat_maxm[2, 1, :], c=(0, 0, 0), alpha=1)
    plt.title('(d) Atm. vs ocean location - NH', fontsize=13, y=1.02)
    plt.xlabel('Atmos. trans. position [degrees of latitude]', fontsize=11)
    plt.ylabel('Oceanic trans. position [degrees of latitude]', fontsize=11)
    plt.grid()
    plt.savefig(pdir + '/{}_scatpeak.png'.format(model))
    plt.close(fig)


def plot_1m_transp(lats, yval, ylim, strings):
    """Plot a meridional section of enthalpy transport for one model.

    This function plots total, atmospheric and oceanic meridional enthalpy
    transports on the same panel.

    Arguments:
    - lats: the latitudinal dimension as a 1D array;
    - yval: the meridional enthalpy transports as a 2D array (3,lat), where
    row 1 is the total, row 2 the atmospheric, row 3 the oceanic transport;
    - ylim: a range for the y-axis;
    - strings: a list of strings containing the title of the figure, the names
    of the x and y axes;
    """
    plt.subplot(111)
    plt.plot(lats, yval)
    plt.title(strings[0], fontsize=10)
    plt.xlabel(strings[1], fontsize=10)
    plt.ylabel(strings[2])
    plt.tight_layout()
    plt.ylim(ylim)
    plt.xlim(-90, 90)


def plot_mm_ebscatter(pdir, eb_list):
    """Plot multi-model scatter plots of EB mean values vs. their variability.

    The function produces a plot containing 4 scatter plots:
    - (a) TOA mean energy budget vs. its interannual variability;
    - (b) Atmospheric mean energy budget vs. its interannual variability;
    - (c) Surface mean energy budget vs. its interannual variability;
    - (d) Atmospheric vs. surface energy budget with whiskers encompassing the
    1sigma uncertainty range;

    Arguments:
    - pdir: a plots directory;
    - eb_list: a list containing the TOA, atmospheri and surface energy budgets
    as a 2D array (model, 2), with the first column being the mean value and
    the second column being the inter-annual variance;
    """
    toab_all = eb_list[0]
    atmb_all = eb_list[1]
    surb_all = eb_list[2]
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    axi = plt.subplot(221)
    plt.ylim(bottom=0)
    title = '(a) TOA energy budget'
    xlabel = 'R_t [W m-2]'
    ylabel = 'Sigma (R_t) [W m-2]'
    varlist = [toab_all[:, 0], toab_all[:, 1]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(222)
    plt.ylim(bottom=0)
    title = '(b) Atmospheric energy budget'
    xlabel = 'F_a [W m-2]'
    ylabel = 'Sigma (F_a) [W m-2]'
    varlist = [atmb_all[:, 0], atmb_all[:, 1]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(223)
    plt.ylim(bottom=0)
    title = '(b) Surface energy budget'
    xlabel = 'F_s [W m-2]'
    ylabel = 'Sigma (F_s) [W m-2]'
    varlist = [surb_all[:, 0], surb_all[:, 1]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(224)
    axi.set_figsize = (50, 50)
    plt.errorbar(x=atmb_all[:, 0],
                 y=surb_all[:, 0],
                 xerr=atmb_all[:, 1],
                 yerr=surb_all[:, 1],
                 fmt='none',
                 ecolor=(0, 0, 0))
    title = '(b) Atmospheric vs. Surface budget'
    xlabel = 'F_a [W m-2]'
    ylabel = 'F_s [W m-2]'
    varlist = [atmb_all[:, 0], surb_all[:, 0]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    plt.savefig(pdir + '/scatters_variability.png')
    plt.close(fig)


def plot_mm_scatter(axi, varlist, title, xlabel, ylabel):
    """Plot a multi-model scatter plot.

    The function produces a scatter plot of a multi-model ensemble, with an
    ellipse encompassing the 1sigma uncertainty around the multi-model mean.

    Arguments:
    - axi: an axis identifier;
    - varlist: a list containing the array for the x and y values (they have to
    be the same length);
    - title: a string containing the title of the plot;
    - xlabel: a string containing the x-axis label;
    - ylabel: a string containing the y-axis label;
    """
    xval = varlist[0]
    yval = varlist[1]
    modnum = len(xval)
    axi.set_figsize = (50, 50)
    plt.scatter(xval, yval, c=(0, 0, 0), alpha=1)
    plt.scatter(np.nanmean(xval), np.nanmean(yval), c='red')
    s_l, _, _, _, _ = stats.linregress(xval, yval)
    semimaj = np.max([np.nanstd(xval), np.nanstd(yval)])
    semimin = np.min([np.nanstd(xval), np.nanstd(yval)])
    plot_ellipse(semimaj,
                 semimin,
                 phi=np.arctan(s_l),
                 x_cent=np.nanmean(xval),
                 y_cent=np.nanmean(yval),
                 a_x=axi)
    plt.title(title, fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    d_x = 0.01 * (max(xval) - min(xval))
    d_y = 0.01 * (max(yval) - min(yval))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (xval[i_m], yval[i_m]),
                     xytext=(xval[i_m] + d_x, yval[i_m] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()


def plot_mm_scatter_spec(axi, varlist, title, xlabel, ylabel):
    """Plot a multi-model scatter plot ("special version").

    The function produces a scatter plot of a multi-model ensemble, with dashed
    diagonal lines containing the sum of the x and y values, an ellipse
    encompassing the 1sigma uncertainty around the multi-model mean.

    Arguments:
    - axi: an axis identifier;
    - varlist: a list containing the array for the x and y values (they have to
    be the same length);
    - title: a string containing the title of the plot;
    - xlabel: a string containing the x-axis label;
    - ylabel: a string containing the y-axis label;
    """
    xval = varlist[0]
    yval = varlist[1]
    axi.set_figsize = (50, 50)
    xrang = abs(max(xval) - min(xval))
    yrang = abs(max(yval) - min(yval))
    plt.xlim(min(xval) - 0.1 * xrang, max(xval) + 0.1 * xrang)
    plt.ylim(min(yval) - 0.1 * yrang, max(yval) + 0.1 * yrang)
    x_x = np.linspace(min(xval) - 0.1 * xrang, max(xval) + 0.1 * xrang, 10)
    y_y = np.linspace(min(yval) - 0.1 * yrang, max(yval) + 0.1 * yrang, 10)
    x_m, y_m = np.meshgrid(x_x, y_y)
    z_m = x_m + y_m
    c_p = plt.contour(x_m,
                      y_m,
                      z_m,
                      colors='black',
                      linestyles='dashed',
                      linewidths=1.)
    plt.clabel(c_p, inline=True, inline_spacing=-4, fontsize=8)
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)


def plot_mm_summaryscat(pdir, summary_varlist):
    """Plot multi-model scatter plots of some key quantities.

    The function produces a plot containing 6 scatter plots:
    - (a) TOA vs. atmospheric energy budget;
    - (b) Baroclinic efficiency vs. Intensity of LEC;
    - (c) Vertical vs. horizontal component;
    - (d) Indirect vs. direct method;
    - (e) Indirect vs. emission temperature;
    - (f) Baroclinic efficiency vs. emission temperature;

    Arguments:
    - pdir: a plots directory;
    - summary_varlist: a list containing the quantities to be plotted as a 1D
    (model) array, or a 2D array (model, 2), with the first column being the
    mean value and the second column being the inter-annual variance;
    """
    atmb_all = summary_varlist[0]
    baroceff_all = summary_varlist[1]
    horzentr_all = summary_varlist[2]
    lec_all = summary_varlist[3]
    matentr_all = summary_varlist[4]
    te_all = summary_varlist[5]
    toab_all = summary_varlist[6]
    vertentr_all = summary_varlist[7]
    indentr_all = horzentr_all[:, 0] + vertentr_all[:, 0]
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    axi = plt.subplot(321)
    title = '(a) TOA vs. atmospheric energy budget'
    xlabel = 'R_t [W m-2]'
    ylabel = 'F_a [W m-2]'
    varlist = [toab_all[:, 0], atmb_all[:, 0]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(322)
    title = '(b) Baroclinic efficiency vs. Intensity of LEC'
    xlabel = 'Eta'
    ylabel = 'W [W/m2]'
    varlist = [baroceff_all, lec_all[:, 0]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(323)
    title = '(c) Vertical vs. horizontal component'
    xlabel = 'S_hor [W m-2 K-1]'
    ylabel = 'S_ver [W m-2 K-1]'
    varlist = [horzentr_all[:, 0], vertentr_all[:, 0]]
    plot_mm_scatter_spec(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(324)
    title = '(d) Indirect vs. direct method'
    xlabel = 'S_ind [W m-2 K-1]'
    ylabel = 'S_dir [W m-2 K-1]'
    varlist = [indentr_all, matentr_all[:, 0]]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(325)
    title = '(e) Indirect vs. emission temperature'
    xlabel = 'T_E [K]'
    ylabel = 'S_mat [W m-2 K-1]'
    varlist = [te_all, indentr_all]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    axi = plt.subplot(326)
    title = '(f) Baroclinic efficiency vs. emission temperature'
    xlabel = 'T_E [K]'
    ylabel = 'Eta'
    varlist = [te_all, indentr_all]
    plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    oname = pdir + '/scatters_summary.png'
    plt.savefig(oname)
    plt.subplots_adjust(hspace=.3)


def plot_mm_transp(model_names, wdir, pdir):
    """Plot multi-model meridional enthalpy transports.

    The function plots in three panels the total, atmospheric and oceanic
    enthalpy transports, respectively.

    Arguments:
    - model_names: a list of model names contained in the ensemble;
    - wdir: a working directory;
    - pdir: a plots directory;
    """
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    axi = plt.subplot(311)
    yrange = [-6.25E15, 6.25E15]
    plot_mm_transp_panel(model_names, wdir, axi, 'total', yrange)
    axi = plt.subplot(312)
    plot_mm_transp_panel(model_names, wdir, axi, 'atmos', yrange)
    axi = plt.subplot(313)
    yrange = [-3E15, 3E15]
    plot_mm_transp_panel(model_names, wdir, axi, 'ocean', yrange)
    oname = pdir + '/meridional_transp.png'
    plt.savefig(oname)
    plt.close(fig)


def plot_mm_transp_panel(model_names, wdir, axi, domn, yrange):
    """Plot a meridional section of enthalpy transport from a model ensemble.

    Arguments:
    - model_names: a list of model names contained in the ensemble;
    - wdir: a working directory;
    - axis: the axis of the pllot;
    - domn: the domain (total, atmospheric or oceanic);
    - yrange: a range for the y-axis;
    """
    axi.set_figsize = (50, 50)
    for model in model_names:
        tot_transp_file = (wdir + '/{}_transp_mean_{}.nc'.format(domn, model))
        name = '{}_{}'.format(domn, model)
        with Dataset(tot_transp_file) as dataset:
            toat = dataset.variables[name][:]
            lats = dataset.variables['lat_{}'.format(model)][:]
        plt.plot(np.array(lats), np.array(toat), color='black', linewidth=1.)
    plt.title('(a) {} heat transports'.format(domn), fontsize=18)
    plt.xlabel('Latitude [deg]', fontsize=14)
    plt.ylabel('[W]', fontsize=14)
    plt.tight_layout()
    plt.ylim(yrange)
    plt.xlim(-90, 90)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.grid()


def pr_output(varout, filep, nc_f, nameout, latn):
    """Print processed ta field to NetCDF file.

    Save fields to NetCDF, retrieving information from an existing
    NetCDF file. Metadata are transferred from the existing file to the
    new one.

    Arguments:
        - varout: the field to be stored, with shape (time,level,lat,lon);
        - filep: the existing dataset, from where the metadata are
          retrieved. Coordinates time,level, lat and lon have to be the
          same dimension as the fields to be saved to the new files;
        - nc_f: the name of the output file;
        - nameout: the name of the variable to be saved;
        - latn: the name of the latitude dimension;

    PROGRAMMER(S)
        Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    fourc = fourier_coefficients
    nc_fid = Dataset(filep, 'r')
    w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
    w_nc_fid.description = ("Total, atmospheric and oceanic annual ",
                            "mean meridional heat transports")
    fourc.extr_lat(nc_fid, w_nc_fid, latn)
    w_nc_var = w_nc_fid.createVariable(nameout, 'f8', (latn))
    varatts(w_nc_var, nameout)
    w_nc_fid.variables[nameout][:] = varout
    w_nc_fid.close()
    nc_fid.close()


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def transport(zmean, gmean, lat):
    """Integrate the energy/water mass budgets to obtain meridional transp.

    Arguments:
    - zmean: zonal mean input fields;
    - gmean: the global mean of the input fields;
    - lat: a latitudinal array (in degrees of latitude);
    """
    p_i = math.pi
    dlat = np.zeros(len(lat))
    for i in range(len(lat) - 1):
        dlat[i] = abs(lat[i + 1] - lat[i])
    dlat[len(lat) - 1] = dlat[len(lat) - 2]
    zmn_ub = np.zeros((np.shape(zmean)[0], np.shape(zmean)[1]))
    for index, value in enumerate(gmean):
        for j_l in range(np.shape(zmean)[1]):
            zmn_ub[index, j_l] = zmean[index, j_l] - value
    zmn_ub[np.isnan(zmn_ub)] = 0
    cumb = np.zeros((np.shape(zmean)[0], np.shape(zmean)[1]))
    transp = np.zeros((np.shape(zmean)[0], np.shape(zmean)[1]))
    for j_l in range(len(lat) - 1):
        cumb[:, j_l] = (-2 * np.nansum(
            latwgt(lat[j_l:len(lat)], zmn_ub[:, j_l:len(lat)]), axis=1))
    r_earth = 6.371 * 10**6
    transp = 2 * p_i * cumb * r_earth * r_earth
    return [zmn_ub, transp]


def transp_max(lat, transp, lim):
    """Obtain transport peak magnitude and location from interpolation.

    Arguments:
    - lat: a latitudinal array;
    - transp: the meridional transport a 1D array (lat);
    - lim: limits to constrain the peak search in
    (necessary for ocean transp.)
    """
    deriv = np.gradient(transp)
    x_c = zerocross1d(lat, deriv)
    y_i = np.zeros(2)
    xc_cut = np.zeros(2)
    j_p = 0
    for value in x_c:
        if abs(value) <= lim:
            xc_cut[j_p] = value
            y_i[j_p] = interpolate.interp1d(lat, transp, kind='cubic')(value)
            j_p = j_p + 1
            if j_p == 2:
                break
    return [xc_cut, y_i]


def transports_preproc(lats, yrs, lim, transp):
    """Compute the peaks magnitude and locations of a meridional transport.

    This function computes the peaks magnitudes and locations recursively at
    each time through the function transp_max and stores them in a list.

    Arguments:
    - lats: a latitudinal array;
    - yrs: the number of years through which iterating;
    - lim: the range (-lim,lim) in which the function transp_max has to search
    for the peaks;
    - transp: the array containing the transport;
    """
    transpp = transp[1]
    transp_mean = np.nanmean(transpp, axis=0)
    yr_ext = []
    lat_maxm = np.zeros([2, yrs])
    tr_maxm = np.zeros([2, yrs])
    lat_max = list()
    tr_max = list()
    for t_t in np.arange(int(yrs)):
        yr_ext = transp_max(lats, transpp[t_t, :], lim)
        lat_max.append(yr_ext[0])
        tr_max.append(yr_ext[1])
    for t_t in np.arange(int(yrs)):
        lat_maxm[:, t_t] = lat_max[t_t]
        tr_maxm[:, t_t] = tr_max[t_t]
    list_peak = [lat_maxm, tr_maxm]
    return transp_mean, list_peak


def varatts(w_nc_var, varname):
    """Add attibutes to the variables, depending on name and time res.

    Arguments:
    - w_nc_var: a variable object;
    - varname: the name of the variable, among total, atmos, ocean, wmb,
    latent;
    """
    if varname == 'total':
        w_nc_var.setncatts({
            'long_name': "Total merid. heat transport",
            'units': "W",
            'level_desc': 'TOA'
        })
    elif varname == 'atmos':
        w_nc_var.setncatts({
            'long_name': "Atmos. merid. heat transport",
            'units': "W",
            'level_desc': 'Vertically integrated'
        })
    elif varname == 'ocean':
        w_nc_var.setncatts({
            'long_name': "Ocean. merid. heat transport",
            'units': "W",
            'level_desc': 'sfc'
        })
    elif varname == 'wmb':
        w_nc_var.setncatts({
            'long_name': "Merid. water mass transport",
            'units': "Kg*s-1",
            'level_desc': 'sfc'
        })
    elif varname == 'latent':
        w_nc_var.setncatts({
            'long_name': "Merid. latent heat transport",
            'units': "W",
            'level_desc': 'sfc'
        })


def zerocross1d(x_x, y_y):
    """Find the zero crossing points in 1d data.

    Find the zero crossing events in a discrete data set. Linear interpolation
    is used to determine the actual locations of the zero crossing between
    two data points showing a change in sign. Data point which are zero
    are counted in as zero crossings if a sign change occurs across them.
    Note that the first and last data point will not be considered whether
    or not they are zero.

    Arguments:
    x_x, y_y : arrays. Ordinate and abscissa data values.

    Credits:
    The PyA group (https://github.com/sczesla/PyAstronomy).
    Modified by Valerio Lembo (valerio.lembo@uni-hamburg.de).

    License:
    Copyright (c) 2011, PyA group.
    """
    indi = np.where(y_y[1:] * y_y[0:-1] < 0.0)[0]
    d_x = x_x[indi + 1] - x_x[indi]
    d_y = y_y[indi + 1] - y_y[indi]
    z_c = -y_y[indi] * (d_x / d_y) + x_x[indi]
    z_i = np.where(y_y == 0.0)[0]
    z_i = z_i[np.where((z_i > 0) & (z_i < x_x.size - 1))]
    z_i = z_i[np.where(y_y[z_i - 1] * y_y[z_i + 1] < 0.0)]
    zzindi = np.concatenate((indi, z_i))
    z_z = np.concatenate((z_c, x_x[z_i]))
    sind = np.argsort(z_z)
    z_z, zzindi = z_z[sind], zzindi[sind]
    return z_z
