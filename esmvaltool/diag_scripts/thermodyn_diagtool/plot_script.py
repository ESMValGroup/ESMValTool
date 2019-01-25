"""Module for plot scripts in Thermodyn_diagtool.

The module provides plots for a single model of:
- climatological mean maps of TOA, atmospheric and surface energy budgets;
- annual mean time series of TOA, atmospheric and surface energy budgets anom.;
- climatological mean maps of latent energy and water mass budgets;
- annual mean time series of latent energy and water mass budget anom.;
- meridional section of meridional enthalpy transports;
- meridional section of meridional water mass transports;
- scatter plots of atmospheric vs. oceani peak magnitudes in the two hem.;
- climatological mean maps of every component of the entropy budget.

CONTENT
    - latwgt: compute weighted average over latitudes;
    - hemean: compute hemispheric averages;
    - transport: compute meridional transports computation;
    - transp_max: compute meridional transport peak magnitudes and locations;
    - balances: produce plots/maps/scatter of energy and water mass budgets;
    - entropy: produce maps of material entropy production budget components;
    - plot_ellipse: plot an ellipse in a scatter plot;
    - pr_output: print fields to NetCDF file;
    - varatts: retrieve attributes from a NetCDF file;
    - removeif: remove file if it exists;

@author: Valerio Lembo, Meteorologisches Institut, University of Hamburg, 2018.
"""
from __future__ import division
import os
from shutil import move
import math
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import cartopy.crs as ccrs
from cdo import Cdo
from esmvaltool.diag_scripts.thermodyn_diagtool import fourier_coefficients


# pylint: disable-msg=R0914
# Sixtyone is reasonable in this case.
# pylint: disable-msg=R0915
# Two hundreds and sixteen is reasonable in this case.
# pylint: disable=too-many-arguments
# Fourteen is reasonable in this case.
# from plot_script import PlotScript
def latwgt(lat, t_r):
    """Compute weighted average over latitudes.

    Arguments:
    - lat: latitude (in degrees);
    - tr: the field to be averaged (time,lat);

    @author: Valerio Lembo, 2018.
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


def hemean(cls, hem, lat, inp):
    """Compute hemispheric averages.

    Arguments:
    - hem: a parameter for the choice of the hemisphere (1 stands for SH);
    - lat: latitude (in degrees);
    - inp: input field;

    @author: Valerio Lembo, 2018.
    """
    j_end = np.shape(inp)[1]
    zmn = latwgt(lat, inp)
    hmean = []
    if hem == 1:
        if j_end % 2 == 0:
            hmean = np.nansum(zmn[:, j_end / 2 + 1:j_end], axis=1)
        else:
            hmean = (np.nansum(zmn[:, (j_end + 3) / 2:j_end], axis=1)
                     + 0.5 * zmn[:, (j_end + 3) / 2 - 1])
    else:
        if j_end % 2 == 0:
            hmean = np.nansum(zmn[:, 1:j_end / 2], axis=1)
        else:
            hmean = (np.nansum(zmn[:, 1:(j_end - 1) / 2], axis=1)
                     + 0.5 * zmn[:, (j_end - 1) / 2 + 1])
    return hmean


def transport(zmean, gmean, lat):
    """Integrate the energy/water mass budgets to obtain meridional transp.

    Arguments:
    - zmean: zonal mean input fields;
    - gmean: the global mean of the input fields;
    - lat: a latitudinal array (in degrees of latitude);

    @author: Valerio Lembo, 2018.
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
        cumb[:, j_l] = (-2 *
                        np.nansum(latwgt(lat[j_l:len(lat)],
                                         zmn_ub[:, j_l:len(lat)]),
                                  axis=1))
    r_earth = 6.371 * 10 ** 6
    transp = 2 * p_i * cumb * r_earth * r_earth
    return [zmn_ub, transp]


def transp_max(lat, transp, lim):
    """Obtain transport peak magnitude and location from interpolation.

    Arguments:
    - lat: a latitudinal array;
    - transp: the meridional transport a 1D array (lat);
    - lim: limits to constrain the peak search in
    (necessary for ocean transp.)

    @author: Valerio Lembo, 2018.
    """
    deriv = np.gradient(transp)
    x_c = zerocross1d(lat, deriv)
    y_i = np.zeros(2)
    xc_cut = np.zeros(2)
    j_p = 0
    for value in x_c:
        if abs(value) <= lim:
            xc_cut[j_p] = value
            y_i[j_p] = interpolate.interp1d(lat, transp,
                                            kind='cubic')(value)
            j_p = j_p + 1
            if j_p == 2:
                break
        else:
            pass
    return [xc_cut, y_i]


def balances(workdir, plotpath, filena, name, model_name):
    """Method for various plots related to energy and water mass budgets.

    This method provides climatological annal mean maps of TOA, atmospheric
    and surface energy budgets, time series of annual mean anomalies in the
    two hemispheres and meridional sections of meridional enthalpy
    transports. Scatter plots of oceanic vs. atmospheric meridional
    enthalpy transports are also provided.

    Arguments:
    - workdir: the working directory;
    - plotpath: the path where the plot has to be saved;
    - filena: the file containing input fields;
    - name: the name of the variable associated with the input field;
    - model_name: the name of the model to be analysed;

    @author: Valerio Lembo, 2018.
    """
    cdo = Cdo()
    nsub = len(filena)
    sep = '.nc'
    path = plotpath
    model = model_name
    timesery = np.zeros([nsub, 2])
    # Import files
    filena[0] = filena[0].split(sep, 1)[0]
    filename = filena[0] + '.nc'
    dataset = Dataset(filename)
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]
    time = dataset.variables['time'][:]
    nlats = len(lats)
    nlons = len(lons)
    ntime = len(time)
    yr_0 = len(time) / 12
    timey = np.linspace(0, yr_0 - 1, num=yr_0)
    var = np.zeros([nsub, ntime, nlats, nlons])
    for i in np.arange(nsub):
        filena[i] = filena[i].split(sep, 1)[0]
        filename = filena[i] + '.nc'
        dataset = Dataset(filename)
        var[i, :, :, :] = dataset.variables[name[i]][:, :, :]
    # Compute annual mean values
    var_r = np.reshape(var, (nsub, np.shape(var)[1] / 12, 12,
                             nlats, nlons))
    vary = np.nanmean(var_r, axis=2)
    # Compute the zonal mean
    zmean = np.nanmean(vary, axis=3)
    # Compute the climatological mean map
    tmean = np.nanmean(vary, axis=1)
    # Compute global and hemispheric means as function of years
    transp_mean = np.zeros([nsub, nlats])
    lat_maxm = np.zeros([nsub, 2, len(timey)])
    tr_maxm = np.zeros([nsub, 2, len(timey)])
    lim = [55, 55, 25]
    for i_f in np.arange(nsub):
        zmean_w = latwgt(lats, zmean[i_f, :, :])
        gmean = np.nansum(zmean_w, axis=1)
        shmean = hemean(0, lats, zmean[i_f, :, :])
        nhmean = hemean(1, lats, zmean[i_f, :, :])
        timeser = np.column_stack((gmean, shmean, nhmean))
        # Compute transports
        transp = transport(zmean[i_f, :, :], gmean, lats)
        transpp = transp[1]
        transp_mean[i_f, :] = np.nanmean(transpp, axis=0)
        yr_ext = []
        lat_max = list()
        tr_max = list()
        for t_t in range(len(timey)):
            yr_ext = transp_max(lats, transpp[t_t, :], lim[i_f])
            lat_max.append(yr_ext[0])
            tr_max.append(yr_ext[1])
        for t_t in range(len(timey)):
            lat_maxm[i_f, :, t_t] = lat_max[t_t]
            tr_maxm[i_f, :, t_t] = tr_max[t_t]
    c_m = 'bwr'
    if nsub == 3:
        ext_name = ['TOA Energy Budget', 'Atmospheric Energy Budget',
                    'Surface Energy Budget']
        timesery[0, :] = (-2, 2)
        rangect = [-100, 100]
        transpty = (-6E15, 6E15)
        timesery[1, :] = (-1, 1)
        timesery[2, :] = (-3, 3)
        fig = plt.figure(figsize=(12, 22))
        axi = plt.subplot(311, projection=ccrs.PlateCarree())
        axi.coastlines()
        plt.contourf(lons, lats, tmean[0, :, :], 60,
                     transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean[0, :, :], vmin=rangect[0],
                   vmax=rangect[1], cmap=c_m, antialiaseds='True')
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name[0]))
        plt.grid()
        axi = plt.subplot(312, projection=ccrs.PlateCarree())
        axi.coastlines()
        plt.contourf(lons, lats, tmean[1, :, :], 60,
                     transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean[1, :, :], vmin=rangect[0],
                   vmax=rangect[1], cmap=c_m, antialiaseds='True')
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name[1]))
        plt.grid()
        axi = plt.subplot(313, projection=ccrs.PlateCarree())
        axi.coastlines()
        plt.contourf(lons, lats, tmean[2, :, :], 60,
                     transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean[2, :, :], vmin=rangect[0],
                   vmax=rangect[1], cmap=c_m, antialiaseds='True')
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name[2]))
        plt.grid()
        plt.savefig(path + '/{}_energy_climap.png'.format(model))
        plt.close(fig)
        fig = plt.figure()
        axi = plt.subplot(111)
        for i in np.arange(nsub):
            filename = filena[i] + '.nc'
            if name[i] == 'toab':
                nameout = 'total'
            elif name[i] == 'atmb':
                nameout = 'atmos'
            elif name[i] == 'surb':
                nameout = 'ocean'
            nc_f = workdir + '/{}_transp_mean_{}.nc'.format(nameout,
                                                            model_name)
            removeif(nc_f)
            pr_output(transp_mean[i, :], filename, nc_f, nameout)
            name_model = '{}_{}'.format(nameout, model_name)
            lat_model = 'lat_{}'.format(model_name)
            cdo.chname('{},{}'.format(nameout, name_model), input=nc_f,
                       output='aux.nc')
            move('aux.nc', nc_f)
            cdo.chname('lat,{}'.format(lat_model), input=nc_f,
                       output='aux.nc')
            move('aux.nc', nc_f)
            plt.plot(lats, transp_mean[i, :])
        plt.title('Meridional heat transports')
        plt.xlabel('Latitude [deg]', fontsize=10)
        plt.ylabel('[W]', fontsize=10)
        plt.tight_layout()
        plt.ylim(transpty)
        plt.xlim(-90, 90)
        plt.grid()
        plt.savefig(path + '/{}_transp.png'.format(model))
        plt.close(fig)
        colors = (0, 0, 0)
        fig = plt.figure()
        fig.set_size_inches(12, 12)
        axi = plt.subplot(221)
        axi.set_figsize = (50, 50)
        plt.scatter(tr_maxm[1, 0, :], tr_maxm[2, 0, :], c=colors, alpha=1)
        plt.title('(a) Atm. vs ocean magnitude - SH', fontsize=13, y=1.02)
        plt.xlabel('Atmos. trans. [W]', fontsize=11)
        plt.ylabel('Oceanic trans. [W]', fontsize=11)
        plt.grid()
        axi = plt.subplot(222)
        axi.set_figsize = (50, 50)
        plt.scatter(tr_maxm[1, 1, :], tr_maxm[2, 1, :], c=colors, alpha=1)
        plt.title('(b) Atm. vs ocean magnitude - NH', fontsize=13, y=1.02)
        plt.xlabel('Atmos. trans. [W]', fontsize=11)
        plt.ylabel('Oceanic trans. [W]', fontsize=11)
        plt.grid()
        axi = plt.subplot(223)
        axi.set_figsize = (50, 50)
        plt.scatter(lat_maxm[1, 0, :], lat_maxm[2, 0, :], c=colors,
                    alpha=1)
        plt.title('(c) Atm. vs ocean location - SH', fontsize=13, y=1.02)
        plt.xlabel('Atmos. trans. position [degrees of latitude]',
                   fontsize=11)
        plt.ylabel('Oceanic trans. position [degrees of latitude]',
                   fontsize=11)
        plt.grid()
        axi = plt.subplot(224)
        axi.set_figsize = (50, 50)
        plt.scatter(lat_maxm[1, 1, :], lat_maxm[2, 1, :], c=colors, alpha=1)
        plt.title('(d) Atm. vs ocean location - NH', fontsize=13, y=1.02)
        plt.xlabel('Atmos. trans. position [degrees of latitude]',
                   fontsize=11)
        plt.ylabel('Oceanic trans. position [degrees of latitude]',
                   fontsize=11)
        plt.grid()
        plt.savefig(path + '/{}_scatpeak.png'.format(model))
        plt.close(fig)
    elif nsub == 2:
        ext_name = ['Water mass budget', 'Latent heat budget']
        timesery[0, :] = (-3E-6, 3E-6)
        rangecw = [-1E-4, 1E-4]
        transpwy = (-2E9, 2E9)
        timesery[1, :] = (-20, 20)
        rangecl = [-150, 150]
        transply = (-6E15, 6E15)

        fig = plt.figure()
        axi = plt.subplot(111, projection=ccrs.PlateCarree())
        axi.coastlines()
        plt.contourf(lons, lats, tmean[0, :, :], 60,
                     transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean[0, :, :], vmin=rangecw[0],
                   vmax=rangecw[1], cmap=c_m, antialiaseds='True')
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name[0]))
        plt.grid()
        plt.savefig(path + '/{}_{}_climap.png'.format(model, name[0]))
        plt.close(fig)
        fig = plt.figure()
        axi = plt.subplot(111, projection=ccrs.PlateCarree())
        axi.coastlines()
        plt.contourf(lons, lats, tmean[1, :, :], 60,
                     transform=ccrs.PlateCarree())
        plt.pcolor(lons, lats, tmean[1, :, :], vmin=rangecl[0],
                   vmax=rangecl[1], cmap=c_m, antialiaseds='True')
        plt.colorbar()
        plt.title('Climatological Mean {}'.format(ext_name[1]))
        plt.grid()
        plt.savefig(path + '/{}_{}_climap.png'.format(model, name[1]))
        plt.close(fig)
        fig = plt.figure()
        axi = plt.subplot(111)
        nc_f = workdir + '/{}_transp_mean_{}.nc'.format('wmass', model)
        removeif(nc_f)
        pr_output(transp_mean[0, :], filename, nc_f, 'wmass')
        plt.plot(lats, transp_mean[0, :])
        plt.title('Water mass transports', fontsize=10)
        plt.xlabel('Latitude [deg]', fontsize=10)
        plt.ylabel('[W]')
        plt.tight_layout()
        plt.ylim(transpwy)
        plt.xlim(-90, 90)
        plt.grid()
        plt.savefig(path + '/{}_wmass_transp.png'.format(model))
        plt.close(fig)
        fig = plt.figure()
        axi = plt.subplot(111)
        nc_f = workdir + '/{}_transp_mean_{}.nc'.format('latent', model)
        removeif(nc_f)
        pr_output(transp_mean[1, :], filename, nc_f, 'latent')
        plt.plot(lats, transp_mean[1, :])
        plt.title('Latent heat transports', fontsize=10)
        plt.xlabel('Latitude [deg]', fontsize=10)
        plt.ylabel('[W]')
        plt.tight_layout()
        plt.ylim(transply)
        plt.xlim(-90, 90)
        plt.grid()
        plt.savefig(path + '/{}_latent_transp.png'.format(model))
        plt.close(fig)
    else:
        quit()
    for i_f in np.arange(nsub):
        fig = plt.figure()
        axi = plt.subplot(111)
        axi.plot(timey, timeser[:, 0], 'k', label='Global')
        axi.plot(timey, timeser[:, 1], 'r', label='SH')
        axi.plot(timey, timeser[:, 2], 'b', label='NH')
        plt.title('Annual mean {}'.format(ext_name[i_f]))
        plt.xlabel('Years')
        plt.ylabel('[W/m2]')
        axi.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
                   shadow=True, ncol=3)
        plt.tight_layout()
        plt.ylim(timesery[i_f, :])
        plt.grid()
        plt.savefig(path + '/{}_{}_timeser.png'.format(model, name[i_f]))
        plt.close(fig)


# flake8: noqa
def entropy(plotpath, filename, name, ext_name, model_name):
    """Method for plots of annual mean maps of mat. entr. prod.

    Arguments:
    - plotpath: the path where the plot has to be saved;
    - filename: the file containing input fields;
    - name: the name of the variable associated with the input field;
    - ext_name: the long name of the input field
    - model_name: the name of the model to be analysed;

    @author: Valerio Lembo, 2018.
    """
    path = plotpath
    model = model_name
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
    elif ext_name == 'Phase changes ice -> rain entropy production':
        rangec = [0, 0.05]
        c_m = 'YlOrBr'
    elif ext_name == 'Phase changes vapor -> snow entropy production':
        rangec = [0, 0.001]
        c_m = 'YlOrBr'
    elif ext_name == 'Snow melting entropy production':
        rangec = [0, 0.05]
        c_m = 'YlOrBr'
    elif ext_name == 'Potential energy entropy production':
        rangec = [0, 0.1]
        c_m = 'YlOrBr'
    else:
        quit()
    dataset = Dataset(filename)
    var = dataset.variables[name][:, :, :]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]
    # Compute the climatological mean map
    tmean = np.nanmean(var, axis=0)
    fig = plt.figure()
    axi = plt.axes(projection=ccrs.PlateCarree())
    axi.coastlines()
    plt.contourf(lons, lats, tmean, 60, transform=ccrs.PlateCarree())
    plt.pcolor(lons, lats, tmean, vmin=rangec[0], vmax=rangec[1],
               cmap=c_m, antialiaseds='True')
    plt.colorbar()
    plt.title('Climatological Mean {}'.format(ext_name))
    plt.tight_layout()
    plt.grid()
    plt.savefig(path + '/{}_{}_climap.png'.format(model, name))
    plt.close(fig)


def plot_ellipse(semimaj, semimin, phi, x_cent, y_cent, a_x):
    """A simple method for plotting ellipses in Python.

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

    @author: Nicholas Kern, 2016 - revised by Valerio Lembo, 2018
    """
    # Generate data for ellipse structure
    theta = np.linspace(0, 2 * np.pi, 100)
    r_r = 1 / np.sqrt((np.cos(theta)) ** 2 + (np.sin(theta)) ** 2)
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
    # Plot!
    a_x.plot(data[0], data[1], color='b', linestyle='-')


def pr_output(varout, filep, nc_f, nameout):
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

    PROGRAMMER(S)
        Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    fourc = fourier_coefficients()
    nc_fid = Dataset(filep, 'r')
    # Writing NetCDF files
    w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
    w_nc_fid.description = "Total, atmospheric and oceanic annual \
                            mean meridional heat transports"
    fourc.extr_lat(nc_fid, w_nc_fid)
    w_nc_var = w_nc_fid.createVariable(nameout, 'f8', ('lat'))
    varatts(w_nc_var, nameout)
    w_nc_fid.variables[nameout][:] = varout
    w_nc_fid.close()
    nc_fid.close()

def varatts(w_nc_var, varname):
    """Add attibutes to the variables, depending on name and time res.

    Arguments:
    - w_nc_var: a variable object;
    - varname: the name of the variable, among total, atmos, ocean, wmass,
    latent;

    @author: Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    if varname == 'total':
        w_nc_var.setncatts({'long_name': u"Total merid. heat transport",
                            'units': u"W", 'level_desc': 'TOA'})
    elif varname == 'atmos':
        w_nc_var.setncatts({'long_name': u"Atmos. merid. heat transport",
                            'units': u"W",
                            'level_desc': 'Vertically integrated'})
    elif varname == 'ocean':
        w_nc_var.setncatts({'long_name': u"Ocean. merid. heat transport",
                            'units': u"W", 'level_desc': 'sfc'})
    elif varname == 'wmass':
        w_nc_var.setncatts({'long_name': u"Merid. water mass transport",
                            'units': u"W", 'level_desc': 'sfc'})
    elif varname == 'latent':
        w_nc_var.setncatts({'long_name': u"Merid. latent heat transport",
                            'units': u"W", 'level_desc': 'sfc'})

    
def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def zerocross1d(x_x, y_y):
    """Find the zero crossing points in 1d data.

    Find the zero crossing events in a discrete data set.
    Linear interpolation is used to determine the actual
    locations of the zero crossing between two data points
    showing a change in sign. Data point which are zero
    are counted in as zero crossings if a sign change occurs
    across them. Note that the first and last data point will
    not be considered whether or not they are zero.

    Parameters
    ----------
    x_x, y_y : arrays. Ordinate and abscissa data values.

    Returns
    -------
    xvals : array. The locations of the zero crossing events.

    Credits
    -------
    The PyA group (https://github.com/sczesla/PyAstronomy).
    Modified by Valerio Lembo (valerio.lembo@uni-hamburg.de).

    License
    -------
    Copyright (c) 2011, PyA group

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files
    (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge,
    publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    """
    # Indices of points *before* zero-crossing
    indi = np.where(y_y[1:] * y_y[0:-1] < 0.0)[0]
    # Find the zero crossing by linear interpolation
    d_x = x_x[indi + 1] - x_x[indi]
    d_y = y_y[indi + 1] - y_y[indi]
    z_c = - y_y[indi] * (d_x / d_y) + x_x[indi]
    # What about the points, which are actually zero
    z_i = np.where(y_y == 0.0)[0]
    # Do nothing about the first and last point should they be zero
    z_i = z_i[np.where((z_i > 0) & (z_i < x_x.size - 1))]
    # Select those point, where zero is crossed
    # (sign change across the point)
    z_i = z_i[np.where(y_y[z_i - 1] * y_y[z_i + 1] < 0.0)]
    # Concatenate indices
    zzindi = np.concatenate((indi, z_i))
    # Concatenate zc and locations corresponding to zi
    z_z = np.concatenate((z_c, x_x[z_i]))
    # Sort by x-value
    sind = np.argsort(z_z)
    z_z, zzindi = z_z[sind], zzindi[sind]
    return z_z
