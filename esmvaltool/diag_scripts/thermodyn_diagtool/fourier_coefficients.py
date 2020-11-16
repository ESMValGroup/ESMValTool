"""Module retrieving Fourier coefficients computation from lonlat grid.

Computation of the Fourier coefficients from lonlat grids
on pressure levels at every timestep.

The spectral truncation is determined by the number of longitudinal
gridsteps. The outputs are given as (time,level,wave,lat) where wave stands
for the zonal wavenumber. In the context of the thermodynamic diagnostic tool,
this is used for the computation of the Lorenz Energy Cycle.

@author: valerio.lembo@uni-hamburg.de, Valerio Lembo, Hamburg University, 2018.
"""

import numpy as np
from netCDF4 import Dataset

GP_RES = np.array([16, 32, 48, 64, 96, 128, 256, 384, 512, 1024, 2048, 4096])
FC_RES = np.array([5, 10, 15, 21, 31, 43, 85, 127, 171, 341, 683, 1365])
G_0 = 9.81  # Gravity acceleration
GAM = 0.0065  # Standard atmosphere lapse rate
GAS_CON = 287.0  # Gas constant
P_0 = 10000  # Reference tropospheric pressure


def fourier_coeff(tadiagfile, outfile, ta_input, tas_input):
    """Compute Fourier coefficients in lon direction.

    Arguments:
    ---------
    - tadiagfile: the name of a file to store modified t fields;
    - outfile: the name of a file to store the Fourier coefficients;
    - ta_input: the name of a file containing t,u,v,w fields;
    - tas_input: the name of a file containing t2m field.
    """
    with Dataset(ta_input) as dataset:
        lon = dataset.variables['lon'][:]
        lat = dataset.variables['lat'][:]
        lev = dataset.variables['plev'][:]
        time = dataset.variables['time'][:]
        t_a = dataset.variables['ta'][:, :, :, :]
        u_a = dataset.variables['ua'][:, :, :, :]
        v_a = dataset.variables['va'][:, :, :, :]
        wap = dataset.variables['wap'][:, :, :, :]
    nlon = len(lon)
    nlat = len(lat)
    nlev = len(lev)
    ntime = len(time)
    i = np.min(np.where(2 * nlat <= GP_RES))
    trunc = FC_RES[i] + 1
    wave2 = np.linspace(0, trunc - 1, trunc)
    with Dataset(tas_input) as dataset:
        tas = dataset.variables['tas'][:, :, :]
    tas = tas[:, ::-1, :]
    ta1_fx = np.array(t_a)
    deltat = np.zeros([ntime, nlev, nlat, nlon])
    p_s = np.full([ntime, nlat, nlon], P_0)
    for i in np.arange(nlev - 1, 0, -1):
        h_1 = np.ma.masked_where(ta1_fx[:, i, :, :] != 0, ta1_fx[:, i, :, :])
        if np.any(h_1.mask > 0):
            deltat[:, i - 1, :, :] = np.where(ta1_fx[:, i - 1, :, :] != 0,
                                              deltat[:, i - 1, :, :],
                                              (ta1_fx[:, i, :, :] - tas))
            deltat[:, i - 1, :, :] = ((1 * np.array(h_1.mask)) *
                                      np.array(deltat[:, i - 1, :, :]))
            d_p = -((P_0 * G_0 /
                     (GAM * GAS_CON)) * deltat[:, i - 1, :, :] / tas)
            p_s = np.where(ta1_fx[:, i - 1, :, :] != 0, p_s, lev[i - 1] + d_p)
            for k in np.arange(0, nlev - i - 1, 1):
                h_3 = np.ma.masked_where(ta1_fx[:, i + k, :, :] != 0,
                                         ta1_fx[:, i + k, :, :])
                if np.any(h_3.mask > 0):
                    deltat[:, i - 1, :, :] = np.where(
                        ta1_fx[:, i + k, :, :] != 0, deltat[:, i - 1, :, :],
                        (ta1_fx[:, i + k + 1, :, :] - tas))
                    d_p = -((P_0 * G_0 /
                             (GAM * GAS_CON)) * deltat[:, i - 1, :, :] / tas)
                    p_s = np.where(ta1_fx[:, i + k, :, :] != 0, p_s,
                                   lev[i + k] + d_p)
    ta2_fx = np.array(t_a)
    mask = np.zeros([nlev, ntime, nlat, nlon])
    dat = np.zeros([nlev, ntime, nlat, nlon])
    tafr_bar = np.zeros([nlev, ntime, nlat, nlon])
    deltap = np.zeros([ntime, nlev, nlat, nlon])
    for i in np.arange(nlev):
        deltap[:, i, :, :] = p_s - lev[i]
        h_2 = np.ma.masked_where(ta2_fx[:, i, :, :] == 0, ta2_fx[:, i, :, :])
        mask[i, :, :, :] = np.array(h_2.mask)
        tafr_bar[i, :, :, :] = (1 * np.array(mask[i, :, :, :]) *
                                (tas - GAM * GAS_CON /
                                 (G_0 * p_s) * deltap[:, i, :, :] * tas))
        dat[i, :, :, :] = (ta2_fx[:, i, :, :] *
                           (1 - 1 * np.array(mask[i, :, :, :])))
        t_a[:, i, :, :] = dat[i, :, :, :] + tafr_bar[i, :, :, :]
    pr_output_diag(t_a, ta_input, tadiagfile, 'ta')
    tafft_p = np.fft.fft(t_a, axis=3)[:, :, :, :int(trunc / 2)] / (nlon)
    uafft_p = np.fft.fft(u_a, axis=3)[:, :, :, :int(trunc / 2)] / (nlon)
    vafft_p = np.fft.fft(v_a, axis=3)[:, :, :, :int(trunc / 2)] / (nlon)
    wapfft_p = np.fft.fft(wap, axis=3)[:, :, :, :int(trunc / 2)] / (nlon)
    tafft = np.zeros([ntime, nlev, nlat, trunc])
    uafft = np.zeros([ntime, nlev, nlat, trunc])
    vafft = np.zeros([ntime, nlev, nlat, trunc])
    wapfft = np.zeros([ntime, nlev, nlat, trunc])
    tafft[:, :, :, 0::2] = np.real(tafft_p)
    tafft[:, :, :, 1::2] = np.imag(tafft_p)
    uafft[:, :, :, 0::2] = np.real(uafft_p)
    uafft[:, :, :, 1::2] = np.imag(uafft_p)
    vafft[:, :, :, 0::2] = np.real(vafft_p)
    vafft[:, :, :, 1::2] = np.imag(vafft_p)
    wapfft[:, :, :, 0::2] = np.real(wapfft_p)
    wapfft[:, :, :, 1::2] = np.imag(wapfft_p)
    dict_v = {'ta': tafft, 'ua': uafft, 'va': vafft, 'wap': wapfft}
    file_desc = 'Fourier coefficients'
    pr_output(dict_v, ta_input, outfile, file_desc, wave2)


def pr_output(dict_v, nc_f, fileo, file_desc, wave2):
    """Print outputs to NetCDF.

    Save fields to NetCDF, retrieving information from an existing
    NetCDF file. Metadata are transferred from the existing file to the
    new one.

    Arguments:
    ---------
        - var1, var2, var3, var4: the fields to be stored, with shape
          (time,level,wave,lon);
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Coordinates time,level and lon have to be the same
          dimension as the fields to be saved to the new files;
        - fileo: the name of the output file;
        - wave2: an array containing the zonal wavenumbers;
        - name1, name2, name3, name4: the name of the variables to be
          saved;

    PROGRAMMER(S)
        Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    # Writing NetCDF files
    with Dataset(fileo, 'w', format='NETCDF4') as var_nc_fid:
        var_nc_fid.description = file_desc
        with Dataset(nc_f, 'r') as nc_fid:
            extr_time(nc_fid, var_nc_fid)
            extr_lat(nc_fid, var_nc_fid, 'lat')
            extr_plev(nc_fid, var_nc_fid)
            # Write the wave dimension
            var_nc_fid.createDimension('wave', len(wave2))
            var_nc_fid.createVariable('wave', nc_fid.variables['plev'].dtype,
                                      ('wave', ))
        var_nc_fid.variables['wave'][:] = wave2
        for key in dict_v:
            value = dict_v[key]
            var1_nc_var = var_nc_fid.createVariable(
                key, 'f8', ('time', 'plev', 'lat', 'wave'))
            varatts(var1_nc_var, key)
            var_nc_fid.variables[key][:, :, :, :] = value


def pr_output_diag(var1, nc_f, fileo, name1):
    """Print processed ta field to NetCDF file.

    Save fields to NetCDF, retrieving information from an existing
    NetCDF file. Metadata are transferred from the existing file to the
    new one.

    Arguments:
    ---------
        - var1: the field to be stored, with shape (time,level,lat,lon);
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Coordinates time,level, lat and lon have to be the
          same dimension as the fields to be saved to the new files;
        - fileo: the name of the output file;
        - name1: the name of the variable to be saved;

    PROGRAMMER(S)
        Chris Slocum (2014), modified by Valerio Lembo (2018).
    """
    with Dataset(fileo, 'w', format='NETCDF4') as var_nc_fid:
        var_nc_fid.description = "Fourier coefficients"
        with Dataset(nc_f, 'r') as nc_fid:
            # Extract data from NetCDF file nad write them to the new file
            extr_time(nc_fid, var_nc_fid)
            extr_lat(nc_fid, var_nc_fid, 'lat')
            extr_lon(nc_fid, var_nc_fid)
            extr_plev(nc_fid, var_nc_fid)
        var1_nc_var = var_nc_fid.createVariable(name1, 'f8',
                                                ('time', 'plev', 'lat', 'lon'))
        varatts(var1_nc_var, name1)
        var_nc_fid.variables[name1][:, :, :, :] = var1


def extr_lat(nc_fid, var_nc_fid, latn):
    """Extract lat coord. from NC files and save them to a new NC file.

    Arguments:
    ---------
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Time,level and lon dimensions
          are retrieved;
        - var_nc_fid: the id of the new NC dataset previously created;
        - latn: the name of the latitude dimension;
    """
    # Extract coordinates from NetCDF file
    lats = nc_fid.variables['lat'][:]
    var_nc_fid.createDimension(latn, len(lats))
    var_nc_dim = var_nc_fid.createVariable(latn, nc_fid.variables['lat'].dtype,
                                           (latn, ))
    for ncattr in nc_fid.variables['lat'].ncattrs():
        var_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
    var_nc_fid.variables[latn][:] = lats


def extr_lon(nc_fid, var_nc_fid):
    """Extract lat coord. from NC files and save them to a new NC file.

    Arguments:
    ---------
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Time,level and lon dimensions
          are retrieved;
        - var_nc_fid: the id of the new NC dataset previously created;
    """
    # Extract coordinates from NetCDF file
    lons = nc_fid.variables['lon'][:]
    var_nc_fid.createDimension('lon', len(lons))
    var_nc_dim = var_nc_fid.createVariable('lon',
                                           nc_fid.variables['lon'].dtype,
                                           ('lon', ))
    for ncattr in nc_fid.variables['lon'].ncattrs():
        var_nc_dim.setncattr(ncattr, nc_fid.variables['lon'].getncattr(ncattr))
    var_nc_fid.variables['lon'][:] = lons


def extr_plev(nc_fid, var_nc_fid):
    """Extract plev coord. from NC files and save them to a new NC file.

    Arguments:
    ---------
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Time,level and lon dimensions
          are retrieved;
        - var_nc_fid: the id of the new NC dataset previously created;
    """
    plev = nc_fid.variables['plev'][:]
    var_nc_fid.createDimension('plev', len(plev))
    var_nc_dim = var_nc_fid.createVariable('plev',
                                           nc_fid.variables['plev'].dtype,
                                           ('plev', ))
    for ncattr in nc_fid.variables['plev'].ncattrs():
        var_nc_dim.setncattr(ncattr,
                             nc_fid.variables['plev'].getncattr(ncattr))
    var_nc_fid.variables['plev'][:] = plev


def extr_time(nc_fid, var_nc_fid):
    """Extract time coord. from NC files and save them to a new NC file.

    Arguments:
    ---------
        - nc_f: the existing dataset, from where the metadata are
          retrieved. Time,level and lon dimensions
          are retrieved;
        - var_nc_fid: the id of the new NC dataset previously created;
    """
    # Extract coordinates from NetCDF file
    time = nc_fid.variables['time'][:]
    # Using our previous dimension info, we can create the new dimensions.
    var_nc_fid.createDimension('time', len(time))
    var_nc_dim = var_nc_fid.createVariable('time',
                                           nc_fid.variables['time'].dtype,
                                           ('time', ))
    for ncattr in nc_fid.variables['time'].ncattrs():
        var_nc_dim.setncattr(ncattr,
                             nc_fid.variables['time'].getncattr(ncattr))
    var_nc_fid.variables['time'][:] = time


def varatts(w_nc_var, varname):
    """Add attibutes to the variables, depending on their name.

    Arguments:
    ---------
    - w_nc_var: a variable object;
    - varname: the name of the variable, among ta, ua, va and wap.
    """
    if varname == 'ta':
        w_nc_var.setncatts({
            'long_name': "Air temperature",
            'units': "K",
            'level_desc': 'pressure levels'
        })
    elif varname == 'ua':
        w_nc_var.setncatts({
            'long_name': "Eastward wind",
            'units': "m s-1",
            'level_desc': 'pressure levels'
        })
    elif varname == 'va':
        w_nc_var.setncatts({
            'long_name': "Northward wind",
            'units': "m s-1",
            'level_desc': 'pressure levels'
        })
    elif varname == 'wap':
        w_nc_var.setncatts({
            'long_name': 'Lagrangian tendency of '
                         'air pressure',
            'units': "Pa s-1",
            'level_desc': 'pressure levels'
        })
