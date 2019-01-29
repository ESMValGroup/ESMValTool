"""Module for computation of the auxiliary variables needed by the tool.

Here the thermodynamic diagnostic tool script computes
some auxiliary variables.

It computes equivalent potential temperatures and temperatures representative
of the sensible and latent heat exchanges in the lower layers of the
troposphere. Estimates of the boundary layer height and lifting condensation
level are also provided.

It ingests monthly mean fields of:
- specific humidity (near-surface or 3D) (hus);
- skin temperature (ts);
- surface pressure (ps);
- near-surface horizontal velocity (uas and vas);
- surface turbulent sensible heat fluxes (hfss);
- emission temperature (te).

Authors: Frank Lunkeit and Valerio Lembo (University of Hamburg)

Created on Fri Jun 15 10:06:30 2018
"""
import os
from netCDF4 import Dataset
import numpy as np
from cdo import Cdo
from esmvaltool.diag_scripts.thermodyn_diagtool import fourier_coefficients

ALV = 2.5008e6    # Latent heat of vaporization
G_0 = 9.81        # Gravity acceleration
P_0 = 100000.     # reference pressure
RV = 461.51       # Gas constant for water vapour
T_MELT = 273.15   # freezing temp.
AKAP = 0.286      # Kappa (Poisson constant R/Cp)
GAS_CON = 287.0   # Gas constant
RA_1 = 610.78     # Parameter for Magnus-Teten-Formula
H_S = 300.        # stable boundary layer height (m)
H_U = 1000.       # unstable boundary layer height (m)
RIC_RS = 0.39     # Critical Richardson number for stable layer
RIC_RU = 0.28     # Critical Richardson number for unstable layer


# pylint: disable-msg=R0914
# Fortynine is reasonable in this case.
# pylint: disable-msg=R0915
# One hundred and twentythree is reasonable in this case.
# flake8: noqa
def mkthe_main(wdir, file_list, modelname):
    """The main script in the module for computation of aux. variables.

    Arguments:
    - wdir: the working directory path;
    - file_list: the list of file containing ts, hus,
    ps, uas, vas, hfss, te;
    - modelname: the name of the model from which the fields are;
    """
    hfss, huss, p_s, t_e, t_s, vv_hor = input_data(wdir, file_list)
    ricr = RIC_RU
    h_bl = H_U
    ricr = np.where(hfss >= 0.75, ricr, RIC_RS)
    h_bl = np.where(hfss >= 0.75, h_bl, H_S)
    ev_p = huss * p_s / (huss + GAS_CON / RV)      # Water vapour pressure
    td_inv = (1 / T_MELT) - (RV / ALV) * np.log(ev_p / RA_1)  # Dewpoint t.
    t_d = 1 / td_inv
    hlcl = 125. * (t_s - t_d)  # Empirical formula for LCL height
    #  Negative heights are replaced by the height of the stable
    #  boundary layer (lower constraint to the height of the cloud layer)
    hlcl = np.where(hlcl >= 0., hlcl, h_bl)
    cp_d = GAS_CON / AKAP
    ztlcl = t_s - (G_0 / cp_d) * hlcl
    # Compute the pseudo-adiabatic lapse rate to obtain the height of cloud
    # top knowing emission temperature.
    gw_pa = (G_0 / cp_d) * (1 + ((ALV * huss) / (GAS_CON * ztlcl)) /
                            (1 + ((ALV ** 2 * huss * 0.622) /
                                  (cp_d * GAS_CON * ztlcl ** 2))))
    htop = - (t_e - ztlcl) / gw_pa + hlcl
    #  Use potential temperature and critical Richardson number to compute
    #  temperature and height of the boundary layer top
    ths = t_s * (P_0 / p_s) ** AKAP
    thz = ths + 0.03 * ricr * (vv_hor) ** 2 / h_bl
    p_z = p_s * np.exp((- G_0 * h_bl) / (GAS_CON * t_s))  # Barometric eq.
    t_z = thz * (P_0 / p_z) ** (-AKAP)
    outlist = [ztlcl, t_z, htop]
    htop_file, tabl_file, tlcl_file = write_output(wdir, modelname, file_list,
                                                   outlist)
    return htop_file, tabl_file, tlcl_file


def input_data(wdir, file_list):
    """Manipulate input fields and read datasets.

    Arguments:
    - wdir: the working directory path;
    - file_list: the list of file containing ts, hus,
    ps, uas, vas, hfss, te;

    Author:
    Valerio Lembo, University of Hamburg, 2019
    """
    cdo = Cdo()
    ts_miss_file = wdir + '/ts.nc'
    removeif(ts_miss_file)
    cdo.setctomiss('0', input=file_list[0], output=ts_miss_file)
    hus_miss_file = wdir + '/hus.nc'
    removeif(hus_miss_file)
    cdo.setctomiss('0', input=file_list[1], output=hus_miss_file)
    ps_miss_file = wdir + '/ps.nc'
    removeif(ps_miss_file)
    cdo.setctomiss('0', input=file_list[2], output=ps_miss_file)
    vv_missfile = wdir + '/V.nc'
    removeif(vv_missfile)
    vv_file = wdir + '/V_miss.nc'
    removeif(vv_file)
    cdo.sqrt(input='-add -sqr {} -sqr {}'.format(file_list[3],
                                                 file_list[4]),
             options='-b F32', output=vv_file)
    cdo.setctomiss('0', input=vv_file, output=vv_missfile)
    hfss_miss_file = wdir + '/hfss.nc'
    removeif(hfss_miss_file)
    cdo.setctomiss('0', input=file_list[5], output=hfss_miss_file)
    te_miss_file = wdir + '/te.nc'
    removeif(te_miss_file)
    cdo.setctomiss('0', input=file_list[6], output=te_miss_file)
    dataset = Dataset(ts_miss_file)
    t_s = dataset.variables['ts'][:, :, :]
    dataset = Dataset(hus_miss_file)
    hus = dataset.variables['hus'][:, :, :, :]
    lev = dataset.variables['plev'][:]
    dataset = Dataset(ps_miss_file)
    p_s = dataset.variables['ps'][:, :, :]
    dataset = Dataset(vv_missfile)
    vv_hor = dataset.variables['uas'][:, :, :]
    dataset = Dataset(hfss_miss_file)
    hfss = dataset.variables['hfss'][:, :, :]
    dataset = Dataset(te_miss_file)
    t_e = dataset.variables['rlut'][:, :, :]
    huss = hus[:, 0, :, :]
    huss = np.where(lev[0] >= p_s, huss, 0.)
    nlev = len(lev)
    for l_l in range(nlev):
        aux = hus[:, l_l, :, :]
        aux = np.where((p_s >= lev[l_l]), aux, 0.)
        huss = huss + aux
    return hfss, huss, p_s, t_e, t_s, vv_hor


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def write_output(wdir, model, file_list, varlist):
    """Write auxiliary variables to new NC files, write new attributes.

    Arguments:
    - wdir: the work directory where the outputs are stored;
    - model: the name of the model;
    - file_list: the list containing the input fields;
    - varlist: a list containing the variables to be written to NC files, i.e.
      tlcl (the temperature at the LCL), t_z (the temperature at the boundary
      layer top), htop (the height of the boundary layer top); their dimensions
      are as (time, lat, lon);

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    fourc = fourier_coefficients
    dataset = Dataset(file_list[0])
    ztlcl = varlist[0]
    t_z = varlist[1]
    htop = varlist[2]
    tlcl_temp = wdir + '/tlcl.nc'
    removeif(tlcl_temp)
    w_nc_fid = Dataset(tlcl_temp, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean LCL temperature from {} model. \
                            Calculated by Thermodynamics model diagnostics\
                            in ESMValTool. Author Valerio Lembo, \
                            Meteorologisches Institut, Universitaet \
                            Hamburg.".format(model)
    fourc.extr_time(dataset, w_nc_fid)
    fourc.extr_lat(dataset, w_nc_fid)
    fourc.extr_lon(dataset, w_nc_fid)
    w_nc_var = w_nc_fid.createVariable('tlcl', 'f8',
                                       ('time', 'lat', 'lon'))
    w_nc_var.setncatts({'long_name': u"LCL Temperature",
                        'units': u"K", 'level_desc': u"surface",
                        'var_desc': u"LCL temperature from LCL \
                        height (Magnus formulas and dry adiabatic \
                        lapse ratio", 'statistic': 'monthly mean'})
    w_nc_fid.variables['tlcl'][:] = ztlcl
    w_nc_fid.close()  # close the new file
    tabl_temp = wdir + '/tabl.nc'
    removeif(tabl_temp)
    w_nc_fid = Dataset(tabl_temp, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean temperature at BL top for {} \
                            model. Calculated by Thermodynamics model \
                            diagnostics in ESMValTool. Author Valerio \
                            Lembo, Meteorologisches Institut, \
                            Universitaet Hamburg.".format(model)
    fourc.extr_time(dataset, w_nc_fid)
    fourc.extr_lat(dataset, w_nc_fid)
    fourc.extr_lon(dataset, w_nc_fid)
    w_nc_var = w_nc_fid.createVariable('tabl', 'f8',
                                       ('time', 'lat', 'lon'))
    w_nc_var.setncatts({'long_name': u"Temperature at BL top",
                        'units': u"K", 'level_desc': u"surface",
                        'var_desc': u"Temperature at the Boundary Layer \
                        top, from boundary layer thickness and barometric \
                        equation", 'statistic': u'monthly mean'})
    w_nc_fid.variables['tabl'][:] = t_z
    w_nc_fid.close()  # close the new file
    htop_temp = wdir + '/htop.nc'
    removeif(htop_temp)
    w_nc_fid = Dataset(htop_temp, 'w', format='NETCDF4')
    w_nc_fid.description = "Monthly mean height of the BL top for {} \
                            model. Calculated by Thermodynamics model \
                            diagnostics in ESMValTool. Author Valerio \
                            Lembo, Meteorologisches Institut, \
                            Universitaet Hamburg.".format(model)
    fourc.extr_time(dataset, w_nc_fid)
    fourc.extr_lat(dataset, w_nc_fid)
    fourc.extr_lon(dataset, w_nc_fid)
    w_nc_var = w_nc_fid.createVariable('htop', 'f8',
                                       ('time', 'lat', 'lon'))
    w_nc_var.setncatts({'long_name': u"Height at BL top",
                        'units': u"m", 'level_desc': u"surface",
                        'var_desc': u"Height at the Boundary Layer top, \
                        from boundary layer thickness and barometric \
                        equation", 'statistic': u'monthly mean'})
    w_nc_fid.variables['htop'][:] = htop
    w_nc_fid.close()  # close the new file
    tlcl_file = wdir + '/{}_tlcl.nc'.format(model)
    cdo.setrtomiss('400,1e36', input=tlcl_temp, output=tlcl_file)
    tabl_temp = wdir + '/tabl.nc'
    tabl_file = wdir + '/{}_tabl.nc'.format(model)
    cdo.setrtomiss('400,1e36', input=tabl_temp, output=tabl_file)
    htop_temp = wdir + '/htop.nc'
    htop_file = wdir + '/{}_htop.nc'.format(model)
    cdo.setrtomiss('12000,1e36', input=htop_temp, output=htop_file)
    return htop_file, tabl_file, tlcl_file
