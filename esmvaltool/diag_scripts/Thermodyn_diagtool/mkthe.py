"""Module for computation of some auxiliary variables.

Module needed by the main thermodynamic diagnosticl tool script for
computation of some auxiliary variables.

It computes equivalent potential temperatures and temperatures representative
of the sensible and latent heat exchanges in the lower layers of the
troposphere. Estimates of the boundary layer height and lifting condensation
level are also provided.

Authors: Frank Lunkeit and Valerio Lembo (University of Hamburg)

Created on Fri Jun 15 10:06:30 2018
"""
import os
from netCDF4 import Dataset
import numpy as np
from cdo import Cdo

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


class Mkthe():
    """The auxiliary variables module.

    A class to compute LCL temperature, boundary layer top temperature,
    boundary layer thickness.

    It ingests monthly mean fields of:
    - specific humidity (near-surface or 3D) (hus);
    - skin temperature (ts);
    - surface pressure (ps);
    - near-surface horizontal velocity (uas and vas);
    - surface turbulent sensible heat fluxes (hfss);
    - emission temperature (te).
    """

    from mkthe import Mkthe

    def mkthe_main(self, wdir, ts_file, hus_file, ps_file, uas_file,
                   vas_file, hfss_file, te_file, modelname):
        """The main script in the module for computation of aux. variables.

        Arguments:
        - wdir: the working directory path;
        - ts_file: the path to the NetCDF containing ts;
        - hus_file: the path to the NetCDF containing hus;
        - ps_file: the path to the NetCDF containing ps;
        - uas_file: the path to the NetCDF containing uas;
        - vas_file: the path to the NetCDF containing vas;
        - hfss_file: the path to the NetCDF containing hfss;
        - te_file: the path to the NetCDF containing te;
        - modelname: the name of the model from which the fields are;
        """
        cdo = Cdo()
        mkthe = Mkthe()
        ts_miss_file = wdir + '/ts.nc'
        mkthe.removeif(ts_miss_file)
        cdo.setctomiss('0', input=ts_file, output=ts_miss_file)
        hus_miss_file = wdir + '/hus.nc'
        mkthe.removeif(hus_miss_file)
        cdo.setctomiss('0', input=hus_file, output=hus_miss_file)
        ps_miss_file = wdir + '/ps.nc'
        mkthe.removeif(ps_miss_file)
        cdo.setctomiss('0', input=ps_file, output=ps_miss_file)
        vv_missfile = wdir + '/V.nc'
        mkthe.removeif(vv_missfile)
        vv_file = wdir + '/{}_V.nc'.format(modelname)
        mkthe.removeif(vv_file)
        cdo.sqrt(input='-add -sqr {} -sqr {}'.format(uas_file, vas_file),
                 options='-b F32', output=vv_file)
        cdo.setctomiss('0', input=vv_file, output=vv_missfile)
        hfss_miss_file = wdir+'/hfss.nc'
        mkthe.removeif(hfss_miss_file)
        cdo.setctomiss('0', input=hfss_file, output=hfss_miss_file)
        te_miss_file = wdir + '/te.nc'
        mkthe.removeif(te_miss_file)
        cdo.setctomiss('0', input=te_file, output=te_miss_file)
        dataset0 = Dataset(ts_miss_file)
        t_s = dataset0.variables['ts'][:, :, :]
        lats = dataset0.variables['lat'][:]
        lons = dataset0.variables['lon'][:]
        time = dataset0.variables['time'][:]
        dataset = Dataset(hus_miss_file)
        hus = dataset.variables['hus'][:, :, :, :]
        lev = dataset.variables['plev'][:]
        nlev = len(lev)
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
        for l_l in range(nlev):
            aux = hus[:, l_l, :, :]
            aux = np.where((p_s >= lev[l_l]), aux, 0.)
            huss = huss + aux
        ricr = RIC_RU
        h_bl = H_U
        ricr = np.where(hfss >= 0.75, ricr, RIC_RS)
        h_bl = np.where(hfss >= 0.75, h_bl, H_S)
        ev_p = huss * p_s / (huss + GAS_CON / RV)        # Water vapour pressure
        td_inv = (1 / T_MELT) - (RV / ALV) * np.log(ev_p / RA_1)  # Dewpoint t.
        t_d = 1 / td_inv
        hlcl = 125. * (t_s - t_d)  # Empirical formula for LCL height
    #
    #  Negative heights are replaced by the height of the stable
    #  boundary layer (lower constraint to the height of the cloud layer)
    #
        hlcl = np.where(hlcl >= 0., hlcl, h_bl)
        cp_d = GAS_CON / AKAP
        ztlcl = t_s - (G_0 / cp_d) * hlcl
    #
    # Compute the pseudo-adiabatic lapse rate to obtain the height of cloud
    # top knowing emission temperature.
    #
        gw_pa = (G_0 / cp_d) * (1 + ((ALV * huss) / (GAS_CON * ztlcl)) /
                                (1 + ((ALV ** 2 * huss * 0.622) /
                                      (cp_d * GAS_CON * ztlcl ** 2))))
        htop = - (t_e - ztlcl) / gw_pa + hlcl
    #
    #  Compute equivalent potential temperature (optional output)
    #
    #   thes=ths*np.exp((Alv*huss)/(cp*ztlcl))
    #
    #  Use potential temperature and critical Richardson number to compute
    #  temperature and height of the boundary layer top
    #
        ths = t_s * (P_0 / p_s) ** AKAP
        thz = ths + 0.03 * ricr * (vv_hor) ** 2 / h_bl
        p_z = p_s * np.exp((- G_0 * h_bl) / (GAS_CON * t_s))  # Barometric eq.
        t_z = thz * (P_0 / p_z) ** (-AKAP)
        nc_attrs, nc_dims = mkthe.ncdump(dataset0, 'ts', True)
        nc_f = wdir + '/tlcl.nc'
        mkthe.removeif(nc_f)
        w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "Monthly mean LCL temperature from {} model. \
                                Calculated by Thermodynamics model diagnostics\
                                in ESMValTool. Author Valerio Lembo, \
                                Meteorologisches Institut, Universität \
                                Hamburg.".format(modelname)
        w_nc_fid.createDimension('time', None)
        w_nc_dim = w_nc_fid.createVariable('time',
                                           dataset0.variables['time'].dtype,
                                           ('time',))
        for ncattr in dataset0.variables['time'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['time'].getncattr(ncattr))
        # Assign the dimension data to the new NetCDF file.
        w_nc_fid.variables['time'][:] = time
        w_nc_fid.createDimension('lat', len(lats))
        w_nc_dim = w_nc_fid.createVariable('lat',
                                           dataset0.variables['lat'].dtype,
                                           ('lat',))
        for ncattr in dataset0.variables['lat'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lat'].getncattr(ncattr))
        w_nc_fid.variables['lat'][:] = lats
        w_nc_fid.createDimension('lon', len(lons))
        w_nc_dim = w_nc_fid.createVariable('lon',
                                           dataset0.variables['lon'].dtype,
                                           ('lon',))
        for ncattr in dataset0.variables['lon'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lon'].getncattr(ncattr))
        w_nc_fid.variables['lon'][:] = lons
        w_nc_var = w_nc_fid.createVariable('tlcl', 'f8',
                                           ('time', 'lat', 'lon'))
        w_nc_var.setncatts({'long_name': u"LCL Temperature",
                            'units': u"K", 'level_desc': u"surface",
                            'var_desc': u"LCL temperature from LCL \
                            height (Magnus formulas and dry adiabatic \
                            lapse ratio", 'statistic': 'monthly mean'})
        w_nc_fid.variables['tlcl'][:] = ztlcl
        w_nc_fid.close()  # close the new file
        nc_f = wdir + '/tabl.nc'
        mkthe.removeif(nc_f)
        w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "Monthly mean temperature at BL top for {} \
                                model. Calculated by Thermodynamics model \
                                diagnostics in ESMValTool. Author Valerio \
                                Lembo, Meteorologisches Institut, Universität \
                                Hamburg.".format(modelname)
        w_nc_fid.createDimension('time', None)
        w_nc_dim = w_nc_fid.createVariable('time',
                                           dataset0.variables['time'].dtype,
                                           ('time',))
        for ncattr in dataset0.variables['time'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['time'].getncattr(ncattr))
        # Assign the dimension data to the new NetCDF file.
        w_nc_fid.variables['time'][:] = time
        w_nc_fid.createDimension('lat', len(lats))
        w_nc_dim = w_nc_fid.createVariable('lat',
                                           dataset0.variables['lat'].dtype,
                                           ('lat',))
        for ncattr in dataset0.variables['lat'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lat'].getncattr(ncattr))
        w_nc_fid.variables['lat'][:] = lats
        w_nc_fid.createDimension('lon', len(lons))
        w_nc_dim = w_nc_fid.createVariable('lon',
                                           dataset0.variables['lon'].dtype,
                                           ('lon',))
        for ncattr in dataset0.variables['lon'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lon'].getncattr(ncattr))
        w_nc_fid.variables['lon'][:] = lons
        w_nc_var = w_nc_fid.createVariable('tabl', 'f8',
                                           ('time', 'lat', 'lon'))
        w_nc_var.setncatts({'long_name': u"Temperature at BL top",
                            'units': u"K", 'level_desc': u"surface",
                            'var_desc': u"Temperature at the Boundary Layer \
                            top, from boundary layer thickness and barometric \
                            equation", 'statistic': u'monthly mean'})
        w_nc_fid.variables['tabl'][:] = t_z
        w_nc_fid.close()  # close the new file
        nc_f = wdir + '/htop.nc'
        mkthe.removeif(nc_f)
        w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "Monthly mean height of the BL top for {} \
                                model. Calculated by Thermodynamics model \
                                diagnostics in ESMValTool. Author Valerio \
                                Lembo, Meteorologisches Institut, Universität \
                                Hamburg.".format(modelname)
        w_nc_fid.createDimension('time', None)
        w_nc_dim = w_nc_fid.createVariable('time',
                                           dataset0.variables['time'].dtype,
                                           ('time',))
        for ncattr in dataset0.variables['time'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['time'].getncattr(ncattr))
        # Assign the dimension data to the new NetCDF file.
        w_nc_fid.variables['time'][:] = time
        w_nc_fid.createDimension('lat', len(lats))
        w_nc_dim = w_nc_fid.createVariable('lat',
                                           dataset0.variables['lat'].dtype,
                                           ('lat',))
        for ncattr in dataset0.variables['lat'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lat'].getncattr(ncattr))
        w_nc_fid.variables['lat'][:] = lats
        w_nc_fid.createDimension('lon', len(lons))
        w_nc_dim = w_nc_fid.createVariable('lon',
                                           dataset0.variables['lon'].dtype,
                                           ('lon',))
        for ncattr in dataset0.variables['lon'].ncattrs():
            w_nc_dim.setncattr(ncattr,
                               dataset0.variables['lon'].getncattr(ncattr))
        w_nc_fid.variables['lon'][:] = lons
        w_nc_var = w_nc_fid.createVariable('htop', 'f8',
                                           ('time', 'lat', 'lon'))
        w_nc_var.setncatts({'long_name': u"Height at BL top",
                            'units': u"m", 'level_desc': u"surface",
                            'var_desc': u"Height at the Boundary Layer top, \
                            from boundary layer thickness and barometric \
                            equation", 'statistic': u'monthly mean'})
        w_nc_fid.variables['htop'][:] = htop
        w_nc_fid.close()  # close the new file

    def ncdump(self, nc_fid, key):
        """Print the NetCDF file attributes for a given key.

        Arguments:
        - nc_fid: the ID of a NetCDF file containing variable 'key';
        - key: the name of a variable to obtain the attributes from.
        """
        nc_attrs = nc_fid.ncattrs()
        nc_dims = [dim for dim in nc_fid.dimensions] # list of nc dimensions
        return nc_attrs, nc_dims

    def removeif(self, filename):
        """Remove filename if it exists."""
        try:
            os.remove(filename)
        except OSError:
            pass
