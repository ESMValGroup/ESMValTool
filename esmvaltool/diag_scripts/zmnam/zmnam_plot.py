"""
Zonal-mean annular mode plot routine.

Author: Federico Serva (ISAC-CNR & ISMAR-CNR, Italy)
Copernicus C3S 34a lot 2 (MAGIC)

"""

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from cartopy.util import add_cyclic_point


def zmnam_plot(file_gh_mo, datafolder, figfolder, src_props, fig_fmt,
               write_plots):
    """Plotting of timeseries and maps for zmnam diagnostics."""
    plot_files = []
    # Open daily and monthly PCs
    file_name = '_'.join(src_props) + '_pc_da.nc'
    # print(datafolder + file_name)
    with netCDF4.Dataset(datafolder + file_name, "r") as in_file:
        lev = np.array(in_file.variables['plev'][:], dtype='d')
        pc_da = np.array(in_file.variables['PC_da'][:], dtype='d')

    file_name = '_'.join(src_props) + '_pc_mo.nc'
    # print(datafolder + file_name)
    with netCDF4.Dataset(datafolder + file_name, "r") as in_file:
        time_mo = np.array(in_file.variables['time'][:], dtype='d')
        time_mo_uni = in_file.variables['time'].units
        time_mo_cal = in_file.variables['time'].calendar

        pc_mo = np.array(in_file.variables['PC_mo'][:], dtype='d')

    # Open monthly gh field
    file_name = file_gh_mo
    # print(datafolder + file_name)
    with netCDF4.Dataset(file_name, "r") as in_file:
        dims = list(in_file.dimensions.keys())[::-1]  # py3
        print('mo full dims', dims)

        # Double check on lat/lon names, possibly redundant
        if 'latitude' in dims:
            latn = 'latitude'
        if 'lat' in dims:
            latn = 'lat'
        if 'longitude' in dims:
            lonn = 'longitude'
        if 'lon' in dims:
            lonn = 'lon'
        lat = np.array(in_file.variables[latn][:])
        lon = np.array(in_file.variables[lonn][:])

        zg_mo = np.array(in_file.variables['zg'][:])

        # Record attributes for output netCDFs
        time_lnam = getattr(in_file.variables['time'], 'long_name', '')
        time_snam = getattr(in_file.variables['time'], 'standard_name', '')
        time_uni = in_file.variables['time'].units
        time_cal = in_file.variables['time'].calendar
        lev_lnam = getattr(in_file.variables['plev'], 'long_name', '')
        lev_snam = getattr(in_file.variables['plev'], 'standard_name', '')
        lev_uni = in_file.variables['plev'].units
        lev_pos = in_file.variables['plev'].positive
        lev_axi = in_file.variables['plev'].axis

        lat_uni = in_file.variables[latn].units
        lat_axi = in_file.variables[latn].axis

        lon_uni = in_file.variables[lonn].units
        lon_axi = in_file.variables[lonn].axis

    # Save dates for timeseries
    date_list = []
    for i_date in np.arange(len(time_mo)):
        yydate = netCDF4.num2date(time_mo, time_mo_uni,
                                  time_mo_cal)[i_date].year
        mmdate = netCDF4.num2date(time_mo, time_mo_uni,
                                  time_mo_cal)[i_date].month
        date_list.append(str(yydate) + '-' + str(mmdate))

    # Prepare array for outputting regression maps (lev/lat/lon)
    regr_arr = np.zeros((len(lev), len(lat), len(lon)), dtype='f')

    for i_lev in np.arange(len(lev)):

        # Plot monthly PCs
        plt.figure()
        plt.plot(time_mo, pc_mo[:, i_lev])

        # Make only a few ticks
        plt.xticks(time_mo[0:len(time_mo) + 1:60],
                   date_list[0:len(time_mo) + 1:60])
        plt.title(
            str(int(lev[i_lev])) + ' Pa  ' + src_props[1] + ' ' + src_props[2])
        plt.xlabel('Time')
        plt.ylabel('Zonal mean NAM')

        if write_plots:
            fname = (figfolder + '_'.join(src_props) + '_' +
                     str(int(lev[i_lev])) + 'Pa_mo_ts.' + fig_fmt)
            plt.savefig(fname, format=fig_fmt)
            plot_files.append(fname)

        plt.figure()

        # PDF of the daily PC
        plt.figure()
        min_var = -5
        max_var = 5
        n_bars = 50

        _, bins, _ = plt.hist(pc_da[:, i_lev],
                              n_bars,
                              density=True,
                              range=(min_var, max_var),
                              facecolor='b',
                              alpha=0.75)

        # Reference normal Gaussian
        plt.plot(bins,
                 1. / (np.sqrt(2 * np.pi)) * np.exp(-bins**2 / 2.),
                 linewidth=2,
                 color='k',
                 linestyle='--')

        plt.xlim(min_var, max_var)
        plt.title('Daily PDF ' + str(int(lev[i_lev])) + ' Pa  ' +
                  src_props[1] + ' ' + src_props[2])
        plt.xlabel('Zonal mean NAM')
        plt.ylabel('Normalized probability')
        plt.tight_layout()

        if write_plots:
            fname = (figfolder + '_'.join(src_props) + '_' +
                     str(int(lev[i_lev])) + 'Pa_da_pdf.' + fig_fmt)
            plt.savefig(fname, format=fig_fmt)
            plot_files.append(fname)

        plt.close('all')

        # Regression of 3D zg field onto monthly PC
        slope = np.zeros((len(lat), len(lon)), dtype='d')

        for j_lat in np.arange(len(lat)):

            for k_lon in np.arange(len(lon)):

                # Following BT09, the maps are Z_m^l*PC_m^l/|PC_m^l|^2
                slope[j_lat,
                      k_lon] = np.dot(zg_mo[:, i_lev, j_lat,
                                            k_lon], (pc_mo[:, i_lev]) /
                                      np.dot(pc_mo[:, i_lev], pc_mo[:, i_lev]))

        # Plots of regression maps
        plt.figure()

        # Fixed contour levels. May be improved somehow.
        regr_levs = -1000 + np.arange(201) * 10

        # Create the projections
        ortho = ccrs.Orthographic(central_longitude=0, central_latitude=90)
        ccrs.Geodetic()

        # Create the geoaxes for an orthographic projection
        axis = plt.axes(projection=ortho)

        # Add wrap-around point in longitude.
        slopew, lonw = add_cyclic_point(slope, lon)

        # lons, lats = np.meshgrid(lonw, lat)

        plt.contourf(lonw,
                     lat,
                     slopew,
                     colors=('#cccccc', '#ffffff'),
                     levels=[-10000, 0, 10000],
                     transform=ccrs.PlateCarree())

        # Switch temporarily to solid negative lines
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contour(lonw,
                    lat,
                    slopew,
                    levels=regr_levs,
                    colors='k',
                    transform=ccrs.PlateCarree(),
                    zorder=1)

        # Invisible contours, only for labels.
        # Change zorder for cartopy/matplotlib label issue, as of June 2019
        inv_map = plt.contour(lonw,
                              lat,
                              slopew,
                              levels=regr_levs,
                              colors='k',
                              transform=ccrs.PlateCarree(),
                              zorder=15)

        mpl.rcParams['contour.negative_linestyle'] = 'dashed'

        for cmap in inv_map.collections:
            cmap.set_visible(False)

        # Add contour labels over white boxes
        kwargs = {'fontsize': 8, 'fmt': '%1.0f'}
        if mpl.__version__.split('.') >= ['3', '3']:
            kwargs['zorder'] = 30  # new in matplotlib version 3.3
        plt.clabel(inv_map, **kwargs)
        # work around https://github.com/SciTools/cartopy/issues/1554
        # in cartopy 0.18
        clabs = inv_map.labelTextsList
        bbox_dict = dict(boxstyle='square,pad=0',
                         edgecolor='none',
                         fc='white',
                         zorder=25)
        clabs = [txt.set_bbox(bbox_dict) for txt in clabs]

        axis.coastlines()
        axis.set_global()

        plt.text(0.20,
                 0.80,
                 str(int(lev[i_lev])) + ' Pa',
                 fontsize=12,
                 transform=plt.gcf().transFigure)
        plt.text(0.75,
                 0.80,
                 src_props[1],
                 fontsize=12,
                 transform=plt.gcf().transFigure)
        plt.text(0.75,
                 0.75,
                 src_props[2],
                 fontsize=12,
                 transform=plt.gcf().transFigure)

        if write_plots:
            fname = (figfolder + '_'.join(src_props) + '_' +
                     str(int(lev[i_lev])) + 'Pa_mo_reg.' + fig_fmt)
            plt.savefig(fname, format=fig_fmt)
            plot_files.append(fname)

        plt.close('all')

        # Save regression results in array
        regr_arr[i_lev, :, :] = slope

    # Save 3D regression results in output netCDF
    with netCDF4.Dataset(datafolder + '_'.join(src_props) + '_regr_map.nc',
                         mode='w') as file_out:
        file_out.title = 'Zonal mean annular mode (4)'
        file_out.contact = 'F. Serva (federico.serva@artov.ismar.cnr.it); \
        C. Cagnazzo (chiara.cagnazzo@cnr.it)'

        #
        file_out.createDimension('time', None)
        file_out.createDimension('plev', np.size(lev))
        file_out.createDimension('lat', np.size(lat))
        file_out.createDimension('lon', np.size(lon))
        #
        time_var = file_out.createVariable('time', 'd', ('time', ))
        if time_lnam:
            time_var.setncattr('long_name', time_lnam)
        if time_snam:
            time_var.setncattr('standard_name', time_snam)
        time_var.setncattr('units', time_uni)
        time_var.setncattr('calendar', time_cal)
        time_var[:] = 0  # singleton
        #
        lev_var = file_out.createVariable('plev', 'd', ('plev', ))
        if lev_lnam:
            lev_var.setncattr('long_name', lev_lnam)
        if lev_snam:
            lev_var.setncattr('standard_name', lev_snam)
        lev_var.setncattr('units', lev_uni)
        lev_var.setncattr('positive', lev_pos)
        lev_var.setncattr('axis', lev_axi)
        lev_var[:] = lev[:]
        #
        lat_var = file_out.createVariable('lat', 'd', ('lat', ))
        lat_var.setncattr('units', lat_uni)
        lev_var.setncattr('axis', lat_axi)
        lat_var[:] = lat[:]
        #
        lon_var = file_out.createVariable('lon', 'd', ('lon', ))
        lon_var.setncattr('units', lon_uni)
        lon_var.setncattr('axis', lon_axi)
        lon_var[:] = lon[:]
        #
        regr_var = file_out.createVariable('regr', 'f', ('plev', 'lat', 'lon'))
        regr_var.setncattr('long_name',
                           'Zonal mean annular mode regression map')
        regr_var.setncattr(
            'comment',
            'Reference: Baldwin and Thompson ' + '(2009), doi:10.1002/qj.479')
        regr_var[:] = regr_arr[:, :, :]

    return plot_files
