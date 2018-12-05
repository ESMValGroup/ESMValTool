import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import sys
from scipy.stats import linregress
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point


def zmnam_plot(datafolder, figfolder, src_props,fig_fmt):

    figs_path = figfolder
    fig_num = 1

    in_folder = datafolder

    # Open daily and monthly PCs
    file_name = '_'.join(src_props) + '_pc_da.nc'
    print(in_folder + file_name)
    in_file = nc4.Dataset(in_folder + file_name, "r")
    #
    time_da = in_file.variables['time'][:]
    time_da_uni = in_file.variables['time'].units
    time_da_cal = in_file.variables['time'].calendar
    #
    lev = np.array(in_file.variables['plev'][:], dtype='d')
    lev_units = in_file.variables['plev'].units
    #
    pc_da = np.array(in_file.variables['PC_da'][:], dtype='d')
    in_file.close()

    file_name = '_'.join(src_props) + '_pc_mo.nc'
    print(in_folder + file_name)
    in_file = nc4.Dataset(in_folder + file_name, "r")
    #
    time_mo = np.array(in_file.variables['time'][:], dtype='d')
    time_mo_uni = in_file.variables['time'].units
    time_mo_cal = in_file.variables['time'].calendar
    #
    pc_mo = np.array(in_file.variables['PC_mo'][:], dtype='d')
    in_file.close()

    # Open monthly gh field
    file_name = 'tmp_gh_mo_an_hem.nc'
    print(in_folder + file_name)
    in_file = nc4.Dataset(in_folder + file_name, "r")
    dims = list(in_file.dimensions.keys())[::-1]  # py3
    print('mo full dims', dims)

    if 'latitude' in dims: latn = 'latitude'
    if 'lat' in dims: latn = 'lat'
    if 'longitude' in dims: lonn = 'longitude'
    if 'lon' in dims: lonn = 'lon'
    lat = np.array(in_file.variables[latn][:])
    lon = np.array(in_file.variables[lonn][:])

    zg_mo = np.array(in_file.variables['zg'][:])
    in_file.close()

    # Open regression files

    date_list = []

    for i_date in np.arange(len(time_mo)):
        yy = nc4.num2date(time_mo, time_mo_uni, time_mo_cal)[i_date].year
        mm = nc4.num2date(time_mo, time_mo_uni, time_mo_cal)[i_date].month
        date_list.append(str(yy) + '-' + str(mm))

    for i_lev in np.arange(len(lev)):

        # Plot monthly PCs
        plt.figure()
        plt.plot(time_mo, pc_mo[:, i_lev])

        # Make only a few ticks
        plt.xticks(time_mo[0:len(time_mo) + 1:60],
                   date_list[0:len(time_mo) + 1:60])
        plt.title(str(int(lev[i_lev]))+' Pa  '+\
        src_props[1]+' '+src_props[2])
        plt.xlabel('Time')
        plt.ylabel('Zonal mean NAM')
        plt.savefig(figfolder+'_'.join(src_props)+'_'+\
        str(int(lev[i_lev]))+'Pa_mo_ts.'+fig_fmt,format=fig_fmt)

        plt.figure()
        
        # PDF of the daily PC
        plt.figure()
        min_var = -5
        max_var = 5
        n_bars = 50

        n, bins, patches = plt.hist(pc_da[:,i_lev], n_bars, density=True, \
        range=(min_var,max_var), facecolor='b', alpha=0.75)

        # Reference normal Gaussian
        mu = 0.
        sigma = 1.
        plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * \
                 np.exp( - (bins - mu)**2 / (2 * sigma**2) ),\
                 linewidth=2, color='k',linestyle='--')

        plt.xlim(min_var, max_var)
        plt.title('Daily PDF ' + str(int(lev[i_lev]))+\
        ' Pa  '+src_props[1]+' '+src_props[2])
        plt.xlabel('Zonal mean NAM')
        plt.ylabel('Normalized probability')
        plt.tight_layout()
        plt.savefig(figfolder+'_'.join(src_props)+'_'+\
        str(int(lev[i_lev]))+'Pa_da_pdf.'+fig_fmt,format=fig_fmt)

        # Regression of 3D zg field onto monthly PC
        slope = np.zeros((len(lat), len(lon)), dtype='d')

        for j_lat in np.arange(len(lat)):

            for k_lon in np.arange(len(lon)):

                # Following BT09, the maps are Z_m^l*PC_m^l/|PC_m^l|^2
                slope[j_lat,k_lon] = np.dot(zg_mo[:,i_lev,j_lat,k_lon],pc_mo[:,i_lev])/\
                np.dot(pc_mo[:,i_lev],pc_mo[:,i_lev])

        plt.figure()

        # Fixed contour levels. May be improved somehow.
        regr_levs = -1000 + np.arange(201) * 10

        # Create the projections
        ortho = ccrs.Orthographic(central_longitude=0, central_latitude=90)
        geo = ccrs.Geodetic()

        # Create the geoaxes for an orthographic projection
        ax = plt.axes(projection=ortho)

        # Add wrap-around point in longitude.
        slopew, lonw = add_cyclic_point(slope, lon)

        lons, lats = np.meshgrid(lonw, lat)

        bkg_plot = plt.contourf(lonw,lat,slopew,colors=('#cccccc','#ffffff'),\
        levels=[-10000,0,10000],transform=ccrs.PlateCarree())

        # Switch temporarily to solid negative lines
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        regr_map = plt.contour(lonw,lat,slopew,levels=regr_levs,\
        colors='k',transform=ccrs.PlateCarree(),zorder=5)

        # Invisible contours, only for labels.
        # Workaround for cartopy issue, as of Dec 18
        inv_map = plt.contour(lonw,lat,slopew,levels=regr_levs,\
        colors='k',transform=ccrs.PlateCarree(),zorder=10)

        mpl.rcParams['contour.negative_linestyle'] = 'dashed'

        for c in inv_map.collections:
            c.set_visible(False)

        plt.clabel(inv_map, fontsize=8, fmt='%1.0f', zorder=15)

        ax.coastlines()
        ax.set_global()

        plt.text(0.20, 0.80, str(int(lev[i_lev]))+ ' Pa', \
        fontsize=12, transform=plt.gcf().transFigure)
        plt.text(0.75, 0.80, src_props[1],\
        fontsize=12, transform=plt.gcf().transFigure)
        plt.text(0.75, 0.75, src_props[2],\
        fontsize=12, transform=plt.gcf().transFigure)
        plt.savefig(figfolder+'_'.join(src_props)+'_'+\
        str(int(lev[i_lev]))+'Pa_mo_reg.'+fig_fmt,format=fig_fmt)

        # TODO Save regression in netCDF for further postprocessing
            # with file name etc


    return
