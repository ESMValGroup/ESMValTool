import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import sys
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib as mpl

#plt.style.use('classic')

def zmnam_plot(datafolder,figfolder,src_props):

    print('plot: ',src_props)
    # TODO use src props for naming plots

    figs_path = figfolder
    fig_num = 1

    in_folder = datafolder

    # See below
    #lev_list = [100000,25000,5000] #,10] 
    #lev_list = [1000,300,30,3]  # BT 09
    #lev_list = [100000,25000,5000]  # test

    # Open daily and monthly PCs
    #file_name = 'daily_pc.nc'
    file_name = '_'.join(src_props)+'_pc_da.nc'
    print(in_folder+file_name)
    in_file = nc4.Dataset(in_folder+file_name,"r")
    #
    time_da = in_file.variables['time'][:]
    time_da_uni = in_file.variables['time'].units
    time_da_cal = in_file.variables['time'].calendar
    #
    lev     = np.array(in_file.variables['plev'][:],dtype='d') 
    lev_units = in_file.variables['plev'].units
    #
    pc_da   = np.array(in_file.variables['PC_da'][:],dtype='d')
    in_file.close()

    #file_name = 'monthly_pc.nc'
    file_name = '_'.join(src_props)+'_pc_mo.nc'
    print(in_folder+file_name)
    in_file = nc4.Dataset(in_folder+file_name,"r")
    #
    time_mo = np.array(in_file.variables['time'][:],dtype='d')
    time_mo_uni = in_file.variables['time'].units 
    time_mo_cal =in_file.variables['time'].calendar
    #
    pc_mo = np.array(in_file.variables['PC_mo'][:],dtype='d')
    in_file.close()

    # Open monthly gh field
    file_name = 'tmp_gh_mo_an_hem.nc'
    print(in_folder+file_name)
    in_file = nc4.Dataset(in_folder+file_name,"r")
    dims =  list(in_file.dimensions.keys())[::-1] # py3
    print('mo full dims', dims)

    if 'latitude' in dims: latn = 'latitude'
    if 'lat' in dims: latn = 'lat'
    if 'longitude' in dims : lonn = 'longitude'
    if 'lon' in dims: lonn = 'lon' 
    lat = np.array(in_file.variables[latn][:])
    lon = np.array(in_file.variables[lonn][:])

    zg_mo = np.array(in_file.variables['zg'][:])
    in_file.close()

    # Open regression files


    date_list = []

    
    if lev_units == 'Pa':
        lev_list = [100000,25000,5000]
        lev_fac = 100. 

    if lev_units == 'hPa' or lev_units == 'hpa' or lev_units == 'mbar' or lev_units == 'millibar':
        lev_list = [1000,250,50] 
        lev_fac = 1.

    for i_date in np.arange(len(time_mo)):
        yy = nc4.num2date(time_mo,time_mo_uni,time_mo_cal)[i_date].year
        mm = nc4.num2date(time_mo,time_mo_uni,time_mo_cal)[i_date].month
        date_list.append(str(yy)+'-'+str(mm))

    for i_lev in np.arange(len(lev)):
        if lev[i_lev] in lev_list:
     
            # Plot monthly PCs
            plt.figure()
            plt.plot(time_mo,pc_mo[:,i_lev])
            plt.xticks(time_mo[0:len(time_mo)+1:60],date_list[0:len(time_mo)+1:60])
            plt.title(str(int(lev[i_lev]/lev_fac))+' hPa  '+\
            src_props[1]+' '+src_props[2])
            plt.xlabel('Time')
            plt.ylabel('Zonal mean NAM')
            plt.savefig(figfolder+'_'.join(src_props)+'_'+\
            str(int(lev[i_lev]/lev_fac))+'hPa_mo_ts.png',format='png')

            plt.figure()
            # PDF of the daily PC
            plt.figure()
            min_var = -5 ; max_var = 5; n_bars=50

            """
            n, bins, patches = plt.hist(pc_da[:,i_lev], bins=n_bars,range=(min_var,max_var),\
                           histtype='bar',normed=True,alpha=.5,color='white')
            """
            n, bins, patches = plt.hist(pc_da[:,i_lev], n_bars, density=True, \
            range=(min_var,max_var), facecolor='b', alpha=0.75)


            # Reference normal Gaussian 
            mu = 0.
            sigma = 1.
            plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * \
                 np.exp( - (bins - mu)**2 / (2 * sigma**2) ),\
                 linewidth=2, color='k',linestyle='--')#,alpha=0.5)

            plt.xlim(min_var,max_var)
            plt.title('Daily PDF ' + str(int(lev[i_lev]/lev_fac))+\
            ' hPa  '+src_props[1]+' '+src_props[2])
            plt.xlabel('Zonal mean NAM')
            plt.ylabel('Normalized probability')  
            plt.tight_layout()
            plt.savefig(figfolder+'_'.join(src_props)+'_'+\
            str(int(lev[i_lev]/lev_fac))+'hPa_da_pdf.png',format='png')



            # Regression of monthly PC onto gh field
            slope = np.zeros((len(lat),len(lon)),dtype='d')
        
            for j_lat in np.arange(len(lat)):
                
                for k_lon in np.arange(len(lon)):
                
                    # Following BT09, the maps are Z_m^l*PC_m^l/|PC_m^l|^2
                    slope[j_lat,k_lon] = np.dot(zg_mo[:,i_lev,j_lat,k_lon],pc_mo[:,i_lev])/\
                    np.dot(pc_mo[:,i_lev],pc_mo[:,i_lev])
         

            plt.figure()
            """
            regr_levs = 0.
            if lev[i_lev] == 100000 : regr_levs = -1005 + np.arange(403)*5 
            if lev[i_lev] == 30000  : regr_levs = -1065 + np.arange(143)*15
            if lev[i_lev] == 3000   : regr_levs = -1000 + np.arange(41)*50 
            if lev[i_lev] == 300    : regr_levs = -1000 + np.arange(41)*50 
            """
            regr_levs = -1000 + np.arange(201)*10 

            bkg_map = Basemap(projection='npstere', boundinglat=20, lon_0=0, \
                        lat_0=90, resolution='c',round=True)
                  
  
            # add wrap-around point in longitude.
            slopew, lonw = addcyclic(slope, lon)

            lons, lats = np.meshgrid(lonw,lat)
            x,y = bkg_map(lons,lats)

            bkg_plot = plt.contourf(x,y,slopew,colors=('#cccccc','#ffffff'),levels=[-10000,0,10000],latlon=True)

            # Switch temporarily to solid negative lines
            mpl.rcParams['contour.negative_linestyle'] = 'solid'
            #regr_map = bkg_map.contour(x,y,slopew,levels=regr_levs,colors='k')#,latlon=True)
            regr_map = plt.contour(x,y,slopew,levels=regr_levs,colors='k',latlon=True,zorder=10)
            #plt.clabel(regr_map,regr_levs[1::2],fontsize=8,fmt='%1.0f') 
            plt.clabel(regr_map,fontsize=8,fmt='%1.0f') 
            mpl.rcParams['contour.negative_linestyle'] = 'dashed'
            
            # coastlines and title
            bkg_map.drawcoastlines(linewidth=0.5,zorder=10)
            bkg_map.fillcontinents(color='None',lake_color='None')

            # No border for the map
            plt.axis('off')
            
            # Write current level - trouble above 100 Pa ...
            plt.text(0.20, 0.80, str(int(lev[i_lev]/lev_fac))+ ' hPa', \
            fontsize=12, transform=plt.gcf().transFigure)
            plt.text(0.75, 0.80, src_props[1],\
            fontsize=12, transform=plt.gcf().transFigure)
            plt.text(0.75, 0.75, src_props[2],\
            fontsize=12, transform=plt.gcf().transFigure)
            #plt.savefig(figfolder+'test'+str(int(lev[i_lev]/lev_fac))+'.png',format='png')
            plt.savefig(figfolder+'_'.join(src_props)+'_'+\
            str(int(lev[i_lev]/lev_fac))+'hPa_mo_reg.png',format='png')
        else: continue

    return 

