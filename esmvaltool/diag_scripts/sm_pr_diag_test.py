"""
;;#############################################################################
;; PRECIPITATION DEPENDANCE ON SOIL MOISTURE DIAGNOSTIC
;; Authors:     Belen Gallego-Elvira (CEH, UK, belgal@nerc.ac.uk)
;;              Chris Taylor (CEH, UK, cmt@ceh.ac.uk)
;;              Luis Garcia-Carreras (University of Leeds, 
;;                                    L.Garcia-Carreras@leeds.ac.uk)
;; EMBRACE project 
;;#############################################################################
;;
;; Description
;;    This script computes and plot the diagnostic "preference for afternoon 
;;    precipitation over soil moisture anomalies" as in Fig.3 of
;;    Taylor et al. 2012, doi:10.1038/nature11377
;;    
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    att1: short description
;;          keep the indentation if more lines are needed
;;    att2: short description
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    att1: short description
;;    att2: short description
;;
;; Required variable_info attributes (variable specific)
;;    att1: short description
;;    att2: short description
;;
;; Optional variable_info attributes (variable specific)
;;    att1: short description
;;    att2: short description
;;
;; Caveats
;;    List possible caveats or limitations of this diagnostic
;;    Features to-be-implemented shall also be mentioned here
;;
;; Modification history
;;    YYYYMMDD-A_xxxx_yy: extended...
;;    YYYYMMDD-A_xxxx_yy: bug-fixed...
;;    YYYYMMDD-A_xxxx_yy: adapted to...
;;    YYYYMMDD-A_xxxx_yy: written.
;;
;;#############################################################################
"""


# Basic Python packages
import os
import subprocess
import glob
import numpy as np
import random as rd


# Common Python packages
import netCDF4 as nc4
import iris
from iris.coord_categorisation import add_month, add_day_of_month, add_year

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm 
import matplotlib.colors as col
from mpl_toolkits.basemap import Basemap


# Fortran routines
import global_rain_sm as grs
import sample_events as se



# ESMValTool defined Python packages
# sys.path.append("./interface_scripts")
# from esmval_lib import ESMValProject
# from auxiliary import info


def main():

    """ 
    ;; Description
    ;;    Main fuction
    ;;    Call all callable fuctions to
    ;;    read CMIP5 data, compute and plot diagnostic 
    """

    file_in = '/localscratch/wllf012/belgal/test_case_data/' 

    fileout = '/localscratch/wllf012/belgal/test_output/' 
    
    years_in = range(1999,2009) 

    # --------------------------------------
    # Get input data to compute diagnostic
    # --------------------------------------

    # Read cmip5 model data
        
    pr, sm, topo, lon, lat, time = read_pr_sm_topo(file_in, years_in) 
    
    # --------------------------------------
    # Get sm monthly climatology 
    # at selected local time: 6:00 am
    # --------------------------------------

    smclim = get_smclim(sm, lon, time)

    # -------------------------------
    # Compute diagnostic per month
    # -------------------------------

    samplefileout = fileout + 'sample_events/'

    for mn in np.arange(1,13):

        # -------------------------------------------------
        # Create montly arrays required by fortran routines
        # -------------------------------------------------

        prbef, smbef, praft, smaft, \
            monthlypr, monthlysm = get_monthly_input(mn, time, lon, lat,  \
                                                     pr, sm, fileout, samplefileout)

        # -----------------------
        # Run fortran routines 
        # -----------------------

        print 'Executing global_rain_sm for month '+str(mn)

        grs.global_rain_sm( np.asfortranarray(monthlypr), 
                                np.asfortranarray(prbef),
                                np.asfortranarray(praft),
                                np.asfortranarray(monthlysm), 
                                np.asfortranarray(smbef),
                                np.asfortranarray(smaft), 
                                np.asfortranarray(smclim[mn-1,:,:]),
                                np.asfortranarray(topo), 
                                np.asfortranarray(lon), 
                                np.asfortranarray(mn),
                                fileout )

        print 'Executing sample_events for month '+str(mn)

        se.sample_events( np.asfortranarray(monthlysm),
                              np.asfortranarray(smbef),
                              np.asfortranarray(smaft),
                              np.asfortranarray(lon),
                              np.asfortranarray(lat),
                              np.asfortranarray(mn),
                              fileout, samplefileout)

    
    # ---------------------------------------------------
    # Compute p_values (as in Fig. 3, Taylor et al 2012)
    # --------------------------------------------------

    print 'Computing diagnostic'

    xs, ys, p_vals = get_p_val(samplefileout)

    # --------------------------------------------------
    # Save diagnostic to netCDF file and plot 
    # --------------------------------------------------
   
    write_nc(fileout, xs, ys, p_vals)
    
    plot_diagnostic(fileout) 

   
"""
;;#############################################################################
;;  CALLABLE FUNCTIONS (alphabetical order)
;;#############################################################################
"""


def coord_change(cubelist):

    """ 
    ;; Arguments
    ;;    cubelist: list
    ;;          list of iris cubes
    ;;
    ;; Return cubelist
    ;;    list of rolled iris cubes
    ;;
    ;; Description
    ;;    Roll iris cubes with longitude 0-360 to -180-180.
    ;;
    """

    if cubelist[0].coord('longitude').points[0] >= 0:

        newcubelist = []

        for cube in cubelist:

            lon = cube.coord('longitude') 
            londim = cube.coord_dims(lon)[0] 
            cube.data = np.roll(cube.data, len(lon.points) // 2, axis=londim)
            newlon = lon.points - 180.
            lon.points = newlon 

            #if lon.bounds is not None:
            #    lon.bounds -= 180.  

            newcubelist.append(cube) 

        return newcubelist 

    else:

        return cubelist 
 

        
def get_monthly_input(mn, time, lon, lat,  \
                      pr, sm, fileout, samplefileout):

    """ 
    ;; Arguments
    ;;    mn: int
    ;;          month, values from 1 to 12
    ;;    time: iris cube coords
    ;;          time info of cube
    ;;    lon: array [lon]
    ;;          longitude
    ;;    lat: array [lat]
    ;;          latitude
    ;;    pr: iris cube [time, lat, lon]
    ;;          3-hourly precipitation time series
    ;;    sm: iris cube [time, lat, lon]
    ;;          3-hourly soil moisture time series
    ;;    fileout: dir
    ;;          output directory
    ;;    samplefileout: dir
    ;;          temporary outpout directory used by fortran routines
    ;;
    ;; Return 
    ;;    prbef: array [year, day time steps (=8), lat, lon]
    ;;          3-hourly precipitation in the previous day
    ;;    smbef: array [year, day time steps (=8), lat, lon]
    ;;          3-hourly soil moisture in the previous day
    ;;    praf: array [year, day time steps (=8), lat, lon]
    ;;          3-hourly precipitation in the following day
    ;;    smaft: array [year, day time steps (=8), lat, lon]
    ;;          3-hourly soil moisture in the following day
    ;;    monthlypr: array [year, days in month * day time steps, lat, lon]
    ;;          3-hourly precipitation in month mn for the whole analysis period 
    ;;    monthlysm: array [year, days in month * day time steps, lat, lon]
    ;;          3-hourly soil moisture in month mn for the whole analysis period  
    ;;
    ;; Description
    ;;    Prepare monthly input data for fortran routines
    ;;
    """

    # --------------------------------------
    # Get grid and calendar info 
    # --------------------------------------

    utimedate = time[0].units

    first_year = utimedate.num2date(time[0].points[0]).year
    last_year = utimedate.num2date(time[0].points[-1]).year
    nyr = last_year - first_year + 1
    years = np.arange(nyr)+first_year

    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    calendar = utimedate.calendar

    nx = len(lon)
    ny = len(lat)

    days_permonth_360 = [30 for i in range(0,12)]
    days_permonth_noleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]    
    days_permonth_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if calendar == '360_day':
        days_per_year = 360
        days_permonth = days_permonth_360
        nts = [8*dpm for dpm in days_permonth_360]
    elif any([calendar == '365_day', calendar == 'noleap']):
        days_per_year = 365
        days_permonth = days_permonth_noleap
        nts = [8*dpm for dpm in days_permonth_noleap]
    elif any([calendar == 'gregorian', calendar == 'standard']): 
        days_permonth = days_permonth_noleap
        nts = [8*dpm for dpm in days_permonth_noleap]
        # Leap days are considered by get_monthly_input()
    else:
        print 'Error: missing calendar info'

    # --------------------------------------
    # Create pr, sm  before and after arrays
    # --------------------------------------
      
    prbef = np.zeros( (nyr, 8, ny, nx), dtype = 'f4' )
    praft = np.zeros( (nyr, 8, ny, nx), dtype = 'f4' )
    
    smbef = np.zeros( (nyr, 8, ny, nx), dtype = 'f4' )
    smaft = np.zeros( (nyr, 8, ny, nx), dtype = 'f4' )

    string = 'mkdir -p ' + str(samplefileout) + '5x5G_mon' + str(mn).zfill(2)
    subprocess.call(string, shell=True)
    print string

    string = 'mkdir -p ' + str(fileout) + 'event_output/mon' + str(mn).zfill(2)
    subprocess.call(string, shell=True)
    print string

    nt = nts[mn-1]

    monthlypr_list = [np.zeros((nt, ny, nx), dtype = 'f4') for y in years]
    monthlysm_list = [np.zeros((nt, ny, nx), dtype = 'f4') for y in years]

    for yr, year in enumerate(years):

        print mn, year

        if all([mn == 2, any([calendar == 'gregorian', calendar == 'standard'])]):

            monthlysm_list[yr][:,:,:] = (
                    sm.extract(iris.Constraint(month=months[mn-1],year=year, dom = range(1,29) )).data)

            monthlypr_list[yr][:,:,:] = (
                    pr.extract(iris.Constraint(month=months[mn-1],year=year, dom = range(1,29) )).data)
        else:

            monthlysm_list[yr][:,:,:] = (
                    sm.extract(iris.Constraint(month=months[mn-1],year=year)).data)

            monthlypr_list[yr][:,:,:] = (
                    pr.extract(iris.Constraint(month=months[mn-1],year=year)).data)

        # last day of previous month

        if all([mn == 1, year == first_year]):

            prbef[yr,:,:,:] = np.zeros( (8, ny, nx) ) - 999.
            smbef[yr,:,:,:] = np.zeros( (8, ny, nx) ) - 999.

        elif (mn == 1):

            prbef[yr,:,:,:] = (
                pr.extract(iris.Constraint(year=year-1,month='Dec',dom=days_permonth[-1])).data)

            smbef[yr,:,:,:] = (
                sm.extract(iris.Constraint(year=year-1,month='Dec',dom=days_permonth[-1])).data)

        else:

            if any([calendar == '360_day', calendar == '365_day', calendar == 'noleap']):

                prbef[yr,:,:,:] = (
                     pr.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth[mn-2])).data)

                smbef[yr,:,:,:] = (
                     sm.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth[mn-2])).data)

            elif any([calendar == 'gregorian', calendar == 'standard']): 

                if isleap(year):

                    prbef[yr,:,:,:] = (
                         pr.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth_leap[mn-2])).data)

                    smbef[yr,:,:,:] = (
                         sm.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth_leap[mn-2])).data)

                else:
                    prbef[yr,:,:,:] = (
                         pr.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth[mn-2])).data)

                    smbef[yr,:,:,:] = (
                         sm.extract(iris.Constraint(year=year,month=months[mn-2],dom=days_permonth[mn-2])).data)

        # first day of following month

        if (mn == 12 and year == last_year):

            praft[yr,:,:,:] = np.zeros( (8, ny, nx) ) - 999.
            smaft[yr,:,:,:] = np.zeros( (8, ny, nx) ) - 999.

        elif (mn == 12):

            praft[yr,:,:,:] = (
                pr.extract(iris.Constraint(year=year+1,month='Jan',dom=1)).data)

            smaft[yr,:,:,:] = (
                sm.extract(iris.Constraint(year=year+1,month='Jan',dom=1)).data) 

        else:

            if any([calendar == '360_day', calendar == '365_day', calendar == 'noleap']):

                praft[yr,:,:,:] = (
                    pr.extract(iris.Constraint(year=year,month=months[mn],dom=1)).data)

                smaft[yr,:,:,:] = (
                    sm.extract(iris.Constraint(year=year,month=months[mn],dom=1)).data)

            elif any([calendar == 'gregorian', calendar == 'standard']): 

                if all([isleap(year),mn == 2]):

                    praft[yr,:,:,:] = (
                        pr.extract(iris.Constraint(year=year,month=months[1],dom=29)).data)

                    smaft[yr,:,:,:] = (
                        sm.extract(iris.Constraint(year=year,month=months[1],dom=29)).data)     
              
                else:
                    praft[yr,:,:,:] = (
                        pr.extract(iris.Constraint(year=year,month=months[mn],dom=1)).data)

                    smaft[yr,:,:,:] = (
                        sm.extract(iris.Constraint(year=year,month=months[mn],dom=1)).data)


    monthlypr = np.vstack(tuple(monthlypr_list))
    monthlysm = np.vstack(tuple(monthlysm_list))
    monthlypr  = monthlypr.reshape(nyr, nt, ny, nx)
    monthlysm = monthlysm.reshape(nyr, nt, ny, nx)

    return prbef, smbef, praft, smaft, monthlypr, monthlysm




def get_p_val(in_dir):

    """ 
    ;; Arguments
    ;;    in_dir: dir
    ;;          directory with intermediary files from fortran routines
    ;;
    ;; Return 
    ;;    xs: array [lon]
    ;;          regridding coordinates
    ;;    ys: array [lat]
    ;;          regridding coordinates
    ;;    p_vals: list
    ;;          p_values of 5x5 deg grid-boxes
    ;;
    ;; Description
    ;;    Computes percentiles (p_values) of "preference for afternoon 
    ;;    precipitation over soil moisture anomalies" as in Fig.3 of
    ;;    Taylor et al. 2012, doi:10.1038/nature11377
    ;;
    """

    # Find gridboxes (xs, ys) with events

    f_ev = glob.glob(in_dir +'*/*.4')
    xss = []
    yss = []
    for f in f_ev:
        xss.append(int(f[-12:-10]))
        yss.append(int(f[-9:-7]))   
    xs = np.unique(xss) 
    ys = np.unique(yss)
    
    p_vals = [] 

    for x in xs: 
        for y in ys:
            
            # Read events

            ev_f = glob.glob(in_dir +'*/i'+str(x).zfill(2) +'j'+str(y).zfill(2)+'_fort.4')

            levents = []

            for ef in ev_f:

                levents.append((np.genfromtxt(ef,dtype=(float), \
                               skip_header=0, delimiter=' ')).tolist())

            if levents != []: 
                events = []
                for l in levents:
                    if isinstance(l, list):
                        for e in l:
                            events.append(e)
                    else:
                        events.append(l)
            else:
               events = levents

            # Minumun number of event to consider the gridbox = 25   
            if len(events) > 24:

                # Read non events

                nev_f = glob.glob(in_dir +'*/i'+str(x).zfill(2) +'j'+ str(y).zfill(2) +'_fort.3')

                lnon_events = []

                for ef in nev_f:

                    lnon_events.append((np.genfromtxt(ef,dtype=(float), \
                                   skip_header=0, delimiter=' ')).tolist())

                if lnon_events != []: 
                    non_events = []                     
                    for l in lnon_events:
                        if isinstance(l, list):
                            for e in l:
                                non_events.append(e)
                        else:
                            non_events.append(l) 
                else:
                    non_events = lnon_events 

                # delta = events_average - non_events_average
                delta = np.mean(events) - np.mean(non_events)

                #-----------------------------
                # Merge events and non events
                # Shuffle 1000 times
                # Get delta of shuffled events
                #-----------------------------

                all_ev_nev = events + non_events

                sdeltas = []

                for i in range(0,1000):

                    rd.shuffle(all_ev_nev)
                    sev = all_ev_nev[0:len(events)]
                    snev = all_ev_nev[len(events)::]
                    sdelta = np.mean(sev) - np.mean(snev)
                    sdeltas.append(sdelta)
                #-----------------------------    
                # Sort (min to max) shuffled deltas
                # Find p_val
                #-----------------------------

                sdeltas.sort()

                closet_val = min(sdeltas, key=lambda x:abs(x-delta))

                p_val = sdeltas.index(closet_val)/1000.

                p_vals.append(p_val)

            else:

                p_vals.append(-999.)  
 
    return xs, ys, p_vals




def get_smclim(sm, lon, time):

    """ 
    ;; Arguments
    ;;    sm: iris cube [time, lat, lon]
    ;;          3-hourly soil moisture time series
    ;;    lon: array [lon]
    ;;          longitude in degrees east
    ;;    time: iris cube coords
    ;;          time info of cube
    ;;
    ;; Return 
    ;;    var_2d : array[month, lat, lon]
    ;;          monthly sm climalogy
    ;;          NOTE: 
    ;;          Month 0 for climatologies = month in which files start
    ;;          ex. start on 199812010300, month 0 = Dec 
    ;; Description
    ;;    Compute monthly soil moisture climatology at local solar time 6:00 am
    ;;
    """

    data_start_time = (time[0].points[0] - int(time[0].points[0]))*24
    utimedate = time[0].units

    local_time = 6.0
    time_loc_sol = local_time + 3 - data_start_time

    steps_perday = 8 # 3h data 
    steps_window = steps_perday * 3
    glob_degrees = 360.0
    hours_perday = 24.0

    #-----------------------------------------------
    # Convert longitudes to the range [-180,180] regardless of grid type.
    # Calculate the fractional indices of [0,23] that correspond to the
    # interpolation points for each longitude.
    #-----------------------------------------------

    lon = np.array([l-360.0 if l >= 180.0 else l for l in lon])
     
    utc_time = time_loc_sol - (lon*hours_perday/glob_degrees)
    kindex = utc_time*steps_perday/hours_perday + steps_perday

    # Calculate the interpolation weights.  NB For each longitude no more than
    # two weights should be non-zero and their sum should be 1.
    weight_lon = np.zeros((steps_window,len(lon)),dtype=lon.dtype)
    for i, (kf, kw) in enumerate(zip(*np.modf(kindex))):
        weight_lon[int(kw),i] = 1.0 - kf
        weight_lon[int(kw)+1,i] = kf

    #-----------------------------------------------
    # Istead of interpolating, get closest value to 6:00    
    # Get the closest UTC file to time_loc_sol (= 6:00)
    #-----------------------------------------------

    weights_lon = np.where((weight_lon>0), 1, weight_lon )
    weights_lon = weights_lon - np.concatenate((np.zeros((1,\
                                               weights_lon.shape[1])),\
                                               weights_lon[:-1]))
    weights_lon = np.where((weights_lon<0), 0, weights_lon)

    #-----------------------------------------------
    # Find longitude at which local day changes
    # Convert weights from local day to UTC day
    #-----------------------------------------------

    lon_day_change = (time_loc_sol/3.0)*45
    indlist = []
    for j,i in enumerate((lon>=lon_day_change)*(lon<=lon_day_change+90)):
         if i:
             indlist.append(j)
    for i in indlist:
        weights_lon[:,i] = np.roll(weights_lon[:,i],8)
    indlist = []
    for j,i in enumerate(lon.tolist()):
        if i in range(-180,180,45):
            indlist.append(j)
    for i in indlist:
        weights_lon[:,i] = np.roll(weights_lon[:,i],-1)

    #-----------------------------------------------
    # get variable values at local time ( = 6:00 am)
    #-----------------------------------------------

    iris_cube_var = sm 
    temp = np.sum(iris_cube_var.data[:, :, :], axis =0)# fill_value = 0
    land_indexes = np.int64(np.where(temp !=0)) # land gridboxes

    weights = np.tile(weights_lon.reshape(weights_lon.shape[0],1,\
                      weights_lon.shape[1]),(iris_cube_var.data.shape[1],1))
    weights_lp =  weights[:, land_indexes[0], land_indexes[1]]

    var_3h = iris_cube_var.data[:, land_indexes[0], land_indexes[1]]

    varday = []
    for i in range(0,len(var_3h)-16,8):
        threedays_3h = var_3h[i:i+24,:]
        varday.append(np.sum(threedays_3h*weights_lp, axis=0))

    varday.insert(0,varday[0])
    varday.append(varday[-1])

    var_array =  np.zeros((len(varday), varday[0].shape[0]))

    for i, l in enumerate(varday):
        var_array[i,:] = l

    #-----------------------------------------------
    # get monthly climatology
    # if calendar is gregorian, leap days are removed
    #-----------------------------------------------

    calendar = utimedate.calendar

    if any([calendar == 'gregorian', calendar == 'standard']) : 
        var_array = removeleap(var_array, utimedate)   
        calendar = '365_day'
    
    if calendar == '360_day':

        days_peryear = 360
        days_month = 30
        var_array = var_array.reshape(var_array.shape[0]/days_peryear,\
                                          12, days_month, var_array.shape[1])
        mvar_array = np.ma.masked_equal(var_array, -999.0)
        var_mean = np.mean(mvar_array, axis = 2)
        var_mean = np.mean(var_mean, axis =  0)

    elif any([calendar == '365_day', calendar == 'noleap']):

        days_peryear = 365
        var_array = var_array.reshape(var_array.shape[0]/days_peryear,\
                                       days_peryear, var_array.shape[1])
        mvar_array = np.ma.masked_equal(var_array, 1e20)
        days_permonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        ind = 0
        var_mean = np.zeros((12, var_array.shape[2]), dtype = 'float32')
        for m, mdays in enumerate(days_permonth):
            temp = np.mean(mvar_array[:, ind:ind+mdays, :], axis = 0)
            var_mean[m] = np.mean(temp, axis = 0)
            ind = mdays+ind

    var_2d = np.zeros((12, iris_cube_var.data.shape[1], iris_cube_var.data.shape[2]), dtype = 'float32')
    var_2d[:,:] = -999.0
    var_2d[:, land_indexes[0],land_indexes[1]] = var_mean

    return var_2d



def plot_diagnostic(fileout):

    """ 
    ;; Arguments
    ;;    fileout: dir
    ;;          directory to save the plot
    ;;
    ;; Description
    ;;    Plot diagnostic and save .png plot
    ;;
    """
 
    # Read code output in netCDF format
    ncf = nc4.Dataset(fileout + 'Taylor2012_diagnostic.nc','r')
    p_val =  ncf.variables['T12_diag'][:,:]
    ncf.close()

    mp_val = np.ma.masked_equal(p_val, -999)

    # Gridding info: Global diagnostic 5x5 deg
    LON = np.arange(-180, 185, 5)
    LAT = np.arange(-60, 65, 5)

    # Define figure

    F, ax = plt.subplots(nrows = 1,ncols = 1,**dict(figsize = [15,8]))

    # Cmap
    cmap  = cm.get_cmap(name = 'bwr_r', lut = 7)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [0., 1, 5, 10, 90, 95, 99, 100]
    norm = col.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad('GainsBoro')

    # Basemap
    map_ax = Basemap(ax = ax, projection = 'cyl', resolution = 'l', \
                                    llcrnrlat = -60, urcrnrlat = 60, \
                                    llcrnrlon = -180, urcrnrlon = 180,\
                                    fix_aspect = False)
    map_ax.drawcoastlines(color='k', linewidth = .7)

    # Plot p_val on basemap
    I = map_ax.pcolormesh(LON, LAT, mp_val, cmap = cmap, norm = norm)

    # Colorbar
    cax = F.add_axes([0.94, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(I, cax = cax, cmap = cmap)

    plt.suptitle('Preference for afternoon precipitation over soil moisture anomalies',
                  fontsize = 14)

    # Save figure to fileout
    plt.savefig(fileout + 'sm_pr_diag_plot')




def read_pr_sm_topo(filedir, years):

    """ 
    ;; Arguments
    ;;    filedir: dir
    ;;          directory with input data
    ;;    years: list of int
    ;;          list of years for the analysis
    ;;
    ;; Return 
    ;;    pr: iris cube [time, lat, lon]
    ;;          precipitation time series
    ;;    sm: iris cube [time, lat, lon]
    ;;          soil moisture time series
    ;;    topo: array [lat, lon]
    ;;          topography
    ;;    lon: array [lon]
    ;;          longitude
    ;;    lat: array [lat]
    ;;          latitude
    ;;    time: iris cube coords
    ;;          time info of cube
    ;;    
    ;;
    ;; Description
    ;;    Read cmip5 input data for computing the diagnostic
    ;;
    """
      
    #-------------------------
    # Input data directories
    #-------------------------
   
    Patt = filedir+'pr_3hr_inmcm4_amip_r1i1p1_{}010101-{}123122.nc'
    pr_files = [Patt.format(y,y) for y in years]
    
    Patt = filedir+'mrsos_3hr_inmcm4_amip_r1i1p1_{}010100-{}123121.nc'
    sm_files = [Patt.format(y,y) for y in years]

    #----------------------
    # Read in precipitation
    #----------------------

    pr_list = []

    for pr_file in pr_files:

        print 'Reading precipitation from '+ pr_file

        pr = iris.load(pr_file)[0]

        for at_k in pr.attributes.keys():
            pr.attributes.pop(at_k)

        pr_list.append(pr)

    pr = iris.cube.CubeList(pr_list)
    pr = pr.concatenate()[0]
    
    # Convert longitude from 0_360 to -180_180

    pr = coord_change([pr])[0]
   
    # Add metadata: day, month, year

    add_month(pr,'time')
    add_day_of_month(pr,'time',name='dom')
    add_year(pr,'time')

    # Convert units to kg m-2 hr-1

    pr.convert_units('kg m-2 hr-1')

    #-----------------------
    # Read in soil moisture
    #-----------------------

    sm_list = []

    for sm_file in sm_files:

        print 'Reading soil moisture from '+ sm_file

        sm = iris.load(sm_file)[0]

        for at_k in sm.attributes.keys():
            sm.attributes.pop(at_k)

        sm_list.append(sm)

    sm = iris.cube.CubeList(sm_list)
    sm = sm.concatenate()[0]

    # Convert longitude from 0_360 to -180_180

    sm = coord_change([sm])[0]

    # Add metadata: day, month, year

    add_month(sm,'time')
    add_day_of_month(sm,'time',name='dom')
    add_year(sm,'time')

    #----------------------------------------------
    # Constrain pr and sm data to latitude 60S_60N
    #----------------------------------------------
    
    latconstraint = iris.Constraint(latitude = lambda cell: -59.0<=cell<=59.0)

    pr = pr.extract(latconstraint)
    sm = sm.extract(latconstraint)

    #---------------------------------------------------
    # Read in grid info: latitude, longitude, timestamp
    #---------------------------------------------------

    lon = sm.coords('longitude')[0].points
    lat = sm.coords('latitude')[0].points
    time = sm.coords('time')

    # --------------------------------------
    # Convert missing data (if any) to -999.
    # --------------------------------------

    try: 
        sm.data.set_fill_value(-999)
        sm.data.data[sm.data.mask]=-999. 

    except:
        print 'no missing data conversion'

    #----------------------
    # Read in topography 
    #----------------------
    
    # Topography map specs:
    # latitude 60S_60N
    # longitude 180W_180E
    # model resolution 

    ftopo = filedir + 'topo_var_5x5_inmcm4.gra'
    dt = '>f4'
    topo = (np.fromfile(ftopo, dtype=dt)).reshape(len(lat), len(lon))

    #-----------------------------------------------
    # Return input data to compute sm_pr diagnostic 
    #-----------------------------------------------
    return pr, sm, topo, lon, lat, time



def write_nc(fileout, xs, ys, p_vals):

    """ 
    ;; Arguments
    ;;    fileout: dir
    ;;          directory to save output
    ;;    xs: array [lon]
    ;;          regridding coordinates
    ;;    ys: array [lat]
    ;;          regridding coordinates
    ;;    p_vals: list
    ;;          p_values of 5x5 deg grid-boxes
    ;;
    ;; Description
    ;;    Save netCDF file with diagnostic in a regular 5x5 deg grid 
    ;;
    """
       
    #------------------------------
    # Transform p_vals in 2D array
    #------------------------------
    p_val_2d = np.zeros((360/5, 120/5), dtype = 'f4')
    p_val_2d[:,:] = -999.
    i = 0    
    for x in xs:
        for y in ys:
            if p_vals[i]>-999.:
                p_val_2d[x-1, y-1] = p_vals[i]*100
            i = i+1

    #------------------------------
    # Write nc file
    #------------------------------

    root_grp = nc4.Dataset(fileout+'Taylor2012_diagnostic.nc', 'w', format='NETCDF4')

    root_grp.description = 'Diagnostic Taylor2012: Precipitation dependance on soil moisture'

    root_grp.fillvalue = -999.0

    root_grp.createDimension('lon',360/5)
    root_grp.createDimension('lat',120/5)

    lat = root_grp.createVariable('latitude', 'f4', ('lat',))
    lon = root_grp.createVariable('longitude', 'f4', ('lon',))
    pval = root_grp.createVariable('T12_diag', 'f4', ('lat','lon'))

    lat[:] = np.arange(-60+2.5, 60, 5)
    lon[:] = np.arange(-180+2.5, 180, 5)
    pval[:,:] = np.transpose(p_val_2d)

    root_grp.close()




if __name__ == "__main__": main()

