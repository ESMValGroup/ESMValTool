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
;;    
;;    Long run times for long time series or/and high resolution models
;;    The code assumes that the time series starts in January
;;
;;    Features to-be-implemented shall also be mentioned here
;;    A better performance version, faster and more flexible
;;    is currently being implemented
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
import re
import glob
import errno
import numpy as np
import random as rd
import ConfigParser
import pdb
import stat
import sys
import shutil


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


# ESMValTool defined Python packages
sys.path.append("./interface_scripts")
from esmval_lib import ESMValProject
from auxiliary import info
from auxiliary import error


# Fortran routines
import global_rain_sm as grs
import sample_events as se


def main(project_info):
    """
    ;; Description
    ;;    Main fuction
    ;;    Call all callable fuctions to
    ;;    read CMIP5 data, compute and plot diagnostic
    """

    E = ESMValProject(project_info)
    plot_dir = E.get_plot_dir()
    work_dir = E.get_work_dir()
    verbosity = E.get_verbosity()
    fileout = work_dir
    
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    for model in project_info['MODELS']:
        info(model, verbosity, required_verbosity=1) 

        if not os.path.exists(work_dir+'/sample_events'):
            os.makedirs(work_dir+'/sample_events')

        if not os.path.exists(work_dir+'/event_output'):
            os.makedirs(work_dir+'/event_output')


        # --------------------------------------
        # Get input data to compute diagnostic
        # --------------------------------------

        # Read cmip5 model data
        pr, sm, topo, lon, lat, time, time_bnds_1 = read_pr_sm_topo(project_info, model)

        # --------------------------------------
        # Get sm monthly climatology
        # at selected local time: 6:00 am
        # --------------------------------------

        smclim = get_smclim(sm, lon, time)

        # -------------------------------
        # Compute diagnostic per month
        # -------------------------------

        samplefileout = fileout + 'sample_events/'

        for mn in np.arange(1, 13):

            # -------------------------------------------------
            # Create montly arrays required by fortran routines
            # -------------------------------------------------

            prbef, smbef, praft, smaft, \
                monthlypr, monthlysm, days_per_year = get_monthly_input(project_info, mn, time, 
                                                         lon, lat, time_bnds_1, pr, sm, fileout,
                                                         samplefileout, model,
                                                         verbosity)

            # -----------------------
            # Run fortran routines
            # -----------------------

            info('Executing global_rain_sm for month ' + str(mn), verbosity, required_verbosity=1)

            grs.global_rain_sm(np.asfortranarray(monthlypr),
                               np.asfortranarray(prbef),
                               np.asfortranarray(praft),
                               np.asfortranarray(monthlysm),
                               np.asfortranarray(smbef),
                               np.asfortranarray(smaft),
                               np.asfortranarray(smclim[mn - 1, :, :]),
                               np.asfortranarray(topo),
                               np.asfortranarray(lon),
                               np.asfortranarray(mn),
                               fileout, days_per_year)

            info('Executing sample_events for month ' + str(mn), verbosity, required_verbosity=1)

            se.sample_events(np.asfortranarray(monthlysm),
                             np.asfortranarray(smbef),
                             np.asfortranarray(smaft),
                             np.asfortranarray(lon),
                             np.asfortranarray(lat),
                             np.asfortranarray(mn),
                             fileout, days_per_year, samplefileout)

        # ---------------------------------------------------
        # Compute p_values (as in Fig. 3, Taylor et al 2012)
        # --------------------------------------------------
        info('Computing diagnostic', verbosity, required_verbosity=1)

        xs, ys, p_vals = get_p_val(samplefileout)

        # --------------------------------------------------
        # Save diagnostic to netCDF file and plot
        # --------------------------------------------------
        write_nc(fileout, xs, ys, p_vals, project_info, model)
        
        plot_diagnostic(fileout, plot_dir, project_info, model)

        # --------------------------------------------------
        # Remove temporary folders
        # --------------------------------------------------
        shutil.rmtree(str(fileout) + 'event_output')
        shutil.rmtree(str(fileout) + 'sample_events') 
    


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


def get_monthly_input(project_info, mn, time, lon, lat,
                      time_bnds_1, pr, sm, fileout,
                      samplefileout, model,
                      verbosity):

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
    ;;    time_bnds_1: float
    ;;          first time_bnd of time series
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
    ;;    days_per_year: int
    ;;          number of days per year 
    ;;
    ;; Description
    ;;    Prepare monthly input data for fortran routines
    ;;
    """

    import projects
    E = ESMValProject(project_info)
    verbosity = E.get_verbosity()

    #-------------------------
    # Read model info
    #-------------------------
    
    model = project_info['MODELS'][0]
    currProject = getattr(vars()['projects'], model.split_entries()[0])()

    model_info = model.split_entries()



    # --------------------------------------
    # Get grid and calendar info
    # --------------------------------------

    utimedate = time[0].units

    first_year = int(model_info[5])
    last_year = int(model_info[6])

    nyr = last_year - first_year + 1
    years = np.arange(nyr) + first_year

    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    calendar = utimedate.calendar

    nx = len(lon)
    ny = len(lat)

    days_permonth_360 = [30 for i in range(0, 12)]
    days_permonth_noleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days_permonth_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if calendar == '360_day':
        days_per_year = 360
        days_permonth = days_permonth_360
        nts = [8 * dpm for dpm in days_permonth_360]
    elif any([calendar == '365_day', calendar == 'noleap']):
        days_per_year = 365
        days_permonth = days_permonth_noleap
        nts = [8 * dpm for dpm in days_permonth_noleap]
    elif any([calendar == 'gregorian', calendar == 'standard']):
        days_per_year = 365 
        days_permonth = days_permonth_noleap
        nts = [8 * dpm for dpm in days_permonth_noleap]
        # Leap days are considered by get_monthly_input()
    else:
        error('Missing calendar info')

    # --------------------------------------
    # Create pr, sm  before and after arrays
    # --------------------------------------

    prbef = np.zeros((nyr, 8, ny, nx), dtype='f4')
    praft = np.zeros((nyr, 8, ny, nx), dtype='f4')

    smbef = np.zeros((nyr, 8, ny, nx), dtype='f4')
    smaft = np.zeros((nyr, 8, ny, nx), dtype='f4')

    try:
        os.mkdir(os.path.join(samplefileout, "5x5G_mon" + str(mn).zfill(2)), stat.S_IWUSR | stat.S_IRUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

    try:
        os.mkdir(os.path.join(fileout, "event_output", "mon" + str(mn).zfill(2)), stat.S_IWUSR | stat.S_IRUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

    nt = nts[mn - 1]

    monthlypr_list = [np.zeros((nt, ny, nx), dtype='f4') for y in years]
    monthlysm_list = [np.zeros((nt, ny, nx), dtype='f4') for y in years]

    for yr, year in enumerate(years):

        info('month, year: ' + str(mn) + ", " + str(year), verbosity, required_verbosity=1)

        if all([mn == 2, any([calendar == 'gregorian', calendar == 'standard'])]):

            monthlysm_list[yr][:,:,:] = (
                    sm.extract(iris.Constraint(month=months[mn-1],year=year, dom = range(1,29) )).data)

            monthlypr_list[yr][:,:,:] = (
                    pr.extract(iris.Constraint(month=months[mn-1],year=year, dom = range(1,29) )).data)

        else: 
            try:
                monthlypr_list[yr][:,:,:] = (
                      pr.extract(iris.Constraint(month=months[mn-1],year=year)).data)  
            except:
                try:
                     monthlypr_list[yr][0,:,:] = -999.
                     monthlypr_list[yr][1::,:,:] = (
                          pr.extract(iris.Constraint(month=months[mn-1],year=year)).data)    
                except:                              
                     info('omitted pr: ' + str(mn) + ", " + str(year), verbosity, required_verbosity=1)

            try:
                monthlysm_list[yr][:,:,:] = (
                      sm.extract(iris.Constraint(month=months[mn-1],year=year)).data)  
            except:
                try:
                     monthlysm_list[yr][0,:,:] = -999.
                     monthlysm_list[yr][1::,:,:] = (
                          sm.extract(iris.Constraint(month=months[mn-1],year=year)).data)    
                except:                              
                     info('omitted sm: ' + str(mn) + ", " + str(year), verbosity, required_verbosity=1)


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

                if cal.isleap(year):

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

                if all([cal.isleap(year),mn == 2]):

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
    
    if time_bnds_1 == 0.0625:
        monthlypr[:, 0:-1, :, :]  = monthlypr[:, 1::, :, :]
        monthlypr[:, -1, :, :]  = -9999.
        prbef[:, 0:-1, :, :] = prbef[:, 1::, :, :]
        praft[:, 0:-1, :, :] = praft[:, 1::, :, :]
        prbef[:, -1, :, :] = -9999.
        praft[:, -1, :, :] = -9999.

    return prbef, smbef, praft, smaft, monthlypr, monthlysm, days_per_year




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
    ;; NOTE:
    ;; Change shuffle_times for a quicker run
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
                # Shuffle n times (n=500)
                # Get delta of shuffled events
                #-----------------------------

                all_ev_nev = events + non_events

                sdeltas = []

                shuffle_times = 500.0
                for i in range(0, int(shuffle_times)):

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

                p_val = sdeltas.index(closet_val)/shuffle_times

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

    data_start_time = (time[0].points[0] - int(time[0].points[0])) * 24
    utimedate = time[0].units

    local_time = 6.0
    time_loc_sol = local_time + 3 - data_start_time

    steps_perday = 8  # 3h data
    steps_window = steps_perday * 3
    glob_degrees = 360.0
    hours_perday = 24.0

    #-----------------------------------------------
    # Convert longitudes to the range [-180, 180] regardless of grid type.
    # Calculate the fractional indices of [0, 23] that correspond to the
    # interpolation points for each longitude.
    #-----------------------------------------------

    lon = np.array([l - 360.0 if l >= 180.0 else l for l in lon])

    utc_time = time_loc_sol - (lon * hours_perday / glob_degrees)
    kindex = utc_time * steps_perday / hours_perday + steps_perday

    # Calculate the interpolation weights.  NB For each longitude no more than
    # two weights should be non-zero and their sum should be 1.
    weight_lon = np.zeros((steps_window, len(lon)), dtype=lon.dtype)
    for i, (kf, kw) in enumerate(zip(*np.modf(kindex))):
        weight_lon[int(kw), i] = 1.0 - kf
        weight_lon[int(kw) + 1, i] = kf

    #-----------------------------------------------
    # Istead of interpolating, get closest value to 6:00
    # Get the closest UTC file to time_loc_sol (= 6:00)
    #-----------------------------------------------

    weights_lon = np.where((weight_lon > 0), 1, weight_lon)
    weights_lon = weights_lon - np.concatenate((np.zeros((1,
                                               weights_lon.shape[1])),
                                               weights_lon[:-1]))
    weights_lon = np.where((weights_lon < 0), 0, weights_lon)

    #-----------------------------------------------
    # Find longitude at which local day changes
    # Convert weights from local day to UTC day
    #-----------------------------------------------

    lon_day_change = (time_loc_sol / 3.0) * 45
    indlist = []
    for j, i in enumerate((lon >= lon_day_change) * (lon <= lon_day_change + 90)):
        if i:
            indlist.append(j)
    for i in indlist:
        weights_lon[:, i] = np.roll(weights_lon[:, i], 8)
    indlist = []
    for j, i in enumerate(lon.tolist()):
        if i in range(-180, 180, 45):
            indlist.append(j)
    for i in indlist:
        weights_lon[:, i] = np.roll(weights_lon[:, i], -1)

    #-----------------------------------------------
    # get variable values at local time ( = 6:00 am)
    #-----------------------------------------------

    iris_cube_var = sm
    temp = np.sum(iris_cube_var.data[:, :, :], axis=0)  # fill_value = 0
    land_indexes = np.int64(np.where(temp != 0))  # land gridboxes

    weights = np.tile(weights_lon.reshape(weights_lon.shape[0], 1,
                      weights_lon.shape[1]), (iris_cube_var.data.shape[1], 1))
    weights_lp = weights[:, land_indexes[0], land_indexes[1]]

    var_3h = iris_cube_var.data[:, land_indexes[0], land_indexes[1]]

    varday = []
    for i in range(0, len(var_3h) - 16, 8):
        threedays_3h = var_3h[i:i + 24, :]
        varday.append(np.sum(threedays_3h * weights_lp, axis=0))

    varday.insert(0, varday[0])
    varday.append(varday[-1])

    var_array = np.zeros((len(varday), varday[0].shape[0]))

    for i, l in enumerate(varday):
        var_array[i, :] = l

    #-----------------------------------------------
    # get monthly climatology
    # if calendar is gregorian, leap days are removed
    #-----------------------------------------------

    calendar = utimedate.calendar

    if any([calendar == 'gregorian', calendar == 'standard']):
        var_array = removeleap(var_array, utimedate)
        calendar = '365_day'

    if calendar == '360_day':

        days_peryear = 360
        days_month = 30
        var_array = var_array.reshape(var_array.shape[0] / days_peryear,
                                      12, days_month, var_array.shape[1])
        mvar_array = np.ma.masked_equal(var_array, -999.0)
        var_mean = np.mean(mvar_array, axis=2)
        var_mean = np.mean(var_mean, axis=0)

    elif any([calendar == '365_day', calendar == 'noleap']):

        days_peryear = 365
        var_array = var_array.reshape(var_array.shape[0] / days_peryear,
                                      days_peryear, var_array.shape[1])
        mvar_array = np.ma.masked_equal(var_array, 1e20)
        days_permonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        ind = 0
        var_mean = np.zeros((12, var_array.shape[2]), dtype='float32')
        for m, mdays in enumerate(days_permonth):
            temp = np.mean(mvar_array[:, ind:ind + mdays, :], axis=0)
            var_mean[m] = np.mean(temp, axis=0)
            ind = mdays + ind

    var_2d = np.zeros((12, iris_cube_var.data.shape[1], iris_cube_var.data.shape[2]), dtype='float32')
    var_2d[:, :] = -999.0
    var_2d[:, land_indexes[0], land_indexes[1]] = var_mean

    return var_2d


def get_topo(project_info, longi, lati, model):

    """ 
    ;; Arguments
    ;;    in_dir: dir
    ;;          directory with input file "topo_var_5x5.gra"
    ;;    longi: array [lon]
    ;;          longitude in degrees east
    ;;    lati: array [lat]
    ;;          latitude
    ;; Return 
    ;;    topo: array [lat, lon]
    ;;          topography ranges in model grid
    ;;
    ;; Description
    ;;    Computes topography information
    ;;
    """
    import projects
    E = ESMValProject(project_info)
    verbosity = E.get_verbosity()

    #-------------------------
    # Read model info
    #-------------------------
    
    currProject = getattr(vars()['projects'], model.split_entries()[0])()

    model_info = model.split_entries()


    # load topo max-min topography(m)
    # within 1.25x1.25 area

    nx = 1440
    ny = 480

    ftopo = currProject.get_cf_fx_file(project_info, model)
    
    info('topo file: ' + ftopo , verbosity, required_verbosity=1)

    dt = '>f4'
    topo = (np.fromfile(ftopo, dtype=dt)).reshape(4, ny, nx)
    topo_range = np.transpose(topo[3,:,:])

    # Average onto model grid

    nxo = len(longi)
    nyo = len(lati)
    lon1 = longi[0]
    lat1 = lati[0]
    dlon = np.abs(longi[0]-longi[1])
    dlat = np.abs(lati[0]-lati[1])

    topo_model = np.zeros((nxo, nyo), dtype = 'f4')
    n = np.zeros((nxo, nyo), dtype = 'f4')

    for j in range(0, ny):
        for i in range(0, nx):

            lon = i*0.25-179.875
            lat = j*0.25-59.875

            jj = int(round((lat-lat1)/dlat))
            ii = int(round((lon-lon1)/dlon))

            if(ii==nxo):
                ii = 0
            if(jj==nyo):
                jj = nyo -1
            if(topo_range[i,j]!=-999.):
                topo_model[ii,jj] = topo_model[ii,jj] + topo_range[i,j]
                n[ii,jj] = n[ii,jj] + 1.

    topo_model[n>0] = topo_model[n>0]/n[n>0]
    topo_model[n==0] = -999.
    topo_model[topo_model==-999] = 0

    return np.transpose(topo_model)


def plot_diagnostic(fileout, plot_dir, project_info, model):

    """
    ;; Arguments
    ;;    fileout: dir
    ;;          directory to save the plot
    ;;
    ;; Description
    ;;    Plot diagnostic and save .png plot
    ;;
    """
    import projects
    E = ESMValProject(project_info)
    verbosity = E.get_verbosity()

    #-------------------------
    # Read model info
    #-------------------------
    
    currProject = getattr(vars()['projects'], model.split_entries()[0])()

    start_year = currProject.get_model_start_year(model)
    end_year = currProject.get_model_end_year(model)
    model_info = model.split_entries()


    # Read code output in netCDF format
    diag_name = 'sm-pr-diag-plot'
    netcdf_file = E.get_plot_output_filename(diag_name=diag_name,
                                             variable='pr-mrsos',
                                             model=model_info[1],
                                             specifier='Taylor2012-diagnostic')
    suffix = E.get_graphic_format() + "$"
    netcdf_file = re.sub(suffix, "nc", netcdf_file)
    ncf = nc4.Dataset(os.path.join(fileout, netcdf_file), 'r')
    p_val =  ncf.variables['T12_diag'][:,:]
    ncf.close()

    mp_val = np.ma.masked_equal(p_val, -999)

    # Gridding info: Global diagnostic 5x5 deg
    LON = np.arange(-180, 185, 5)
    LAT = np.arange(-60, 65, 5)

    # Define figure

    F, ax = plt.subplots(nrows=1, ncols=1, **dict(figsize=[15, 8]))

    # Cmap
    cmap = cm.get_cmap(name='bwr_r', lut=7)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = [0., 1, 5, 10, 90, 95, 99, 100]
    norm = col.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad('GainsBoro')

    # Basemap
    map_ax = Basemap(ax=ax, projection='cyl', resolution='l',
                     llcrnrlat=-60, urcrnrlat=60,
                     llcrnrlon=-180, urcrnrlon=180,
                     fix_aspect=False)
    map_ax.drawcoastlines(color='k', linewidth=.7)

    # Plot p_val on basemap
    I = map_ax.pcolormesh(LON, LAT, mp_val, cmap=cmap, norm=norm)

    # Colorbar
    cax = F.add_axes([0.94, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(I, cax=cax, cmap=cmap)

    plt.suptitle('Preference for afternoon precipitation over soil moisture anomalies, ' + model_info[1] + " (" + start_year + "-" + end_year + ")",
                  fontsize = 14)

    diag_name = 'sm-pr-diag-plot'
    figure_filename = E.get_plot_output_filename(diag_name=diag_name,
                                                 variable='pr-mrsos',
                                                 model=model_info[1])

    # Save figure to fileout
    plt.savefig(os.path.join(plot_dir, figure_filename))



def read_pr_sm_topo(project_info, model):

    """
    ;; Arguments
    ;;    project_info: dictionary
    ;;          all info from namelist
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
    ;;    time_bnds_1: float
    ;;          first time_bnd of time series
    ;;
    ;;
    ;; Description
    ;;    Read cmip5 input data for computing the diagnostic
    ;;
    """
    
    import projects
    E = ESMValProject(project_info)
    verbosity = E.get_verbosity()
    #-------------------------
    # Read model info
    #-------------------------

    currProject = getattr(vars()['projects'], model.split_entries()[0])()

    model_info = model.split_entries()

    mip = currProject.get_model_mip(model)
    exp = currProject.get_model_exp(model)
    start_year = currProject.get_model_start_year(model)
    end_year = currProject.get_model_end_year(model)

    years = range(int(start_year), int(end_year) + 1)
    
    '''
    #-------------------------
    # Read model info
    #-------------------------

    model_name = model_info[1]
    time_step = model_info[2]
    exp_fam = model_info[3]
    model_run = model_info[4]
    year_start = model_info[5]
    year_end = model_info[6]
    filedir = model_info[7]

    years = range(int(year_start), int(year_end)+1)
    '''

    
    #-------------------------
    # Input data directories
    #-------------------------
    currDiag = project_info['RUNTIME']['currDiag']

    pr_index = currDiag.get_variables().index('pr')
    pr_field = currDiag.get_field_types()[pr_index]

    sm_index = currDiag.get_variables().index('mrsos')
    sm_field = currDiag.get_field_types()[sm_index]

    indir = currProject.get_cf_outpath(project_info, model)
    in_file = currProject.get_cf_outfile(project_info, model, pr_field, 'pr', mip, exp)
    pr_files = [os.path.join(indir, in_file)]

    in_file = currProject.get_cf_outfile(project_info, model, sm_field, 'mrsos', mip, exp)
    sm_files = [os.path.join(indir, in_file)]
    
    '''
    #-------------------------
    # Input data directories
    #-------------------------
    pr_files = []
    sm_files = []

    for yy in years:

        Patt = filedir+'pr_'+time_step+'_'+model_name+'_'+exp_fam+'_'+\
               model_run+'_'+str(yy)+'*.nc'
        pr_files.append(glob.glob(Patt))

        Patt = filedir+'mrsos_'+time_step+'_'+model_name+'_'+exp_fam+'_'+\
                model_run+'_'+str(yy)+'*.nc'
        sm_files.append(glob.glob(Patt))

    pr_files = [l[0] for l in pr_files if len(l)>0]
    pr_files = sorted(pr_files)

    sm_files = [l[0] for l in sm_files if len(l)>0]
    sm_files = sorted(sm_files)
    '''

    #----------------------
    # Read in precipitation
    #----------------------

    pr_list = []

    for pr_file in pr_files:

        info('Reading precipitation from ' + pr_file, verbosity, required_verbosity=1)

        pr = iris.load(pr_file)[0]

        for at_k in pr.attributes.keys():
            pr.attributes.pop(at_k)

        pr_list.append(pr)

    pr = iris.cube.CubeList(pr_list)
    pr = pr.concatenate()[0]

    # Convert longitude from 0_360 to -180_180

    pr = coord_change([pr])[0]

    # Add metadata: day, month, year

    add_month(pr, 'time')
    add_day_of_month(pr, 'time', name='dom')
    add_year(pr, 'time')

    # Convert units to kg m-2 hr-1

    pr.convert_units('kg m-2 hr-1')

    #-----------------------
    # Read in soil moisture
    #-----------------------

    sm_list = []

    for sm_file in sm_files:

        info('Reading soil moisture from ' + sm_file, verbosity, required_verbosity=1)

        sm = iris.load(sm_file)[0]

        for at_k in sm.attributes.keys():
            sm.attributes.pop(at_k)

        sm_list.append(sm)

    sm = iris.cube.CubeList(sm_list)
    sm = sm.concatenate()[0]

    # Convert longitude from 0_360 to -180_180

    sm = coord_change([sm])[0]

    # Add metadata: day, month, year

    add_month(sm, 'time')
    add_day_of_month(sm, 'time', name='dom')
    add_year(sm, 'time')

    #----------------------------------------------
    # Constrain pr and sm data to latitude 60S_60N
    #----------------------------------------------

    latconstraint = iris.Constraint(latitude=lambda cell: -59.0 <= cell <= 59.0)

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
        sm.data.data[sm.data.mask] = -999.

    except:
        info('no missing data conversion', verbosity, required_verbosity=1)

    #----------------------
    # Read in topography
    #----------------------

    # Topography map specs:
    # latitude 60S_60N
    # longitude 180W_180E
    # model resolution

    #ftopo = currProject.get_cf_fx_file(project_info, model)

    #dt = '>f4'
    #topo = (np.fromfile(ftopo, dtype=dt)).reshape(len(lat), len(lon))

    topo = get_topo(project_info, lon, lat, model)

    #----------------------
    # Read in time bounds
    #----------------------

    indir, infiles = currProject.get_cf_infile(project_info, model, pr_field, 'pr', mip, exp)
    Patt = os.path.join(indir, infiles)
    pr_files = sorted(glob.glob(Patt))

    ncf = nc4.Dataset(pr_files[0])
    time_bnds_1 = ncf.variables['time_bnds'][0][0]
    time_bnds_1 = time_bnds_1 - int(time_bnds_1)
    ncf.close()

    #-----------------------------------------------
    # Return input data to compute sm_pr diagnostic
    #-----------------------------------------------
    return pr, sm, topo, lon, lat, time, time_bnds_1


def removeleap(var_array, utimedate):
    """ 
    ;; Arguments
    ;;    var_array : ndarray
    ;;
    ;;    utimedate : netcdftime object     
    ;;
    ;; Description
    ;;    Removes 29-feb from leap years
    ;;
    """

    datesrange = range(0, var_array.shape[0])
    months, days = [], []
    for date in datesrange:
        t = utimedate.num2date(date)
        months.append(t.month)
        days.append(t.day)
    leapindex = []
    for i, (j,z) in enumerate(zip(days, months)):
         if all([ item in [j,z] for item in [29,2]]):
             leapindex.append(i)
    var_array = np.delete(var_array, tuple(leapindex), axis = 0)
    
    return var_array



def write_nc(fileout, xs, ys, p_vals, project_info, model):

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
    import projects
    E = ESMValProject(project_info)
    verbosity = E.get_verbosity()

    #-------------------------
    # Read model info
    #-------------------------
    
    currProject = getattr(vars()['projects'], model.split_entries()[0])()

    model_info = model.split_entries()

       
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
    diag_name = 'sm-pr-diag-plot'
    netcdf_file = E.get_plot_output_filename(diag_name=diag_name,
                                             variable='pr-mrsos',
                                             model=model_info[1],
                                             specifier='Taylor2012-diagnostic')

    suffix = E.get_graphic_format() + "$"
    netcdf_file = re.sub(suffix, "nc", netcdf_file)
    root_grp = nc4.Dataset(os.path.join(fileout, netcdf_file), 'w', format='NETCDF4')

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



if __name__ == "__main__":
    main()
