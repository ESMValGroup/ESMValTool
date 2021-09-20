"""
This is an independent one off script to calculated the observations
for the Ascension Islands plots.
It will save the file in a easy to read output in a given aux folder.
This saves it from being calculated every time.

Ascensionb island region is:
    central_longitude = -14.25  +/-3 # Northh -11.25
    central_latitude = -7.56 +/3 # East

"""
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import numpy as np

from shelve import open as shopen
import os
from glob import glob
from netCDF4 import Dataset, num2date
from matplotlib import pyplot

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools


def get_shelve_path(field, pane='timeseries'):
    shelve_path = diagtools.folder(['aux', 'obs_shelves'])
    shelve_path += '_'.join([field, pane]) +'.shelve'
    return shelve_path

def load_times_from_nc(nc):
    """
    Converts nc into datetimes
    """
    datetimes = num2date(nc.variables['time'],  units=nc.variables['time'].units, calendar=nc.variables['time'].calendar)
    return datetimes


def nc_time_to_float(nc):
    """
    Converts nc time to float time in units of decimal years.
    """
    datetimes = load_times_from_nc(nc)
    if nc.variables['time'].calendar == 'gregorian':
        daysperyear = 365.25
    else:
        assert 0

    times = []

    for dtime in datetimes:
        try:
            dayofyr = dtime.dayofyr
        except AttributeError:
            time = datetime(dtime.year, dtime.month, dtime.day)
            time0 = datetime(dtime.year, 1, 1, 0, 0)
            dayofyr = (time - time0).days

        floattime = dtime.year + dayofyr / daysperyear + dtime.hour / (
                24. * daysperyear)
        times.append(floattime)
    return times


def time_series(field='tos', pane='timeseries'):
    """
    Calculates the time series.
    """
    shelve_path = get_shelve_path(field, pane='ts')
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    if glob(shelve_path+'*'):
        print('shelve exists:', shelve_path+'*')
        sh = shopen(shelve_path)
        times, data = sh['times'], sh['data']
        annual_times, annual_data = sh['annual_times'], sh['annual_data']
        clim = sh['clim']
        sh.close()

        if pane == 'monthly_timeseries':
            return times, data
        if pane == 'timeseries':
            return annual_times, annual_data
        if pane == 'clim':
            return months, clim

    if field == 'tos':
        files = sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
    times, data = [], []
    annual_times, annual_data = [], []
    months_dat = {m:[] for m in months}

    central_longitude = -14.25 # +/-3 # West -11.25
    central_latitude = -7.56 # +/3 # North
    for nc_path in sorted(files):
        print('loading:', nc_path)
        nc = Dataset(nc_path, 'r')
        nctimes = nc_time_to_float(nc)
       
        # lons = nc.variables['lon']
        # np.argmin(np.abs(lons[:] -(360-central_longitude -3)))
        # np.argmin(np.abs(lons[:] -(360-central_longitude +3)))

        # lons = nc.variables['lon']
        # np.argmin(np.abs(lons[:] -(360-central_longitude -3)))
        # np.argmin(np.abs(lons[:] -(360-central_longitude +3)))
        # np.argmin(np.abs(lats[:] -(central_latitude -3)))
        # np.argmin(np.abs(lats[:] -(central_latitude +3)))
        ncdata =  nc.variables['tos'][:, 106:114+1, 457:465+1].mean(axis=(1,2))
        times.extend(nctimes)
        data.extend(ncdata)

        annual_times.append(np.mean(nctimes))
        annual_data.append(np.mean(ncdata))
        print('data:', np.mean(nctimes), np.mean(ncdata))

        #calculate clim
        if np.min(nctimes) < 2000.: continue
        if np.max(nctimes) > 2010.: continue
        for m, d in zip(months, ncdata[:]):
            print('clim:', m, d)
            months_dat[m].append(d)

    clim = [np.mean(months_dat[m]) for m in months]

    sh = shopen(shelve_path)
    sh['times'], sh['data'] = times, data
    sh['annual_times'], sh['annual_data'] = annual_times, annual_data
    sh['clim'] = clim

    sh['header'] = 'TOS calculated from ERA-Interim, monthly data'
    sh['files'] = files
    sh.close()

    if pane == 'monthly_timeseries':
        return times, data
    if pane == 'timeseries':
        return annual_times, annual_data
    if pane == 'clim':
        return months, clim



def make_figure(field):

    path = diagtools.folder('images/obs/timeseries')
    path +='_'.join([field, 'ts'])+'.png'

    if os.path.exists(path):
        print('Already exists:', path)
        return

    fig = pyplot.figure()
    ax = fig.add_subplot(311)
    time, data = time_series(field=field, pane='timeseries')
    pyplot.plot(time,data)
    pyplot.title('Annual ' + field)

    ax = fig.add_subplot(312)
    time, data = time_series(field=field, pane='monthly_timeseries')
    pyplot.plot(time,data)
    pyplot.title('Monthly ' + field)

    ax = fig.add_subplot(313)
    months, data = time_series(field=field, pane='clim')
    print(months, data)
    pyplot.plot([t for t,d in enumerate(data)], data)
    ax.set_xticks([t for t,d in enumerate(data)])
    ax.set_xticklabels(months)
    pyplot.title('Climatological ' + field)

    print('saving figure:', path)
    pyplot.savefig(path)
    pyplot.close()


#
# def load_WOA_data(cfg, short_name, plot, grid='1'):
#     """
#     plan to calculate the WOA data here.
#     """
#
#
#     # If shelve exists, then return that data
#     if os.path.exists(shelve_path):
#         sh = shopen(shelve_path)
#         times, data = sh['times'], sh['data']
#         sh.close()
#         return times,data
#
#
#
#
#     # Look for the original netcdf data
#     files = sorted(glob.glob()'/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
#     # load the data
#
#     # Extract the relevant data region
#     #106:115 is array([-10.5 ,  -9.75,  -9.  ,  -8.25,  -7.5 ,  -6.75,  -6.  ,  -5.25, -4.5 ]) # Latitude
#     #
#
#      cube = cube[:,106:115, ]
#     # apply some kind of masking.
#
#
#     sh = shopen(shelve_path)
#     sh['times'], sh['data'] = times, data
#     sh.close()

def main():
     make_figure('tos')

if __name__ == '__main__':
    main()
