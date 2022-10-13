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

import iris

import numpy as np

from shelve import open as shopen
import os
from glob import glob
from netCDF4 import Dataset, num2date
from matplotlib import pyplot

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvalcore.preprocessor._time import extract_time, annual_statistics, regrid_time
from esmvalcore.preprocessor._regrid import extract_levels, regrid
from esmvalcore.preprocessor._area import extract_region

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

def regrid_intersect(cube, region='global'):
    central_longitude = -14.25 #W #-160.+3.5
    central_latitude = -7.56
    cube = regrid(cube, '1x1', 'linear')
    #cube = regrid_to_1x1(cube)
    if region=='global':
        cube = cube.intersection(longitude=(central_longitude-180., central_longitude+180.))
    if region=='midatlantic':
        lat_bnd = 20.
        lon_bnd = 30.
        cube = cube.intersection(longitude=(central_longitude-lon_bnd, central_longitude+lon_bnd),
                                 latitude=(central_latitude-lat_bnd, central_latitude+lat_bnd), )
    return cube


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


def get_obsdata_paths_ts(field, data_type='3d'):
    """
    Time series obs files:
    """
    t2 = '/gws/nopw/j04/esmeval/obsdata-v2/Tier2/WOA/'
    # temperature
    if field in ['tos', 'thetao'] and data_type=='surface':
    #    return [t2 + 'OBS6_WOA_clim_2018_Omon_tos_200007-200007.nc']
        #ob('/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
    #if field in ['tos', 'thetao'] and data_type=='3d':
        #return [t2+'OBS6_WOA_clim_2018_Omon_thetao_200001-200012.nc', ]
        return sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
        #return [t2 + 'OBS6_WOA_clim_2018_Omon_thetao_200007-200007.nc', ]
    # sal
    if field in ['sos', 'so', 'sal'] and data_type=='surface':
        return [t2 + 'OBS6_WOA_clim_2018_Omon_sos_200007-200007.nc', ]
    if field in ['sos', 'so', 'sal'] and data_type=='3d':
        return [t2 + 'OBS6_WOA_clim_2018_Omon_so_200007-200007.nc', ]

    # nutrients
    if field in ['no3', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_no3_200001-200012.nc', ]

    if field in ['o2', ]:
       return ['/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o01_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o02_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o03_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o04_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o05_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o06_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o07_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o08_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o09_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o10_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o11_01.nc', '/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o12_01.nc', ]

       #return sorted(glob('/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/WOA_oxygen/woa18_all_o{01,02,03,04,05,06,07,08,09,10,11,12}_01.nc'))

    if field in ['po4', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_po4_200001-200012.nc', ]
    if field in ['si', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_si_200001-200012.nc', ]
    if field in ['pH', 'ph']:
        return ["/gws/nopw/j04/esmeval/obsdata-v2/Tier2/GLODAP/OBS6_GLODAP_clim_v2.2016b_Oyr_ph_200001-200012.nc", ]

        assert 0
        # return ["/gws/nopw/j04/esmeval/obsdata-v2/Tier2/GLODAP/OBS6_GLODAP_clim_v2.2016b_Oyr_ph_200001-200012.nc", ]

    if field in ['intpp', ]:
        return ['/gws/nopw/j04/esmeval/obsdata-v2/Tier2/Eppley-VGPM-MODIS/OBS_Eppley-VGPM-MODIS_sat_R2018_Omon_intpp_200207-201903.nc', ]

    if field in ['chl', ]:
        return sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier2/ESACCI-OC/OBS6_ESACCI-OC_sat_fv5.0_Omon_chl_199709-202012.nc'))

    if field in ['mld', ]:
        return ['/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/MLD/mld_DR003_c1m_reg2.0.nc', ]



    print('Unable to find:', field, data_type)
    assert 0



def get_obsdata_paths(field, data_type='3d'):

    #files = sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier2/WOA/OBS6_WOA_clim_2018_Omon_thetao_*'))
    #OBS6_WOA_clim_2018_Omon_so_200007-200007.nc      OBS6_WOA_clim_2018_Oyr_no3_200001-200012.nc  OBS_WOA_clim_2013v2_Omon_so_200001-200012.nc      OBS_WOA_clim_2013v2_Oyr_no3_200001-200012.nc
    #OBS6_WOA_clim_2018_Omon_sos_200007-200007.nc     OBS6_WOA_clim_2018_Oyr_o2_200001-200012.nc   OBS_WOA_clim_2013v2_Omon_sos_200007-200007.nc     OBS_WOA_clim_2013v2_Oyr_o2_200001-200012.nc
    #OBS6_WOA_clim_2018_Omon_thetao_200007-200007.nc  OBS6_WOA_clim_2018_Oyr_po4_200001-200012.nc  OBS_WOA_clim_2013v2_Omon_thetao_200001-200012.nc  OBS_WOA_clim_2013v2_Oyr_po4_200001-200012.nc
    #OBS6_WOA_clim_2018_Omon_tos_200007-200007.nc     OBS6_WOA_clim_2018_Oyr_si_200001-200012.nc   OBS_WOA_clim_2013v2_Omon_tos_200007-200007.nc     OBS_WOA_clim_2013v2_Oyr_si_200001-200012.nc

    t2 = '/gws/nopw/j04/esmeval/obsdata-v2/Tier2/WOA/'
    # temperature
    if field in ['tos', 'thetao'] and data_type=='surface':
        return [t2 + 'OBS6_WOA_clim_2018_Omon_tos_200007-200007.nc']
        #ob('/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
    if field in ['tos', 'thetao'] and data_type=='3d':
        #return [t2+'OBS6_WOA_clim_2018_Omon_thetao_200001-200012.nc', ]
        return [t2 + 'OBS6_WOA_clim_2018_Omon_thetao_200007-200007.nc', ]
    # sal
    if field in ['sos', 'so', 'sal'] and data_type=='surface':
        return [t2 + 'OBS6_WOA_clim_2018_Omon_sos_200007-200007.nc', ]
    if field in ['sos', 'so', 'sal'] and data_type=='3d':
        return [t2 + 'OBS6_WOA_clim_2018_Omon_so_200007-200007.nc', ]

    # nutrients
    if field in ['no3', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_no3_200001-200012.nc', ]
    if field in ['o2', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_o2_200001-200012.nc', ]
    if field in ['po4', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_po4_200001-200012.nc', ]
    if field in ['si', ]:
       return [t2 + 'OBS6_WOA_clim_2018_Oyr_si_200001-200012.nc', ]
    if field in ['pH', 'ph']:
        if data_type in ['3d', 'profile', 'map']:
            return ["/gws/nopw/j04/esmeval/obsdata-v2/Tier2/GLODAP/OBS6_GLODAP_clim_v2.2016b_Oyr_ph_200001-200012.nc", ]

    if field in ['chl', ]:
        return sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier2/ESACCI-OC/OBS6_ESACCI-OC_sat_fv5.0_Omon_chl_199709-202012.nc'))

    if field in ['intpp', ]:
        return ['/gws/nopw/j04/esmeval/obsdata-v2/Tier2/Eppley-VGPM-MODIS/OBS_Eppley-VGPM-MODIS_sat_R2018_Omon_intpp_200207-201903.nc', ]

    if field in ['mld', ]:
        return ['/home/users/ldemora/workspace/ESMValTool/run-mpas/auxiliary_data/MLD/mld_DR003_c1m_reg2.0.nc', ]

    print('Unable to find:', field, data_type)
    assert 0


def time_series(field='tos', pane='timeseries', overwrite=False):
    """
    Calculates the time series.
    """
    shelve_path = get_shelve_path(field, pane='ts')
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    if not overwrite and glob(shelve_path+'*'):
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

    files = get_obsdata_paths_ts(field, data_type='surface')
    if not len(files):
        print('No files found:', files, field)
        assert 0
    annual_times, annual_data = [], []
    months_dat = {m:[] for m in months}

    datas = {}
    central_longitude = -14.25 # +/-3 # West -11.25
    central_latitude = -7.56 # +/3 # North

    for nc_path in sorted(files):
        print('loading:', nc_path)
        if field in ['mld', ]:
            print('loading:', nc_path)
            cube = iris.load_raw(nc_path)[2]
            cube.data = np.ma.masked_where(cube.data == 1000000000., cube.data)
            nctimes = diagtools.cube_time_to_float(cube)
            #print(cube)
            #assert 0
        elif field in ['o2', ]:
    #        float o_an(time, depth, lat, lon) ;
    #            o_an:standard_name = "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water" ;
    #            o_an:long_name = "Objectively analyzed mean fields for mole_concentration_of_dissolved_molecular_oxygen_in_sea_water at standard depth levels." ;
    #            o_an:coordinates = "time lat lon depth" ;
    #            o_an:cell_methods = "area: mean depth: mean time: mean within years time: mean over years" ;
    #            o_an:grid_mapping = "crs" ;
    #            o_an:units = "micromoles_per_kilogram" ;
    #            o_an:_FillValue = 9.96921e+36f ;


            cube = iris.load_raw(nc_path, 'o_an')[0]
            print(cube.units)
            #assert 0

#            for c in cubes:
 
            cube.data = np.ma.masked_where(cube.data.mask+(cube.data==0.)+(cube.data>1e10), cube.data)
            #cube = diagtools.bgc_units(cube, field)


            # woa18_all_o01_01.nc
            nctimes = [(float(nc_path[-8:-6])-1.)/12. + 2000., ]
#           nctimes = np.arange(12)/12.+2000.
            #cube = diagtools.bgc_units(cube, field)
            print(nc_path, nctimes, cube.data.max())
            #assert 0

        else:
            cube = iris.load_cube(nc_path)
            cube = diagtools.bgc_units(cube, field)

            nctimes = diagtools.cube_time_to_float(cube)
        if field == 'chl':
            cube.data = cube.data *1000.
        nc = Dataset(nc_path)

        if field == 'o2':
            print('depth',  np.abs(cube.coord(axis='Z').points)*-1.)
            #assert 0

            #cube = extract_levels(cube,
            #    scheme='linear',
            #    levels =  [550., 600., 650.,700., 750.,800.,850., 900., 950., ],
            #)
            #cube.data = np.ma.masked_where(cube.data.mask+(cube.data==0.), cube.data)
            #cube = cube.collapsed(['depth',], iris.analysis.MEAN) # can do this as layers are equal.

            print('ncdata 0a :range:  min', cube.data.min(), cube.data.mean(), cube.data.max())


            cube = extract_levels(cube,
                scheme='linear',
                levels =  [500.,]
            )
            print('ncdata 0b :range:  min', cube.data.min(), cube.data.mean(), cube.data.max())

            cube.data = np.ma.masked_where(cube.data.mask+(cube.data==0.)+(cube.data>1e10), cube.data)
            print('ncdata 0c :range:  min', cube.data.min(), cube.data.mean(), cube.data.max())

            print('file:', nc_path)
            print('times:', nctimes)
            print('ncdata:range:  min', cube.data.min(), cube.data.mean(), cube.data.max())
            print('file:', nc_path)
#            assertlater = True
            #assert 0

            #cube = cube.collapsed(['depth',], iris.analysis.MEAN) # can do this as layers are equal.

            #cube.data = np.ma.masked_where(cube.data.mask+(cube.data==0.), cube.data)

        lons = nc.variables['lon']
        lats = nc.variables['lat']

        lon_min = np.argmin(np.abs(lons[:] -(central_longitude +360.-3.)))
        lon_max = np.argmin(np.abs(lons[:] -(central_longitude +360.+3.)))
        print('lon_min', lon_min, lons[lon_min])
        print('lon_max', lon_max, lons[lon_max])
        print('source lon range', np.min(lons), np.max(lons))
        lon_mins = [np.argmin(np.abs(lons[:] -(360.+central_longitude -3.))), np.argmin(np.abs(lons[:] -(360.+central_longitude +3.)))]
        lat_maxs = [np.argmin(np.abs(lats[:] -(central_latitude -3.))), np.argmin(np.abs(lats[:] -(central_latitude +3.)))]
        print('lonsmin', lon_mins, lons[lon_mins[0]], lons[lon_mins[1]])
        print('latsmin', lat_maxs, lats[lat_maxs[0]], lats[lat_maxs[1]])

        print('ncdata:1  range:  min', cube.data.min(), cube.data.mean(), cube.data.max())

 #       print(cube.dimensions)
        if cube.ndim == 3:
            latdim = 1
            londim = 2
        elif cube.ndim == 4:
            latdim = 2
            londim = 3
        else: assert 0

        if cube.dim_coords[latdim].var_name == 'lat' and cube.dim_coords[londim].var_name == 'lon':
            pass
        else: assert 0
        print('ncdata:2  range:  min', cube.data.min(), cube.data.mean(), cube.data.max())


        if cube.ndim == 3:
            ncdata = cube[...,lat_maxs[0]:lat_maxs[1]+1, lon_mins[0]:lon_mins[1]+1]
        elif cube.ndim == 4:
            ncdata = cube[:,0,lat_maxs[0]:lat_maxs[1]+1, lon_mins[0]:lon_mins[1]+1]
        else: assert 0

        print('ncdata:3  range:  min', ncdata.data.min(), ncdata.data.mean(), ncdata.data.max())

        if len(files) == 1 and len(nctimes) ==1:

            print('ONLY one time point:')
            print(cube)
            print('times:', nctimes)
            print('ncdata 4a :range:  min', ncdata.data.min(), ncdata.data.mean(), ncdata.data.max())
            print('file:', nc_path)
            assert 0
        else:
            print('file:', nc_path)
            print('times:', nctimes)
            print('ncdata 4b:range:  min', ncdata.data.min(), ncdata.data.mean(), ncdata.data.max())
            print('file:', nc_path)
#            assertlater = True
#            assert 0

        try:
            ncdata = ncdata.collapsed(['lat', 'lon'], iris.analysis.MEAN).data
        except:
            ncdata = ncdata.collapsed(['latitude', 'longitude'], iris.analysis.MEAN).data

        print('ncdata 5 :range:  min', ncdata.min(), ncdata.mean(), ncdata.max())


#       if field == 'tos':
#           ncdata =  cube[:, 106:114+1, 457:465+1].collapsed(['lat', 'lon'], iris.analysis.MEAN).data #.mean(axis=(1,2))
#       elif field =='chl':
#           ncdata =  cube[:, 0, 317:341+1,1439:1439+1].collapsed(['latitude', 'longitude'], iris.analysis.MEAN).data

        for t,d in zip(nctimes, ncdata):
            print('adding data:', field, t, d, cube.units)
            #if 08
            datas[t] = d
        #times.extend(nctimes)
        #data.extend(ncdata)
    #print(files)
    #assert 0
    times = [t for t in sorted(datas.keys())]
    data = [datas[t] for t in times]
    years = {int(t):[] for t in times}

    # calculate annual and climatological data.
    for t, d in zip(times, data):
        years[int(t)].append(d)
        print('calculationg clim', int(t), t, d)

    for yr in sorted(years.keys()):
        annual_times.append(yr + 0.5)
        annual_data.append(np.mean(years[yr]))

#        if yr < 2000.: continue
#        if yr > 2010.: continue
        for m, d in zip(months, years[yr]):
            print('clim:', yr, m, d)
            months_dat[m].append(d)

    clim = [np.mean(months_dat[m]) for m in months]

    sh = shopen(shelve_path)
    sh['times'], sh['data'] = times, data
    sh['annual_times'], sh['annual_data'] = annual_times, annual_data
    sh['clim'] = clim

    if field == 'tos':
        sh['header'] = 'TOS calculated from ERA-Interim, monthly data'
    elif field =='chl':
        sh['header'] = 'chl calculated from ESACCI-OC_sat, monthly data'
    else:
        sh['header'] = 'WOA 2018 ' +field
    sh['files'] = files
    sh.close()

    if pane == 'monthly_timeseries':
        return times, data
    if pane == 'timeseries':
        return annual_times, annual_data
    if pane == 'clim':
        return months, clim


def load_map_netcdf(field='tos', pane='map'):
    """
    Time range is 2000-2010.
    """
    path = diagtools.folder(['aux', 'obs_ncs'])
    path += '_'.join([field, pane]) +'.nc'

    if os.path.exists(path):
        print('netcdf exists:', path)
        return iris.load_cube(path)

    files = get_obsdata_paths(field, data_type='surface')

#    if field=='tos':
#        files = sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier3/ERA-Interim/OBS6_ERA-Interim_reanaly_1_Omon_tos_*.nc'))
#    elif field == 'chl':
#        files = sorted(glob('/gws/nopw/j04/esmeval/obsdata-v2/Tier2/ESACCI-OC/OBS6_ESACCI-OC_sat_fv5.0_Omon_chl_199709-202012.nc'))
#    else: assert 0

    cube_list = []
    for fn in files:
        if field in ['mld', ]:
            cube = iris.load_raw(fn)[2]
            cube.data = np.ma.masked_where(cube.data == 1000000000., cube.data)
        else:
            cube = iris.load_cube(fn)

            cube = diagtools.bgc_units(cube, field)
            times = diagtools.cube_time_to_float(cube)
#           if cube.ndim == 4:
#               cube = cube[:,0]
            cube = extract_time(cube, 2000, 1, 1, 2010, 1, 1)
            if field == 'o2':
                cube = extract_levels(cube,
                scheme='linear',
                levels =  [500., ]
                )
            elif cube.ndim == 4:
                cube = cube[:,0]


#        if np.min(times) < 2000.: continue
#        if np.max(times) > 2010.: continue

        # assumes tium
#       if field in ['tos',]:
#           if np.min(times) < 2000.: continue
#           if np.max(times) > 2010.: continue
#        if field in ['chl',]:
#            cube = cube[:,0] # extract surface layer
#           if np.min(times) > 2010.: continue
#           if np.max(times) < 2000.: continue
#            cube = extract_time(cube, 2000, 1, 1, 2010, 1, 1)

        print('loaded:', fn)
        new_cube = cube.collapsed('time', iris.analysis.MEAN)
        cube_list.append(new_cube.copy())

    outcube  = diagtools.make_mean_of_cube_list_notime(cube_list)

    print('saving netcdf:', path)
    iris.save(outcube, path)
    return outcube


def make_ts_figure(field, overwrite=True):
    """
    3 Pane Time series figure:
        Monthly, annual and climatology panes.
    """

    path = diagtools.folder('images/obs/timeseries')
    path +='_'.join([field, 'ts'])+'.png'

    if not overwrite and os.path.exists(path):
        print('Already exists:', path)
        return

    fig = pyplot.figure()
    ax = fig.add_subplot(311)
    time, data = time_series(field=field, pane='timeseries',overwrite=overwrite)
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


def make_map_figure(field):
    """
    Make a map of the surface in the historical period.
    """
    path = diagtools.folder('images/obs/timeseries')
    path +='_'.join([field, 'map'])+'.png'

    if os.path.exists(path):
        print('Already exists:', path)
        return

    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    cube = load_map_netcdf(field=field, pane='map')
    if field == 'chl':
        cube.data = np.clip(cube.data, 0., 5.)

    qplot = iris.plot.contourf(
        cube,
        linewidth=0,
        )

    pyplot.colorbar()
    pyplot.title(' '.join([field, ',', str(cube.units)]))
    print('saving figure:', path)
    pyplot.savefig(path)
    pyplot.close()


def load_profile_netcdf(field='tos', pane='profile'):
    """
    Perform the calculation of the depth

    """

    path = diagtools.folder(['aux', 'obs_ncs'])
    path += '_'.join([field, pane]) +'.nc'

    if os.path.exists(path):
        print('netcdf exists:', path)
        return iris.load_cube(path)

    files = get_obsdata_paths(field, data_type='3d')
    cube_list = []
    central_longitude = -14.25 #W #-160.+3.5
    central_latitude = -7.56

    for fn in files:
        print('load_profile_netcdf:', fn)
        cube = iris.load_cube(fn)
        cube = diagtools.bgc_units(cube, field)
        times = diagtools.cube_time_to_float(cube)
        if np.min(times) < 2000.: continue
        if np.max(times) > 2010.: continue
        cube = extract_region(cube, central_longitude -3., central_longitude+3.,
            central_latitude-3., central_latitude+3.)

        new_cube = cube.collapsed('time', iris.analysis.MEAN)
        #new_cube = new_cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
        cube_list.append(new_cube.copy())

    if len(cube_list) == 1:
        outcube = cube_list[0]
    else:
        outcube  = diagtools.make_mean_of_cube_list_notimecube_list(cube_list)

    outcube = regrid_intersect(outcube, region='midatlantic')
    outcube = outcube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

    outcube = extract_levels(outcube,
        scheme='linear',
        levels =  [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 150., 200.0, 250., 300.0, 350., 400.0, 450., 500.0,
                   600.0, 650., 700.0, 750., 800.0, 850., 900.0, 950., 999.0,
                   1001., 1250., 1500.0, 1750., 2000.0, 2250., 2500.0, 2750., 3000.0, 3250., 3500.0, 3750.,
                   4000.0, 4250., 4500.0, 4750., 5000.0]
        )

    print('saving netcdf:', path)
    iris.save(outcube, path)
    return outcube



def make_profile_figure(field):
    """
    Make a map of the surface in the historical period.
    """
    path = diagtools.folder('images/obs/timeseries')
    path +='_'.join([field, 'profile'])+'.png'

    if os.path.exists(path):
        print('Already exists:', path)
        return

    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    cube = load_profile_netcdf(field=field, pane='profile')
    depths = np.abs(cube.coord(axis='Z').points)*-1.
    ax.plot(cube.data, depths,
        lw=2,
        ls='-',
        c='k',
        label=field)
    pyplot.title(field+ ' profile')

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
    twodfields = []#ld', ] #'intpp', 'chl', ]# 'mld' ]
    threedfields = ['po4', ] #'no3','po4', ]#'tos',] #'o2',] #'no3', ] #'ph', ]#'so', ]# 'o2',] #'tos', ] #'o2', 'so','ph',  'tos',]#  'no3', 'si',]
    for field in twodfields:
#        make_map_figure(field)
        make_ts_figure(field)

    for field in threedfields:
        make_profile_figure(field)
        make_map_figure(field)
        make_ts_figure(field)

if __name__ == '__main__':
    main()
