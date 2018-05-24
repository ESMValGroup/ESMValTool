'''
This file contains several modules which are needed for the calculation
of various different monsoon metrics
'''

import os

import numpy as np

import iris
import iris.fileformats.pp
import iris.coord_categorisation
import iris.analysis.cartography

from ..utility.area_utils import area_average
from ..auto_assess_deprecated.loaddata import load_run_ss



OBS_TO_META = {
    'precip': ('m01s05i216', 'precipitation_flux', 'kg m-2 day-1'),
    'temp': ('m01s03i236', 'air_temperature', 'Celsius'),
    'mslp': ('m01s16i222', 'air_pressure_at_sea_level', 'Pa'),
    'Uwind850': ('m01s30i201', 'x_wind', 'm s-1'),
    'Uwind300': ('m01s30i201', 'x_wind', 'm s-1'),
    'Uwind200': ('m01s30i201', 'x_wind', 'm s-1'),
    'Vwind850': ('m01s30i202', 'y_wind', 'm s-1'),
    'Vwind200': ('m01s30i202', 'y_wind', 'm s-1'),
    'geopotential': ('m01s30i207', 'geopotential_height', 'm'),
    '10mUwind': ('m01s03i225', 'x_wind', 'm s-1'),
    }

OBS_TO_FILES = {
    'precip': 'GPCP/monthly_v2.2/gpcp_v2.2.nc',
    'temp': 'CRU/cru_ts3.23.1901.2014.tmp.dat.nc',
    'mslp': 'ERA-Interim/MonthlyPP/ERAI_19892011_16222.pp',
    'Uwind850': 'ERA-Interim/MonthlyPP/ERAI_19892011_30201_850.pp',
    'Uwind300': 'ERA-Interim/MonthlyPP/ERAI_19892011_30201_300.pp',
    'Uwind200': 'ERA-Interim/MonthlyPP/ERAI_19892011_30201_200.pp',
    'Vwind850': 'ERA-Interim/MonthlyPP/ERAI_19892011_30202_850.pp',
    'Vwind200': 'ERA-Interim/MonthlyPP/ERAI_19892011_30202_200.pp',
    'geopotential': 'ERA-Interim/MonthlyPP/ERAI_19892011_30207_500.pp',
    '10mUwind': 'ERA-Interim/MonthlyPP/ERAI_19892011_3225.pp',
    }

STASH_TO_KEYS = {
    'm01s05i216': ('m01s05i216', dict(lbtim=[121, 122])),
    30207500: ('m01s30i207', dict(lbproc=128, lblev=500)),
    30201850: ('m01s30i201', dict(lbproc=128, lblev=850)),
    30201300: ('m01s30i201', dict(lbproc=128, lblev=300)),
    30201200: ('m01s30i201', dict(lbproc=128, lblev=200)),
    30202850: ('m01s30i202', dict(lbproc=128, lblev=850)),
    30202200: ('m01s30i202', dict(lbproc=128, lblev=200)),
    }


def regll(x1, x2, y1, y2):
    '''Create dictionary that describes regional cutout'''
    return dict(longitude=(x1, x2), latitude=(y1, y2))


def _data_convert(cube):
    '''Converts units of data in cube'''

    stash = str(cube.attributes['STASH'])
    if stash == 'm01s05i216':     # Precip
        cube.convert_units('kg m-2 day-1')
    elif stash == 'm01s03i236':   # 1.5m Temperature
        cube.convert_units('Celsius')
    elif stash == 'm01s16i222':   # MSLP
        cube.convert_units('hPa')


def extract_model_month(run, stash, months):
    '''This routine processes monthly data'''

    # Get proper stash code and search keys
    (stash, keys) = STASH_TO_KEYS.get(stash, (stash, {}))

    # Extract monthly data and create a jjas mean
    data = load_run_ss(run, 'monthly', stash, lbmon=months, **keys)

    # Convert units of data
    _data_convert(data)

    # Return stats
    return data


def extract_model_season(run, stash, season):
    '''This routine processes seasonal data'''

    # Get proper stash code and search keys
    (stash, keys) = STASH_TO_KEYS.get(stash, (stash, {}))

    # Extract seasonal data
    data = load_run_ss(run, 'seasonal', stash, **keys)

    # Extract season
    iris.coord_categorisation.add_season(data, 'time', name='clim_season')
    data_seas = data.extract(iris.Constraint(clim_season=season))

    # Convert units of data
    _data_convert(data_seas)

    # Return stats
    return data_seas


def mag(cubeA, cubeB, longitude, latitude):
    '''
    Returns the magnitude sum of the two input cubes
    eg use U and V wind to calculate wind speed
    '''
    regll = dict(longitude=longitude, latitude=latitude)
    cubeA_reg = area_average(cubeA, weighted=True, **regll)
    cubeB_reg = area_average(cubeB, weighted=True, **regll)
    magnitude = ((cubeA_reg * cubeA_reg) + (cubeB_reg * cubeB_reg)) ** 0.5
    return magnitude


def calc_pattern_corr(run, obs, season, region):
    '''
    Routine to calculate pattern correlation

    rcorr = np.corrcoef(obs_cube.data.flatten(),
                        model_cube.data.flatten())

    Result is 2x2 matrix where we want the [0,1] entry
    '''

    # Extract obs and model data
    (obs_data, mod_data, mask_data) = _extract_obs_model(run, obs, season)

    # Extract region to compare
    obs_reg = obs_data.intersection(**region)
    mod_reg = mod_data.intersection(**region)

    # Calculate correlation coefficients
    rcorr = np.corrcoef(obs_reg.data.flatten(), mod_reg.data.flatten())

    # Note that result above is a 2x2 matrix where diagonal elements are 1,
    # and [0, 1] and [1, 0] should have the same value between -1 and 1.
    return float(rcorr[0, 1])


def calc_rms_error(run, obs, season, region, mask=False):
    '''Routine to calculate RMS Errors'''
    # TODO (use iris analysis function):
    # error = data - obs
    # rmse = error.collapsed('longitude', iris.analysis.RMS)

    # Extract data
    (obs_data, mod_data, mask_data) = \
        _extract_obs_model(run, obs, season, mask=mask)

    # Calculate square error
    sqerr = (mod_data - obs_data) ** 2

    # Calculate MSE
    if mask:
        mse = area_average(sqerr, mask=mask_data, **region)
    else:
        mse = area_average(sqerr, **region)

    # Return RMSE
    return float((mse ** 0.5).data)


def _extract_obs_model(run, obs, season, mask=False):
    '''Routine to extract matching obs and model data, regridded to obs grid'''
    jjas = [6, 7, 8, 9]

    # Extract obs data
    data = extract_observations(run, obs, season)
    obs_data = data.collapsed('time', iris.analysis.MEAN)

    # Extract model data
    if season == 'jjas':
        mod = extract_model_month(run, OBS_TO_META[obs][0], jjas)
    else:
        mod = extract_model_season(run, OBS_TO_META[obs][0], season)
    mod_data = mod.collapsed('time', iris.analysis.MEAN)

    # Regrid model data to observational grid
    regrid_scheme = iris.analysis.AreaWeighted()
    guess_bounds_for_regrid(obs_data)
    guess_bounds_for_regrid(mod_data)
    rmod_data = mod_data.regrid(obs_data, regrid_scheme)

    # Extract regridded model mask if required
    if mask:
        if season == 'jjas':
            mask = extract_model_month(run, 'm01s03i395', jjas)
        else:
            mask = extract_model_season(run, 'm01s03i395', season)
        mask.data = np.asarray(mask.data, dtype=np.float64)
        mask_avg = mask.collapsed('time', iris.analysis.MEAN)
        mask_avg.data = np.asarray(mask_avg.data, dtype=np.float32)
        guess_bounds_for_regrid(mask_avg)
        rmask_data = mask_avg.regrid(obs_data, regrid_scheme)
    else:
        rmask_data = None

    return (obs_data, rmod_data, rmask_data)


def guess_bounds_for_regrid(cube):
    '''Checks the horizontal coords of a cube, and adds bounds to them:'''

    # check what horizontal coords are used:
    coord_list = [coord.name() for coord in cube.coords()]
    if 'latitude' in coord_list:
        latitude_coord_name = 'latitude'
        longitude_coord_name = 'longitude'
    elif 'grid_latitude' in [coord.name() for coord in cube.coords()]:
        latitude_coord_name = 'grid_latitude'
        longitude_coord_name = 'grid_longitude'
    else:
        mesg = 'Cannot deal with coordinates "{:s}"'.format(coord_list)
        raise NotImplementedError(mesg)

    # add bounds to those coords:
    if cube.coord(longitude_coord_name).has_bounds():
        cube.coord(longitude_coord_name).bounds = None
    cube.coord(longitude_coord_name).guess_bounds()

    if cube.coord(latitude_coord_name).has_bounds():
        cube.coord(latitude_coord_name).bounds = None
    cube.coord(latitude_coord_name).guess_bounds()


def extract_observations(run, stash, season):
    '''Routine to load observations'''
    season_numbers = {'djf': [12, 1, 2],
                      'mam': [3, 4, 5],
                      'jja': [6, 7, 8],
                      'son': [9, 10, 11],
                      'jjas': [6, 7, 8, 9],
                      'as': [8, 9]}

    # Generate callback
    (scode, lname, units) = OBS_TO_META[stash]
    callback = generate_callback(scode, lname, units)

    # Load obs
    fname = os.path.join(run['clim_root'], OBS_TO_FILES[stash])
    if stash == 'precip':
        has_var_name = lambda cube: cube.var_name == 'precip'
        obs = iris.load_cube(fname,
                             iris.Constraint(cube_func=has_var_name),
                             callback=callback)
    else:
        obs = iris.load_cube(fname, callback=callback)

    iris.coord_categorisation.add_month_number(obs, 'time',
                                               name='month_number')
    iris.coord_categorisation.add_season_year(obs, 'time',
                                              name='season_year')
    lamfn = lambda year: 1989 <= year <= 2008
    obs = obs.extract(iris.Constraint(month_number=season_numbers[season],
                                      season_year=lamfn))

    # Convert units
    _data_convert(obs)

    return obs


def generate_callback(scode, stdname, units):
    '''Routine to generate a callback function for iris.load'''
    def callback(cube, field, filename):
        cube.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi(scode)
        if cube.var_name:
            varname = cube.var_name
        else:
            varname = None
        cube.rename(stdname)
        if varname:
            cube.var_name = varname
        cube.units = units
        coordsys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        if cube.coords('latitude'):
            if cube.coord('latitude').coord_system is None:
                cube.coord('latitude').coord_system = coordsys
        if cube.coords('longitude'):
            if cube.coord('longitude').coord_system is None:
                cube.coord('longitude').coord_system = coordsys
    return callback

