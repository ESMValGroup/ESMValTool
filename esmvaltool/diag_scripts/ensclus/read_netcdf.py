"""Read netCDF file of 3D field."""

# Standard packages
import os

import numpy as np
from netCDF4 import Dataset, num2date

import iris


def read_iris(ifile):
    """Read netCDF file of 3D field using iris.

    USAGE: var, lat, lon, dates = read_3d_ncfield_iris(filename)
    """
    cube = iris.load_cube(ifile)
    variabs = [coord.name() for coord in cube.coords()]

    if 'lat' in variabs:
        lat = cube.coord('lat').points
    elif 'latitude' in variabs:
        lat = cube.coord('latitude').points
    if 'lon' in variabs:
        lon = cube.coord('lon').points
    elif 'longitude' in variabs:
        lon = cube.coord('longitude').points
    time = cube.coord('time')
    time_units = str(cube.coord('time').units)
    dates = time.units.num2date(time.points)
    var_units = str(cube.units)
    var = cube.data
    if isinstance(var, np.ma.masked_array):
        var = var.filled(fill_value=np.nan)

    return var, var_units, lat, lon, dates, time_units


def read_3d_ncfield(ifile):
    """Read netCDF file of 3Dfield.

    USAGE: var, lat, lon, dates = read_3d_ncfield(filename)
    """
    fileh = Dataset(ifile, mode='r')
    variabs = []
    for variab in fileh.variables:
        variabs.append(variab)
    print('The variables in the nc file are: ', variabs)

    if 'lat' in variabs:
        lat = fileh.variables['lat'][:]
    elif 'latitude' in variabs:
        lat = fileh.variables['latitude'][:]
    if 'lon' in variabs:
        lon = fileh.variables['lon'][:]
    elif 'longitude' in variabs:
        lon = fileh.variables['longitude'][:]
    time = fileh.variables['time'][:]
    time_units = fileh.variables['time'].units
    var_units = fileh.variables[variabs[0]].units
    var = fileh.variables[variabs[0]][:, :, :]
    dates = num2date(time, time_units)
    fileh.close()

    return var, var_units, lat, lon, dates, time_units


def save_n_2d_fields(lats, lons, variab, varname, varunits, ofile):
    """Save var in ofile netCDF file.

    Save a number N of 2D fields [lat x lon]
    """
    try:
        os.remove(ofile)  # Remove the outputfile
    except OSError:
        pass
    dataset = Dataset(ofile, 'w', format='NETCDF4_CLASSIC')
    # print(dataset.file_format)

    num = dataset.createDimension('num', variab.shape[0])
    lat = dataset.createDimension('lat', variab.shape[1])
    lon = dataset.createDimension('lon', variab.shape[2])

    # Create coordinate variables for 3-dimensions
    num = dataset.createVariable('num', np.int32, ('num',))
    lat = dataset.createVariable('lat', np.float32, ('lat',))
    lon = dataset.createVariable('lon', np.float32, ('lon',))
    # Create the actual 3-d variable
    var = dataset.createVariable(varname, np.float64, ('num', 'lat', 'lon'))

    # print('variable:', dataset.variables[varname])
    # for varn in dataset.variables.keys():
    #    print(varn)
    # Variable Attributes
    lat.units = 'degree_north'
    lon.units = 'degree_east'
    var.units = varunits

    num[:] = np.arange(variab.shape[0])
    lat[:] = lats
    lon[:] = lons
    var[:, :, :] = variab

    dataset.close()

    # -----------------------------------------------------------------------
    print('The {0} 2D fields [num x lat x lon] are saved as \n{1}'
          .format(variab.shape[0], ofile))
    print('__________________________________________________________')
    # -----------------------------------------------------------------------


def read_n_2d_fields(ifile):
    """Read a number N of 2D fields [latxlon].

    USAGE: var, lat, lon, dates = read_n_2d_fields(filename)
    """
    fileh = Dataset(ifile, mode='r')
    variabs = []
    for variab in fileh.variables:
        variabs.append(variab)
    # print('The variables in the nc file are: ', variabs)

    # num = fh.variables['num'][:]
    if 'lat' in variabs:
        lat = fileh.variables['lat'][:]
    elif 'latitude' in variabs:
        lat = fileh.variables['latitude'][:]
    if 'lon' in variabs:
        lon = fileh.variables['lon'][:]
    elif 'longitude' in variabs:
        lon = fileh.variables['longitude'][:]
    var = fileh.variables[variabs[3]][:, :, :]
    var_units = fileh.variables[variabs[3]].units
    # print(fh.variables)
    fileh.close()

    # print('\n'+txt)

    return var, var_units, lat, lon
