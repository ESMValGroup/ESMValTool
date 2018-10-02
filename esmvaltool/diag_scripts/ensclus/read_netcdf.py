# Standard packages
from netCDF4 import Dataset, num2date, date2num, datetime
import numpy as np
import os

def read3Dncfield(ifile):
    '''
    GOAL
        Read netCDF file of 3Dfield
    USAGE
        var, lat, lon, dates = read3Dncfield(filename)
    '''
    #----------------------------------------------------------------------------------------
    #print('__________________________________________________________')
    #print('Reading the 3D field [time x lat x lon]: \n{0}'.format(ifile))
    #----------------------------------------------------------------------------------------
    fh = Dataset(ifile, mode='r')
    variabs=[]
    for variab in fh.variables:
        variabs.append(variab)
    print('The variables in the nc file are: ', variabs)
    
    if ('lat' in variabs):
        lat     = fh.variables['lat'][:]
    elif ('latitude' in variabs):
        lat     = fh.variables['latitude'][:]
    if ('lon' in variabs):
        lon     = fh.variables['lon'][:]
    elif ('longitude' in variabs):
        lon     = fh.variables['longitude'][:]
    time        = fh.variables['time'][:]
    time_units  = fh.variables['time'].units
    var_units   = fh.variables[variabs[0]].units
    var         = fh.variables[variabs[0]][:,:,:]
    txt='{0} dimension [time x lat x lon]: {1}'.format(variabs[0],var.shape)
    #print(fh.variables)
    dates=num2date(time,time_units)
    fh.close()
       
    #print('\n'+txt)
       
    return var, var_units, lat, lon, dates, time_units


def save_N_2Dfields(lats,lons,variab,varname,varunits,ofile):
    '''
    GOAL
        Save var in ofile netCDF file
    USAGE
        save a number N of 2D fields [latxlon]
    '''
    try:
        os.remove(ofile) # Remove the outputfile
    except OSError:
        pass
    dataset = Dataset(ofile, 'w', format='NETCDF4_CLASSIC')
    #print(dataset.file_format)
    
    num = dataset.createDimension('num', variab.shape[0])
    lat = dataset.createDimension('lat', variab.shape[1])
    lon = dataset.createDimension('lon', variab.shape[2])
    
    # Create coordinate variables for 3-dimensions
    num = dataset.createVariable('num', np.int32, ('num',))
    lat = dataset.createVariable('lat', np.float32, ('lat',))
    lon = dataset.createVariable('lon', np.float32, ('lon',))
    # Create the actual 3-d variable
    var = dataset.createVariable(varname, np.float64,('num','lat','lon'))
    
    #print('variable:', dataset.variables[varname])
    
    #for varn in dataset.variables.keys():
    #    print(varn)
    # Variable Attributes
    lat.units='degree_north'
    lon.units='degree_east'
    var.units = varunits
    
    num[:]=np.arange(variab.shape[0])    
    lat[:]=lats
    lon[:]=lons
    var[:,:,:]=variab

    dataset.close()

    #----------------------------------------------------------------------------------------
    print('The {0} 2D fields [num x lat x lon] are saved as \n{1}'.format(variab.shape[0], ofile))
    print('__________________________________________________________')
    #----------------------------------------------------------------------------------------


def read_N_2Dfields(ifile):
    '''
    GOAL
        read a number N of 2D fields [latxlon]
    USAGE
        var, lat, lon, dates = read_N_2Dfields(filename)
    '''
    fh = Dataset(ifile, mode='r')
    variabs=[]
    for variab in fh.variables:
        variabs.append(variab)
    #print('The variables in the nc file are: ', variabs)
    
    num         = fh.variables['num'][:]
    if ('lat' in variabs):
        lat     = fh.variables['lat'][:]
    elif ('latitude' in variabs):
        lat     = fh.variables['latitude'][:]
    if ('lon' in variabs):
        lon     = fh.variables['lon'][:]
    elif ('longitude' in variabs):
        lon     = fh.variables['longitude'][:]
    var         = fh.variables[variabs[3]][:,:,:]
    var_units   = fh.variables[variabs[3]].units
    txt='{0} dimension [num x lat x lon]: {1}'.format(variabs[3],var.shape)
    #print(fh.variables)
    fh.close()
       
    #print('\n'+txt)
       
    return var, var_units, lat, lon

