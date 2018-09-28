"""Diagnostic to select grid points within a shapefile."""
from copy import deepcopy
import logging
import os
import sys
import numpy as np
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import iris
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

def main(cfg):
    """Calculate drought indices."""
    for filename, attributes in cfg['input_data'].items():
        logger.info("Processing variable %s from dataset %s",
                    attributes['standard_name'], attributes['dataset'])
        logger.debug("Loading %s", filename)
        cube = iris.load_cube(filename)
        drymaxcube, fqthcube = droughtindex(cube, cfg)
        name = os.path.splitext(os.path.basename(filename))[0]
        # Write to file
        if cfg['write_netcdf']:
            path = os.path.join(cfg['work_dir'],
                                name + '_drymax.nc',)
            write_netcdf(path, drymaxcube, cfg)
            #iris.save(drymaxcube, target = path)
            path = os.path.join(cfg['work_dir'],
                                name + '_dryfreq.nc',)
            write_netcdf(path, fqthcube, cfg)
            #iris.save(fqthcube, target = path)
        if cfg['write_plots'] and cfg.get('quickplot'):
            path = os.path.join(
                cfg['plot_dir'],
                name + '.' + cfg['output_file_type'],
            )
            logger.debug("Plotting analysis results to %s", path)
            quickplot(drymaxcube, filename=path, **cfg['quickplot'])


def droughtindex(cube, cfg):
    """Calculates drought stats."""
    if cfg['dryindex'] == 'cdd':
        plim = float(cfg['plim']) / 86400. # units of kg m-2 s-1
        frlim = float(cfg['frlim'])
        pr = deepcopy(cube.data)
        pr[cube.data<plim] = 1
        pr[cube.data>=plim] = 0
        cube.data[0,:,:] = pr[0,:,:]
        for t in range(1,cube.data.shape[0]):
            cube.data[t,:,:] = ((pr[t,:,:] + cube.data[t-1,:,:]) *
                                pr[t,:,:])
        dif = cube.data[0:-1,:,:]-cube.data[1:cube.data.shape[0],:,:]
        wh = np.where(dif!=cube.data[0:-1])
        cube.data[wh] = 0
        # Longest consecutive period
        drymaxcube = cube.collapsed('time',iris.analysis.MAX)
        #drymaxcube.standard_name = 'consecutive_dry_days'
        drymaxcube.long_name = 'Consecutive dry days is the greatest number of consecutive days per time period with daily precipitation amount below ' + str(cfg['plim']) + 'mm.'
        #drymaxcube.units.origin = '1'
        whth = np.where(cube.data>frlim)
        cube.data = cube.data*0
        cube.data[whth] = 1
        fqthcube = cube.collapsed('time', iris.analysis.SUM)
        print(fqthcube.shape)

    return drymaxcube, fqthcube

def write_netcdf(path, cube, cfg):
    """Write results to a netcdf file."""
    # wgtmet = cfg['wgtmet']
    ncout = Dataset(path, mode='w')
    ncout.createDimension('time', cube.coord('time').shape[0])
    ncout.createDimension('longitude', cube.coord('longitude').shape[0])
    ncout.createDimension('latitude', cube.coord('latitude').shape[0])
    times = ncout.createVariable('time', 'f8', ('time'), zlib=True)
    times.setncattr_string('standard_name', cube.coord('time').standard_name)
    times.setncattr_string('long_name', cube.coord('time').long_name)
    times.setncattr_string('calendar', cube.coord('time').units.calendar)
    times.setncattr_string('units', cube.coord('time').units.origin)
    lon = ncout.createVariable('longitude', 'f8', ('longitude'), zlib=True)
    lon.setncattr_string('standard_name', cube.coord('longitude').standard_name)
    lon.setncattr_string('long_name', cube.coord('longitude').long_name)
    lon.setncattr_string('units', cube.coord('longitude').units.origin)
    lat = ncout.createVariable('latitude', 'f8', ('latitude'), zlib=True)
    lat.setncattr_string('standard_name', cube.coord('latitude').standard_name)
    lat.setncattr_string('long_name', cube.coord('latitude').long_name)
    lat.setncattr_string('units', cube.coord('latitude').units.origin)
    data = ncout.createVariable(cube.var_name, 'f4',
                                ('time', 'latitude', 'longitude'),
                                zlib=True)
    data.setncattr_string('standard_name', cube.standard_name)
    data.setncattr_string('long_name', cube.long_name)
    data.setncattr_string('units', "1") #cube.units.origin)
    for key, val in cube.metadata[-2].items():
        ncout.setncattr_string(key, val)
    times[:] = cube.coord('time').points[0]
    lon[:] = cube.coord('longitude').points
    lat[:] = cube.coord('latitude').points
    print(cube.data.shape)
    data[:] = cube.data
    ncout.close()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
