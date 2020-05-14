# pylint: disable=invalid-name
"""ESMValTool CMORizer for PHC data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://psc.apl.washington.edu/nonwp_projects/PHC/Data3.html

Last access
   20190131

Go to `DOWNLOAD DATA (NetCDF)` and download the `ANNUAL` fields
for both `TEMPERATURE` and `SALINITY`.

"""
import logging
import os
from collections import OrderedDict
import iris
import xarray as xr
import numpy as np
import seawater as sw

from esmvaltool.cmorizers.obs.utilities import (
    fix_coords, fix_var_metadata, save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def save_fx_variable(cube, var, out_dir, attrs):
    """Saver function for fx variable."""
    file_name = '_'.join([
        attrs['project_id'], attrs['dataset_id'], attrs['modeling_realm'],
        attrs['version'], attrs['mip'], var
    ]) + '.nc'
    file_path = os.path.join(out_dir, file_name)
    iris.save(cube, file_path, fill_value=1e20)


def _fix_fx_areacello(xr_time, var):
    """Specific data fix for areacello."""
    cube = xr_time.salt.to_iris()
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    grid_areas_xr = xr.DataArray(grid_areas[0, 0, :, :],
                                 coords={
                                     'lat': xr_time.temp.coords['lat'],
                                     'lon': xr_time.temp.coords['lon']
    },
        dims=['lat', 'lon'],
        name=var)
    grid_areas_xr.attrs = OrderedDict([('cell_area', 'Ocean Grid-Cell Area'),
                                       ('units', 'm2')])
    cube = grid_areas_xr.to_iris()
    return cube.copy()


def _fix_data(xr_time, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    if var == 'thetao':
        depth3d = np.zeros(xr_time.temp.shape[1:])
        for i in range(xr_time.depth.shape[0]):
            depth3d[i, :, :] = xr_time.depth[i]
        ptemp = sw.ptmp(xr_time.salt[0, :], xr_time.temp[0, :], depth3d)
        ptemp = np.expand_dims(ptemp, axis=0)
        temp_new = xr.DataArray(ptemp + 273.15,
                                coords={
                                    'time': xr_time.temp.coords['time'],
                                    'depth': xr_time.temp.coords['depth'],
                                    'lat': xr_time.temp.coords['lat'],
                                    'lon': xr_time.temp.coords['lon'],
                                },
                                dims=['time', 'depth', 'lat', 'lon'])

        temp_new.attrs = OrderedDict([('standard_name',
                                       'sea_water_potential_temperature'),
                                      ('units', 'K')])
        cube = temp_new.to_iris()
        return cube.copy()
    elif var == 'so':
        cube = xr_time.salt.to_iris()
        return cube.copy()
    elif var == 'areacello':
        cube = _fix_fx_areacello(xr_time, var)
        return cube.copy()
    return None


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    xr_file = xr.open_dataset(raw_info['file'])
    xr_time = xr_file.expand_dims('time')
    xr_time = xr_time.assign_coords(time=[1])
    xr_time.time.attrs = OrderedDict(
        [('standard_name', 'time'), ('units', 'days since 1950-1-1 00:00:00')])

    cube = _fix_data(xr_time, var)
    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    set_global_atts(cube, attrs)
    print(out_dir)
    if var != "areacello":
        save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])
    else:
        save_fx_variable(cube, var, out_dir, attrs)


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        inpfile = os.path.join(in_dir, vals['file'])
        logger.info("CMORizing var %s from file %s", var, inpfile)
        var_info = cmor_table.get_variable(vals['mip'], var)
        raw_info = {'name': vals['raw'], 'file': inpfile}
        glob_attrs['mip'] = vals['mip']
        extract_variable(var_info, raw_info, out_dir, glob_attrs)
