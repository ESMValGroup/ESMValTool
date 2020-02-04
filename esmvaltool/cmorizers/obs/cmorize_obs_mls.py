"""ESMValTool CMORizer for MLS data.

Tier
    Tier 3: restricted dataset.

Source
    https://disc.gsfc.nasa.gov/datasets/ML2RHI_004/summary
    https://disc.gsfc.nasa.gov/datasets/ML2T_004/summary

Last access
    20200203

Download and processing instructions
    Select "Data Access" -> "Subset/Get Data" -> "Get Data" and follow the
    "Instructions for downloading". All *.he5 files need to be saved in the
    $RAWOBS/Tier3/MLS directory, where $RAWOBS refers to the RAWOBS directory
    defined in the user configuration file. The temperature fields are
    necessary for quality control of the RHI data (see Data Quality Document
    for MLS for more information).
    A registration is required for downloading the data.

"""

import glob
import logging
import os
from datetime import datetime

import iris
import iris.coord_categorisation
import netCDF4 as nc
import numpy as np
import pandas as pd
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


ALL_LATS = np.linspace(-90.0, 90.0, 181)
ALL_LONS = np.linspace(-180.0, 180.0, 361)
LAT_COORD = iris.coords.DimCoord(ALL_LATS,
                                 var_name='lat',
                                 standard_name='latitude',
                                 long_name='Latitude',
                                 units='degrees')
LON_COORD = iris.coords.DimCoord(ALL_LONS,
                                 var_name='lon',
                                 standard_name='longitude',
                                 long_name='Longitude',
                                 units='degrees')
TIME_UNITS = Unit('days since 1850-01-01 00:00:00', calendar='standard')


def _get_cube(variable, nc_rhi, nc_loc, mask, time):
    """Get :class:`iris.cube.Cube` with correct data."""
    data = np.ma.array(nc_rhi.variables[variable][:], mask=mask)
    pressure = nc_loc.variables['Pressure'][:]
    lat = nc_loc.variables['Latitude'][:]
    lon = nc_loc.variables['Longitude'][:]

    # Place on 0.5x0.5 degree grid
    lat = _round(lat)
    lon = _round(lon)

    # Iterate over pressure levels
    gridded_data = []
    for (p_idx, _) in enumerate(pressure):
        p_data = data[:, p_idx]
        data_frame = pd.DataFrame({'lat': lat, 'lon': lon, 'data': p_data})

        # Create daily-mean gridded data using pivot table
        data_frame = pd.pivot_table(data_frame,
                                    values='data',
                                    index='lat',
                                    columns='lon',
                                    aggfunc=np.mean,
                                    dropna=False)
        data_frame = data_frame.reindex(index=ALL_LATS, columns=ALL_LONS)
        gridded_data.append(data_frame.values)

    # Create cube
    time_coord = iris.coords.DimCoord(TIME_UNITS.date2num(time),
                                      var_name='time',
                                      standard_name='time',
                                      long_name='Time',
                                      units=TIME_UNITS)
    pressure_coord = iris.coords.DimCoord(pressure,
                                          var_name='plev',
                                          standard_name='air_pressure',
                                          long_name='Pressure',
                                          units='hPa')
    coord_spec = [
        (time_coord, 0),
        (pressure_coord, 1),
        (LAT_COORD, 2),
        (LON_COORD, 3),
    ]
    gridded_data = np.expand_dims(np.array(gridded_data), 0)
    cube = iris.cube.Cube(gridded_data, dim_coords_and_dims=coord_spec,
                          units='%')
    return cube


def _get_file_attributes(filename):
    """Get global file attributes."""
    dataset = nc.Dataset(filename, mode='r')
    add_info = dataset.groups['HDFEOS'].groups['ADDITIONAL']
    attrs = add_info.groups['FILE_ATTRIBUTES']
    return {key: attrs.getncattr(key) for key in attrs.ncattrs()}


def _get_files(var_config, in_dir, cfg):
    """Get all files for a given variable."""
    logger.info("Searching files")
    ext = cfg['extension']
    filename_rhi = var_config['file_pattern'].format(var='RHI')
    filename_t = var_config['file_pattern'].format(var='Temperature')
    file_pattern_rhi = f'{filename_rhi}*.{ext}'
    # TODO: remove :20
    files_rhi = glob.glob(os.path.join(in_dir, file_pattern_rhi))[:20]
    prefix = os.path.join(in_dir, filename_rhi)
    files_t = []
    for file_ in files_rhi:
        suffix = file_.replace(prefix, '')
        suffix = suffix.replace(f'.{ext}', '')
        suffix = suffix.split('_')[1]
        file_pattern_t = os.path.join(in_dir, filename_t)
        file_pattern_t += f'*_{suffix}.{ext}'
        file_t = glob.glob(file_pattern_t)
        if len(file_t) != 1:
            raise ValueError(
                f"Expected exactly one corresponding temperature file for RHI "
                f"file {file_}, found {len(file_t):d} ({file_t})")
        files_t.append(file_t[0])
    files = zip(files_rhi, files_t)
    logger.info("Found %d files", len(files_rhi))
    return files


def _get_mask(nc_rhi, nc_t, nc_loc):
    """Remove invalid data (see Data Quality Document of MLS)."""
    mask = np.full(nc_rhi.variables['L2gpValue'][:].shape, False)

    # Pressure range (accept only 320 hPa and above)
    new_mask = np.where(nc_loc.variables['Pressure'][:] < 320, False, True)
    new_mask = np.broadcast_to(new_mask, mask.shape)
    mask |= new_mask

    # Status (accept only even status flags)
    status = np.expand_dims(nc_rhi.variables['Status'][:], -1)
    status = np.broadcast_to(status, mask.shape)
    mask |= np.where(status % 2, True, False)

    # Precision of RHI (accept only positive numbers)
    precision = nc_rhi.variables['L2gpPrecision'][:]
    mask |= np.where(precision > 0, False, True)

    # Quality of RHI (accept only values greater than 1.45)
    quality_rhi = np.expand_dims(nc_rhi.variables['Quality'][:], -1)
    quality_rhi = np.broadcast_to(quality_rhi, mask.shape)
    mask |= np.where(quality_rhi > 1.45, False, True)

    # Quality of Temperature (accept only values greater than 0.2/0.9)
    pressure_greater_90 = np.where(
        nc_loc.variables['Pressure'][:] > 90, True, False)
    quality_t = np.expand_dims(nc_t.variables['Quality'][:], -1)
    quality_t = np.broadcast_to(quality_t, mask.shape)
    new_mask = np.full(mask.shape, False)
    new_mask[:, pressure_greater_90] = np.where(
        quality_t[:, pressure_greater_90] > 0.9, False, True)
    new_mask[:, ~pressure_greater_90] = np.where(
        quality_t[:, ~pressure_greater_90] > 0.2, False, True)
    mask |= new_mask

    # Convergence of RHI (accept only values smaller than 2.0)
    convergence_rhi = np.expand_dims(nc_rhi.variables['Convergence'][:], -1)
    convergence_rhi = np.broadcast_to(convergence_rhi, mask.shape)
    mask |= np.where(convergence_rhi < 2.0, False, True)

    # Convergence of Temperature (accept only values smaller than 1.03)
    convergence_t = np.expand_dims(nc_t.variables['Convergence'][:], -1)
    convergence_t = np.broadcast_to(convergence_t, mask.shape)
    mask |= np.where(convergence_t < 1.03, False, True)

    return mask


def _open_nc_file(filename, variable):
    """Open :class:`netCDF4.Dataset`."""
    dataset = nc.Dataset(filename, mode='r')
    swaths = dataset.groups['HDFEOS'].groups['SWATHS']
    var = swaths.groups[variable]
    return (var.groups['Data Fields'], var.groups['Geolocation Fields'])


def _save_cube(cube, cmor_info, attrs, out_dir):
    """Save :class:`iris.cube.Cube`."""
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _round(array):
    """Round floating point array to next 0.5."""
    array *= 2.0
    array = np.around(array)
    return array / 2.0


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        files = _get_files(var_info, in_dir, cfg)

        # Extract cubes
        cubes = iris.cube.CubeList()
        for (filename_rhi, filename_t) in files:
            logger.info("Processing %s", filename_rhi)
            attributes = _get_file_attributes(filename_rhi)

            # Read files
            (nc_rhi, nc_loc) = _open_nc_file(filename_rhi, 'RHI')
            (nc_t, _) = _open_nc_file(filename_t, 'Temperature')

            # Get mask to remove invalid values
            mask = _get_mask(nc_rhi, nc_t, nc_loc)

            # Extract data
            time = datetime(year=attributes['GranuleYear'],
                            month=attributes['GranuleMonth'],
                            day=attributes['GranuleDay'],
                            hour=12)
            cube = _get_cube(var_info['raw_var'], nc_rhi, nc_loc, mask, time)
            cubes.append(cube)

        # Create final cube and save it
        final_cube = cubes.concatenate_cube()
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        _save_cube(final_cube, cmor_info, glob_attrs, out_dir)
