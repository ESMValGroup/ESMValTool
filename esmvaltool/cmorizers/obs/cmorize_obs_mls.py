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

ALL_LATS = np.linspace(-90.0, 90.0, 91)
ALL_LONS = np.linspace(-180.0, 180.0, 181)
LAT_COORD = iris.coords.DimCoord(ALL_LATS,
                                 var_name='lat',
                                 standard_name='latitude',
                                 long_name='latitude',
                                 units='degrees')
LON_COORD = iris.coords.DimCoord(ALL_LONS,
                                 var_name='lon',
                                 standard_name='longitude',
                                 long_name='Longitude',
                                 units='degrees')
TIME_UNITS = Unit('days since 1850-01-01 00:00:00', calendar='standard')


def _cut_cube(cube, var_info):
    """Cut cube if desired."""
    if 'cut_levels_outside' in var_info:
        lims = var_info['cut_levels_outside']
        constraint = iris.Constraint(
            air_pressure=lambda cell: lims[0] < cell < lims[1])
        cube = cube.extract(constraint)
    return cube


def _extract_cubes(cfg, files):
    """Extract cubes from files."""
    cubes_dict = {}
    for (var, _) in cfg['variables'].items():
        logger.info("Found variable '%s'", var)
        cubes_dict[var] = iris.cube.CubeList()
    for (filename_rhi, filename_t) in files:
        logger.info("Processing %s", filename_rhi)

        # Read files
        (nc_rhi, nc_loc) = _open_nc_file(filename_rhi, 'RHI')
        (nc_t, _) = _open_nc_file(filename_t, 'Temperature')

        # Get cubes for all desired variables
        for (var, var_info) in cfg['variables'].items():
            (gridded_data, time,
             pressure) = _get_gridded_data(var_info['raw_var'], nc_rhi, nc_loc,
                                           nc_t, filename_rhi)
            cubes_dict[var].append(_get_cube(gridded_data, time, pressure))

    # Create final cube and return it
    for (var, cubes) in cubes_dict.items():
        var_info = cfg['variables'][var]
        cubes_dict[var] = cubes.concatenate_cube()
        cubes_dict[var] = _cut_cube(cubes_dict[var], var_info)

    return cubes_dict


def _get_cube(gridded_data, time, pressure):
    """Get :class:`iris.cube.Cube` with correct data."""
    time_coord = iris.coords.DimCoord(TIME_UNITS.date2num(time),
                                      var_name='time',
                                      standard_name='time',
                                      long_name='time',
                                      units=TIME_UNITS)
    pressure_coord = iris.coords.DimCoord(pressure,
                                          var_name='plev',
                                          standard_name='air_pressure',
                                          long_name='pressure',
                                          units='hPa')
    coord_spec = [
        (time_coord, 0),
        (pressure_coord, 1),
        (LAT_COORD, 2),
        (LON_COORD, 3),
    ]
    cube = iris.cube.Cube(gridded_data,
                          dim_coords_and_dims=coord_spec,
                          units='%')
    return cube


def _get_file_attributes(filename):
    """Get global file attributes."""
    dataset = nc.Dataset(filename, mode='r')
    add_info = dataset.groups['HDFEOS'].groups['ADDITIONAL']
    attrs = add_info.groups['FILE_ATTRIBUTES']
    return {key: attrs.getncattr(key) for key in attrs.ncattrs()}


def _get_files(in_dir, cfg):
    """Get all files for a given variable."""
    logger.info("Searching files")
    if 'year' in cfg:
        year = cfg['year']
        logger.info("Only considering year %d", year)
    else:
        year = None
    ext = cfg['extension']
    filename_rhi = cfg['file_pattern'].format(var='RHI')
    filename_t = cfg['file_pattern'].format(var='Temperature')
    if year is None:
        file_pattern_rhi = f'{filename_rhi}*.{ext}'
    else:
        file_pattern_rhi = f'{filename_rhi}*_{year:d}*.{ext}'
    files_rhi = glob.glob(os.path.join(in_dir, file_pattern_rhi))
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


def _get_gridded_data(variable, nc_rhi, nc_loc, nc_t, filename):
    """Get gridded data."""
    file_attrs = _get_file_attributes(filename)

    # Extract coords
    time = datetime(year=file_attrs['GranuleYear'],
                    month=file_attrs['GranuleMonth'],
                    day=file_attrs['GranuleDay'],
                    hour=12)
    pressure = nc_loc.variables['Pressure'][:]
    lat = nc_loc.variables['Latitude'][:]
    lon = nc_loc.variables['Longitude'][:]

    # Extract data
    data = np.ma.array(nc_rhi.variables[variable][:],
                       mask=_get_mask(nc_rhi, nc_t, nc_loc))

    # For version 4.20, remove last four profiles (see Data Quality Document)
    if file_attrs['PGEVersion'] == 'V04-20':
        data = data[:-4]
        lat = lat[:-4]
        lon = lon[:-4]

    # Place on 1x1 degree grid
    lat = np.around(lat)
    lon = np.around(lon)

    # Iterate over pressure levels
    gridded_data = []
    for (p_idx, _) in enumerate(pressure):
        data_frame = pd.DataFrame({
            'lat': lat,
            'lon': lon,
            'data': data[:, p_idx],
        })

        # Create daily-mean gridded data using pivot table
        data_frame = pd.pivot_table(data_frame,
                                    values='data',
                                    index='lat',
                                    columns='lon',
                                    aggfunc=np.mean,
                                    dropna=False)
        data_frame = data_frame.reindex(index=ALL_LATS, columns=ALL_LONS)
        gridded_data.append(data_frame.values)
    gridded_data = np.expand_dims(np.array(gridded_data), 0)

    return (gridded_data, time, pressure)


def _get_mask(nc_rhi, nc_t, nc_loc):
    """Remove invalid data (see Data Quality Document of MLS)."""
    mask = np.full(nc_rhi.variables['L2gpValue'][:].shape, False)

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
    pressure_greater_90 = np.where(nc_loc.variables['Pressure'][:] > 90, True,
                                   False)
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
    cube.coord('air_pressure').convert_units('Pa')
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    files = _get_files(in_dir, cfg)

    # Run the cmorization
    cubes_dict = _extract_cubes(cfg, files)

    # Save data
    for (var, cube) in cubes_dict.items():
        logger.info("Saving variable '%s'", var)
        var_info = cfg['variables'][var]
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        _save_cube(cube, cmor_info, glob_attrs, out_dir)
