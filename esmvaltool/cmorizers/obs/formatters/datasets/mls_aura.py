"""ESMValTool CMORizer for MLS-AURA data.

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
    $RAWOBS/Tier3/MLS-AURA directory, where $RAWOBS refers to the RAWOBS
    directory defined in the user configuration file. Apply this procedure to
    both links provided above. The temperature fields are necessary for quality
    control of the RHI data (see Data Quality Document for MLS-AURA for more
    information).
    A registration is required for downloading the data.

"""

import glob
import logging
import os
from datetime import datetime

import iris
import iris.coord_categorisation
import netCDF4
import numpy as np
import pandas as pd
from cf_units import Unit

from esmvaltool.cmorizers.obs import utilities as utils

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


def _extract_cubes(files_dict, cfg):
    """Extract cubes from files."""
    cubes_dict = _get_cubes_dict(files_dict, cfg)

    # Create final cubes and return it
    cube_dict = {}
    for (var, cubes) in cubes_dict.items():
        var_info = cfg['variables'][var]
        cube = cubes.concatenate_cube()
        cube = _cut_cube(cube, var_info)

        # Calculate monthly mean if desired
        if 'mon' in cfg['mip']:
            logger.info("Calculating monthly mean")
            iris.coord_categorisation.add_month_number(cube, 'time')
            iris.coord_categorisation.add_year(cube, 'time')
            cube = cube.aggregated_by(['month_number', 'year'],
                                      iris.analysis.MEAN)
            cube.remove_coord('month_number')
            cube.remove_coord('year')

        # Save cube
        cube_dict[var] = cube

    return cube_dict


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


def _get_cubes_dict(files_dict, cfg):
    """Get :obj:`dict` of :class:`iris.cube.CubeList`."""
    cubes_dict = {var: iris.cube.CubeList() for var in cfg['variables']}

    # Process files
    file_idx = 1
    for (filename_rhi, filename_t) in files_dict.values():
        logger.info("Processing file %5d/%5d [%s]", file_idx, len(files_dict),
                    filename_rhi)

        # Read files
        (nc_rhi, nc_loc) = _open_nc_file(filename_rhi, 'RHI')
        (nc_t, _) = _open_nc_file(filename_t, 'Temperature')

        # Get cubes for all desired variables
        for (var, var_info) in cfg['variables'].items():
            (gridded_data, time,
             pressure) = _get_gridded_data(var_info['raw_var'], nc_rhi, nc_loc,
                                           nc_t, filename_rhi)
            cubes_dict[var].append(_get_cube(gridded_data, time, pressure))
        file_idx += 1

    return cubes_dict


def _get_date(filename, variable, cfg):
    """Extract date from a filename."""
    file_pattern = cfg['file_pattern'].format(var=variable)
    filename = os.path.basename(filename)
    filename = os.path.splitext(filename)[0]
    filename = filename.replace(file_pattern, '')
    date = filename.split('_')[1]
    return date


def _get_file_attributes(filename):
    """Get global file attributes."""
    dataset = netCDF4.Dataset(filename, mode='r')
    add_info = dataset.groups['HDFEOS'].groups['ADDITIONAL']
    attrs = add_info.groups['FILE_ATTRIBUTES']
    return {key: attrs.getncattr(key) for key in attrs.ncattrs()}


def _get_files_single_var(variable, in_dir, cfg):
    """Get files for a single variable."""
    filename = cfg['file_pattern'].format(var=variable)
    ext = cfg['extension']
    file_pattern = f'{filename}*.{ext}'

    # Get all files
    files = glob.glob(os.path.join(in_dir, file_pattern))

    # Only accept certain years if desired
    if 'start_year' in cfg:
        start_year = cfg['start_year']
        logger.info("Only considering year %d and above", start_year)
    else:
        start_year = -np.inf
    if 'end_year' in cfg:
        end_year = cfg['end_year']
        logger.info("Only considering year %d and below", end_year)
    else:
        end_year = np.inf
    files_dict = {}
    for file_ in files:
        date = _get_date(file_, variable, cfg)
        year = int(date[:4])
        if start_year <= year <= end_year:
            files_dict[date] = file_

    return files_dict


def _get_files(in_dir, cfg):
    """Get all files for a given variable."""
    logger.info("Searching files")

    # Get file dictionaries
    files_dict_rhi = _get_files_single_var('RHI', in_dir, cfg)
    files_dict_t = _get_files_single_var('Temperature', in_dir, cfg)

    # Check if all files are available
    all_files = {}
    for (date, filename_rhi) in files_dict_rhi.items():
        if date not in files_dict_t:
            raise ValueError(f"No corresponding temperature file for RHI file "
                             f"{filename_rhi} found")
        all_files[date] = (filename_rhi, files_dict_t[date])
    logger.info("Found %d files", len(all_files))
    return all_files


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
            'data': data[:, p_idx].filled(np.nan),
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
    gridded_data = np.ma.masked_invalid(gridded_data)

    return (gridded_data, time, pressure)


def _get_mask(nc_rhi, nc_t, nc_loc):
    """Remove invalid data (see Data Quality Document of MLS-AURA)."""
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
    dataset = netCDF4.Dataset(filename, mode='r')
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
    glob_attrs['mip'] = cfg['mip']
    cmor_table = cfg['cmor_table']
    files_dict = _get_files(in_dir, cfg)

    # Run the cmorization
    cube_dict = _extract_cubes(files_dict, cfg)

    # Save data
    for (var, cube) in cube_dict.items():
        logger.info("Saving variable '%s'", var)
        var_info = cfg['variables'][var]
        if 'mip' in var_info:
            glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(glob_attrs['mip'], var)
        _save_cube(cube, cmor_info, glob_attrs, out_dir)
