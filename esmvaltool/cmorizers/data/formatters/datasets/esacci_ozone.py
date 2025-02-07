import logging
import re
from pathlib import Path
import iris
import numpy as np
import os
import glob
from datetime import datetime

from ...utilities import (
    fix_var_metadata,
    fix_dim_coordnames,
    fix_bounds,
    fix_coords,
    save_variable,
    set_global_atts
)

logger = logging.getLogger(__name__)


def _calculate(cube, short_name):
    """Perform variable-specific calculations."""
    if short_name == 'o3':  # Ozone mole fraction
        r = 8.31446261815324  # Ideal gas constant (J mol-1 K-1)
        try:
            t_cube = cube.extract('air_temperature')
            p_cube = cube.extract('air_pressure')

            if t_cube and p_cube and t_cube.data.size > 0 and p_cube.data.size > 0:  # Check if cubes exist AND contain data
                t = t_cube.data[0]
                p = p_cube.data[0]
            else:
                raise iris.exceptions.ConstraintNotFoundError()  # Raise the exception to use the default values
        except iris.exceptions.ConstraintNotFoundError:
            logger.warning("Temperature or pressure variables not found or empty. Using default values.")
            t = 273.15  # Default temperature (K)
            p = 101325.0  # Default pressure (Pa)

        air_density = p / (r * t)  # mol m-3
        o3_concentration = cube.core_data()  # mol m-3
        mixing_ratio = o3_concentration / air_density
        mole_fraction = mixing_ratio / (1.0 + mixing_ratio)
        cube.data = mole_fraction
        cube.units = 'mol mol-1'

    elif short_name == 'toz':  # Total ozone column (DU or m)
        # Check if the units are in Dobson Units (DU)
        if cube.units == 'DU':
            # Convert DU to m (1 DU = 2.24115E-5 m)
            cube.data = cube.data * 2.24115e-5
            cube.units = 'm'
        elif cube.units == 'molecules cm-2':
            cube.data = cube.data * 1.e-4 #molecules m-2
            cube.units = 'molecules m-2'
        elif cube.units == 'mol m-2':
            pass #no conversion is needed
        else:
            raise ValueError(f"Unexpected units for total ozone column: {cube.units}")

    return cube

def _extract_variable(short_name, var, cfg, filename, out_dir):
    """Extract variable, add time coordinate, and scalar longitude."""

    mip = var['mip']
    cmor_info = cfg['cmor_table'].get_variable(mip, short_name)

    cube = iris.load_cube(filename, var['raw'])

    if not cube:
        raise ValueError(f"Variable '{short_name}' not found in file {filename}")

    cube = _calculate(cube, short_name)

    logger.info(f"Checking CMOR info for {short_name}: {cmor_info}")
    if cmor_info is None:
        raise ValueError(f"CMOR info for {short_name} in MIP {mip} not found!")

# Check if the CMOR unit is meters and the current unit is mol m-2
    if cmor_info.units == 'm' and cube.units == 'mol m-2':
# Convert mol m-2 to m (assuming total ozone column)
        cube.data = cube.data * 1.0 # Replace 1.0 with the correct factor.
        cube.units = 'm'
        logger.warning("Converting mol m-2 to m for total ozone column.")

    cube.convert_units(cmor_info.units)  # Convert to CMOR units
    fix_dim_coordnames(cube)

    match = re.search(r'(\d{4})(\d{2})-', filename)  # Capture YYYYMM before the hyphen
    if match:
        year = int(match.group(1))
        month = int(match.group(2))
        day = 15  # Mid-month
        time = datetime(year, month, day, 12)  # Midday
        print(f"Filename: {filename}, Extracted Year: {year}, Month: {month}")  # Debug print
    else:
        raise ValueError(f"Could not extract date (YYYYMM) from filename: {filename}")

# Time coordinate handling:
    if 'time' not in (coord.standard_name for coord in cube.coords()):
        time_coord = iris.coords.DimCoord(np.array([time.timestamp()]),
                                           standard_name='time', units='seconds since 1970-01-01 00:00:00 UTC')
        #print (cube)
        #cube.add_dim_coord(time_coord, None)
    else:
        time_coord = cube.coord('time')
        time_coord.points = np.array([time.timestamp()])
        time_coord.units = 'seconds since 1970-01-01 00:00:00 UTC'
        if time_coord.bounds is not None:
            time_coord.bounds = None
        logger.warning("Time coordinate already exists. Updating its value.")

# Add scalar longitude coordinate (ALTERNATIVE SOLUTION - NO new_axis):
    if 'longitude' not in (coord.standard_name for coord in cube.coords()):
        # Create the longitude coordinate:
        lon_coord = iris.coords.DimCoord(180, standard_name='longitude', units='degrees_east')
        lon_coord.bounds = [0, 360]

        # Create the new cube's dim_coords_and_dims list:
        new_dim_coords_and_dims = [(lon_coord, 0)]  # Longitude at index 0

        # Add the existing dimensions and coordinates:
        for i, dim in enumerate(cube.dimensions):
            coord = cube.coord(dim)  # Get the coordinate for this dimension
            new_dim_coords_and_dims.append((coord, i + 1))  # Add to the list, incrementing index

    else:
        lon_coord = cube.coord('longitude')
        lon_coord.points = 180
        lon_coord.bounds = [0, 360]
        logger.warning("Longitude coordinate already exists. Updating its value.")

    # Handle lat/lon and pressure coordinates:
    for coord in cube.coords():
        if coord.var_name == 'air_pressure':
            coord.standard_name = 'air_pressure'
            coord.var_name = 'plev'
            fix_bounds(cube, coord)
        elif 'lat' in coord.var_name.lower():
            cube.coord(coord.var_name).standard_name = 'latitude'
            fix_bounds(cube, coord)
        elif 'lon' in coord.var_name.lower() and coord.standard_name != 'longitude':  # Avoid overwriting scalar lon
            cube.coord(coord.var_name).standard_name = 'longitude'
            fix_bounds(cube, coord)

    fix_var_metadata(cube, cfg)
    set_global_atts(cube, cfg['attributes'])
    save_variable(cube, short_name, out_dir, cfg['attributes'], unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization process."""

    if not start_date:
        start_date = datetime(1984, 1, 1)
        end_date = datetime(2023, 12, 31)

    for var_name, var in cfg['variables'].items():
        if start_date.year > var['max_year']:
            continue
        if end_date.year < var['min_year']:
            continue

        start_year = max(start_date.year, var['min_year'])
        end_year = min(end_date.year, var['max_year'])

        for year in range(start_year, end_year + 1):
            for month in range(1, 13):
                date_str = f"{year}{month:02}"  # YYYYMM format

                # Create glob pattern, handling potential prefixes:
                fname_pattern = f"*{date_str}-{var['filename']}"  # * before and after date
                fname = glob.glob(os.path.join(in_dir, fname_pattern))

                if not fname:
                    # Try with the old pattern if the new one fails:
                    fname_pattern = var['filename'].replace("_DATE_", date_str)  # Original pattern
                    fname = glob.glob(os.path.join(in_dir, fname_pattern))
                    if not fname:
                        logger.warning(f"No file found for {var_name} in {year}-{month:02}")
                        continue

                filename = fname[0]  # Take the first matching file (if any)
                logger.info(f"CMORizing variable '{var_name}' from file '{filename}'")
                _extract_variable(var_name, var, cfg, filename, out_dir)
