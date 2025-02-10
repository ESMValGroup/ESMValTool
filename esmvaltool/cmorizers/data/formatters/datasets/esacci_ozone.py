"""ESMValTool CMORizer for ESACCI-OZONE data.

Tier
    Tier 2: other freely-available dataset.

Source

Last access
    20250107

Download and processing instructions
    Download the data from:

    Put all files under a single directory (no subdirectories with years).
    in ${RAWOBS}/Tier2/ESACCI-OZONE

"""
import glob
import logging
from datetime import datetime
import os
import iris
import iris.util
from cf_units import Unit
from esmvalcore.preprocessor import (
    concatenate
)

from ...utilities import (
    fix_var_metadata,
    fix_coords,
    save_variable,
    set_global_atts
)

logger = logging.getLogger(__name__)


def _convert_units(cubes, short_name, var):
    """Perform variable-specific calculations."""
    cube = cubes.extract_cube(var['raw'])
    if short_name == 'o3':  # Ozone mole fraction
        r = 8.31446261815324  # Ideal gas constant (J mol-1 K-1)
        t_cube = cubes.extract_cube('air_temperature')
        p_cube = cubes.extract_cube('air_pressure')

        air_mol_concentration = p_cube / (r * t_cube)  # mol m-3
        mixing_ratio = cube / air_mol_concentration
        cube = mixing_ratio / (1.0 + mixing_ratio)
        cube.units = 'mol mol-1'

    elif short_name == 'toz':  # Total ozone column (DU or m)
        cube = cube * 2241.399
        cube.units = 'DU'
        cube = cube * 2.24115e-5
        cube.units = 'm'
    return cube


def _extract_variable(short_name, var, cfg, filename, year, month, out_dir):
    """Extract variable, add time coordinate, and scalar longitude."""

    mip = var['mip']
    cmor_info = cfg['cmor_table'].get_variable(mip, short_name)

    cubes = iris.load(filename)

    cube = _convert_units(cubes, short_name, var)

    logger.info("Checking CMOR info for %s: %s", short_name, cmor_info)
    if cmor_info is None:
        raise ValueError(f"CMOR info for {short_name} in MIP {mip} not found!")

    year = year
    month = month
    day = 15  # Mid-month
    time_units = Unit("days since 1950-01-01")
    time_points = time_units.date2num(datetime(year, month, day))
    print(f"Filename: {filename}, Extracted Year: {year}, Month: {month}")

    # Add time coordinate to cube.
    time_coord = iris.coords.DimCoord(
        time_points,
        var_name='time',
        standard_name='time',
        long_name='time',
        units=time_units)
    time_coord.guess_bounds(monthly=True)
    cube.add_aux_coord(time_coord, ())
    cube = iris.util.new_axis(cube, time_coord)

    # Add longitude coordinate to cube only for o3.
    if short_name == 'o3':
        lon_coord = iris.coords.DimCoord(
            [180.0],
            bounds=[[0.0, 360.0]],
            var_name='lon',
            standard_name='longitude',
            long_name='longitude',
            units='degrees_east')
        cube.add_aux_coord(lon_coord, ())
        cube = iris.util.new_axis(cube, lon_coord)
        cube.transpose([1, 3, 2, 0])
    fix_var_metadata(cube, cmor_info)
    cube = fix_coords(cube)
    set_global_atts(cube, cfg['attributes'])
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization process."""
    glob_attrs = cfg['attributes']
    if not start_date:
        start_date = datetime(1984, 1, 1)
    if not end_date:
        end_date = datetime(2023, 12, 31)

    for var_name, var in cfg['variables'].items():

        all_data_cubes = []
        glob_attrs['mip'] = var['mip']
        for year in range(start_date.year, end_date.year + 1):
            for month in range(1, 13):
                date_str = f"{year}{month:02}"  # YYYYMM format

                # Create glob pattern, handling prefixes:
                fname_pattern = f"*{date_str}-{var['filename']}"
                fname = glob.glob(os.path.join(in_dir, fname_pattern))

                if not fname:
                    fname = glob.glob(os.path.join(in_dir, fname_pattern))
                    if not fname:
                        raise ValueError(
                            "No file found for %s in %s-%02d", var_name, year,
                            month
                        )

                filename = fname[0]  # Take the first matching file (if any)
                logger.info("CMORizing variable '%s' from file '%s'", var_name,
                            filename)
                cube = _extract_variable(var_name, var, cfg, filename, year,
                                         month, out_dir)
                all_data_cubes.append(cube)
        final_cube = concatenate(all_data_cubes)
        save_variable(
            final_cube, var_name, out_dir, glob_attrs,
            unlimited_dimensions=['time'])
