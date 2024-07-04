"""ESMValTool CMORizer for ESACCI-SOILMOISTURE data.

Tier
    Tier 2: other freely-available dataset.

Source
    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/soil_moisture/data/

Last access
    20240626

Download and processing instructions
    Download the data from:
      daily_files/COMBINED/v08.1/
      ancillary/v08.1/
    Put all files under a single directory (no subdirectories with years).
    in ${RAWOBS}/Tier2/ESACCI-SOILMOISTURE

"""

import glob
import logging
import os
from datetime import datetime
import iris
from esmvalcore.cmor.fixes import get_time_bounds
from esmvalcore.preprocessor import concatenate, monthly_statistics
from cf_units import Unit

from ...utilities import (
    convert_timeunits,
    fix_var_metadata,
    fix_dim_coordnames,
    fix_bounds,
    save_variable,
    set_global_atts,
    roll_cube_data
)

logger = logging.getLogger(__name__)


def fix_coords_esacci_soilmoisture(cube,
                                   overwrite_time_bounds=True,
                                   overwrite_lon_bounds=True,
                                   overwrite_lat_bounds=True):
    """Fix coordinates to CMOR standards.

    Fixes coordinates eg time to have correct units, bounds etc;
    longitude to be CMOR-compliant 0-360deg; fixes some attributes
    and bounds - the user can avert bounds fixing by using supplied
    arguments; if bounds are None they will be fixed regardless.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube with coordinates to be fixed.

    overwrite_time_bounds: bool (optional)
        set to False not to overwrite time bounds.

    overwrite_lon_bounds: bool (optional)
        set to False not to overwrite longitude bounds.

    overwrite_lat_bounds: bool (optional)
        set to False not to overwrite latitude bounds.

    Returns
    -------
    cube: iris.cube.Cube
        data cube with fixed coordinates.
    """
    # first fix any completely missing coord var names
    fix_dim_coordnames(cube)
    # fix individual coords
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == 'time':
            logger.info("Fixing time...")
            cube.coord('time').convert_units(
                Unit('days since 1970-01-01T00:00:00+00:00',
                     calendar='proleptic_gregorian'))

        # fix longitude
        if cube_coord.var_name == 'lon':
            logger.info("Fixing longitude...")
            if cube_coord.ndim == 1:
                if cube_coord.points[0] < 0. and \
                        cube_coord.points[-1] < 181.:
                    cube_coord.points = \
                        cube_coord.points + 180.
                    cube.attributes['geospatial_lon_min'] = 0.
                    cube.attributes['geospatial_lon_max'] = 360.
                    nlon = len(cube_coord.points)
                    roll_cube_data(cube, nlon // 2, -1)
        fix_bounds(cube, cube_coord)

    return cube


def extract_variable(var_info, raw_info, attrs, year):
    """Extract variables."""
    rawvar = raw_info['name']
    constraint = iris.Constraint(name=rawvar)
    if rawvar == 'sm_uncertainty':
        sm_cube = iris.load_cube(raw_info['file'],
                                 iris.NameConstraint(var_name='sm'))
        ancillary_var = sm_cube.ancillary_variable(
            'Volumetric Soil Moisture Uncertainty'
        )
        cube = sm_cube.copy(ancillary_var.core_data())
    else:
        cube = iris.load_cube(raw_info['file'], constraint)

    # Continue with processing the cube if it was successfully loaded
    fix_var_metadata(cube, var_info)
    convert_timeunits(cube, year)
    fix_coords_esacci_soilmoisture(cube, overwrite_time_bounds=False)
    set_global_atts(cube, attrs)

    # Remove dysfunctional ancillary data without standard names
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)

    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var_name, vals in cfg['variables'].items():
        all_data_cubes = []
        if not isinstance(vals, dict):  # Ensure vals is a dictionary
            raise ValueError(
                f"Invalid format for variable {var_name}: {type(vals)}"
            )
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var_name)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw']}
        inpfile_pattern = os.path.join(in_dir, vals['filename'])
        logger.info("CMORizing var %s from file type %s",
                    var_name, inpfile_pattern)

        for year in range(vals['start_year'], vals['end_year'] + 1):
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info['file'] = inpfile
                cube = extract_variable(var_info, raw_info, glob_attrs, year)
                all_data_cubes.append(cube)
        final_cube = concatenate(all_data_cubes)
        time = final_cube.coord('time')
        time.bounds = get_time_bounds(time, vals['frequency'])

        # Save daily data with Eday mip
        if var_name == 'sm':
            final_cube.attributes['mip'] = 'Eday'
            final_cube.attributes["history"] = f"Created on {datetime.now()}"
            glob_attrs['mip'] = 'Eday'
            save_variable(final_cube, var_name, out_dir, glob_attrs,
                          unlimited_dimensions=['time'])

            # Calculate and save monthly means with Lmon mip
            monthly_mean_cube = monthly_statistics(final_cube, 'mean')
            monthly_mean_cube.var_name = var_name
            monthly_mean_cube.attributes["history"] = (
                f"Created on {datetime.now()}")
            glob_attrs['mip'] = 'Lmon'
            monthly_mean_cube.attributes.update(glob_attrs)
            save_variable(monthly_mean_cube, var_name, out_dir, glob_attrs,
                          unlimited_dimensions=['time'])
        else:
            # Save the smStderr data
            final_cube.attributes['mip'] = 'Eday'
            final_cube.attributes["history"] = f"Created on {datetime.now()}"
            glob_attrs['mip'] = 'Eday'
            save_variable(final_cube, var_name, out_dir, glob_attrs,
                          unlimited_dimensions=['time'])
