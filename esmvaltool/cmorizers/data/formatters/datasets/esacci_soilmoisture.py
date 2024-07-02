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
Modification history
    20240626-cammarano_diego: written.
"""


import glob
import logging
import os
import iris
import yaml
from esmvalcore.cmor.fixes import get_time_bounds
from esmvalcore.preprocessor import concatenate
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
            if overwrite_time_bounds or not cube.coord('time').has_bounds():
                fix_bounds(cube, cube.coord('time'))

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
            if overwrite_lon_bounds or not cube_coord.has_bounds():
                fix_bounds(cube, cube_coord)

        # fix latitude
        if cube_coord.var_name == 'lat':
            logger.info("Fixing latitude...")
            if overwrite_lat_bounds or not cube.coord('latitude').has_bounds():
                fix_bounds(cube, cube.coord('latitude'))

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
        if isinstance(vals, dict):  # Ensure vals is a dictionary
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
                    logger.info("CMORizing var %s from file type %s",
                                var_name, raw_info['file'])
                    cube = extract_variable(var_info, raw_info, glob_attrs,
                                            year)
                    all_data_cubes.append(cube)
            final_cube = concatenate(all_data_cubes)
            time = final_cube.coord('time')
            time.bounds = get_time_bounds(time, vals['frequency'])
            save_variable(final_cube, var_name, out_dir, glob_attrs,
                          unlimited_dimensions=['time'])

        else:
            raise ValueError(
                f"Invalid format for variable {var_name}: {type(vals)}"
            )


# if __name__ == "__main__":
#     # RAWOBS dir
#     in_dir = '/scratch/b/b309265/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'
#     # in_dir= '/work/bd0854/DATA/ESMValTool2/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'

#     # OBS dir CMOR-compliant
#     out_dir = '/scratch/b/b309265/esmvaltool_output/'
#     # out_dir = '/work/bd0854/DATA/ESMValTool2/OBS/'

#     # Configuration file ESACCI-SOILMOISTURE
#     cfg_path = '../cmor_config/ESACCI-SOILMOISTURE.yml'

#     with open(cfg_path, 'r') as f:
#         cfg = yaml.safe_load(f)

#     start_date = cfg.get('start_date')
#     end_date = cfg.get('end_date')

#     if end_date is None:
#         raise ValueError("end_date cannot be None in configuration."
#                          "Please provide a valid end date.")

#     logger.info("Files in input directory:")
#     logger.info(os.listdir(in_dir))

#     cmorization(in_dir, out_dir, cfg, None, start_date, end_date)
