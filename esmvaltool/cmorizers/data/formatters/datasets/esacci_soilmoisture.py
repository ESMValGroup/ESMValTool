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
)

logger = logging.getLogger(__name__)

def fix_coords_esacci_soilmoisture(cube,
               overwrite_time_bounds=True,
               overwrite_lon_bounds=True,
               overwrite_lat_bounds=True,
               overwrite_lev_bounds=True,
               overwrite_airpres_bounds=True):
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

    overwrite_lev_bounds: bool (optional)
        set to False not to overwrite depth bounds.

    overwrite_airpres_bounds: bool (optional)
        set to False not to overwrite air pressure bounds.

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
                Unit('days since 1970-01-01T00:00:00+00:00', calendar='proleptic_gregorian'))
            if overwrite_time_bounds or not cube.coord('time').has_bounds():
                fix_bounds(cube, cube.coord('time'))
        
    try:
        time_coord = cube.coord('time')
        
        if time_coord is None:
            raise ValueError("Time coordinate 'time' not found in cube.")
        
        # Ensure time bounds are strictly monotonic
        if time_coord.bounds is None:
            logger.warning("Time bounds are not available for 'time' coordinate.")
        else:
            bounds = time_coord.bounds
            if not (np.all(bounds[:, 0] <= bounds[:, 1]) or np.all(bounds[:, 0] >= bounds[:, 1])):
                logger.warning("Time bounds are not monotonic. Adjusting...")
                
                # Sort the bounds array to ensure monotonicity
                time_coord.bounds = bounds[np.argsort(bounds[:, 0])]
                
                logger.info("Adjusted time bounds to ensure monotonicity.")
        
    except Exception as e:
        logger.error(f"Error fixing coordinates: {e}")
        raise  # Propagate the error up if necessary
    
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

        # fix depth
        if cube_coord.var_name == 'lev':
            logger.info("Fixing depth...")
            if overwrite_lev_bounds or not cube.coord('depth').has_bounds():
                fix_bounds(cube, cube.coord('depth'))

        # fix air_pressure
        if cube_coord.var_name == 'air_pressure':
            logger.info("Fixing air pressure...")
            if overwrite_airpres_bounds \
                    or not cube.coord('air_pressure').has_bounds():
                fix_bounds(cube, cube.coord('air_pressure'))

    # remove CS
    cube.coord('latitude').coord_system = None
    cube.coord('longitude').coord_system = None

    return cube

def extract_variable(var_info, raw_info, attrs, year):
    """Extract variables."""
    rawvar = raw_info['name']
    constraint = iris.Constraint(name=rawvar)

    if rawvar == 'sm_uncertainty':
        sm_cube = iris.load_cube(raw_info['file'], iris.NameConstraint(var_name='sm'))
        ancillary_var = sm_cube.ancillary_variable('Volumetric Soil Moisture Uncertainty')
        cube = sm_cube.copy(ancillary_var.core_data())

    try:
        cube = iris.load_cube(raw_info['file'], constraint)
    except iris.exceptions.ConstraintMismatchError:
        logger.warning(f"No data found for variable {rawvar} in file {raw_info['file']} and year {year}")
        return None  # Return None or handle appropriately when data is not found
    except Exception as e:
        logger.error(f"Error loading data for variable {rawvar} and year {year}: {e}")
        return None  # Handle other exceptions gracefully

    if cube is None:
        logger.warning(f"Loaded cube is None for variable {rawvar} and year {year}. Skipping concatenation.")
        return None

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

    all_data_cubes = []
    all_sm_uncertainty_cubes = []

    # run the cmorization
    for var_name, vals in cfg['variables'].items():
        if isinstance(vals, dict):  # Ensure vals is a dictionary
            try:
                var_info = cfg['cmor_table'].get_variable(vals['mip'], vals['raw'])
            except KeyError as e:
                raise ValueError(f"Missing key in variable configuration: {e}")

            glob_attrs['mip'] = vals['mip']
            raw_info = {'name': vals['raw']}
            inpfile_pattern = os.path.join(in_dir, vals['filename'])
            logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)

            for year in range(vals['start_year'], vals['end_year'] + 1):
                year_inpfile_pattern = inpfile_pattern.format(year=year)
                inpfiles = sorted(glob.glob(year_inpfile_pattern))
                for inpfile in inpfiles:
                    raw_info['file'] = inpfile
                    logger.info("CMORizing var %s from file type %s", var_name, raw_info['file'])
                    cube = extract_variable(var_info, raw_info, glob_attrs, year)
                    if cube is not None:
                        all_data_cubes.append(cube)

                    if 'sm_uncertainty' in cfg['variables']:
                        raw_info['name'] = cfg['variables']['sm_uncertainty']['raw']
                        logger.info("CMORizing var sm_uncertainty from file type %s", raw_info['file'])
                        sm_uncertainty_cube = extract_variable(var_info, raw_info, glob_attrs, year)
                        if sm_uncertainty_cube is not None:
                            all_sm_uncertainty_cubes.append(sm_uncertainty_cube)
                        else:
                            logger.info(f"No sm_uncertainty data found for year {year} and file {raw_info['file']}")

            # Process the accumulated data after all years have been processed
            if all_data_cubes:
                final_cube = concatenate(all_data_cubes)
                time = final_cube.coord('time')
                time.bounds = get_time_bounds(time, vals['frequency'])
                save_variable(final_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
            else:
                logger.warning(f"No valid data found for {var_name} in the period {vals['start_year']}-{vals['end_year']}. Skipping.")

            if all_sm_uncertainty_cubes:
                final_sm_uncertainty_cube = concatenate(all_sm_uncertainty_cubes)
                final_sm_uncertainty_cube.rename('smStderr')
                save_variable(final_sm_uncertainty_cube, 'smStderr', out_dir, glob_attrs, unlimited_dimensions=['time'])
            else:
                logger.warning(f"No valid data found for sm_uncertainty in the period {vals['start_year']}-{vals['end_year']}. Skipping smStderr.")
        else:
            raise ValueError(f"Invalid format for variable {var_name}: {type(vals)}")


if __name__ == "__main__":
    # Example usage
    in_dir = '/scratch/b/b309265/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'
    out_dir = '/scratch/b/b309265/esmvaltool_output/'
    cfg_path = '/home/b/b309265/ESMValTool/esmvaltool/cmorizers/data/cmor_config/ESACCI-SOILMOISTURE.yml'

    with open(cfg_path, 'r') as f:
        cfg = yaml.safe_load(f)

    start_date = cfg.get('start_date')
    end_date = cfg.get('end_date')

    if end_date is None:
        raise ValueError("end_date cannot be None in configuration. Please provide a valid end date.")

    logger.info("Files in input directory:")
    logger.info(os.listdir(in_dir))

    cmorization(in_dir, out_dir, cfg, None, start_date, end_date)
