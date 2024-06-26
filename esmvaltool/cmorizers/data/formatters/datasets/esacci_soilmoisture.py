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
import yaml

import iris
from esmvalcore.cmor.fixes import get_time_bounds
from esmvalcore.preprocessor import concatenate

from ...utilities import (
    convert_timeunits,
    fix_coords_esacci_soilmoisture,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)

def extract_variable(var_info, raw_info, attrs, year):
    """Extract variables."""
    rawvar = raw_info['name']
    constraint = iris.Constraint(name=rawvar)

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

    all_data_cubes = {}
    all_sm_uncertainty_cubes = []

    # run the cmorization
    for var_name, vals in cfg['variables'].items():
        if isinstance(vals, dict):  # Ensure vals is a dictionary
            try:
                var_info = cfg['cmor_table'].get_variable(vals['mip'], vals['short_name'])
            except KeyError as e:
                raise ValueError(f"Missing key in variable configuration: {e}")

            glob_attrs['mip'] = vals['mip']
            raw_info = {'name': vals['raw']}
            inpfile_pattern = os.path.join(in_dir, vals['filename'])
            logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)
            
            all_data_cubes[var_name] = []

            for year in range(vals['start_year'], vals['end_year'] + 1):
                year_inpfile_pattern = inpfile_pattern.format(year=year)
                inpfiles = sorted(glob.glob(year_inpfile_pattern))
                for inpfile in inpfiles:
                    raw_info['file'] = inpfile
                    logger.info("CMORizing var %s from file type %s", var_name, raw_info['file'])
                    cube = extract_variable(var_info, raw_info, glob_attrs, year)
                    if cube is not None:
                        all_data_cubes[var_name].append(cube)

                    # if 'sm_uncertainty' in cfg['variables']:
                    #     raw_info['name'] = cfg['variables']['sm_uncertainty']['raw']
                    #     logger.info("CMORizing var sm_uncertainty from file type %s", raw_info['file'])
                    #     sm_uncertainty_cube = extract_variable(var_info, raw_info, glob_attrs, year)
                    #     if sm_uncertainty_cube is not None:
                    #         all_sm_uncertainty_cubes.append(sm_uncertainty_cube)

            # Process the accumulated data after all years have been processed
            if all_data_cubes[var_name]:
                final_cube = concatenate(all_data_cubes[var_name])
                time = final_cube.coord('time')
                time.bounds = get_time_bounds(time, vals['frequency'])
                save_variable(final_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
            else:
                logger.warning(f"No valid data found for {var_name} in the period {vals['start_year']}-{vals['end_year']}. Skipping.")

    # if all_sm_uncertainty_cubes:
    #     final_sm_uncertainty_cube = concatenate(all_sm_uncertainty_cubes)
    #     final_sm_uncertainty_cube.rename('smStderr')
    #     save_variable(final_sm_uncertainty_cube, 'smStderr', out_dir, glob_attrs, unlimited_dimensions=['time'])
    # else:
    #     logger.warning(f"No valid data found for sm_uncertainty in the period {vals['start_year']}-{vals['end_year']}. Skipping smStderr.")


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





# import glob
# import logging
# import os
# import yaml

# import iris
# from esmvalcore.cmor.fixes import get_time_bounds
# from esmvalcore.preprocessor import concatenate

# from ...utilities import (
#     convert_timeunits,
#     fix_coords_esacci_soilmoisture,
#     fix_var_metadata,
#     save_variable,
#     set_global_atts,
# )

# logger = logging.getLogger(__name__)

# def extract_variable(var_info, raw_info, attrs, year):
#     """Extract variables."""
#     rawvar = raw_info['name']
#     constraint = iris.Constraint(name=rawvar)

#     try:
#         cube = iris.load_cube(raw_info['file'], constraint)
#     except iris.exceptions.ConstraintMismatchError:
#         logger.warning(f"No data found for variable {rawvar} in file {raw_info['file']} and year {year}")
#         return None  # Return None or handle appropriately when data is not found
#     except Exception as e:
#         logger.error(f"Error loading data for variable {rawvar} and year {year}: {e}")
#         return None  # Handle other exceptions gracefully

#     if cube is None:
#         logger.warning(f"Loaded cube is None for variable {rawvar} and year {year}. Skipping concatenation.")
#         return None

#     # Continue with processing the cube if it was successfully loaded
#     fix_var_metadata(cube, var_info)
#     convert_timeunits(cube, year)
#     fix_coords_esacci_soilmoisture(cube, overwrite_time_bounds=False)
#     set_global_atts(cube, attrs)

#     # Remove dysfunctional ancillary data without standard names
#     for ancillary_variable in cube.ancillary_variables():
#         cube.remove_ancillary_variable(ancillary_variable)

#     return cube


# def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
#     """Cmorize data."""
#     glob_attrs = cfg['attributes']

#     all_data_cubes = []
#     all_sm_uncertainty_cubes = []

#     # run the cmorization
#     for var_name, vals in cfg['variables'].items():
#         if isinstance(vals, dict):  # Ensure vals is a dictionary
#             try:
#                 var_info = cfg['cmor_table'].get_variable(vals['mip'], vals['short_name'])
#             except KeyError as e:
#                 raise ValueError(f"Missing key in variable configuration: {e}")

#             glob_attrs['mip'] = vals['mip']
#             raw_info = {'name': vals['raw']}
#             inpfile_pattern = os.path.join(in_dir, vals['filename'])
#             logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)
#             for year in range(vals['start_year'], vals['end_year'] + 1):
#                 year_inpfile_pattern = inpfile_pattern.format(year=year)
#                 inpfiles = sorted(glob.glob(year_inpfile_pattern))
#                 for inpfile in inpfiles:
#                     raw_info['file'] = inpfile
#                     logger.info("CMORizing var %s from file type %s", var_name, raw_info['file'])
#                     cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                     if cube is not None:
#                         all_data_cubes.append(cube)

#                     if 'sm_uncertainty' in cfg['variables']:
#                         raw_info['name'] = cfg['variables']['sm_uncertainty']['raw']
#                         logger.info("CMORizing var sm_uncertainty from file type %s", raw_info['file'])
#                         sm_uncertainty_cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                         if sm_uncertainty_cube is not None:
#                             all_sm_uncertainty_cubes.append(sm_uncertainty_cube)

#             # Process the accumulated data after all years have been processed
#             if all_data_cubes:
#                 final_cube = concatenate(all_data_cubes)
#                 time = final_cube.coord('time')
#                 time.bounds = get_time_bounds(time, vals['frequency'])
#                 save_variable(final_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
#             else:
#                 logger.warning(f"No valid data found for {var_name} in the period {vals['start_year']}-{vals['end_year']}. Skipping.")

#             if all_sm_uncertainty_cubes:
#                 final_sm_uncertainty_cube = concatenate(all_sm_uncertainty_cubes)
#                 final_sm_uncertainty_cube.rename('smStderr')
#                 save_variable(final_sm_uncertainty_cube, 'smStderr', out_dir, glob_attrs, unlimited_dimensions=['time'])
#             else:
#                 logger.warning(f"No valid data found for sm_uncertainty in the period {vals['start_year']}-{vals['end_year']}. Skipping smStderr.")
#         else:
#             raise ValueError(f"Invalid format for variable {var_name}: {type(vals)}")


# if __name__ == "__main__":
#     # Example usage
#     in_dir = '/scratch/b/b309265/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'
#     out_dir = '/scratch/b/b309265/esmvaltool_output/'
#     cfg_path = '/home/b/b309265/ESMValTool/esmvaltool/cmorizers/data/cmor_config/ESACCI-SOILMOISTURE.yml'

#     with open(cfg_path, 'r') as f:
#         cfg = yaml.safe_load(f)

#     start_date = cfg.get('start_date')
#     end_date = cfg.get('end_date')

#     if end_date is None:
#         raise ValueError("end_date cannot be None in configuration. Please provide a valid end date.")

#     logger.info("Files in input directory:")
#     logger.info(os.listdir(in_dir))

#     cmorization(in_dir, out_dir, cfg, None, start_date, end_date)


# import glob
# import logging
# import os
# import yaml

# import iris
# from esmvalcore.cmor.fixes import get_time_bounds
# from esmvalcore.preprocessor import concatenate

# from ...utilities import (
#     convert_timeunits,
#     fix_coords_esacci_soilmoisture,
#     fix_var_metadata,
#     save_variable,
#     set_global_atts,
# )

# logger = logging.getLogger(__name__)

# def extract_variable(var_info, raw_info, attrs, year):
#     """Extract variables."""
#     rawvar = raw_info['name']
#     constraint = iris.Constraint(name=rawvar)

#     try:
#         cube = iris.load_cube(raw_info['file'], constraint)
#     except iris.exceptions.ConstraintMismatchError:
#         logger.warning(f"No data found for variable {rawvar} in file {raw_info['file']} and year {year}")
#         return None  # Return None or handle appropriately when data is not found
#     except Exception as e:
#         logger.error(f"Error loading data for variable {rawvar} and year {year}: {e}")
#         return None  # Handle other exceptions gracefully

#     if cube is None:
#         logger.warning(f"Loaded cube is None for variable {rawvar} and year {year}. Skipping concatenation.")
#         return None

#     # Continue with processing the cube if it was successfully loaded
#     fix_var_metadata(cube, var_info)
#     convert_timeunits(cube, year)
#     fix_coords_esacci_soilmoisture(cube, overwrite_time_bounds=False)
#     set_global_atts(cube, attrs)

#     # Remove dysfunctional ancillary data without standard names
#     for ancillary_variable in cube.ancillary_variables():
#         cube.remove_ancillary_variable(ancillary_variable)

#     return cube


# def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
#     """Cmorize data."""
#     glob_attrs = cfg['attributes']

#     all_data_cubes = []
#     all_sm_uncertainty_cubes = []

#     # run the cmorization
#     for var_name, vals in cfg['variables'].items():
#         if isinstance(vals, dict):  # Ensure vals is a dictionary
#             try:
#                 var_info = cfg['cmor_table'].get_variable(vals['mip'], vals['short_name'])
#             except KeyError as e:
#                 raise ValueError(f"Missing key in variable configuration: {e}")

#             glob_attrs['mip'] = vals['mip']
#             raw_info = {'name': vals['raw']}
#             inpfile_pattern = os.path.join(in_dir, vals['filename'])
#             logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)
#             for year in range(vals['start_year'], vals['end_year'] + 1):
#                 year_inpfile_pattern = inpfile_pattern.format(year=year)
#                 inpfiles = sorted(glob.glob(year_inpfile_pattern))
#                 for inpfile in inpfiles:
#                     raw_info['file'] = inpfile
#                     logger.info("CMORizing var %s from file type %s", var_name, raw_info['file'])
#                     cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                     if cube is not None:
#                         all_data_cubes.append(cube)

#                     if 'sm_uncertainty' in cfg['variables']:
#                         raw_info['name'] = cfg['variables']['sm_uncertainty']['raw']
#                         logger.info("CMORizing var sm_uncertainty from file type %s", raw_info['file'])
#                         sm_uncertainty_cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                         if sm_uncertainty_cube is not None:
#                             all_sm_uncertainty_cubes.append(sm_uncertainty_cube)

#             # Process the accumulated data after all years have been processed
#             if all_data_cubes:
#                 final_cube = concatenate(all_data_cubes)
#                 time = final_cube.coord('time')
#                 time.bounds = get_time_bounds(time, vals['frequency'])
#                 save_variable(final_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
#             else:
#                 logger.warning(f"No valid data found for {var_name} in the period {vals['start_year']}-{vals['end_year']}. Skipping.")

#             if all_sm_uncertainty_cubes:
#                 final_sm_uncertainty_cube = concatenate(all_sm_uncertainty_cubes)
#                 final_sm_uncertainty_cube.rename('smStderr')
#                 save_variable(final_sm_uncertainty_cube, 'smStderr', out_dir, glob_attrs, unlimited_dimensions=['time'])
#             else:
#                 logger.warning(f"No valid data found for sm_uncertainty in the period {vals['start_year']}-{vals['end_year']}. Skipping smStderr.")
#         else:
#             raise ValueError(f"Invalid format for variable {var_name}: {type(vals)}")


# if __name__ == "__main__":
#     # Example usage
#     in_dir = '/scratch/b/b309265/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'
#     out_dir = '/scratch/b/b309265/esmvaltool_output/'
#     cfg_path = '/home/b/b309265/ESMValTool/esmvaltool/cmorizers/data/cmor_config/ESACCI-SOILMOISTURE.yml'

#     with open(cfg_path, 'r') as f:
#         cfg = yaml.safe_load(f)

#     start_date = cfg.get('start_date')
#     end_date = cfg.get('end_date')

#     if end_date is None:
#         raise ValueError("end_date cannot be None in configuration. Please provide a valid end date.")

#     logger.info("Files in input directory:")
#     logger.info(os.listdir(in_dir))

#     cmorization(in_dir, out_dir, cfg, None, start_date, end_date)


# import glob
# import logging
# import os
# import yaml

# import iris
# from esmvalcore.cmor.fixes import get_time_bounds
# from esmvalcore.preprocessor import concatenate

# from ...utilities import (
#     convert_timeunits,
#     fix_coords_esacci_soilmoisture,
#     fix_var_metadata,
#     save_variable,
#     set_global_atts,
# )

# logger = logging.getLogger(__name__)

# def extract_variable(var_info, raw_info, attrs, year):
#     """Extract variables."""
#     rawvar = raw_info['name']
#     constraint = iris.Constraint(name=rawvar)

#     try:
#         cube = iris.load_cube(raw_info['file'], constraint)
#     except iris.exceptions.ConstraintMismatchError:
#         logger.warning(f"No data found for variable {rawvar} in file {raw_info['file']} and year {year}")
#         return None  # Return None or handle appropriately when data is not found
#     except Exception as e:
#         logger.error(f"Error loading data for variable {rawvar} and year {year}: {e}")
#         return None  # Handle other exceptions gracefully

#     if cube is None:
#         logger.warning(f"Loaded cube is None for variable {rawvar} and year {year}. Skipping concatenation.")
#         return None

#     # Continue with processing the cube if it was successfully loaded
#     fix_var_metadata(cube, var_info)
#     convert_timeunits(cube, year)
#     fix_coords_esacci_soilmoisture(cube, overwrite_time_bounds=False)
#     set_global_atts(cube, attrs)

#     # Remove dysfunctional ancillary data without standard names
#     for ancillary_variable in cube.ancillary_variables():
#         cube.remove_ancillary_variable(ancillary_variable)

#     return cube


# def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
#     """Cmorize data."""
#     glob_attrs = cfg['attributes']

#     all_data_cubes = []
#     all_sm_uncertainty_cubes = []

#     # run the cmorization
#     for var_name, vals in cfg['variables'].items():
#         if isinstance(vals, dict):  # Ensure vals is a dictionary
#             try:
#                 var_info = cfg['cmor_table'].get_variable(vals['mip'], vals['short_name'])
#             except KeyError as e:
#                 raise ValueError(f"Missing key in variable configuration: {e}")

#             glob_attrs['mip'] = vals['mip']
#             raw_info = {'name': vals['raw']}
#             inpfile_pattern = os.path.join(in_dir, vals['filename'])
#             logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)
#             for year in range(vals['start_year'], vals['end_year'] + 1):
#                 year_inpfile_pattern = inpfile_pattern.format(year=year)
#                 inpfiles = sorted(glob.glob(year_inpfile_pattern))
#                 for inpfile in inpfiles:
#                     raw_info['file'] = inpfile
#                     logger.info("CMORizing var %s from file type %s", var_name, raw_info['file'])
#                     cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                     if cube is not None:
#                         all_data_cubes.append(cube)

#                     # if 'sm_uncertainty' in cfg['variables']:
#                     #     raw_info['name'] = cfg['variables']['sm_uncertainty']['raw']
#                     #     logger.info("CMORizing var sm_uncertainty from file type %s", raw_info['file'])
#                     #     sm_uncertainty_cube = extract_variable(var_info, raw_info, glob_attrs, year)
#                     #     if sm_uncertainty_cube is not None:
#                     #         all_sm_uncertainty_cubes.append(sm_uncertainty_cube)

#             # Process the accumulated data after all years have been processed
#             if all_data_cubes:
#                 final_cube = concatenate(all_data_cubes)
#                 time = final_cube.coord('time')
#                 time.bounds = get_time_bounds(time, vals['frequency'])
#                 save_variable(final_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
#             else:
#                 logger.warning(f"No valid data found for {var_name} in the period {vals['start_year']}-{vals['end_year']}. Skipping.")

#             # if all_sm_uncertainty_cubes:
#             #     final_sm_uncertainty_cube = concatenate(all_sm_uncertainty_cubes)
#             #     final_sm_uncertainty_cube.rename('smStderr')
#             #     save_variable(final_sm_uncertainty_cube, 'smStderr', out_dir, glob_attrs, unlimited_dimensions=['time'])
#             # else:
#             #     logger.warning(f"No valid data found for sm_uncertainty in the period {vals['start_year']}-{vals['end_year']}. Skipping smStderr.")
#         else:
#             raise ValueError(f"Invalid format for variable {var_name}: {type(vals)}")


# if __name__ == "__main__":
#     # Example usage
#     in_dir = '/scratch/b/b309265/RAWOBS/Tier2/ESACCI-SOILMOISTURE/'
#     out_dir = '/scratch/b/b309265/esmvaltool_output/'
#     cfg_path = '/home/b/b309265/ESMValTool/esmvaltool/cmorizers/data/cmor_config/ESACCI-SOILMOISTURE.yml'

#     with open(cfg_path, 'r') as f:
#         cfg = yaml.safe_load(f)

#     start_date = cfg.get('start_date')
#     end_date = cfg.get('end_date')

#     if end_date is None:
#         raise ValueError("end_date cannot be None in configuration. Please provide a valid end date.")

#     logger.info("Files in input directory:")
#     logger.info(os.listdir(in_dir))

#     cmorization(in_dir, out_dir, cfg, None, start_date, end_date)


