import iris
import os
import glob
import numpy as np
import logging
import gc
from datetime import datetime
from esmvalcore.preprocessor import regrid

from ...utilities import (
    fix_dim_coordnames,
    fix_bounds,
    #save_variable,
    fix_dtype,
    set_global_atts,
    fix_coords
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set the future configuration option to enable split attributes handling mode
#iris.FUTURE.save_split_attrs = True

def extract_variable(raw_info, year):
    """Extract the variable from the raw data file."""
    logger.info(f"Extracting variable for year {year} from {raw_info['file']}")
    cube_list = iris.load(raw_info['file'], raw_info['name'])
    if not cube_list:
        logger.warning(f"No cubes found for {raw_info['name']} in file {raw_info['file']}")
    return cube_list

def regrid_tile(tile, target_cube):
    """Regrid a single tile."""
    return tile.regrid(target_cube, iris.analysis.Linear())
    #return tile.regrid(target_cube, iris.analysis.AreaWeighted())


def calculate_bounds(points):
    """Calculate bounds for a set of points."""
    bounds = np.zeros((len(points), 2))
    bounds[1:, 0] = (points[:-1] + points[1:]) / 2
    bounds[:-1, 1] = bounds[1:, 0]
    bounds[0, 0] = points[0] - (bounds[1, 0] - points[0])
    bounds[-1, 1] = points[-1] + (points[-1] - bounds[-2, 1])
    return bounds

def regrid_iris(cube):
    """Regrid the cubes using Iris."""
    logger.info("Regridding using Iris")
    
    target_lats = np.linspace(-90, 90, int(180 / 0.25) + 1)
    target_lons = np.linspace(-180, 180, int(360 / 0.25) + 1)
    
    lat_bounds = calculate_bounds(target_lats)
    lon_bounds = calculate_bounds(target_lons)
    
    target_cube = iris.cube.Cube(np.zeros((len(target_lats), len(target_lons))),
                                 dim_coords_and_dims=[(iris.coords.DimCoord(target_lats, standard_name='latitude', units='degrees', bounds=lat_bounds), 0),
                                                      (iris.coords.DimCoord(target_lons, standard_name='longitude', units='degrees', bounds=lon_bounds), 1)])

    # Fix bounds for the source cube coordinates
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
    
    # Ensure bounds are valid and log for debugging
    logger.debug(f"Latitude bounds: {cube.coord('latitude').bounds}")
    logger.debug(f"Longitude bounds: {cube.coord('longitude').bounds}")

    combined_data = np.zeros((cube.shape[0], len(target_lats), len(target_lons)))

    tile_size = 1012  # Smaller tile size to reduce memory usage
    for time_idx in range(cube.shape[0]):
        for lat_start in range(0, cube.shape[1], tile_size):
            for lon_start in range(0, cube.shape[2], tile_size):
                lat_end = min(lat_start + tile_size, cube.shape[1])
                lon_end = min(lon_start + tile_size, cube.shape[2])
                
                # Extract sub-cube data
                sub_cube_data = cube.data[time_idx, lat_start:lat_end, lon_start:lon_end].astype(np.float32)
                
                # Log sub-cube data to debug
                logger.debug(f"Sub-cube data (time_idx={time_idx}, lat_start={lat_start}, lon_start={lon_start}): "
                             f"{sub_cube_data}")
                
                sub_cube = iris.cube.Cube(sub_cube_data,
                                          dim_coords_and_dims=[(cube.coord('latitude')[lat_start:lat_end], 0),
                                                               (cube.coord('longitude')[lon_start:lon_end], 1)])
                
                # Log sub-cube before regridding
                logger.debug(f"Sub-cube before regridding: {sub_cube.data}")

                try:
                    regridded_sub_cube = regrid_tile(sub_cube, target_cube)
                    
                    # Log regridded sub-cube data to debug
                    logger.debug(f"Regridded sub-cube data (time_idx={time_idx}, lat_start={lat_start}, lon_start={lon_start}): "
                                 f"{regridded_sub_cube.data}")
                    
                except Exception as e:
                    logger.error(f"Regridding error at time_idx={time_idx}, lat_start={lat_start}, lon_start={lon_start}: {e}")
                    continue
                
                lat_indices = [np.abs(target_lats - lat).argmin() for lat in regridded_sub_cube.coord('latitude').points]
                lon_indices = [np.abs(target_lons - lon).argmin() for lon in regridded_sub_cube.coord('longitude').points]
                
                # Log indices to verify correct assignment
                logger.debug(f"Lat indices: {lat_indices}, Lon indices: {lon_indices}")
                
                # Ensure indices are within bounds
                if lat_indices and lon_indices:
                    combined_data[time_idx, lat_indices[0]:lat_indices[-1] + 1, lon_indices[0]:lon_indices[-1] + 1] = regridded_sub_cube.data

    combined_cube = iris.cube.Cube(combined_data,
                                   dim_coords_and_dims=[(cube.coord('time'), 0),
                                                        (iris.coords.DimCoord(target_lats, standard_name='latitude', units='degrees', bounds=lat_bounds), 1),
                                                        (iris.coords.DimCoord(target_lons, standard_name='longitude', units='degrees', bounds=lon_bounds), 2)])
    
    # Log combined data to debug
    logger.debug(f"Combined cube data: {combined_data}")
    
    return combined_cube




def save_variable(cube, var, outdir, attrs, **kwargs):
    """Saver function.

    Saves iris cubes (data variables) in CMOR-standard named files.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to be saved.

    var: str
        Variable short_name e.g. ts or tas.

    outdir: str
        root directory where the file will be saved.

    attrs: dict
        dictionary holding cube metadata attributes like
        project_id, version etc.

    **kwargs: kwargs
        Keyword arguments to be passed to iris.save
    """
    fix_dtype(cube)

    # Ensure the variable metadata is correctly set
    cube.var_name = var

    # Only set standard_name if it is a valid CF standard name
    if hasattr(cube, 'standard_name') and cube.standard_name in iris.std_names.STD_NAMES:
        cube.standard_name = cube.standard_name
    else:
        cube.standard_name = None  # or keep it unset

    # Ensure long_name is set
    cube.long_name = cube.long_name or var

    # Debugging metadata before saving
    logger.debug(f"Saving cube with var_name={cube.var_name}, "
                 f"standard_name={cube.standard_name}, long_name={cube.long_name}")

    # CMOR standard
    try:
        time = cube.coord('time')
    except iris.exceptions.CoordinateNotFoundError:
        time_suffix = None
    else:
        if len(time.points) == 1 and "mon" not in cube.attributes.get('mip'):
            year = str(time.cell(0).point.year)
            time_suffix = '-'.join([year + '01', year + '12'])
        else:
            date1 = (
                f"{time.cell(0).point.year:d}{time.cell(0).point.month:02d}"
            )
            date2 = (
                f"{time.cell(-1).point.year:d}{time.cell(-1).point.month:02d}"
            )
            time_suffix = '-'.join([date1, date2])

    name_elements = [
        attrs['project_id'],
        attrs['dataset_id'],
        attrs['modeling_realm'],
        attrs['version'],
        attrs['mip'],
        var,
    ]
    if time_suffix:
        name_elements.append(time_suffix)
    file_name = '_'.join(name_elements) + '.nc'
    file_path = os.path.join(outdir, file_name)
    logger.info('Saving: %s', file_path)
    status = 'lazy' if cube.has_lazy_data() else 'realized'
    logger.info('Cube has %s data [lazy is preferred]', status)
    iris.save(cube, file_path, fill_value=1e20, **kwargs)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date=None, end_date=None):
    """Cmorize data."""
    glob_attrs = cfg['attributes']
    if not start_date:
        start_date = datetime(2000, 1, 1)
    if not end_date:
        end_date = datetime(2000, 12, 31)

    for var_name, vals in cfg['variables'].items():
        if not isinstance(vals, dict):
            raise ValueError(f"Invalid format for variable {var_name}: {type(vals)}")
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var_name)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['long_name']}
        inpfile_pattern = os.path.join(in_dir, vals['filename'])
        logger.info("CMORizing var %s from file type %s", var_name, inpfile_pattern)

        for year in range(start_date.year, end_date.year + 1):
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info['file'] = inpfile
                try:
                    logger.info("Starting variable extraction")
                    cube_list = extract_variable(raw_info, year=year)
                    logger.info("End variable extraction")
                    for cube in cube_list:
                        logger.info(f"Regridding cube for {var_name}")
                        regridded_cube = regrid_iris(cube)
                        logger.info("Regridding done")
                        regridded_cube = fix_coords(regridded_cube)
                        set_global_atts(regridded_cube, glob_attrs)
                        output_filename = f"{var_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.nc"
                        output_filepath = os.path.join(out_dir, output_filename)
                        logger.info(f"Saving: {output_filepath}")
                        save_variable(regridded_cube, var_name, out_dir, glob_attrs, unlimited_dimensions=['time'])
                        logger.info(f"Saved {var_name} to {output_filepath}")
                    del cube_list  # Free memory
                    gc.collect()  # Explicitly call garbage collection
                except Exception as e:
                    logger.error(f"Failed to process file {inpfile}: {e}")
                    continue