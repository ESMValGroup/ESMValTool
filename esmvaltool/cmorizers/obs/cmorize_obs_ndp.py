"""ESMValTool CMORizer for NDP data.

Tier
    Tier 3: restricted dataset.

Source
    https://data.ess-dive.lbl.gov/view/doi:10.3334/CDIAC/LUE.NDP017.2006

Last access
    20191014

Download and processing instructions
    Download the following files:
        ndp017b.tar.gz
    A registration is required for downloading the data.

"""

import gdal
import glob
import gzip
import logging
import os
import shutil
import tarfile
import zipfile
from datetime import datetime
from PIL import Image

import iris
import iris.coord_categorisation
import numpy as np
from cf_units import Unit

from esmvalcore.preprocessor import regrid

from . import utilities as utils

logger = logging.getLogger(__name__)



def _clean(file_dir):
    """Remove unzipped input files."""
    if os.path.isdir(file_dir):
        shutil.rmtree(file_dir)
        logger.info("Removed cached directory %s", file_dir)


def _extract_variable(cmor_info, attrs, var_file, out_dir, cfg):
    """Extract variable."""
    grid_file = gdal.Open(var_file)
    array = grid_file.ReadAsArray()
    print(np.min(array))
    print(np.max(array))
    print(np.mean(array))
    print(array)
    assert False

    driver = ogr.GetDriverByName('E00GRID')
    grid_file = driver.Open(var_file)
    print(grid_file)
    assert False
    nc_files = []
    for year in _get_years(in_dir, cfg):
        cube_path = _get_cube_for_year(year, in_dir, cfg)
        nc_files.append(cube_path)

    # Build final cube
    logger.info("Building final cube")
    cubes = iris.cube.CubeList()
    for nc_file in nc_files:
        cube = iris.load_cube(nc_file)
        cubes.append(cube)
    final_cube = cubes.concatenate_cube()
    utils.fix_var_metadata(final_cube, cmor_info)
    utils.convert_timeunits(final_cube, 1950)
    utils.fix_coords(final_cube)
    if not cfg.get('regrid'):
        utils.flip_dim_coord(final_cube, 'latitude')
    utils.set_global_atts(final_cube, attrs)
    utils.save_variable(final_cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _get_coords(year, filename, cfg):
    """Get correct coordinates for cube."""
    filename = os.path.basename(filename)
    time_units = Unit('days since 1950-1-1 00:00:00', calendar='standard')

    # Extract date from filename
    time_str = filename.replace(cfg['binary_prefix'], '')
    month = MONTHS[time_str[4:7]]
    day = DAYS[time_str[7:8]]
    date = datetime(year, month, day)

    # Build time coordinate
    time_data = [time_units.date2num(date)]
    time_coord = iris.coords.DimCoord(time_data,
                                      standard_name='time',
                                      long_name='time',
                                      var_name='time',
                                      units=time_units)

    # Build latitude/Longitude coordinates
    latitude_data = np.linspace(UPPER_LEFT_LAT, LOWER_RIGHT_LAT, N_LAT)
    longitude_data = np.linspace(UPPER_LEFT_LON, LOWER_RIGHT_LON, N_LON)
    lat_coord = iris.coords.DimCoord(latitude_data,
                                     standard_name='latitude',
                                     long_name='latitude',
                                     var_name='lat',
                                     units='degrees')
    lon_coord = iris.coords.DimCoord(longitude_data,
                                     standard_name='longitude',
                                     long_name='longitude',
                                     var_name='lon',
                                     units='degrees')

    return [(time_coord, 0), (lat_coord, 1), (lon_coord, 2)]


def _get_cube_for_year(year, in_dir, cfg):
    """Exract cube containing one year from raw file."""
    logger.info("Processing year %i", year)
    bin_files = glob.glob(
        os.path.join(in_dir, f"{cfg['binary_prefix']}{year}*.bin"))

    # Read files of one year
    cubes = iris.cube.CubeList()
    for bin_file in bin_files:
        raw_data = np.fromfile(bin_file, DTYPE,
                               N_LAT * N_LON).reshape(1, N_LAT, N_LON)
        raw_data = np.ma.masked_equal(raw_data, MISSING_VALUE)
        raw_data = raw_data.astype(np.float32)
        raw_data /= SCALE_FACTOR

        # Build coordinates and cube, regrid, and append it
        coords = _get_coords(year, bin_file, cfg)
        cube = iris.cube.Cube(raw_data, dim_coords_and_dims=coords)
        if cfg.get('regrid'):
            cube = regrid(cube, cfg['regrid']['target_grid'],
                          cfg['regrid']['scheme'])
        cubes.append(cube)

    # Build cube for single year with monthly data
    # (Raw data has two values per month)
    cube = cubes.concatenate_cube()
    iris.coord_categorisation.add_month_number(cube, 'time')
    cube = cube.aggregated_by('month_number', iris.analysis.MEAN)

    # Cache cube on disk to save memory
    cached_path = os.path.join(in_dir, f'{year}.nc')
    iris.save(cube, cached_path)
    logger.info("Cached %s", cached_path)
    return cached_path


def _get_years(in_dir, cfg):
    """Get all available years from input directory."""
    bin_files = os.listdir(in_dir)
    bin_files = [f.replace(cfg['binary_prefix'], '') for f in bin_files]
    years = {int(f[:4]) for f in bin_files}
    return years


def _extract(filepath, out_dir):
    """Extract `*.tar.gz` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    with tarfile.open(filepath) as tar:
        tar.extractall()
    new_path = os.path.join(out_dir, 'ndp017b', 'revised')
    logger.info("Succesfully extracted files to %s", new_path)
    return new_path


def _unzip(filepath, out_dir):
    """Extract `*.gz` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    new_path = filepath.replace('.gz', '')
    with gzip.open(filepath, 'rb') as zip_file:
        with open(new_path, 'wb') as new_file:
            shutil.copyfileobj(zip_file, new_file)
    logger.info("Succesfully extracted '%s' to '%s'", filepath, new_path)
    return new_path


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = os.path.join(in_dir, cfg['filename'])
    tar_file = os.path.join(in_dir, filepath)
    logger.info("Found input file '%s'", tar_file)
    file_dir = _extract(tar_file, out_dir)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        zip_file = os.path.join(file_dir, var_info['filename'])
        var_file = _unzip(zip_file, out_dir)
        logger.info("Found input file '%s' for variable '%s'", var_file, var)
        _extract_variable(cmor_info, glob_attrs, var_file, out_dir, cfg)
    _clean(file_dir)
