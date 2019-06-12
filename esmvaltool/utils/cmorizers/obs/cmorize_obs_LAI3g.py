"""ESMValTool CMORizer for LAI3g data.

Tier
    Tier 3: restricted dataset.

Source
    http://cliveg.bu.edu/modismisr/lai3g-fpar3g.html

Last access
    20190503

Download and processing instructions
    To obtain the data sets it is necessary to contance Ranga B. Myneni
    (Department of Earth and Environment, Boston University). See link above
    for more information.

"""

import logging
import os
import shutil
import zipfile
from datetime import datetime

import iris
import numpy as np
from cf_units import Unit

import esmvaltool.utils.cmorizers.obs.utilities as utils

logger = logging.getLogger(__name__)

CFG = utils.read_cmor_config('LAI3g.yml')

# Properties of the binary file (cannot be stored in .yml since computations
# are necessary)
DTYPE = '>i2'
N_LAT = 2160
N_LON = 4320
MISSING_VALUE = -32768
SCALE_FACTOR = 1000.0
UPPER_LEFT_LAT = 90.0 - 1.0 / 24.0
UPPER_LEFT_LON = -180.0 + 1.0 / 24.0
LOWER_RIGHT_LAT = -90.0 + 1.0 / 24.0
LOWER_RIGHT_LON = 180.0 - 1.0 / 24.0
MONTHS = {
    'jan': 1,
    'feb': 2,
    'mar': 3,
    'apr': 4,
    'may': 5,
    'jun': 6,
    'jul': 7,
    'aug': 8,
    'sep': 9,
    'oct': 10,
    'nov': 11,
    'dec': 12,
}
DAYS = {
    'a': 8,
    'b': 23,
}


def _clean(file_dir):
    """Remove unzipped input files."""
    if os.path.isdir(file_dir):
        shutil.rmtree(file_dir)
        logger.info("Removed cached directory %s", file_dir)


def _extract_variable(cmor_info, attrs, in_dir, out_dir):
    """Extract variable."""
    nc_files = []
    for bin_file in os.listdir(in_dir):
        filepath = os.path.join(in_dir, bin_file)
        logger.info("Reading %s", filepath)

        # Data
        raw_data = np.fromfile(filepath, DTYPE,
                               N_LAT * N_LON).reshape(1, N_LAT, N_LON)
        raw_data = np.ma.masked_equal(raw_data, MISSING_VALUE)
        raw_data = raw_data.astype(np.float32)
        raw_data /= SCALE_FACTOR

        # Coordinates
        coords = _get_coords(bin_file)

        # Build cube for single time step and cache it on disk to save memory
        cube = iris.cube.Cube(raw_data, dim_coords_and_dims=coords)
        cached_path = filepath.replace('.bin', '.nc')
        iris.save(cube, cached_path)
        logger.info("Cached %s", cached_path)
        nc_files.append(cached_path)

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
    utils.flip_dim_coord(final_cube, 'latitude')
    utils.set_global_atts(final_cube, attrs)
    utils.save_variable(final_cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _get_coords(filename):
    """Get correct coordinates for cube."""
    time_units = Unit('days since 1950-1-1 00:00:00', calendar='standard')

    # Extract date from filename
    time_str = filename.replace(CFG['binary_prefix'], '')
    year = int(time_str[:4])
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


def _unzip(filepath, out_dir):
    """Unzip `*.zip` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    with zipfile.ZipFile(filepath, 'r') as zip_ref:
        zip_ref.extractall(out_dir)
    new_path = os.path.join(out_dir, 'LAI')
    logger.info("Succefully extracted files to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    glob_attrs = CFG['attributes']
    cmor_table = CFG['cmor_table']
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    filepath = os.path.join(in_dir, CFG['filename'])

    # Run the cmorization
    for (var, var_info) in CFG['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        zip_file = os.path.join(in_dir, filepath)
        if not os.path.isfile(zip_file):
            logger.debug("Skipping '%s', file '%s' not found", var, zip_file)
            continue
        logger.info("Found input file '%s'", zip_file)
        file_dir = _unzip(zip_file, out_dir)
        _extract_variable(cmor_info, glob_attrs, file_dir, out_dir)
        _clean(file_dir)
