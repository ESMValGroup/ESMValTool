"""ESMValTool CMORizer for LAI3g data.

Tier
    Tier 3: restricted dataset.

Source
    http://cliveg.bu.edu/modismisr/lai3g-fpar3g.html

Last access
    20190503

Download and processing instructions
    To obtain the data sets it is necessary to contact Ranga B. Myneni
    (Department of Earth and Environment, Boston University). See link above
    for more information.

    By default, this dataset is regridded to a 1°x1° grid (original resoultion
    is 1/12°). If you want to use the original resolution, remove the `regrid`
    section in the configuration file (`LAI3g.yml`). Note that in this case,
    preprocessing the dataset with ESMValTool (i.e. every time you run the
    tool) can take a very long time (> 30 min).

"""

import glob
import logging
import os
import shutil
import zipfile
from datetime import datetime

import iris
import iris.coord_categorisation
import numpy as np
from cf_units import Unit

from esmvalcore.preprocessor import regrid

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)

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


def _extract_variable(cmor_info, attrs, in_dir, out_dir, cfg):
    """Extract variable."""
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


def _unzip(filepath, out_dir):
    """Unzip `*.zip` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    with zipfile.ZipFile(filepath, 'r') as zip_ref:
        zip_ref.extractall(out_dir)
    new_path = os.path.join(out_dir, 'LAI')
    logger.info("Succefully extracted files to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        if cfg.get('regrid'):
            cfg['regrid'].setdefault('target_grid', '1x1')
            cfg['regrid'].setdefault('scheme', 'nearest')
            logger.info(
                "Final dataset will be regridded to %s grid using scheme '%s'",
                cfg['regrid']['target_grid'], cfg['regrid']['scheme'])
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        zip_file = os.path.join(in_dir, filepath)
        if not os.path.isfile(zip_file):
            logger.debug("Skipping '%s', file '%s' not found", var, zip_file)
            continue
        logger.info("Found input file '%s'", zip_file)
        file_dir = _unzip(zip_file, out_dir)
        _extract_variable(cmor_info, glob_attrs, file_dir, out_dir, cfg)
        _clean(file_dir)
