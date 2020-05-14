"""ESMValTool CMORizer for JMA-TRANSCOM data.

Tier
    Tier 3: restricted dataset.

Source
    http://www.globalcarbonatlas.org/en/content/atmospheric-inversions

Last access
    20190702

Download and processing instructions
    To obtain the data sets it is necessary to contact Takashi Maki
    (Department of Atmosphere, Ocean and Earth System Modeling Research,
    Meteorological Research Institute, Tsukuba City, Japan). See link above
    for more information.

"""

import logging
import os
import shutil
import tarfile
from datetime import datetime, timedelta

import iris
import iris.coord_categorisation
import numpy as np
from cf_units import Unit

from esmvalcore.preprocessor import mask_landsea

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _clean(file_dir):
    """Remove unzipped input files."""
    if os.path.isdir(file_dir):
        shutil.rmtree(file_dir)
        logger.info("Removed cached directory %s", file_dir)


def _extract_variable(cmor_info, attrs, in_dir, out_dir, ctl):
    """Extract variable."""
    filepath = os.path.join(in_dir, ctl['binary_prefix'] + '.dat')
    raw_data = np.fromfile(filepath, ctl['dtype'],
                           ctl['t_size'] * ctl['y_size'] *
                           ctl['x_size']).reshape(ctl['t_size'], ctl['y_size'],
                                                  ctl['x_size'])

    # Get coordinates
    coords = _get_coords(ctl)

    # Build cube
    cube = iris.cube.Cube(raw_data, dim_coords_and_dims=coords)

    # Mask appropriate parts
    if cmor_info.short_name == 'nbp':
        cube = mask_landsea(cube, {}, 'sea')
    elif cmor_info.short_name == 'fgco2':
        cube = mask_landsea(cube, {}, 'land')
    else:
        raise NotImplementedError(
            f"CMORizer for '{cmor_info.short_name}' not implemented yet")

    # Fix metadata
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _add_months(date, delta):
    """Add months to a date."""
    add_years = delta // 12
    add_months = delta % 12
    return date.replace(year=date.year + add_years,
                        month=date.month + add_months)


def _get_coords(ctl):
    """Get correct coordinates for cube."""
    # Time
    time_units = Unit('days since 1950-1-1 00:00:00', calendar='standard')
    time_start = datetime.strptime(ctl['t_start'], '%d%b%Y')
    times = [
        _add_months(time_start, d) + timedelta(days=14)
        for d in int(ctl['t_delta'][0]) * np.arange(ctl['t_size'])
    ]
    times = [time_units.date2num(time) for time in times]
    time_coord = iris.coords.DimCoord(times,
                                      standard_name='time',
                                      long_name='time',
                                      var_name='time',
                                      units=time_units)

    # Latitude
    lats = float(
        ctl['y_start']) + (float(ctl['y_delta']) * np.arange(ctl['y_size']))
    lat_coord = iris.coords.DimCoord(lats,
                                     standard_name='latitude',
                                     long_name='latitude',
                                     var_name='lat',
                                     units='degrees_north')

    # Longitude
    lons = float(
        ctl['x_start']) + (float(ctl['x_delta']) * np.arange(ctl['x_size']))
    lon_coord = iris.coords.DimCoord(lons,
                                     standard_name='longitude',
                                     long_name='longitude',
                                     var_name='lon',
                                     units='degrees_east')

    return [(time_coord, 0), (lat_coord, 1), (lon_coord, 2)]


def _extract_tar(filepath, out_dir):
    """Extract `*.tar.gz` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    with tarfile.open(filepath) as tar:
        tar.extractall()
    new_path = os.path.join(out_dir, 'JMA_2018')
    logger.info("Succesfully extracted files to %s", new_path)
    return new_path


def _read_control_file(file_dir, cfg):
    """Read '*.ctl' file."""
    ctl_path = os.path.join(file_dir, cfg['binary_prefix'] + '.ctl')
    with open(ctl_path, mode='r') as ctl_file:
        contents = ctl_file.read()
    contents = contents.split()
    ctl = {}
    ctl['binary_prefix'] = cfg['binary_prefix']
    endian = contents[contents.index('OPTIONS') + 1].lower()
    if endian == 'big_endian':
        ctl['dtype'] = '>f4'
    elif endian == 'little_endian':
        ctl['dtype'] = '<f4'
    else:
        raise ValueError(f"Unknown endian {endian}")
    t_def = contents[contents.index('TDEF') + 1:contents.index('TDEF') + 5]
    x_def = contents[contents.index('XDEF') + 1:contents.index('XDEF') + 5]
    y_def = contents[contents.index('YDEF') + 1:contents.index('YDEF') + 5]
    ctl['t_size'] = int(t_def[0])
    ctl['x_size'] = int(x_def[0])
    ctl['y_size'] = int(y_def[0])
    ctl['t_start'] = t_def[2]
    ctl['x_start'] = float(x_def[2])
    ctl['y_start'] = float(y_def[2])
    ctl['t_delta'] = t_def[3]
    ctl['x_delta'] = float(x_def[3])
    ctl['y_delta'] = float(y_def[3])
    return ctl


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    tar_file = os.path.join(in_dir, cfg['filename'])
    logger.info("Found input file '%s'", tar_file)
    file_dir = _extract_tar(tar_file, out_dir)

    # Read control file
    ctl = _read_control_file(file_dir, cfg)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        glob_attrs['positive'] = var_info['positive']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        _extract_variable(cmor_info, glob_attrs, file_dir, out_dir, ctl)
    _clean(file_dir)
