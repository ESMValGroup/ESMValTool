# pylint: disable=invalid-name
"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid

Last access
   20190510

Download and processing instructions
   Download the desired files from:
   https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/

"""

import logging
import os
import glob

import numpy as np
from cf_units import Unit
import iris
from iris.coords import AuxCoord, DimCoord

from .utilities import (constant_metadata, fix_coords,
                        fix_var_metadata, read_cmor_config, save_variable,
                        set_global_atts)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = read_cmor_config('PIOMAS.yml')


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in CFG['variables'].items():
        grid_path = os.path.join(in_dir, vals['grid_file'])
        file_expression = os.path.join(in_dir, '{0}.H????'.format(vals['raw']))
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        for file_path in glob.glob(file_expression):
            lon, lat, data = read_binary_file(
                file_path, grid_path
            )
            year = int(file_path[-4:])
            cube = create_cube(lon, lat, year, data, var_info)
            set_global_atts(cube, glob_attrs)
            save_variable(cube, var_info.short_name, out_dir, glob_attrs)


def create_cube(lon, lat, year, data, var_info):
    lon_coord = AuxCoord(
        lon,
        standard_name='longitude',
        var_name='lon',
        units='degrees_east'
    )

    lat_coord = AuxCoord(
        lat,
        standard_name='latitude',
        var_name='lat',
        units='degrees_north'
    )

    time_coord = DimCoord(
        np.arange(0, data.shape[0]),
        standard_name='time',
        var_name='time',
        units=Unit('days since {}-01-01'.format(year), calendar='noleap'),
    )

    cube = iris.cube.Cube(
        data,
        standard_name=var_info.standard_name,
        var_name=var_info.short_name,
        units='m',
    )
    cube.add_dim_coord(time_coord, 0)
    cube.add_aux_coord(lon_coord, (1, 2))
    cube.add_aux_coord(lat_coord, (1, 2))
    return cube


def read_binary_file(data_path, grid_path):
    nx = 360
    ny = 120

    grids = np.loadtxt(grid_path)
    grids = grids.reshape(2, ny, nx)
    lon = grids[0, ...]
    lat = grids[1, ...]

    fd_data = open(data_path, 'rb')
    data = np.fromfile(fd_data, dtype=np.dtype('f'), count=-1)
    days = data.shape[0] // nx // ny
    data = data.reshape(days, ny, nx)

    return lon, lat, data
