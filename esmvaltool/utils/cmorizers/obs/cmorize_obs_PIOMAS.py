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

NX = 360
NY = 120


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    scalar_path = os.path.join(in_dir, CFG['custom']['scalar_file'])
    grids = np.loadtxt(scalar_path).reshape(2, NY, NX)
    lat_sca, lon_sca = _create_lat_lon_coords(grids[1, ...], grids[0, ...])

    vector_path = os.path.join(in_dir, CFG['custom']['vector_file'])
    grids = np.loadtxt(vector_path).reshape(7, NY, NX)
    lat_vec, lon_vec = _create_lat_lon_coords(grids[1, ...], grids[0, ...])
    # Area in m2
    area_cello = grids[2, ...] * grids[3, ...] * 1e6
    # run the cmorization
    for var, vals in CFG['variables'].items():

        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        if vals['type'] == 'scalar':
            lat = lat_sca
            lon = lon_sca
        else:
            lat = lat_vec
            lon = lon_vec

        if var == "areacello":
            cube = _create_areacello(lon, lat, area_cello, var_info)
            set_global_atts(cube, glob_attrs)
            save_variable(cube, var_info.short_name, out_dir, glob_attrs)
            continue

        file_expression = os.path.join(in_dir, '{0}.H????'.format(vals['raw']))
        for file_path in glob.glob(file_expression):
            data = _read_binary_file(
                file_path
            )
            year = int(file_path[-4:])
            cube = _create_cube(lon, lat, year, data, var_info)
            set_global_atts(cube, glob_attrs)
            save_variable(cube, var_info.short_name, out_dir, glob_attrs)


def _create_lat_lon_coords(lat, lon):
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
    return lat_coord, lon_coord


def _create_cube(lon, lat, year, data, var_info):
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
    cube.add_aux_coord(lon, (1, 2))
    cube.add_aux_coord(lat, (1, 2))
    return cube


def _create_areacello(lon, lat, data, var_info):
    cube = iris.cube.Cube(
        data,
        standard_name=var_info.standard_name,
        var_name=var_info.short_name,
        units='m2',
    )
    cube.add_aux_coord(lon, (0, 1))
    cube.add_aux_coord(lat, (0, 1))
    return cube


def _read_binary_file(data_path):
    fd_data = open(data_path, 'rb')
    data = np.fromfile(fd_data, dtype=np.dtype('f'), count=-1)
    days = data.shape[0] // NX // NY
    data = data.reshape(days, NY, NX)
    return data
