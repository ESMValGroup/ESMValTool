"""ESMValTool CMORizer for IABP data.

Tier
    Tier 2: restricted dataset.

Source
    http://iabp.apl.washington.edu/

Last access
    20190627

Download and processing instructions
    Download the full http://iabp.apl.washington.edu/Data_Products/D/ folder

"""

import logging
import os
import glob
import calendar

import math
import numpy as np
import numpy.ma as ma
from cf_units import Unit
import iris
from iris.cube import Cube, CubeList
from iris.coords import AuxCoord, DimCoord
from iris.coord_categorisation import add_day_of_year
from iris.analysis import MEAN
from iris.analysis.cartography import DEFAULT_SPHERICAL_EARTH_RADIUS
import datetime

from . import utilities as utils

logger = logging.getLogger(__name__)


def _save_variable(cube, cmor_info, attrs, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(
        cube, var, out_dir, attrs, zlib=True
    )


def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    time_unit = Unit('days since 1850-01-01', 'gregorian')

    reference = iris.load_cube(cfg['reference_file'])
    shape = reference.coord('latitude').shape
    logger.debug(reference)

    cmor_info = cmor_table.get_variable('fx', 'areacello')
    glob_attrs['mip'] = 'fx'
    _save_variable(reference, cmor_info, glob_attrs, out_dir)

    cmor_info = dict()
    for (var, var_info) in cfg['variables'].items():
        glob_attrs['mip'] = var_info['mip']
        cmor_info[var] = cmor_table.get_variable(var_info['mip'], var)

    file_template = 'im_from_buoy_nh_list_{0}0101_{0}1231_{1}.nc'
    for year in range(cfg['start_year'], cfg['end_year'] + 1):
        logger.info('Cmorizing %s', year)

        cubes = iris.load(os.path.join(
            in_dir,
            file_template.format(year, glob_attrs['version']
        )))

        usi = cubes.extract_strict('sea_ice_x_velocity')
        vsi = cubes.extract_strict('sea_ice_y_velocity')
        time = cubes.extract_strict('time')
        x = cubes.extract_strict('x on 25km EASEgrid')
        y = cubes.extract_strict('y on 25km EASEgrid')

        usi_grid = _get_empy_array(year, shape)
        vsi_grid = _get_empy_array(year, shape)

        first_day = int(time.units.date2num(datetime.datetime(year, 1, 1)))

        for i in range(time.shape[0]):
            time_coord = int(round(time[i].data - first_day))
            x_coord = int(x[i].data)
            y_coord = int(y[i].data)
            usi_grid[time_coord, y_coord, x_coord] = usi[i].data / 100.
            vsi_grid[time_coord, y_coord, x_coord] = usi[i].data / 100.

        if calendar.isleap(year):
            time_coord = range(first_day, first_day + 366)
        else:
            time_coord = range(first_day, first_day + 365)

        time_coord = DimCoord(
            time_coord, 'time', 'time', 'time', units=time.units
        )

        glob_attrs['mip'] = cfg['variables']['usi']['mip']
        _create_var_file('usi', usi_grid, cmor_info['usi'], glob_attrs, out_dir, reference, time_coord)

        glob_attrs['mip'] = cfg['variables']['vsi']['mip']
        _create_var_file('vsi', usi_grid, cmor_info['vsi'], glob_attrs, out_dir, reference, time_coord)


def _get_empy_array(year, grid_shape):
    if calendar.isleap(year):
        data = np.full([366, grid_shape[0], grid_shape[1]], np.nan, dtype=np.float32)
    else:
        data = np.full([365, grid_shape[0], grid_shape[1]], np.nan, dtype=np.float32)
    return data


def _create_var_file(var_name, data, cmor_info, glob_attrs, out_dir, reference, time_coord):
    cube = Cube(
        ma.masked_invalid(data),
        cmor_info.standard_name,
        var_name=var_name,
        units='m s-1'
    )

    cube.add_dim_coord(time_coord, 0)
    cube.add_dim_coord(reference.coord('projection_y_coordinate'), 1)
    cube.add_dim_coord(reference.coord('projection_x_coordinate'), 2)

    cube.add_aux_coord(reference.coord('latitude'), (1, 2))
    cube.add_aux_coord(reference.coord('longitude'), (1, 2))

    _save_variable(cube, cmor_info, glob_attrs, out_dir)