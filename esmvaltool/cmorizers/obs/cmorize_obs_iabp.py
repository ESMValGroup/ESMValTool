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

import numpy as np
import numpy.ma as ma
from cf_units import Unit
import iris
from iris.cube import Cube, CubeList
from iris.coords import AuxCoord, DimCoord
from iris.coord_categorisation import add_day_of_year
from iris.analysis import MEAN
import datetime

from . import utilities as utils

logger = logging.getLogger(__name__)


def _save_variable(cube, cmor_info, attrs, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    utils.fix_var_metadata(cube, cmor_info)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(
        cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    time_unit = Unit('days since 1850-01-01', 'gregorian')

    cmor_info = dict()
    for (var, var_info) in cfg['variables'].items():
        glob_attrs['mip'] = var_info['mip']
        cmor_info[var] = cmor_table.get_variable(var_info['mip'], var)

    file_list = glob.glob(os.path.join(in_dir, 'D????')) + \
        glob.glob(os.path.join(in_dir, 'D????.dat'))
    file_list.sort()

    for filepath in file_list:
        logger.info('Cmorizing %s', filepath)
        current_date = None
        with open(filepath) as file_handle:
            usi_cubes = CubeList()
            vsi_cubes = CubeList()
            for line in file_handle.readlines():
                line = line.split(' ')
                line = [string for string in line if string]
                hour = int(line[3])
                date = _get_date(line, hour)

                if current_date != date:
                    if current_date:
                        latitude = AuxCoord(
                            lat, 'latitude', 'latitude', 'lat', 'degrees_north'
                        )
                        longitude = AuxCoord(
                            lon, 'longitude', 'longitude', 'lon', 'degrees_east'
                        )
                        time = AuxCoord(
                            time_unit.date2num(current_date),
                            'time', 'time', 'time', time_unit
                        )
                        index = DimCoord(range(len(lat)),
                                         None, 'Cell index', 'i')
                        usi = ma.array(usi, mask=mask)
                        usi_cube = Cube(usi, var_name='usi', units='m s-1')
                        usi_cube.add_dim_coord(index, 0)
                        usi_cube.add_aux_coord(time)
                        usi_cube.add_aux_coord(latitude, (0, ))
                        usi_cube.add_aux_coord(longitude, (0, ))
                        usi_cubes.append(usi_cube)

                        vsi = ma.array(usi, mask=mask)
                        vsi_cube = Cube(usi, var_name='vsi', units='m s-1')
                        vsi_cube.add_dim_coord(index, 0)
                        vsi_cube.add_aux_coord(time)
                        vsi_cube.add_aux_coord(latitude, (0, ))
                        vsi_cube.add_aux_coord(longitude, (0, ))
                        vsi_cubes.append(vsi_cube)

                    lat = []
                    lon = []
                    usi = []
                    vsi = []
                    mask = []
                    current_date = date

                lat.append(float(line[4]))
                lon.append(float(line[5]))
                usi.append(float(line[6]) / 100)
                vsi.append(float(line[7]) / 100)
                # SIGMA2   is the variance of the interpolation error in
                # velocity, in dimensionless units.
                # No confidence should be placed in interpolated velocities
                # for which SIGMA2 > 0.5.
                sigma2 = float(line[8])
                mask.append(sigma2 > 0.5)


        if current_date:
            year = current_date.year
            logger.info('Creating file for year %s', year)
            current_date = datetime.datetime(year, 1, 1)
            delta = datetime.timedelta(days=1)
            while current_date.year == year:
                if not any(((cube.coord('time').cell(0) == current_date) for cube in usi_cubes)):
                    time = AuxCoord(
                        time_unit.date2num(current_date),
                        'time', 'time', 'time', time_unit
                    )
                    empty_cube = usi_cubes[0].copy(ma.masked_invalid(np.full(len(lat), np.nan)))
                    empty_cube.remove_coord('time')
                    empty_cube.add_aux_coord(time)
                    usi_cubes.append(empty_cube)

                    empty_cube = vsi_cubes[0].copy(ma.masked_invalid(np.full(len(lat), np.nan)))
                    empty_cube.remove_coord('time')
                    empty_cube.add_aux_coord(time)
                    vsi_cubes.append(empty_cube)
                current_date += delta

        usi = _merge_cube(usi_cubes)
        vsi = _merge_cube(vsi_cubes)

        glob_attrs['mip'] = cfg['variables']['usi']['mip']
        _save_variable(usi, cmor_info['usi'], glob_attrs, out_dir)

        glob_attrs['mip'] = cfg['variables']['vsi']['mip']
        _save_variable(vsi, cmor_info['vsi'], glob_attrs, out_dir)


def _get_date(line, hour):
    year = int(line[0])
    if year < 1900:
        year += 1900
    month = int(line[1])
    day = int(line[2])
    date = datetime.datetime(year, month, day, hour)
    return date

def _merge_cube(cube_list):
    cube = cube_list.merge_cube()
    add_day_of_year(cube, 'time')
    cube = cube.aggregated_by('day_of_year', MEAN)
    cube.remove_coord('day_of_year')
    return cube