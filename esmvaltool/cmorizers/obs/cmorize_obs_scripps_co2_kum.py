"""ESMValTool CMORizer for Scripps-CO2-KUM data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://scrippsco2.ucsd.edu/data/atmospheric_co2/kum.html

Last access
    20200422

Download and processing instructions
    Download the following file:
        monthly_flask_co2_kum.csv

"""

import logging
import os
from datetime import datetime

import iris
import numpy as np
import pandas as pd
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


LAT_COORD = iris.coords.DimCoord([19.5], var_name='lat',
                                 standard_name='latitude',
                                 long_name='latitude', units='degrees')
LON_COORD = iris.coords.DimCoord([205.2], var_name='lon',
                                 standard_name='longitude',
                                 long_name='longitude', units='degrees')


def _get_time_coord(year, month):
    """Get time coordinate."""
    point = datetime(year=year, month=month, day=15)
    bound_low = datetime(year=year, month=month, day=1)
    if month == 12:
        month_bound_up = 1
        year_bound_up = year + 1
    else:
        month_bound_up = month + 1
        year_bound_up = year
    bound_up = datetime(year=year_bound_up, month=month_bound_up, day=1)
    time_units = Unit('days since 1950-01-01 00:00:00', calendar='standard')
    time_coord = iris.coords.DimCoord(
        time_units.date2num(point),
        bounds=time_units.date2num([bound_low, bound_up]),
        var_name='time',
        standard_name='time',
        long_name='time',
        units=time_units,
    )
    return time_coord


def _get_cube(row, column_name, fill_value):
    """Create :class:`iris.cube.Cube` from :class:`pandas.Series`."""
    time_coord = _get_time_coord(int(row['Yr']), int(row['Mn']))
    lat_coord = LAT_COORD.copy()
    lon_coord = LON_COORD.copy()
    data = np.ma.masked_equal(row[tuple(column_name)], fill_value)
    cube = iris.cube.Cube(
        data.reshape((1, 1, 1)),
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        units='ppm',
    )
    return cube


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    data_frame = pd.read_csv(filepath, comment='"', header=[0, 1, 2])
    data_frame = data_frame.rename(columns=lambda x: x.strip())

    # Extract cube
    cubes = iris.cube.CubeList()
    for (_, row) in data_frame.iterrows():
        cube = _get_cube(row, var['column_name'], cfg['fill_value'])
        cubes.append(cube)
    cube = cubes.concatenate_cube()
    cube.var_name = short_name

    # Fix metadata
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    filepath = os.path.join(in_dir, cfg['filename'])
    logger.info("Reading file '%s'", filepath)

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, filepath, out_dir)
