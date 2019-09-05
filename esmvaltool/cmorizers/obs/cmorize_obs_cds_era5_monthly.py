"""ESMValTool CMORizer for ERA5 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means

Last access
    20190830

Download and processing instructions
    The data can be downloaded using the CDS API:
        https://cds.climate.copernicus.eu/api-how-to

    The retrieval command is the following and should be embedded in a
    loop over all years and months.

        datarequest = cds.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
                'variable' : 'swvl1,
                'grid' : ['0.25','0.25'],
                'product_type': 'monthly_averaged_reanalysis',
                'year': year,
                'month': month,
                'time': ['00:00'],
                'format':'netcdf'
                }
        )


History
    20190902 crez_ba adapted from cmorize_obs_era5.py
"""

import logging
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import cf_units
import iris
import numpy as np

from esmvalcore.cmor.table import CMOR_TABLES

from . import utilities as utils

logger = logging.getLogger(__name__)


def _guess_bnds_time_monthly(cube):
    """Guess time bounds from time points for monthly data."""
    from dateutil import relativedelta
    coord = cube.coord('time')
    time_as_datetime = cf_units.num2date(coord.points,
                                         coord.units.origin,
                                         coord.units.calendar)
    time_bnds = []
    time_points = []
    for timestep in time_as_datetime:
        # Deal with time points
        midmonth_datetime = datetime(timestep.year, timestep.month, 15)
        midmonth_point = cf_units.date2num(midmonth_datetime,
            coord.units.origin, coord.units.calendar)
        time_points.append(midmonth_point)
        # Deal with time bounds
        a_datetime = datetime(timestep.year, timestep.month, 1)
        b_datetime = a_datetime + relativedelta.relativedelta(months=1)
        bnd_a = cf_units.date2num(a_datetime, coord.units.origin,
                                  coord.units.calendar)
        bnd_b = cf_units.date2num(b_datetime, coord.units.origin,
                                  coord.units.calendar)
        time_bnds.append([bnd_a, bnd_b])
    # Convert them to an array
    time_bnds = np.array(time_bnds)
    time_points = np.array(time_points)
    # Set them as bounds
    coord.bounds = time_bnds
    coord.points = time_points


def _extract_variable(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = iris.load_cube(
        str(in_file),
        constraint=utils.var_name_constraint(var['raw']))

    # Set correct names
    cube.var_name = definition.short_name
    try:
        cube.standard_name = definition.standard_name
    except ValueError:
        logger.error("Failed setting standard_name for variable\
                     short_name %s", cube.var_name)
    cube.long_name = definition.long_name

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Fix lat/lon coordinates
    cube.coord('latitude').var_name = 'lat'
    cube.coord('longitude').var_name = 'lon'
    for coord_name in 'latitude', 'longitude':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float64')
        coord.guess_bounds()

    # Fix time
    cube.coord('time').points = cube.coord('time').\
        core_points().astype('float64')
    # Set time bounds by setting manually to month start; month end
    _guess_bnds_time_monthly(cube)

    # Now convert them to the recommended calendar
    cube.coord('time').convert_units(cf_units.Unit(
        'days since 1950-1-1 00:00:00', calendar='gregorian'))

    # Convert units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        for in_file in sorted(Path(in_dir).glob(var['file'])):
            _extract_variable(in_file, var, cfg, out_dir)
        logger.info("Finished CMORizing variable %s", short_name)
