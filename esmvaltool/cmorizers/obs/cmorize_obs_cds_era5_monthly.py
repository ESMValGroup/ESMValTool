"""ESMValTool CMORizer for ERA5 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means

Last access
    20190830

Download and processing instructions
    TODO

History
    20190902 crez_ba adapted from cmorize_obs_era5.py
"""

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime
from os import cpu_count
from pathlib import Path
from warnings import catch_warnings, filterwarnings

import iris
import numpy as np

from esmvalcore.cmor.table import CMOR_TABLES

from . import utilities as utils

logger = logging.getLogger(__name__)


def _guess_bnds_time_monthly(cube):
    import cf_units
    import calendar
    from dateutil import relativedelta
    coord = cube.coord('time')
    time_as_datetime = cf_units.num2date(coord.points,
        coord.units.origin,coord.units.calendar)
    time_bnds = []
    for timestep in time_as_datetime:
        _,monthend = calendar.monthrange(timestep.year,timestep.month)
        a_datetime = datetime(timestep.year,timestep.month,1,0,0)
        b_datetime = a_datetime + relativedelta.relativedelta(months=1)
        bnd_a = int(cf_units.date2num(a_datetime,coord.units.origin,
                                  coord.units.calendar))
        bnd_b = int(cf_units.date2num(b_datetime,coord.units.origin,
                                  coord.units.calendar))
        time_bnds.append([bnd_a,bnd_b])
    # Convert them to an array
    time_bnds = np.array(time_bnds)
    # Set them as bounds
    coord.bounds = time_bnds

def _extract_variable(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    with catch_warnings():
        filterwarnings(
            action='ignore',
            message="Ignoring netCDF variable 'tcc' invalid units '(0 - 1)'",
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(
            str(in_file),
            constraint=utils.var_name_constraint(var['raw']),
        )

    # Set correct names
    cube.var_name = definition.short_name
    try:
        cube.standard_name = definition.standard_name
    except ValueError:
        logger.error("Failed setting standard_name for variable short_name %s",cube.var_name)
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
    cube.coord('time').points = cube.coord('time').core_points().astype('float64')
    # Set time bounds by setting manually to month start; month end
    _guess_bnds_time_monthly(cube)
    
    # Convert units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    logger.info("Expected output size is %.1fGB",
                np.prod(cube.shape) * 4 / 2**30)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

#    n_workers = int(cpu_count() / 1.5)
#    logger.info("Using at most %s workers", n_workers)
#    futures = {}
#    with ProcessPoolExecutor(max_workers=1) as executor:
    if True:
        for short_name, var in cfg['variables'].items():
            var['short_name'] = short_name
            for in_file in sorted(Path(in_dir).glob(var['file'])):
#                future = executor.submit(_extract_variable, in_file, var, cfg,
#                                         out_dir)
                _extract_variable(in_file,var,cfg,out_dir)
#                futures[future] = in_file

#    for future in as_completed(futures):
#        try:
#            future.result()
#        except:  # noqa
#            logger.error("Failed to CMORize %s", futures[future])
#            raise
#        logger.info("Finished CMORizing %s", futures[future])
