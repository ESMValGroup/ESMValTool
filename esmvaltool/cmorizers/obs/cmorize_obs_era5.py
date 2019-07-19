"""ESMValTool CMORizer for ERA5 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels

Last access
    20190718

Download and processing instructions
    This cmorization script currently supports hourly data of the following
variables:
        2m_temperature
        potential_evaporation
        total_precipitation

    Downloading ERA5 data can either be done via the Climate Data Store (cds)
web form or era5cli:
        $pip install era5cli
        $era5cli hourly --variable total_precipitation --startyear 1990

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
        if cube.var_name == 'tcc':
            # Change cloud cover units from fraction to percentage
            cube.units = definition.units
            cube.data = cube.core_date() * 100.
        if cube.var_name in ['tp', 'pev']:
            # Change units from meters of water to kg of water
            # and add missing 'per hour'
            cube.units = cube.units * 'kg m-3 h-1'
            cube.data = cube.core_data() * 1000.

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Fix coordinates
    cube.coord('latitude').var_name = 'lat'
    cube.coord('longitude').var_name = 'lon'

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float64')
        coord.guess_bounds()

    if 'height2m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 2.)
    if 'height10m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 10.)

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


def cmorization(in_dir, out_dir, cfg):
    """Cmorization func call."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

    n_workers = int(cpu_count() / 1.5)
    logger.info("Using at most %s workers", n_workers)
    futures = {}
    with ProcessPoolExecutor(max_workers=1) as executor:
        for short_name, var in cfg['variables'].items():
            var['short_name'] = short_name
            for in_file in sorted(Path(in_dir).glob(var['file'])):
                future = executor.submit(_extract_variable, in_file, var, cfg,
                                         out_dir)
                futures[future] = in_file

    for future in as_completed(futures):
        try:
            future.result()
        except:  # noqa
            logger.error("Failed to CMORize %s", futures[future])
            raise
        logger.info("Finished CMORizing %s", in_file)
