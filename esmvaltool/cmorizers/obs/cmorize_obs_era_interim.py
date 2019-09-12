"""ESMValToo CMORizer for ERA-Interim data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    http://apps.ecmwf.int/datasets/data/interim-full-moda/

Last access
    20190905

Download and processing instructions
    Select "Era Interim Fields":
        Daily: for daily values
        Invariant: for time invariant variables (like land-sea mask)
        Monthly Means of Daily Means: for monthly values
        Monthly Means of Daily Forecast Accumulation: for accumulated variables
        like precipitation or radiation fluxes
    Select "Type of level" (Surface or Pressure levels)
    Download the data on a single variable and single year basis, and save
    them as ERA-Interim_<var>_<mean>_YYYY.nc, where <var> is the ERA-Interim
    variable name and <mean> is either monthly or daily. Further download
    "land-sea mask" from the "Invariant" data and save it in
    ERA-Interim_lsm.nc.
    It is also possible to download data in an automated way, see:
        https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
        https://confluence.ecmwf.int/display/WEBAPI/Python+ERA-interim+examples
    A registration is required for downloading the data.

Caveats
    Make sure to select the right steps for accumulated fluxes, see:
        https://confluence.ecmwf.int/pages/viewpage.action?pageId=56658233
        https://confluence.ecmwf.int/display/CKB/ERA-Interim%3A+monthly+means
    for a detailed explanation.
    The data are updated regularly: recent years are added, but also the past
    years are sometimes corrected. To have a consistent timeseries, it is
    therefore recommended to download the full timeseries and not just add
    new years to a previous version of the data.

"""
from datetime import datetime
import logging
from concurrent.futures import as_completed, ProcessPoolExecutor
from copy import deepcopy
from os import cpu_count
from pathlib import Path
from warnings import catch_warnings, filterwarnings

import iris
import numpy as np

from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.preprocessor import extract_month, daily_statistics
from iris.coord_categorisation import add_month_number

from . import utilities as utils

logger = logging.getLogger(__name__)

# Acceleration of gravity [m s-2],
# required for surface geopotential height, see https://confluence.ecmwf.int/pages/viewpage.action?pageId=79955800
G = 9.80665


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
        filterwarnings(
            action='ignore',
            message=("Ignoring netCDF variable 'e' invalid units "
                     "'m of water equivalent'"),
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(
            str(in_file),
            constraint=utils.var_name_constraint(var['raw']),
        )

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    if cube.var_name in {'e', 'sf'}:
        # Change evaporation and snowfall units from
        # 'm of water equivalent' to m
        cube.units = 'm'
    if cube.var_name in {'e', 'ro', 'sf', 'tp', 'pev'}:
        # Change units from meters per day of water to kg of water per h
        cube.units = cube.units * 'kg m-3 h-1'
        cube.data = cube.core_data() * 1000. / 24
    if cube.var_name in {'ssr', 'ssrd', 'tisr'}:
        # Add missing 'per hour'
        cube.units = cube.units * 'h-1'
    if cube.var_name in {'lsm', 'tcc'}:
        # Change cloud cover units from fraction to percentage
        cube.units = definition.units
        cube.data = cube.core_data() * 100
    if cube.var_name in {'z'}:
        cube.data = cube.core_data() / G

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

    # era-interim is in 3hr or 6hr or 12hr freq need to convert to daily
    if var['mip'] in {'day', 'Eday'}:
        if cube.var_name == 'tasmax':
            daily_statistics(cube, 'max')
        elif cube.var_name == 'tasmin':
            daily_statistics(cube, 'min')
        elif cube.var_name in {'pr', 'rsds', 'hfds'}:
            # Sum is not available in daily_statistics so call iris directly
            iris.coord_categorisation.add_day_of_year(cube, 'time')
            iris.coord_categorisation.add_year(cube, 'time')
            cube = cube.aggregated_by(['day_of_year', 'year'], iris.analysis.SUM)
        else:
            daily_statistics(cube, 'mean')

    # Add height coordinate to tas variable (required by the new backend)
    if 'height2m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 2.)
    if 'height10m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 10.)

    # Convert units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # For daily data write a netcdf for each month
    if var['mip'] in {'day', 'Eday'}:
        add_month_number(cube, 'time')
        for month_number in range(1, 13):
            month_cube = extract_month(cube, month_number)
            month_cube.remove_coord(cube.coord('month_number'))
            logger.info("Saving cube\n%s", month_cube)
            logger.info("Expected output size is %.1fGB",
                        np.prod(month_cube.shape) * 4 / 2**30)
            utils.save_variable(month_cube, month_cube.var_name, out_dir, attributes)
    elif var['mip'] in {'Amon'}:
        logger.info("Saving cube\n%s", cube)
        logger.info("Expected output size is %.1fGB",
                    np.prod(cube.shape) * 4 / 2 ** 30)
        utils.save_variable(cube, cube.var_name, out_dir, attributes)
    elif var['mip'] == 'fx':
        logger.info("Saving cube\n%s", cube)
        logger.info("Expected output size is %.1fGB",
                    np.prod(cube.shape) * 4 / 2 ** 30)
        utils.save_variable(cube, cube.var_name, out_dir, attributes)


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

    n_workers = int(cpu_count() / 1.5)
    logger.info("Using at most %s workers", n_workers)
    futures = {}
    with ProcessPoolExecutor(max_workers=1) as executor:
        for short_name, var in cfg['variables'].items():
            if 'short_name' not in var:
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
        logger.info("Finished CMORizing %s", futures[future])
