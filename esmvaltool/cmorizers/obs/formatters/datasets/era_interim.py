"""ESMValTool CMORizer for ERA-Interim data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    http://apps.ecmwf.int/datasets/data/interim-full-moda/

Last access
    20190905

Download and processing instructions
    Select "ERA Interim Fields":
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
    It is alo possible to use the script in:
    esmvaltool/cmorizers/obs/download_scripts/download_era-interim.py
    This cmorization script currently supports daily and monthly data of
the following variables:
        10m u component of wind
        10m v component of wind
        2m dewpoint temperature
        2m temperature
        evaporation
        maximum 2m temperature since previous post processing
        mean sea level pressure
        minimum 2m temperature since previous post processing
        skin temperature
        snowfall
        surface net solar radiation
        surface solar radiation downwards
        temperature of snow layer
        toa incident solar radiation
        total cloud cover
        total precipitation
and daily, monthly (not invariant) data of:
        Geopotential

and monthly data of:
        Inst. eastward turbulent surface stress
        Inst. northward turbulent surface stress
        Sea surface temperature
        Surface net thermal radiation
        Surface latent heat flux
        Surface sensible heat flux
        Relative humidity
        Temperature
        U component of wind
        V component of wind
        Vertical velocity
        Specific humidity

Caveats
    Make sure to select the right steps for accumulated fluxes, see:
        https://confluence.ecmwf.int/pages/viewpage.action?pageId=56658233
        https://confluence.ecmwf.int/display/CKB/ERA-Interim%3A+monthly+means
    for a detailed explanation.
    The data are updated regularly: recent years are added, but also the past
    years are sometimes corrected. To have a consistent timeseries, it is
    therefore recommended to download the full timeseries and not just add
    new years to a previous version of the data.

For further details on obtaining daily values from ERA-Interim,
    see:
    https://confluence.ecmwf.int/display/CKB/ERA-Interim
    https://confluence.ecmwf.int/display/CKB/ERA-Interim+documentation#ERA-Interimdocumentation-Monthlymeans
    https://confluence.ecmwf.int/display/CKB/ERA-Interim%3A+How+to+calculate+daily+total+precipitation

"""
import logging
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timedelta
from os import cpu_count
from pathlib import Path
from warnings import catch_warnings, filterwarnings

import iris
import numpy as np

from esmvalcore.cmor.table import CMOR_TABLES
from esmvalcore.preprocessor import daily_statistics, monthly_statistics

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _fix_units(cube, definition):
    """Fix issues with the units."""
    if cube.var_name in {'evspsbl', 'pr', 'prsn'}:
        # Change units from meters of water per day
        # to kg of water per m2 per day
        cube.units = 'm'  # fix invalid units
        cube.units = cube.units * 'kg m-3 day-1'
        cube.data = cube.core_data() * 1000.
    if cube.var_name in {'hfds', 'rss', 'rsds', 'rsdt', 'rlds'}:
        # Add missing 'per day'
        cube.units = cube.units * 'day-1'
        # Radiation fluxes are positive in downward direction
        cube.attributes['positive'] = 'down'
    if cube.var_name in {'tauu', 'tauv'}:
        cube.attributes['positive'] = 'down'
    if cube.var_name in {'sftlf', 'clt'}:
        # Change units from fraction to percentage
        cube.units = definition.units
        cube.data = cube.core_data() * 100.
    if cube.var_name in {'zg', 'orog'}:
        # Divide by acceleration of gravity [m s-2],
        # required for geopotential height, see:
        # https://apps.ecmwf.int/codes/grib/param-db?id=129
        cube.units = cube.units / 'm s-2'
        cube.data = cube.core_data() / 9.80665


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    # Make latitude increasing
    cube = cube[..., ::-1, :]

    # Make pressure_level decreasing
    coord_long_name = [item.long_name for item in cube.coords()]
    if 'pressure_level' in coord_long_name:
        cube = cube[:, ::-1, ...]

    # Add scalar height coordinates
    if 'height2m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 2.)
    if 'height10m' in definition.dimensions:
        utils.add_scalar_height_coord(cube, 10.)
    for coord_def in definition.coordinates.values():
        axis = coord_def.axis
        coord = cube.coord(axis=axis)
        if axis == 'T':
            coord.convert_units('days since 1850-1-1 00:00:00.0')
        if axis == 'Z':
            coord.convert_units(coord_def.units)
        coord.standard_name = coord_def.standard_name
        coord.var_name = coord_def.out_name
        coord.long_name = coord_def.long_name
        coord.points = coord.core_points().astype('float64')
        if len(coord.points) > 1:
            coord.guess_bounds()
        if coord.var_name == 'plev':
            coord.attributes['positive'] = 'down'
    return cube


def _fix_monthly_time_coord(cube):
    """Set the monthly time coordinates to the middle of the month."""
    coord = cube.coord(axis='T')
    end = []
    for cell in coord.cells():
        month = cell.point.month + 1
        year = cell.point.year
        if month == 13:
            month = 1
            year = year + 1
        end.append(cell.point.replace(month=month, year=year))
    end = coord.units.date2num(end)
    start = coord.points
    coord.points = 0.5 * (start + end)
    coord.bounds = np.column_stack([start, end])


def _fix_monthly_time_coord_eiland(cube):
    """Set the monthly time coordinates to the middle of the month."""
    coord = cube.coord(axis='T')
    start = []
    end = []
    for cell in coord.cells():
        # set start to first day 00 UTC
        start.append(cell.point.replace(day=1, hour=0))
        # now deal with the end
        month = cell.point.month + 1
        year = cell.point.year
        if month == 13:
            month = 1
            year = year + 1
        end.append(cell.point.replace(month=month, year=year, day=1, hour=0))
    end = coord.units.date2num(end)
    start = coord.units.date2num(start)
    coord.points = 0.5 * (start + end)
    coord.bounds = np.column_stack([start, end])


def _compute_monthly(cube):
    """Convert various frequencies to daily frequency.

    ERA-Interim-Land is in 6hr freq need to convert to monthly

    """
    cube = monthly_statistics(cube, operator='mean')
    # Remove monthly statistics aux coordinates
    cube.remove_coord(cube.coord('month_number'))
    cube.remove_coord(cube.coord('year'))
    return cube


def _compute_daily(cube):
    """Convert various frequencies to daily frequency.

    ERA-Interim is in 3hr or 6hr or 12hr freq need to convert to daily
    Only variables with step 12 need accounting time 00 AM as time 24 PM

    """
    # Account for time 00 AM as time 24 PM
    if cube.var_name in {
            'tasmax',
            'tasmin',
            'pr',
            'rsds',
            'rlds',
            'hfds',
            'evspsbl',
            'rsdt',
            'rss',
            'prsn',
    }:
        cube.coord('time').points = cube.coord('time').units.date2num([
            cell.point - timedelta(seconds=1)
            for cell in cube.coord('time').cells()
        ])

    if cube.var_name == 'tasmax':
        cube = daily_statistics(cube, 'max')
    elif cube.var_name == 'tasmin':
        cube = daily_statistics(cube, 'min')
    elif cube.var_name in {
            'pr',
            'rsds',
            'rlds',
            'hfds',
            'evspsbl',
            'rsdt',
            'rss',
            'prsn',
    }:
        cube = daily_statistics(cube, 'sum')
    else:
        cube = daily_statistics(cube, 'mean')

    # Correct the time coordinate
    cube.coord('time').points = cube.coord('time').units.date2num([
        cell.point.replace(hour=12, minute=0, second=0, microsecond=0)
        for cell in cube.coord('time').cells()
    ])
    cube.coord('time').bounds = None
    cube.coord('time').guess_bounds()

    return cube


def _load_cube(in_files, var):
    """Load in_files into an iris cube."""
    ignore_warnings = (
        {
            'raw': 'tcc',
            'units': '(0 - 1)',
        },
        {
            'raw': 'lsm',
            'units': '(0 - 1)',
        },
        {
            'raw': 'e',
            'units': 'm of water equivalent',
        },
        {
            'raw': 'sf',
            'units': 'm of water equivalent',
        },
        {
            'raw': 'tp',
            'units': 'm of water equivalent',
        },
    )

    with catch_warnings():
        msg = "Ignoring netCDF variable '{raw}' invalid units '{units}'"
        for warning in ignore_warnings:
            filterwarnings(action='ignore',
                           message=re.escape(msg.format(**warning)),
                           category=UserWarning,
                           module='iris')

        if len(in_files) == 1:
            cube = iris.load_cube(
                in_files[0],
                constraint=utils.var_name_constraint(var['raw']),
            )
        elif var.get('operator', '') == 'sum':
            # Multiple variables case using sum operation
            cube = None
            for raw_name, filename in zip(var['raw'], in_files):
                in_cube = iris.load_cube(
                    filename,
                    constraint=utils.var_name_constraint(raw_name),
                )
                if cube is None:
                    cube = in_cube
                else:
                    cube += in_cube
        else:
            raise ValueError(
                "Multiple input files found, with operator '{}' configured: {}"
                .format(var.get('operator'), ', '.join(in_files)))

    return cube


def _extract_variable(in_files, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input files '%s'",
                var['short_name'], ', '.join(in_files))
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = _load_cube(in_files, var)

    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    _fix_units(cube, definition)

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    cube = _fix_coordinates(cube, definition)

    if attributes['dataset_id'] == 'ERA-Interim':
        if 'mon' in var['mip']:
            _fix_monthly_time_coord(cube)
        if 'day' in var['mip']:
            cube = _compute_daily(cube)
        if 'fx' in var['mip']:
            cube = iris.util.squeeze(cube)
            cube.remove_coord('time')

    # Specific to ERA Interim Land
    elif attributes['dataset_id'] == 'ERA-Interim-Land':
        if 'mon' in var['mip']:
            cube = _compute_monthly(cube)
            _fix_monthly_time_coord_eiland(cube)
        if 'day' in var['mip']:
            cube = _compute_daily(cube)
    else:
        raise ValueError("Unknown dataset_id for this script:\
                         {attributes['dataset_id']}")

    # Convert units if required
    cube.convert_units(definition.units)

    logger.debug("Saving cube\n%s", cube)
    logger.debug("Expected output size is %.1fGB",
                 np.prod(cube.shape) * 4 / 2**30)
    utils.save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        local_keys=['positive'],
    )
    logger.info("Finished CMORizing %s", ', '.join(in_files))


def _get_in_files_by_year(in_dir, var):
    """Find input files by year."""
    if 'file' in var:
        var['files'] = [var.pop('file')]

    in_files = defaultdict(list)
    for pattern in var['files']:
        for filename in Path(in_dir).glob(pattern):
            year = str(filename.stem).split('_')[-1]
            in_files[year].append(str(filename))

    # Check if files are complete
    for year in in_files.copy():
        if len(in_files[year]) != len(var['files']):
            logger.warning(
                "Skipping CMORizing %s for year '%s', %s input files needed, "
                "but found only %s", var['short_name'], year,
                len(var['files']), ', '.join(in_files[year]))
            in_files.pop(year)

    return in_files.values()


def _run(jobs, n_workers):
    """Run CMORization jobs using n_workers."""
    if n_workers == 1:
        for job in jobs:
            _extract_variable(*job)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            for job in jobs:
                future = executor.submit(_extract_variable, *job)
                futures[future] = job[0]

            for future in as_completed(futures):
                try:
                    future.result()
                except:  # noqa
                    logger.error("Failed to CMORize %s",
                                 ', '.join(futures[future]))
                    raise


def cmorization(in_dir, out_dir, cfg, config_user):
    """Run CMORizer for ERA-Interim."""
    cfg['attributes']['comment'] = cfg['attributes']['comment'].strip().format(
        year=datetime.now().year)
    cfg.pop('cmor_table')

    n_workers = config_user.get('max_parallel_tasks')
    if n_workers is None:
        n_workers = int(cpu_count() / 1.5)
    logger.info("Using at most %s workers", n_workers)

    jobs = []
    for short_name, var in cfg['variables'].items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        for in_files in _get_in_files_by_year(in_dir, var):
            jobs.append([in_files, var, cfg, out_dir])

    _run(jobs, n_workers)
