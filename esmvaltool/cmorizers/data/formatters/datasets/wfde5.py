"""ESMValTool CMORizer for WFDE5 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://doi.org/10.24381/cds.20d54e34

Last access
    20210416

Download and processing instructions
    Download the following variables from the cds:
        Near-surface air temperature ("CRU")
        Rainfall flux ("CRU" as well as "CRU and GPCC")
        Snowfall flux ("CRU" as well as "CRU and GPCC")
    unzip the downloaded files
    rename to follow syntax '{raw_name}_WFDE5_{reference}_*_v1.0.nc'
"""

import copy
import logging
import re
from pathlib import Path

import iris
from cf_units import Unit
from esmvalcore.preprocessor import daily_statistics, monthly_statistics

from ... import utilities as utils

logger = logging.getLogger(__name__)


def _fix_time_coord(cube, var):
    """Correct wrong time points."""
    # Fix units
    cube.coord('time').units = Unit(cube.coord('time').units.origin,
                                    calendar=var['calendar'])
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    # time points are XX:00:00, should be XX:30:00
    time = cube.coord('time')
    time.points = time.points + 1 / 48


def _extract_variable(var, cfg, filenames, out_dir):
    """Extract variable."""
    short_name = var['short_name']
    version = cfg['attributes']['version'] + '-' + var['reference']

    cubes = iris.load(filenames)
    cube = cubes.concatenate_cube()
    if short_name == 'pr':
        # Rainf add Snowf
        snow_filenames = [
            fname.replace(var['raw_name'], var['raw_name_snow'])
            for fname in filenames
        ]
        scubes = iris.load(snow_filenames)
        scube = scubes.concatenate_cube()
        cube.data = scube.core_data() + cube.core_data()

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    if 'raw_units' in var:
        cube.units = var['raw_units']
    cube.convert_units(cmor_info.units)

    # Fix time coord
    _fix_time_coord(cube, var)

    # Fix coordinates
    utils.fix_coords(cube)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    attrs['version'] = version
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    logger.info("Building daily means")
    # Calc daily
    cube = daily_statistics(cube)
    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    if var['add_mon']:
        logger.info("Building monthly means")

        # Calc monthly
        cube = monthly_statistics(cube)
        cube.remove_coord('month_number')
        cube.remove_coord('year')

        # Fix metadata
        attrs['mip'] = 'Amon'
        utils.set_global_atts(cube, attrs)

        # Save variable
        utils.save_variable(cube,
                            short_name,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for var in cfg['variables'].values():
        logger.info("CMORizing variable '%s'", var['short_name'])
        file_names = cfg['filename'].format(**var)

        sorted_filenames = {}
        # Sort files according to year
        pattern = re.compile(r'.*_(\d{4})\d{2}_.*\.nc')
        for filename in sorted(Path(in_dir).glob(file_names)):
            filename = str(filename)
            year = int(pattern.match(filename)[1])
            sorted_filenames.setdefault(year, [])
            sorted_filenames[year].append(filename)

        # Run CMORization for each year
        for (year, filenames) in sorted_filenames.items():
            logger.info("Processing year %i", year)
            _extract_variable(var, cfg, filenames, out_dir)
