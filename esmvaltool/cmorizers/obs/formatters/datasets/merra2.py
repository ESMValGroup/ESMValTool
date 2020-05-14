"""ESMValTool CMORizer for MERRA2 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/

Last access
    20191129

Download and processing instructions
    - For download instructions see the download script `download_merra2.sh`.

"""
import glob
import logging
import os
from copy import deepcopy
from datetime import datetime

import cf_units
import iris
from dask import array as da

from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _fix_time_monthly(cube):
    """Fix time by setting it to 15th of month."""
    # Read dataset time unit and calendar from file
    dataset_time_unit = str(cube.coord('time').units)
    dataset_time_calender = cube.coord('time').units.calendar
    # Convert datetime
    time_as_datetime = cf_units.num2date(cube.coord('time').core_points(),
                                         dataset_time_unit,
                                         dataset_time_calender)
    newtime = []
    for timepoint in time_as_datetime:
        midpoint = datetime(timepoint.year, timepoint.month, 15)
        newtime.append(midpoint)

    newtime = cf_units.date2num(newtime,
                                dataset_time_unit,
                                dataset_time_calender)
    # Put them on the file
    cube.coord('time').points = newtime
    cube.coord('time').bounds = None
    return cube


def _load_cube(in_files, var):
    cube_list = iris.load_raw(in_files)
    selected = [c for c in cube_list if c.var_name == var['raw']]
    selected = iris.cube.CubeList(selected)

    drop_attrs = ['History', 'Filename', 'Comment', 'RangeBeginningDate',
                  'RangeEndingDate', 'GranuleID', 'ProductionDateTime',
                  'Source']
    drop_time_attrs = ['begin_date', 'begin_time',
                       'time_increment', 'valid_range', 'vmax', 'vmin']
    for cube in selected:
        for attr in drop_attrs:
            cube.attributes.pop(attr)
        for attr in drop_time_attrs:
            cube.coord('time').attributes.pop(attr)
        cube.coord('time').points = cube.coord(
            'time').core_points().astype('float64')

    iris.util.unify_time_units(selected)
    cube = selected.concatenate_cube()
    return cube


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    axis2def = {'T': 'time', 'X': 'longitude', 'Y': 'latitude'}
    for axis in 'T', 'X', 'Y':
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)
            if axis == 'T':
                coord.convert_units('days since 1850-1-1 00:00:00.0')
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name
            coord.points = coord.core_points().astype('float64')
            if len(coord.points) > 1:
                coord.guess_bounds()
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
    # cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix units
    cube.units = definition.units

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Roll longitude
    cube.coord('longitude').points = cube.coord('longitude').points + 180.
    nlon = len(cube.coord('longitude').points)
    cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=-1)

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)

    cube.coord('latitude').attributes = None
    cube.coord('longitude').attributes = None

    cube = _fix_time_monthly(cube)

    logger.debug("Saving cube\n%s", cube)
    utils.save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes
    )
    logger.info("Finished CMORizing %s", ', '.join(in_files))


def cmorization(in_dir, out_dir, cfg, _):
    """Run CMORizer for MERRA2."""
    cfg.pop('cmor_table')

    for year in range(1980, 2019):
        for short_name, var in cfg['variables'].items():
            if 'short_name' not in var:
                var['short_name'] = short_name
            # Now get list of files
            filepattern = os.path.join(in_dir, var['file'].format(year=year))
            in_files = glob.glob(filepattern)
            if not in_files:
                logger.warning('Year %s data not found', year)
                continue
            _extract_variable(in_files, var, cfg, out_dir)
