"""ESMValTool CMORizer for MERRA2 data.

Tier
    Tier 3: restricted datasets (i.e., dataset which requires a registration
 to be retrieved or provided upon request to the respective contact or PI).

Source
    https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/
    https://daac.gsfc.nasa.gov/datasets/M2IUNPASM_5.12.4/summary?keywords=MERRA2

Last access
    20220913

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

from esmvaltool.cmorizers.data import utilities as utils

# iris spits out a large amount of warnings
# logging.disable('WARNING')  # careful: this deactivates writing
# to log files in the pytest test environment for other cmorizer tests
logger = logging.getLogger(__name__)


def _fix_time_monthly(cube):
    """Fix time by setting it to 15th of month."""
    # Read dataset time unit and calendar from file
    dataset_time_unit = str(cube.coord('time').units)
    dataset_time_calender = cube.coord('time').units.calendar
    # Convert datetime
    time_as_datetime = cf_units.num2date(
        cube.coord('time').core_points(), dataset_time_unit,
        dataset_time_calender)
    newtime = []
    for timepoint in time_as_datetime:
        midpoint = datetime(timepoint.year, timepoint.month, 15)
        newtime.append(midpoint)

    newtime = cf_units.date2num(newtime, dataset_time_unit,
                                dataset_time_calender)
    # Put them on the file
    cube.coord('time').points = newtime
    cube.coord('time').bounds = None
    return cube


def _var_pairs(cube_list, var_parts, oper):
    """Return a selection composed of two variables."""
    selected_1 = [c for c in cube_list if c.var_name == var_parts[0]]
    selected_2 = [c for c in cube_list if c.var_name == var_parts[1]]
    if not selected_1:
        logger.error("Raw variable %s could not be found "
                     "in str(cube_list) - operation can not be performed.",
                     var_parts[0])
        raise ValueError
    if not selected_2:
        logger.error("Raw variable %s could not be found "
                     "in str(cube_list) - operation can not be performed.",
                     var_parts[1])
        raise ValueError
    if oper == "-":
        selected = [
            cube_1 - cube_2 for cube_1, cube_2 in zip(selected_1, selected_2)
        ]
        selected = iris.cube.CubeList(selected)
    else:
        raise NotImplementedError(f"Pairwise variables operation {oper} "
                                  "not implemented yet, you can do it "
                                  "yourself in the MERRA2 cmorizer.")

    return selected


def _load_cube(in_files, var):
    cube_list = iris.load_raw(in_files)
    pairwise_ops = ["+", "-", ":"]
    var_parts = []
    for oper in pairwise_ops:
        split_var = var['raw'].split(oper)
        if len(split_var) == 2:
            var_parts = [split_var[0],  split_var[1]]
            break
        if len(split_var) > 2:
            logger.error("Splitting raw variable %s by "
                         "operation %s results in more than two"
                         " raw variables, this is not yet implemented.",
                         var['raw'], oper)
            raise NotImplementedError
    if not var_parts:
        selected = [c for c in cube_list if c.var_name == var['raw']]
        selected = iris.cube.CubeList(selected)
    else:
        selected = _var_pairs(cube_list, var_parts, oper)

    drop_attrs = [
        'History', 'Filename', 'Comment', 'RangeBeginningDate',
        'RangeEndingDate', 'GranuleID', 'ProductionDateTime', 'Source'
    ]
    drop_time_attrs = [
        'begin_date', 'begin_time', 'time_increment', 'valid_range', 'vmax',
        'vmin'
    ]
    for cube in selected:
        for attr in drop_attrs:
            cube.attributes.pop(attr)
        for attr in drop_time_attrs:
            cube.coord('time').attributes.pop(attr)
        cube.coord('time').points = cube.coord('time').core_points().astype(
            'float64')

    iris.util.unify_time_units(selected)
    cube = selected.concatenate_cube()
    return cube


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    if cube.ndim == 3:
        axis2def = {'T': 'time', 'X': 'longitude', 'Y': 'latitude'}
        axes = ['T', 'X', 'Y']
    elif cube.ndim == 4:
        axis2def = {'T': 'time', 'X': 'longitude',
                    'Y': 'latitude', 'Z': 'plev19'}
        axes = ['T', 'X', 'Y', 'Z']
    for axis in axes:
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
        else:
            # special case for UV
            # variable "uv" (raw: "V") comes with "alevel" instead
            # of "plev19" in the table; "alevel" has empty fields for
            # standard_name, out_name etc. so we need to set them; it's safe
            # to do so since the cmor checker/fixer will convert that during
            # preprocessing at cmor fix stage
            if cube.var_name == "uv" and axis == "Z":
                coord = cube.coord(axis=axis)
                coord_def = definition.coordinates.get('alevel')
                coord.standard_name = "air_pressure"
                coord.var_name = "plev"
                coord.long_name = "pressure"
                coord.points = coord.core_points().astype('float64')
                if len(coord.points) > 1:
                    coord.guess_bounds()

    return cube


def _extract_variable(in_files, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input files '%s'",
                var['short_name'], ', '.join(in_files))
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    attributes['raw'] = var['raw']
    pairwise_ops = ["+", "-", ":"]
    for oper in pairwise_ops:
        if oper in var['raw']:
            components = var['raw'].split(oper)
            if len(components) == 2:
                attributes['component_raw_1'] = components[0]
                attributes['component_raw_2'] = components[1]
                attributes['component_operation'] = oper
                break
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = _load_cube(in_files, var)

    # keep the following raw cube attributes
    attrs_to_keep = [
        "institution", "Institution",
        "institute_id", "VersionID",
        "experiment_id",
        "source", "Source",  # overrides empty string default
        "model_id", "ModelID",
        "contact", "Contact",
        "references",
        "tracking_id",
        "mip_specs",  # described by "mip" already
        "source_id", "SourceID",
        "product", "Product",
        "frequency", "Frequency",
        "creation_date",
        "project_id", "ProjectID",
        "table_id", "TableID",
        "title", "Title",
        "modeling_realm",
        "doi",
        "VersionID",  # described by "version" already
    ]

    attrs_to_keep_exist = [
        att for att in cube.attributes if att in attrs_to_keep
    ]
    for att in attrs_to_keep_exist:
        attributes[att] = cube.attributes[att]

    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    # cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix units (if needed)
    # input variable reports m-3 m-3 instead of m3 m-3
    if cube.var_name == "sm":
        cube.units = definition.units
    # Convert units to CMOR units
    cube.convert_units(definition.units)

    # Add height2m or height10m if needed
    if 'height2m' in definition.dimensions:
        utils.add_height2m(cube)
    elif 'height10m' in definition.dimensions:
        utils.add_height10m(cube)

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
    utils.save_variable(cube, cube.var_name, out_dir, attributes)
    logger.info("Finished CMORizing %s", ', '.join(in_files))


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for MERRA2."""
    cfg.pop('cmor_table')
    if start_date is None:
        start_date = 1980
    else:
        start_date = start_date.year
    if end_date is None:
        end_date = 2022
    else:
        end_date = end_date.year
    for year in range(start_date, end_date + 1):
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
