"""
ESMValTool CMORizer for JRA-55 data.

Tier
    Tier 2: other freely-available dataset.

Source
    Research Data Archive (RDA):
    https://rda.ucar.edu/datasets/ds628.1/

Last access
    20230322

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/jra_55.py
"""

import copy
import glob
import logging
import os

import iris
import numpy as np

from datetime import datetime
from iris_grib.message import GribMessage
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _load_jra55_grib(filenames, var):
    """Load data from GRIB file and return list of cubes."""
    cubelist = []
    for infile in filenames:
        tmp_cubes = iris.load(infile)
        if len(tmp_cubes) > 1:
            start_element = var.get('start_element', 0)
            i = 0
            # create list of files (needed in case infile contains wildcards)
            listing = glob.glob(infile)
            for fname in listing:
                for message in GribMessage.messages_from_filename(fname):
                    day = message.sections[1]['day']
                    month = message.sections[1]['month']
                    year = ((message.sections[1]['centuryOfReferenceTimeOfData']
                         - 1) * 100 + message.sections[1]['yearOfCentury'])

                    point = datetime(year=year, month=month, day=day)
                    time_units = Unit('days since 1950-01-01 00:00:00',
                        calendar='standard')
                    time_coord = iris.coords.DimCoord(
                        time_units.date2num(point),
                        var_name='time',
                        standard_name='time',
                        long_name='time',
                        units=time_units)

                    tmp_cubes[i].add_aux_coord(time_coord)
                    tmp_cubes[i].remove_coord('originating_centre')

                    i = i + 1

            # message.sections[1]['indicatorOfTypeOfLevel'] always gives
            # 'sfc', so a distinction between "surface" and "top of the
            # atmosphere" is not possible. Instead, we simply use every
            # second cube as "surface" and "top of the atmosphere" are
            # alternating in the GRIB file with "surface" being element
            # 'start_element' read from config file
            # (esmvaltool/cmorizers/data/cmor_config/JRA-55.yml).

            cubelist.append(tmp_cubes[start_element::2].merge_cube())
        else:
            #tmp_cubes.remove_coord('originating_centre')
            cubelist = tmp_cubes

    return(cubelist)


def _extract_variable(short_name, var, in_files, cfg, in_dir,
                      out_dir):
    """Extract variable."""
    # load data (returns a list of cubes)

    cubes = _load_jra55_grib(in_files, var)
    
    # apply operators (if any)

    if len(cubes) > 1:
        if var.get('operator', '') == 'sum':
            # Multiple variables case using sum operation
            cube = None
            for in_cube in cubes:
                if cube is None:
                    cube = in_cube
                else:
                    cube += in_cube
        elif var.get('operator', '') == 'diff':
            # two variables case using diff operation
            if len(cubes) != 2:
                errmsg = (f'operator diff selected for variable {short_name} '
                          f'expects exactly two input variables and two input '
                          f'files')
                raise ValueError(errmsg)
            cube = cubes[0]
            cube -= cubes[1]
        else:
            oper = var.get('operator')
            raise ValueError(
                f'multiple input files found for variable {short_name} '
                f'with unknown operator {oper}')
    else:
        cube = cubes[0]

    print(cube)
    print(cubes)

    # Fix metadata
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
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
    for (short_name, var) in cfg['variables'].items():
        short_name = var['short_name']

        if 'file' in var:
            filename = [os.path.join(in_dir, var['file'])]
        elif 'files' in var:
            filename = []
            for file in var['files']:
                filename.append(os.path.join(in_dir, file))
        else:
            raise ValueError(f"No input file(s) specified for variable "
                             f"{short_name}.")

        logger.info("CMORizing variable '%s' from file '%s'", short_name,
                    filename)
        _extract_variable(short_name, var, filename, cfg, in_dir, out_dir)
