"""
ESMValTool CMORizer for ISCCP-H data (3-hourly).

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/international-satellite-cloud-
        climate-project-isccp-h-series-data/access/isccp-basic/hgg

Last access
    20230524

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/isccp_h.py
"""

import copy
import glob
import logging
import os

import iris
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

from datetime import datetime
from dateutil import relativedelta

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, in_files, cfg, in_dir,
                      out_dir):
    """Extract variable."""
    # load data
    raw_var = var.get('raw', short_name)
    rawcubes = iris.load(in_files, NameConstraint(var_name=raw_var))

    drop_global_atts = [
        'id', 'isccp_input_files', 'time_coverage_start', 'time_coverage_end',
        'isccp_gmt', 'isccp_number_of_satellites_contributing', 'platform',
        'instrument', 'date_issued', 'date_created', 'date_modified',
        'date_metadata_modified', 'history'
    ]

    drop_var_atts = [
        'isccp_day', 'isccp_month', 'isccp_percent_empty_cells',
        'isccp_percent_full_cells', 'isccp_year'
    ]

    for cube in rawcubes:
        for att in drop_global_atts:
            if (att in cube.attributes):
                cube.attributes.pop(att)
        for att in drop_var_atts:
            if (att in cube.attributes):
                cube.attributes.pop(att)

    cube = rawcubes.concatenate_cube()

    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    try:
        cube.convert_units(cmor_info.units)
    except Exception:
        # special case: cloud water path
        if cube.units == 'cm' and cmor_info.units == 'kg m-2':
            # The ISCCP-H documentation available at
            # https://www.ncei.noaa.gov/data/international-
            # satellite-cloud-climate-project-isccp-h-series-data/doc/
            # CDRP-ATBD-0872%20Rev%200%20Cloud%20Properties%20-
            # %20ISCCP%20C-ATBD%20(01B-29)%20(DSR-1109).pdf
            # (page 78) reports that cloud water path is given in units
            # of 'g m-2'.
            # Using units of 'g m-2' from the documentation rather
            # that the units 'cm' reported by the netCDF files gives reasonable
            # values for 'cloud water path' while the units of 'cm' does not.
            # ---> convert from 'g m-2' to 'kg m-2' (ignoring netCDF units)
            cube *= 0.001
            cube.units = cmor_info.units

    # Fix metadata
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
    if start_date is None:
        start_date = datetime(1984, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 12, 31)
    loop_date = start_date

    while loop_date <= end_date:
        year = loop_date.year
        month = f'{loop_date.month:0>2}'

        for short_name, var in cfg['variables'].items():
            if 'short_name' not in var:
                var['short_name'] = short_name

            # Now get list of files
            filepattern = os.path.join(in_dir + '/' + str(year),
                                       var['file'].format(year=year,
                                                          month=month))
            print(filepattern)
            in_files = glob.glob(filepattern)
            if not in_files:
                logger.warning('Year %s data not found', year)
                continue
            _extract_variable(short_name, var, in_files, cfg, in_dir, out_dir)

        loop_date += relativedelta.relativedelta(months=1)
