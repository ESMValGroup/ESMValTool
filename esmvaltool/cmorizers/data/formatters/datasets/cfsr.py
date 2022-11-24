"""ESMValTool CMORizer for CFSR data.
Tier
    Tier 2: other freely-available dataset.
Source
    ESGF:
    https://esgf.nccs.nasa.gov/thredds/fileServer/CREATE-IP/
        reanalysis/NOAA-NCEP/CFSR/CFSR/mon/atmos/
Last access
    20221122
Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/cfsr.py
"""

import copy
import logging
import re
from copy import deepcopy
import os
from pathlib import Path
from warnings import catch_warnings, filterwarnings
from cf_units import Unit

import iris
from iris.coords import CellMethod
from esmvalcore.cmor.table import CMOR_TABLES
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, raw_filepath, out_dir):
    """Extract variable."""
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], short_name)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    if cmor_info.positive != '':
        attributes['positive'] = cmor_info.positive

    # load data
    raw_var = var.get('raw', short_name)
    with catch_warnings():
        filterwarnings('ignore',
                       message='Ignoring netCDF variable .* invalid units .*',
                       category=UserWarning,
                       module='iris')
        cube = iris.load_cube(str(raw_filepath),
                              NameConstraint(var_name=raw_var))

    utils.set_global_atts(cube, attributes)

    utils.fix_var_metadata(cube, cmor_info)

    # adjusting the attribute "cell methods"
    cube.cell_methods = ()
    cube.add_cell_method(CellMethod("mean", coords=["time"]))

    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    # Save variable
    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attributes,
        unlimited_dimensions=['time'],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        filename = var['file']
        logger.info("CMORizing variable '%s' from file '%s'", short_name,
                    filename)
        short_name = var['short_name']
        print(short_name)
        raw_filenames = Path(in_dir).rglob('*.nc')
        filenames = []
        for raw_filename in raw_filenames:
            if re.search(var['file'], str(raw_filename)) is not None:
                filenames.append(raw_filename)

        for filename in sorted(filenames):

            _extract_variable(short_name, var, cfg, filename, out_dir)