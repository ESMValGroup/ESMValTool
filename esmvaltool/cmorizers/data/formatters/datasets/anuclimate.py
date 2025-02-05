"""ESMValTool CMORizer for ANU Climate data.

Tier
    Tier 3: restricted dataset.

Source
    https://dx.doi.org/10.25914/60a10aa56dd1b

Last access
    20231121

Download and processing instructions
    Data from NCI project requiring an NCI account and access to GADI
    Processing is done on GADI

"""
import logging
import os
import re
import calendar

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    regex = re.compile(basename)
    return_files = []
    # Search sub folders of raw data directory
    for root, _dir, files in os.walk(in_dir, followlinks=True):

        for filename in files:
            if regex.match(filename):

                return_files.append(os.path.join(root, filename))

    return return_files


def fix_data_var(cube, var):
    """Convert units in cube for the variable."""
    # get month, year from cube
    tcoord = cube.coord('time')
    tdate = tcoord.units.num2date(tcoord.points[0])
    no_ofdays = calendar.monthrange(tdate.year, tdate.month)[1]

    if var == 'pr':

        cube = cube / (no_ofdays * 86400)  # days in month
        cube.units = 'kg m-2 s-1'

    elif var in ['tas', 'tasmin', 'tasmax']:  # other variables in v1
        cube = cube + 273.15
        cube.units = 'K'
        utils.add_height2m(cube)

    else:
        logger.info("Variable %s not converted", var)

    return cube, tdate.year


def _extract_variable(cmor_info, attrs, filepaths, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    logger.info("Var is %s", var)
    cbls_2 = iris.cube.CubeList()
    cbls_1 = iris.cube.CubeList()
    for filepath in filepaths:
        cubes = iris.load(filepath)

        cube, year = fix_data_var(cubes[0], var)

        utils.fix_var_metadata(cube, cmor_info)

        utils.set_global_atts(cube, attrs)

        if year < 2000:  # split for cube save
            cbls_1.append(cube)
        else:
            cbls_2.append(cube)

    for cbls in [cbls_1, cbls_2]:
        iris.util.equalise_attributes(cbls)
        cubesave = cbls.concatenate_cube()
        cubesave = utils.fix_coords(cubesave)

        logger.info("Saving file")
        utils.save_variable(cubesave,
                            var,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    ver = cfg['attributes']['version']
    logger.info(cfg, cfg_user)

    # Run the cmorization, multiple variables
    for (var, var_info) in cfg['variables'].items():

        glob_attrs['mip'] = var_info['mip']

        raw_filename = cfg['filename'].format(version=ver,
                                              raw=var_info['raw'],
                                              freq=var_info['freq'])
        filepaths = _get_filepaths(in_dir, raw_filename)

        if len(filepaths) == 0:
            logger.info("no files for %s pattern: %s", var, raw_filename)
            logger.info("directory: %s", in_dir)
        else:
            logger.info("Found files, count %s", len(filepaths))

            cmor_info = cmor_table.get_variable(var_info['mip'], var)
            _extract_variable(cmor_info, glob_attrs, filepaths, out_dir)
