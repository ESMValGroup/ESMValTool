"""ESMValTool CMORizer for AGCD data.

Tier
    Tier 3: restricted dataset.

Source
    https://dx.doi.org/10.25914/hjqj-0x55

Last access
    20231121

Download and processing instructions
    Data from NCI project requiring an NCI account and access to GADI
    Processing is done on GADI

"""
import logging
import os
import re

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp).
        Search sub folders of raw data directory"""
    regex = re.compile(basename)
    return_files = []
    for root, dir, files in os.walk(in_dir, followlinks=True):

        for filename in files:
            if regex.match(filename):
                return_files.append(os.path.join(root, filename))

    return return_files


def fix_data_var(cube, var):

    if var == 'pr':
        cube = cube / (30 * 86400)
        cube.units = 'kg m-2 s-1'
    elif var.startswith('tas'):
        cube = cube + 273.15
        cube.units = 'K'
        utils.add_height2m(cube)
    else:
        logger.info(f"Variable {var} not converted")

    return cube


def _extract_variable(cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    logger.info("Var is %s", var)
    cubes = iris.load(filepath)
    for cube in cubes:

        cube = fix_data_var(cube, var)

        utils.fix_var_metadata(cube, cmor_info)

        utils.fix_coords(cube)
        utils.set_global_atts(cube, attrs)

        logger.info("Saving file")
        utils.save_variable(cube,
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

    # Run the cmorization #multiple variables
    for (var, var_info) in cfg['variables'].items():

        glob_attrs['mip'] = var_info['mip']
        logger.info(f"CMORizing variable '{var}', {var_info['mip']}")

        raw_filename = cfg['filename'].format(version=ver,
                                              variable=var_info['raw'],
                                              raw_calc=var_info['raw_calc'],
                                              freq=var_info['freq'])
        filepaths = _get_filepaths(in_dir, raw_filename)

        if len(filepaths) == 0:
            logger.info(f"no files for {var}. pattern:{raw_filename}")
            logger.info(f"directory:{in_dir}")
        for inputfile in filepaths:
            logger.info("Found input file '%s'", inputfile)

            cmor_info = cmor_table.get_variable(var_info['mip'], var)
            _extract_variable(cmor_info, glob_attrs, inputfile, out_dir)
