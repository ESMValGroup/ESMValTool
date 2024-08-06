"""ESMValTool CMORizer for AGCD data.

Tier
    Tier 2: other freely available dataset.

Source
    https://dx.doi.org/10.25914/rses-zh67

Last access
    20231121

Download and processing instructions
    Data from NCI (National Computing Infrastructure Australia)
    https://nci.org.au/,
    requiring an NCI account and access to Gadi(Supercomputer in Australia)
      and the dataset project found in
      catalogue record https://dx.doi.org/10.25914/rses-zh67.
    Access can be requested through NCI.
    NCI is an ESGF node: (https://esgf.nci.org.au/projects/esgf-nci/)
    Processing is done on Gadi.

"""
import logging
import os
import re

import iris

from esmvalcore.cmor._fixes.shared import get_time_bounds
from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    regex = re.compile(basename)
    return_files = []
    for root, _, files in os.walk(in_dir, followlinks=True):

        for filename in files:
            if regex.match(filename):
                return_files.append(os.path.join(root, filename))

    return return_files


def fix_data_var(cube, var):
    """Convert units in cube for the variable."""
    monthdays = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30,
                 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
    if var == 'pr':
        newcubels = []
        for i, m_cube in enumerate(cube.slices(['latitude', 'longitude'])):
            m_cube = m_cube / (monthdays[i + 1] * 86400)  # days in month
            newcubels.append(m_cube)

        cube = iris.cube.CubeList(newcubels).merge()[0]
        cube.units = 'kg m-2 s-1'

    elif var in ['tas', 'tasmin', 'tasmax']:  # other variables in v1
        cube = cube + 273.15
        cube.units = 'K'
        utils.add_height2m(cube)

    else:
        logger.info("Variable %s not converted", var)

    return cube


def _extract_variable(cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    logger.info("Var is %s", var)
    cubes = iris.load(filepath)
    for cube in cubes:

        cube = fix_data_var(cube, var)

        utils.fix_var_metadata(cube, cmor_info)

        cube = utils.fix_coords(cube)
        bounds = get_time_bounds(cube.coords('time')[0], 'mon')
        cube.coords('time')[0].bounds = bounds
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

    # Run the cmorization #multiple variables
    for (var, var_info) in cfg['variables'].items():

        glob_attrs['mip'] = var_info['mip']
        logger.info("CMORizing variable '%s', %s", var, var_info['mip'])

        raw_filename = cfg['filename'].format(version=ver,
                                              variable=var_info['raw'],
                                              raw_calc=var_info['raw_calc'],
                                              freq=var_info['freq'])
        filepaths = _get_filepaths(in_dir, raw_filename)

        if not filepaths:
            logger.info("no files for %s. pattern:%s", var, raw_filename)
            logger.info("directory:%s", in_dir)
        for inputfile in filepaths:
            logger.info("Found input file '%s'", inputfile)

            cmor_info = cmor_table.get_variable(var_info['mip'], var)
            _extract_variable(cmor_info, glob_attrs, inputfile, out_dir)
