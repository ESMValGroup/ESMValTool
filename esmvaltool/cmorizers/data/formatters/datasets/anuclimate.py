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

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp).
        Search sub folders of raw data directory"""
    regex = re.compile(basename)
    return_files = {} #[]
    for root, dir, files in os.walk(in_dir, followlinks=True):

        for filename in files:
            if regex.match(filename):
                # return_files.append(os.path.join(root, filename)) ## dir? group by year - dict
                year = root.split('/')[-1] # if list not exists
                if year in return_files.keys():
                    return_files[year].append(os.path.join(root, filename))
                else:
                    return_files[year] = [os.path.join(root, filename)]
                # logger.info("year: %s, filename: %s", year, filename)

    return return_files


def fix_data_var(cube, var, month):
    """Convert units in cube for the variable."""
    monthdays = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30,
                 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
    if var == 'pr': #month separate files

        cube = cube / (monthdays[month] * 86400)  # days in month
        cube.units = 'kg m-2 s-1'

    elif var in ['tas', 'tasmin', 'tasmax']:  # other variables in v1
        cube = cube + 273.15
        cube.units = 'K'
        utils.add_height2m(cube)

    else:
        logger.info("Variable %s not converted", var)

    return cube


def _extract_variable(cmor_info, attrs, filepaths, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    logger.info("Var is %s", var)
    cbls = iris.cube.CubeList()
    for filepath in filepaths:
        cubes = iris.load(filepath) #load by year? list
        month = filepath[-5:-3]
        # logger.info("month number: %s", month)
    # for cube in cubes: ## ls of filepaths? for year

        cube = fix_data_var(cubes[0], var, int(month))

        utils.fix_var_metadata(cube, cmor_info)

        utils.set_global_atts(cube, attrs)
        cbls.append(cube)

    iris.util.equalise_attributes(cbls)
    cubesave = cbls.concatenate_cube()
    utils.fix_coords(cube)
    logger.info("Saving file") #merge cubes
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

    # Run the cmorization #multiple variables
    for (var, var_info) in cfg['variables'].items():

        glob_attrs['mip'] = var_info['mip']
        logger.info(f"CMORizing variable '{var}', {var_info['mip']}")

        raw_filename = cfg['filename'].format(version=ver,
                                              raw=var_info['raw'],
                                              freq=var_info['freq'])
        filepaths = _get_filepaths(in_dir, raw_filename) #dict

        if len(filepaths) == 0:
            logger.info(f"no files for {var}. pattern:{raw_filename}")
            logger.info(f"directory:{in_dir}")
        # for inputfile in filepaths:
        for year,inputfiles in filepaths.items():
            logger.info("Found input file '%s', count %s", year, len(inputfiles))

            cmor_info = cmor_table.get_variable(var_info['mip'], var)
            _extract_variable(cmor_info, glob_attrs, inputfiles, out_dir)
