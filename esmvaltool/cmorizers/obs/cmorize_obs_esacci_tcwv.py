"""ESMValTool CMORizer for ESACCI-TCWV data.

Tier
   Tier 3: currently still restricted because preliminary.

Source
   Marc Schröder, ftp.brockmann-consult.de

Last access
   20212903

Download and processing instructions
   TBD

Modification history
   20210408-weigel_katja: written.

"""

import logging
import os

import iris

from esmvalcore.preprocessor import concatenate
from .utilities import (convert_timeunits, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def extract_variable(var_info, raw_info, attrs, year):
    """Extract to all vars."""
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            convert_timeunits(cube, year)
            fix_coords(cube)
            set_global_atts(cube, attrs)
            # Remove disfunctional ancillary data without sandard name
            for ancillary_variable_, dim in cube._ancillary_variables_and_dims:
                cube.remove_ancillary_variable(ancillary_variable_)
            return cube


def cmorization(in_dir, out_dir, cfg, _):
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}
        inpfile = os.path.join(in_dir, cfg['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile)
        years = range(vals['start_year'], vals['end_year'] + 1)
        months = ["0" + str(mo) for mo in range(1, 10)] + ["10", "11", "12"]
        for year in years:
            monthly_cubes = []
            for month in months:
                raw_info['file'] = inpfile.format(year=year, month=month)
                logger.info("CMORizing var %s from file type %s", var,
                            raw_info['file'])
                cube = extract_variable(var_info, raw_info, glob_attrs, year)
                monthly_cubes.append(cube)
            yearly_cube = concatenate(monthly_cubes)
            save_variable(yearly_cube,
                          var,
                          out_dir,
                          glob_attrs,
                          unlimited_dimensions=['time'])
