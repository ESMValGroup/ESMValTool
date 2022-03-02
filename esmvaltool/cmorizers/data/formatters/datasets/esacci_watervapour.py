"""ESMValTool CMORizer for ESACCI-WATERVAPOUR data.

Tier
   Tier 3: currently still restricted because preliminary.

Source
   Marc Schröder, ftp.brockmann-consult.de

Last access
   20210329

Download and processing instructions
   FTP server: ftp.brockmann-consult.de, access currently restricted
               data/tcwv/dataset3_1/CDR-*/...
   All files need to be in one directory, not in yearly subdirectories.

Modification history
   20210607-weigel_katja: Fix for monthly time bounds.
   20210408-weigel_katja: written.
"""

import logging
import os

import iris
from esmvalcore.cmor.check import _get_time_bounds
from esmvalcore.preprocessor import concatenate

from ...utilities import (
    convert_timeunits,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def extract_variable(var_info, raw_info, attrs, year):
    """Extract to all vars."""
    rawvar = raw_info['name']
    constraint = iris.NameConstraint(var_name=rawvar)
    try:
        cube = iris.load_cube(raw_info['file'], constraint)
    except iris.exceptions.ConstraintMismatchError as constraint_error:
        raise ValueError(f"No data available for variable {rawvar}"
                         f"and year {year}") from constraint_error

    # Fix cube
    fix_var_metadata(cube, var_info)
    convert_timeunits(cube, year)
    fix_coords(cube, overwrite_time_bounds=False)
    set_global_atts(cube, attrs)
    # Remove dysfunctional ancillary data without sandard name
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    # cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}
        inpfile = os.path.join(in_dir, cfg['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile)
        # years = range(vals['start_year'], vals['end_year'] + 1)
        months = ["0" + str(mo) for mo in range(1, 10)] + ["10", "11", "12"]
        for year in range(vals['start_year'], vals['end_year'] + 1):
            monthly_cubes = []
            for month in months:
                raw_info['file'] = inpfile.format(year=year, month=month)
                logger.info("CMORizing var %s from file type %s", var,
                            raw_info['file'])
                monthly_cubes.append(
                    extract_variable(var_info, raw_info, glob_attrs, year))
            yearly_cube = concatenate(monthly_cubes)
            # Fix monthly time bounds
            time = yearly_cube.coord('time')
            time.bounds = _get_time_bounds(time, 'mon')
            save_variable(yearly_cube,
                          var,
                          out_dir,
                          glob_attrs,
                          unlimited_dimensions=['time'])
