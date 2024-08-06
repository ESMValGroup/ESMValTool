"""ESMValTool CMORizer for ESACCI-WATERVAPOUR data.

Tier
   Tier 3: CDR2 requires registration at EUMETSAT CM SAF.

Source
   https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=COMBI_V001

Last access
   20240221

Download and processing instructions
   CDR2 requires registration at EUMETSAT CM SAF, the information on how to
        download the order will be emailed once the order is ready.
   All files need to be in one directory, not in yearly subdirectories.

Modification history
   20240221-malinina_elizaveta: Adjust for daily cmorization and updated
                                filenames, remove CDR1 due to irrelevance.
   20210607-weigel_katja: Fix for monthly time bounds.
   20210408-weigel_katja: written.
"""

import glob
import logging
import os

import iris
from esmvalcore.cmor.fixes import get_time_bounds
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
    cube = fix_coords(cube, overwrite_time_bounds=False)
    set_global_atts(cube, attrs)
    # Remove dysfunctional ancillary data without sandard name
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var_name, vals in cfg['variables'].items():
        var = vals['short_name']
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw']}
        inpfile_pattern = os.path.join(in_dir, vals['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile_pattern)
        for year in range(vals['start_year'], vals['end_year'] + 1):
            data_cubes = []
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info['file'] = inpfile
                logger.info("CMORizing var %s from file type %s", var,
                            raw_info['file'])
                data_cubes.append(
                    extract_variable(var_info, raw_info, glob_attrs, year))
            yearly_cube = concatenate(data_cubes)
            # Fix monthly time bounds
            time = yearly_cube.coord('time')
            time.bounds = get_time_bounds(time, vals['frequency'])
            save_variable(yearly_cube,
                          var,
                          out_dir,
                          glob_attrs,
                          unlimited_dimensions=['time'])
