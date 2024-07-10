"""ESMValTool CMORizer for ESACCI-SST data.

Tier
   Tier 3: need to register at CEDA

Source
   https://catalogue.ceda.ac.uk/uuid/4a9654136a7148e39b7feb56f8bb02d2

Last access
   20240628

Download and processing instructions
   A donwnloader is provided by ESMValTool. First you need
   to register.
   Go to https://services.ceda.ac.uk/cedasite/register/info/
   and create an account at CEDA if needed.

Modification history
   20240618-bock_lisa: update for v3.0
   20201204-roberts_charles: written.
   20201214-predoi_valeriu: approved.
   20201214-lauer_axel: approved.
"""

import glob
import logging
import os

import iris
from esmvalcore.cmor.fixes import get_time_bounds
from esmvalcore.preprocessor import regrid
from esmvaltool.cmorizers.data import utilities as utils
from esmvalcore.preprocessor import concatenate

from ...utilities import (
    convert_timeunits,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def extract_variable(raw_info, year):
    """Extract to all vars."""
    rawvar = raw_info['name']
    constraint = iris.NameConstraint(var_name=rawvar)
    if rawvar == 'analysed_sst_uncertainty':
        tmp_cube = iris.load_cube(raw_info['file'],
                                  iris.NameConstraint(var_name='analysed_sst'))
        ancillary_var = tmp_cube.ancillary_variable('sea_water_temperature'
                                                    ' standard_error')
        cube = tmp_cube.copy(ancillary_var.core_data())
    else:
        try:
            cube = iris.load_cube(raw_info['file'], constraint)
        except iris.exceptions.ConstraintMismatchError as constraint_error:
            raise ValueError(f"No data available for variable {rawvar}"
                             f" and year {year}") from constraint_error

    # Remove ancillary data
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)
    ## regridding from 0.05x0.05 to 0.5x0.5
    #cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')
    return cube


def get_monthly_cube(var, vals, raw_info, var_info, attrs,
                     inpfile_pattern, year, month):
    data_cubes = []
    month_inpfile_pattern = inpfile_pattern.format(
                            year=str(year)+"{:02}".format(month))
    logger.info("Pattern: %s", month_inpfile_pattern)
    inpfiles = sorted(glob.glob(month_inpfile_pattern))
    if inpfiles == []:
        logger.error("Could not find any files with this"
                     " pattern %s", month_inpfile_pattern)
        raise ValueError
    logger.info("Found input files: %s", inpfiles)

    for inpfile in inpfiles:
        raw_info['file'] = inpfile
        logger.info("CMORizing var %s from file type %s", var,
                    raw_info['file'])
        data_cubes.append(
            extract_variable(raw_info, year))

    cube = concatenate(data_cubes)

    # regridding from 0.05x0.05 to 0.5x0.5
    cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')
    # Fix dtype
    utils.fix_dtype(cube)
    # Fix cube
    fix_var_metadata(cube, var_info)
    convert_timeunits(cube, year)
    fix_coords(cube, overwrite_time_bounds=False)
    set_global_atts(cube, attrs)
    # Fix monthly time bounds
    time = cube.coord('time')
    time.bounds = get_time_bounds(time, vals['frequency'])

    return cube
    

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        if not start_date:
            start_date = vals['start_year']
        if not end_date:
            end_date = vals['end_year']
        var_info = cmor_table.get_variable(vals['mip'][0], var)
        glob_attrs['mip'] = vals['mip'][0]
        raw_info = {'name': vals['raw']}
        inpfile_pattern = os.path.join(in_dir, '{year}*' + vals['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile_pattern)
        mon_cubes = []
        for year in range(start_date.year, end_date.year + 1):
            logger.info("Processing year %s", year)
            for month in range(1, 13):
                monthly_cube = get_monthly_cube(var, vals, raw_info, var_info,
                                                glob_attrs, inpfile_pattern,
                                                year, month)
                # Save daily data
                save_variable(monthly_cube,
                              var,
                              out_dir,
                              glob_attrs,
                              unlimited_dimensions=['time'])
                # Calculate monthly mean
                if 'Stderr' not in var:
                    logger.info("Calculating monthly mean")
                    iris.coord_categorisation.add_month_number(monthly_cube,
                                                               'time')
                    iris.coord_categorisation.add_year(monthly_cube, 'time')
                    monthly_cube = monthly_cube.aggregated_by(['month_number',
                                    'year'], iris.analysis.MEAN)
                    monthly_cube.remove_coord('month_number')
                    monthly_cube.remove_coord('year')
                    mon_cubes.append(monthly_cube)
            # Save monthly data
            if 'Stderr' not in var:
                yearly_cube = concatenate(mon_cubes)
                glob_attrs['mip'] = vals['mip'][1]
                save_variable(yearly_cube,
                              var,
                              out_dir,
                              glob_attrs,
                              unlimited_dimensions=['time'])
                mon_cubes.clear()
