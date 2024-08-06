"""ESMValTool CMORizer for NOAA ERSST data, version 5.

   This is the CMORizer script for the NOAA Extended Reconstructed Sea Surface
   Temperature (ERSST) data of version 5.

Tier
    Tier 2: open dataset.

Source
    https://doi.org/10.7289/V5T72FNM

Last access
    20200520

Download and processing instructions
    The data is provided by NOAA at:
    https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/

"""

import logging
import os
import re

import iris
import cf_units

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    regex = re.compile(basename)
    return_files = []
    return_files_gr08 = []
    for file in os.listdir(in_dir):

        if regex.match(file):
            year = file.split('.')[2][:4]  # ersst.v5.$yr$nm.nc
            # return 2 lists as files differ from 2008
            if int(year) < 2008:
                return_files.append(os.path.join(in_dir, file))
            else:
                return_files_gr08.append(os.path.join(in_dir, file))

    return return_files, return_files_gr08


def _fix_time_coord(cube, _, _filename):
    """Set time points to central day of month and standardise time units."""
    t_coord = cube.coord('time')
    _unit = t_coord.units
    new_time = [d.replace(day=15) for d in _unit.num2date(t_coord.points)]
    t_coord.points = _unit.date2num(new_time).astype('float64')
    t_coord.units = cf_units.Unit(t_coord.units.origin, calendar='standard')
    t_coord.long_name = 'Time'


def _extract_variable(raw_var, cmor_info, attrs, filepaths, out_dir):
    """Extract variable and concatenate months."""
    var = cmor_info.short_name

    cubels = iris.load(filepaths, raw_var, _fix_time_coord)
    iris.util.equalise_attributes(cubels)
    iris.util.unify_time_units(cubels)
    cube = cubels.concatenate_cube()
    cube = iris.util.squeeze(cube)

    utils.fix_var_metadata(cube, cmor_info)
    cube = utils.fix_coords(cube)

    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    filepaths = _get_filepaths(in_dir, cfg['filename'])

    if len(filepaths[0]) > 0 or len(filepaths[1]) > 0:
        totalfiles = len(filepaths[0]) + len(filepaths[1])
        logger.info("%d files before 2008", len(filepaths[0]))
        logger.info("Found %d input files in '%s'", totalfiles, in_dir)
    else:
        logger.info("No files found, basename: %s", cfg['filename'])

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs,
                          filepaths[0], out_dir)
        _extract_variable(raw_var, cmor_info, glob_attrs,
                          filepaths[1], out_dir)
