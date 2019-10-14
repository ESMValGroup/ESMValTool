
"""ESMValTool CMORizer for Landschuetzer2016 data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/data/0-data/

Last access
   20190308

Download and processing instructions
   Download the file spco2_1982-2015_MPI_SOM-FFN_v2016.nc

Modification history
   20190227-lovato_tomas: written.

"""

import logging
import os
from warnings import catch_warnings, filterwarnings

import iris
from dask import array as da

from .utilities import (constant_metadata, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube) as metadata:
        if var == 'fgco2':
            # Assume standard year 365_day
            cube *= -12.01 / 1000. / (86400. * 365.)
            metadata.attributes['positive'] = 'down'
        elif var == 'dpco2':
            cube *= -1.0 * 101325. / 1.e06
        elif var == 'spco2':
            cube *= 101325. / 1.e06
    return cube


# pylint: disable=unused-argument
def _fix_fillvalue(cube, field, filename):
    """Create masked array from missing_value."""
    if hasattr(field.cf_data, 'missing_value'):
        # fix for bad missing value definition
        cube.data = da.ma.masked_equal(cube.core_data(),
                                       field.cf_data.missing_value)


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(raw_info['file'], callback=_fix_fillvalue)
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            fix_coords(cube)
            _fix_data(cube, var)
            set_global_atts(cube, attrs)
            save_variable(
                cube,
                var,
                out_dir,
                attrs,
                local_keys=['positive'],
                unlimited_dimensions=['time'],
            )


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        inpfile = os.path.join(in_dir, vals['file'])
        logger.info("CMORizing var %s from file %s", var, inpfile)
        var_info = cmor_table.get_variable(vals['mip'], var)
        raw_info = {'name': vals['raw'], 'file': inpfile}
        glob_attrs['mip'] = vals['mip']
        with catch_warnings():
            filterwarnings(
                action='ignore',
                message=('WARNING: missing_value not used since it\n'
                         'cannot be safely cast to variable data type'),
                category=UserWarning,
                module='iris',
            )
            extract_variable(var_info, raw_info, out_dir, glob_attrs)
