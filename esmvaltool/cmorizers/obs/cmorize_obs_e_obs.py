"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles

Last access
    20200225

Download and processing instructions
    Download the ensemble mean files for:
        TG TN TX RR PP
"""

import logging
import os
import iris
import iris.coord_categorisation

from esmvalcore.preprocessor import monthly_statistics
from . import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, res, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.units = var.get('raw_units', short_name)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    utils.fix_var_metadata(cube, cmor_info)
    attrs = cfg['attributes'].copy()
    attrs['version'] = 'v'+attrs['version']+'_'+str(res)
    attrs.pop('resolution')
    attrs['mip'] = var['mip']
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    #####
    # also derive monthly data
    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Building monthly means")

            # Calc monthly
            cube = monthly_statistics(cube)
            cube.remove_coord('month_number')
            cube.remove_coord('year')

            # Fix metadata
            attrs['mip'] = 'Amon'

            # Fix coordinates
            utils.fix_coords(cube)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    ver = cfg['attributes']['version']
    for res in cfg['attributes']['resolution'].values():
        for (short_name, var) in cfg['variables'].items():
            logger.info("CMORizing variable '%s' on %s°x%s°",
                        short_name, res, res)
            raw_var = var.get('raw', short_name)
            filepath = raw_filepath.format(raw_name=raw_var, resolution=res,
                                           version=ver)
            _extract_variable(short_name, var, res, cfg, filepath, out_dir)
