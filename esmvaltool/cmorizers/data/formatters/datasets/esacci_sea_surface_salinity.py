"""ESMValTool CMORizer for ESACCI-SOS data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_surface_salinity/data

Last access
   20200921

Download and processing instructions
   Download the data from:
       v01.8/30days/
       v02.31/30days/
   Put all files under a single directory (no subdirectories with years).
"""

import logging
import os

import iris
from iris.util import equalise_attributes, unify_time_units

from esmvaltool.cmorizers.data.utilities import (
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize ESACCI-SOS dataset."""
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    for version in glob_attrs['versions']:
        logger.info('Cmorizing version %s', version)
        file_expr = os.path.join(
            in_dir,
            "ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_Monthly_CENTRED_15Day_"
            f"25km-20????15-{version}.nc")

        for var, vals in cfg['variables'].items():
            var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
            logger.info('Cmorizing var %s', var)
            cubes = iris.load_raw(file_expr, vals['raw'])
            equalise_attributes(cubes)
            unify_time_units(cubes)
            cube = cubes.concatenate_cube()
            cube.units = '0.001'
            logger.info(cube)
            glob_attrs['mip'] = vals['mip']
            glob_attrs['version'] = version
            fix_var_metadata(cube, var_info)
            set_global_atts(cube, glob_attrs)
            save_variable(cube, var, out_dir, glob_attrs)
