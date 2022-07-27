"""ESMValTool CMORizer for Kadow2020 data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://users.met.fu-berlin.de/~ChristopherKadow/

Last access
    20220329

Download and processing instructions
    Download the following file:
        [SOURCE]/HadCRUT.5.0.1.0.anomalies.Kadow_et_al_2020_20crAI-infilled
                        .ensemble_mean_185001-202012.nc
"""

import copy
import logging
import os

import iris
from cf_units import Unit
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, version, filename, cfg, in_dir,
                      out_dir):
    """Extract variable."""
    # load data
    filepath = os.path.join(in_dir, filename)
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    # Fix coordinates
    utils.fix_dim_coordnames(cube)
    # fix flipped latitude
    utils.flip_dim_coord(cube, 'latitude')
    utils.fix_dim_coordnames(cube)
    cube_coord = cube.coord('latitude')
    utils.fix_bounds(cube, cube_coord)
    cube_coord = cube.coord('longitude')
    utils.fix_bounds(cube, cube_coord)

    # add heigt2m coordinate
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata and  update version information
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    attrs['version'] = version
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for (version, filename) in cfg['filenames'].items():
            logger.info("CMORizing variable '%s' '%s'", short_name, version)
            _extract_variable(short_name, var, version, filename, cfg, in_dir,
                              out_dir)
