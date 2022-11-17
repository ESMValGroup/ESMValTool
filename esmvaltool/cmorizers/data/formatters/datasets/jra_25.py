"""ESMValTool CMORizer for "NOAAGlobalTemp data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/products/land-based-station/noaa-global-temp

Last access
    20220628

Download and processing instructions
    Download the following files:
        [SOURCE]/v5/access/gridded/
          NOAAGlobalTemp_v5.0.0_gridded_s188001_e202205_c20220608T133245.nc
"""

import copy
import logging
import os

import iris
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, filename, cfg, in_dir,
                      out_dir):
    """Extract variable."""
    # load data
    filepath = os.path.join(in_dir, filename)
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']

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
        filename = var['file']
        logger.info("CMORizing variable '%s' from file '%s'", short_name, filename)
        _extract_variable(short_name, var, filename, cfg, in_dir, out_dir)
