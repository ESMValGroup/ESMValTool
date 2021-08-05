"""ESMValTool CMORizer for RAOBCOREv1.5.1 data.

Tier
    Tier 3: restricted dataset, available from Leo Haimberger

Source
    Leo Haimberger

Last access
    20200924

Download and processing instructions
    Download the ensemble member files for rio and rit
"""

import logging
import os
import iris
import numpy as np
from cf_units import Unit
from . import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, ver, var, cfg, filepath, out_dir):

    # extract variable
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))   

    # fix plev
    cube.coord("pressure level").standard_name = "air_pressure"
    cube.coord("air_pressure").long_name = "air_pressure"
    cube.coord("air_pressure").var_name = "plev"
    cube.coord("air_pressure").attributes = None  
 
    # fix lat
    cube = iris.util.reverse(cube, "latitude")
    cube.coord("latitude").attributes = None

    # fix lon
    cube.coord("longitude").attributes = None

    # fix metadata
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    attrs = cfg['attributes']
    attrs['version'] = ver
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)
    
    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])
 

def cmorization(in_dir, out_dir, cfg, config_user):

    # Cmorization func call
    raw_filepath = os.path.join(in_dir, cfg['filename'])
    
    # Run the cmorization
    for ver in cfg['attributes']['version'].values():
       for (short_name, var) in cfg['variables'].items():
            logger.info("CMORizing variable %s in '%s' ", short_name, ver)
            filepath = raw_filepath.format(version=ver)
            _extract_variable(short_name, ver, var, cfg, filepath, out_dir)
