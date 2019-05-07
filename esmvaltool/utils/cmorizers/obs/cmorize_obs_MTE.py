"""ESMValTool CMORizer for MTE data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://www.bgc-jena.mpg.de/geodb/BGI/Home

Last access
    20190507

Download and processing instructions
    Download the following files:
        EnsembleGPP_GL.nc

"""

import logging
import os

import iris

import esmvaltool.utils.cmorizers.obs.utilities as utils

logger = logging.getLogger(__name__)

CFG = utils.read_cmor_config('MTE.yml')


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    cube = utils.flip_dim_coord(cube, 'latitude')
    utils.save_variable(
        cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    glob_attrs = CFG['attributes']
    filepath = os.path.join(in_dir, glob_attrs['original_filename'])
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # Run the cmorization
    for (var, var_info) in CFG['variables'].items():
        logger.info("CMORizing var %s from file %s", var, filepath)
        glob_attrs['mip'] = var_info['mip']
        cmor_table = (CFG['custom_cmor_table'] if
                      var_info.get('custom_cmor_table') else CFG['cmor_table'])
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepath, out_dir)
