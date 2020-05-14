"""ESMValTool CMORizer for CERES-EBAF data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAF4Selection.jsp

Last access
    20191126

Download and processing instructions
    Select: "TOA Fluxes" ("Shortwave Flux" and "Longwave Flux", "All Sky"
    and "Clear Sky"), "Monthly Mean", "Regional 1x1 global grid".

"""

import logging
import os
import warnings

import iris

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def filter_warnings():
    """Filter certain :mod:`iris` warnings."""
    for msg in ('min', 'max'):
        warnings.filterwarnings(
            'ignore',
            message=f"WARNING: valid_{msg} not used",
            category=UserWarning,
            module='iris',
        )


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    with warnings.catch_warnings():
        filter_warnings()
        cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    with warnings.catch_warnings():
        filter_warnings()
        utils.save_variable(cube,
                            short_name,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, raw_filepath, out_dir)
