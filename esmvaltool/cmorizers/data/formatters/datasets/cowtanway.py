"""ESMValTool CMORizer for CowtanWay.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html

Last access
    20200226

Download and processing instructions
    Download the following files:
        'had4_krig_v1_0_0.nc.gz'
        'had4_uah_v1_0_0.nc.gz'
        'had4_short_krig_v2_0_0.nc.gz'
        'had4_short_uah_v2_0_0.nc.gz'
        'ghcn_short_krig_v2_0_0.nc.gz'
        'ghcn_short_uah_v2_0_0.nc.gz'
        'had4sst4_krig_v2_0_0.nc.gz'
        'had4_krig_v2_0_0.nc.gz'
"""

import logging
import os

import iris
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, vkey, version, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name).copy()
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    attrs = cfg['attributes'].copy()
    attrs['mip'] = var['mip']
    attrs['version'] = version
    baseline = cfg['attributes']['baseline'][vkey]
    attrs['baseline'] = baseline
    attrs['comment'] = attrs['comment'].format(baseline=baseline)
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
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for (vkey, version) in cfg['attributes']['version'].items():
            logger.info("CMORizing variable '%s' version '%s'", short_name,
                        version)
            filepath = raw_filepath.format(version=version)
            _extract_variable(short_name, var, vkey, version, cfg, filepath,
                              out_dir)
