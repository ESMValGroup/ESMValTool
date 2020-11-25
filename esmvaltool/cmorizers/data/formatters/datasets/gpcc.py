"""ESMValTool CMORizer for GPCC data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
    https://opendata.dwd.de/climate_environment/GPCC/full_data_2018/full_data_monthly_v2018_[025 05 10 25].nc.gz # noqa
Last access
    20200225

Download and processing instructions
    Download the following files:
        full_data_monthly_{version}.nc.gz

"""

import logging
import os
from warnings import catch_warnings, filterwarnings
import cf_units

import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _extract_variable(short_name, var, version, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    utils.set_units(cube, var.get('raw_units', short_name))
    # fix calendar type
    cal_time = var.get('calendar', short_name)
    origin_time = cube.coord('time').units.origin
    cube.coord('time').units = cf_units.Unit(origin_time, calendar=cal_time)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes']
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


def cmorization(in_dir, out_dir, cfg, _, __, ___):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for version in cfg['attributes']['version'].values():
        for (short_name, var) in cfg['variables'].items():
            raw_var = var.get('raw', short_name)
            filepath = raw_filepath.format(version=version, raw_name=raw_var)
            logger.info("CMORizing variable '%s'", short_name)
            _extract_variable(short_name, var, version, cfg, filepath, out_dir)
            _clean(filepath)
