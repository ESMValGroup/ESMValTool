"""ESMValTool CMORizer for SWOOSH data.

Tier
    Tier 3: restricted dataset.

Source
    https://csl.noaa.gov/groups/csl8/swoosh/

Last access
    20210826

Download and processing instructions
    First tick the box that you agree to the terms and conditions for using the
    data. Then download at least one of the following files:

    * 10° resolution: swoosh-v02.6-198401-202104-latpress-10deg-L31.nc

    * 5° resolution: swoosh-v02.6-198401-202104-latpress-5deg-L31.nc

    * 2.5° resolution: swoosh-v02.6-198401-202104-latpress-2.5deg-L31.nc

"""

import logging
import re
from pathlib import Path
from pprint import pformat

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)


def _get_input_file_dicts(in_dir, cfg):
    """Get input files."""
    in_path = Path(in_dir).resolve()
    all_files = list(in_path.glob('*.nc'))
    input_files_dicts = []
    for path in all_files:
        match = re.search(rf"{cfg['input_file_pattern']}", str(path))
        if match is not None:
            metadata = dict(match.groupdict())
            metadata['filename'] = path
            input_files_dicts.append(metadata)
    logger.debug("Found input files:\n%s",
                 pformat([str(f['filename']) for f in input_files_dicts]))
    return input_files_dicts


def _extract_variable(short_name, var, cfg, file_dict, out_dir):
    """Extract variable."""
    raw_name = var.get('raw_name', short_name)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    filename = file_dict['filename']

    # Extract data
    filename = str(file_dict['filename'])
    cube = iris.load_cube(filename, utils.var_name_constraint(raw_name))

    # Fix units
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)
    plev_coord = cube.coord('pressure')
    plev_coord.standard_name = 'air_pressure'
    plev_coord.var_name = 'plev'
    plev_coord.convert_units('Pa')

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    attrs['version'] = f"{file_dict['version']}-{file_dict['resolution']}"
    if 'version_suffix' in var:
        attrs['version'] += f"-{var['version_suffix']}"
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    input_files_dicts = _get_input_file_dicts(in_dir, cfg)

    # Run the cmorization for every file
    # Note: The different files correspond to different versions; each file
    # contains all variables.
    for file_dict in input_files_dicts:
        for (short_name, var) in cfg['variables'].items():
            logger.info(
                "CMORizing variable '%s' from file '%s'", short_name,
                str(file_dict['filename']))
            _extract_variable(short_name, var, cfg, file_dict, out_dir)
