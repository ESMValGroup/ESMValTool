"""ESMValTool CMORizer for BSVertOzone data.

Tier
    Tier 3: restricted dataset.

Source
    https://zenodo.org/record/1217184/

Last access
    20210910

Download and processing instructions
    TBA

"""

import logging
from pathlib import Path

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)


def _get_input_file(in_dir, cfg):
    """Get input file."""
    in_path = Path(in_dir).resolve()
    return in_path / cfg['input_file']


def _extract_variable(short_name, var, cfg, input_file, out_dir):
    """Extract variable."""
    mip = var['mip']
    cmor_info = cfg['cmor_table'].get_variable(mip, short_name)
    cube = iris.load_cube(str(input_file))

    # Convert O3 mixing ratios given by BSVertOzone to O3 mole fractions used
    # by CMOR (n: moles of a given substance)
    # Mixing ratio of O3 in air (BSVertOzone): r = n(O3) / n(Air)
    # Mole fraction of O3 in air (CMOR): x = n(O3) / (n(Air) + n(O3))
    # --> x = r / (1 + r)
    cube.convert_units('1')
    cube.data = cube.core_data() / (1.0 + cube.core_data())

    # Fix units
    cube.convert_units(cmor_info.units)

    # Fix coordinates
    cube.coord('air_pressure').var_name = 'plev'
    utils.fix_coords(cube)

    # Switch plev and lat coordinate
    new_data = cube.core_data().swapaxes(1, 2)
    new_coord_spec = [
        (cube.coord('time'), 0),
        (cube.coord('air_pressure'), 1),
        (cube.coord('latitude'), 2),
    ]
    new_metadata = cube.metadata
    cube = iris.cube.Cube(new_data, dim_coords_and_dims=new_coord_spec)
    cube.metadata = new_metadata

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = mip
    utils.set_global_atts(cube, attrs)
    utils.fix_var_metadata(cube, cmor_info)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    input_file = _get_input_file(in_dir, cfg)

    # Run the cmorization for every variable
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, input_file, out_dir)
