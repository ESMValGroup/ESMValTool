"""ESMValTool CMORizer for CFSR data.
Tier
    Tier 2: other freely-available dataset.
Source
    ESGF:
    https://esgf.nccs.nasa.gov/thredds/fileServer/CREATE-IP/
        reanalysis/NOAA-NCEP/CFSR/CFSR/mon/atmos/
Last access
    20221122
Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/cfsr.py
"""

import logging
import re
from copy import deepcopy
from pathlib import Path
from warnings import catch_warnings, filterwarnings
from cf_units import Unit

import iris
from iris.coords import CellMethod
from iris import NameConstraint
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.data import utilities as utils


logger = logging.getLogger(__name__)


def _fix_coordinates(cube, definition, cmor_info):
    # fix flipped latitude
    utils.flip_dim_coord(cube, 'latitude')
    # fix other coordinates
    utils.fix_coords(cube)

    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)
    if 'height10m' in cmor_info.dimensions:
        utils.add_scalar_height_coord(cube, height=10.)

    for coord_def in definition.coordinates.values():
        axis = coord_def.axis
        coord = cube.coord(axis=axis)
        if axis == 'Z':
            coord.convert_units(coord_def.units)
        coord.standard_name = coord_def.standard_name
        coord.var_name = coord_def.out_name
        coord.long_name = coord_def.long_name
        coord.points = coord.core_points().astype('float64')
        if coord.var_name == 'plev':
            coord.attributes['positive'] = 'down'

    return cube


def _extract_variable(short_name, var, cfg, raw_filepath, out_dir):
    """Extract variable."""
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], short_name)
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    if cmor_info.positive != '':
        attributes['positive'] = cmor_info.positive

    # load data
    raw_var = var.get('raw', short_name)
    with catch_warnings():
        filterwarnings('ignore',
                       message='Ignoring netCDF variable .* invalid units .*',
                       category=UserWarning,
                       module='iris')
        cube = iris.load_cube(str(raw_filepath),
                              NameConstraint(var_name=raw_var))

    utils.set_global_atts(cube, attributes)

    utils.fix_var_metadata(cube, cmor_info)

    # adjusting the attribute "cell methods"
    cube.cell_methods = ()
    cube.add_cell_method(CellMethod("mean", coords=["time"]))

    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    cube = _fix_coordinates(cube, definition, cmor_info)

    # Save variable
    utils.save_variable(
        cube,
        short_name,
        out_dir,
        attributes,
        unlimited_dimensions=['time'],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        filename = var['file']
        logger.info("CMORizing variable '%s' from file '%s'", short_name,
                    filename)
        short_name = var['short_name']
        print(short_name)
        raw_filenames = Path(in_dir).rglob('*.nc')
        filenames = []
        for raw_filename in raw_filenames:
            if re.search(var['file'], str(raw_filename)) is not None:
                filenames.append(raw_filename)

        for filename in sorted(filenames):

            _extract_variable(short_name, var, cfg, filename, out_dir)

