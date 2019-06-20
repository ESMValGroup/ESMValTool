"""CMORizer for ERA5."""
import logging
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import iris
import numpy as np

from . import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(in_file, raw_name, definition, attributes, out_dir):

    cube = iris.load_cube(
        str(in_file),
        constraint=utils.var_name_constraint(raw_name),
    )

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Fix coordinates
    cube.coord('latitude').var_name = 'lat'
    cube.coord('longitude').var_name = 'lon'

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float64')
        coord.guess_bounds()

    if 'height2m' in definition.dimensions:
        utils.add_height2m(cube)

    # Fix units if required
    cube.convert_units(definition.units)

    # Make latitude increasing
    cube = cube[:, ::-1, ...]

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    logger.info("Saving cube\n%s", cube)
    logger.info("Expected output size is %.1fGB",
                np.prod(cube.shape) * 4 / 2**30)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cfg = utils.read_cmor_config('ERA5.yml')
    cfg['attributes']['comment'] = cfg['attributes']['comment'].format(
        year=datetime.now().year)

    for short_name, var in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        attributes = deepcopy(cfg['attributes'])
        attributes['mip'] = var['mip']
        definition = cfg['cmor_table'].get_variable(var['mip'], short_name)

        for in_file in Path(in_dir).glob(var['file']):
            logger.info("CMORizing input file '%s'", in_file)
            _extract_variable(in_file, var['raw'], definition, attributes,
                              out_dir)
