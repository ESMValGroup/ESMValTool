"""Common tools to CMORize NSIDC-0116 northern and sothern data."""

import logging
import os
import glob
import iris
from iris.coords import AuxCoord


from .utilities import fix_var_metadata, save_variable, set_global_atts

logger = logging.getLogger(__name__)


def cmorize(cfg, region, in_dir, out_dir):
    """Cmorize NSIDC-0116 dataset."""
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    file_expr = os.path.join(in_dir, 'icemotion_daily_{}_*.nc'.format(region))
    for filepath in glob.glob(file_expr):
        logger.info('Cmorizing file %s', filepath)
        cubes = iris.load(filepath)
        logger.debug(cubes)
        lat_coord = _create_coord(cubes, 'lat', 'latitude')
        lon_coord = _create_coord(cubes, 'lon', 'longitude')
        for var, vals in cfg['variables'].items():
            logger.info('Cmorizing var %s', var)
            cube = cubes.extract_strict(iris.Constraint(vals['raw']))
            cube.add_aux_coord(lat_coord, (1, 2))
            cube.add_aux_coord(lon_coord, (1, 2))
            logger.debug(cube)
            var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
            glob_attrs['mip'] = vals['mip']
            fix_var_metadata(cube, var_info)
            set_global_atts(cube, glob_attrs)
            zlib = vals.get('compress', False)
            if zlib:
                # Realize data to speed-up writing
                # pylint: disable=pointless-statement
                cube.data
            save_variable(cube, var, out_dir, glob_attrs, zlib=zlib)
            cubes.remove(cube)
            del cube


def _create_coord(cubes, var_name, standard_name):
    cube = cubes.extract_strict(standard_name)
    coord = AuxCoord(
        cube.data,
        standard_name=standard_name,
        long_name=cube.long_name,
        var_name=var_name,
        units=cube.units,
    )
    return coord
