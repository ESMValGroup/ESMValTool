"""Common tools to CMORize NSIDC-0116 northern and sothern data."""

import logging
import os
import glob
import numpy as np
import iris
from iris.coords import AuxCoord
from iris.cube import Cube


from esmvaltool.cmorizers.obs.utilities import fix_var_metadata, \
    save_variable, set_global_atts

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
        lon_coord.points[lon_coord.points < 0] += 360

        for var, vals in cfg['variables'].items():
            var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
            logger.info('Cmorizing var %s', var)
            cube = cubes.extract_strict(iris.Constraint(vals['raw']))
            cube.add_aux_coord(lat_coord, (1, 2))
            cube.add_aux_coord(lon_coord, (1, 2))
            cube.convert_units(var_info.units)
            logger.debug(cube)
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

    _create_areacello(cfg, cube, glob_attrs, out_dir)


def _create_areacello(cfg, sample_cube, glob_attrs, out_dir):
    if not cfg['custom'].get('create_areacello', False):
        return
    var_info = cfg['cmor_table'].get_variable('fx', 'areacello')
    glob_attrs['mip'] = 'fx'
    lat_coord = sample_cube.coord('latitude')
    cube = Cube(
        np.full(lat_coord.shape, cfg['custom']['grid_cell_size'], np.float32),
        standard_name=var_info.standard_name,
        long_name=var_info.long_name,
        var_name=var_info.short_name,
        units='m2',
    )
    cube.add_aux_coord(lat_coord, (0, 1))
    cube.add_aux_coord(sample_cube.coord('longitude'), (0, 1))
    cube.add_dim_coord(sample_cube.coord('projection_y_coordinate'), 0)
    cube.add_dim_coord(sample_cube.coord('projection_x_coordinate'), 1)
    fix_var_metadata(cube, var_info)
    set_global_atts(cube, glob_attrs)
    save_variable(
        cube, var_info.short_name, out_dir, glob_attrs, zlib=True
    )


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
