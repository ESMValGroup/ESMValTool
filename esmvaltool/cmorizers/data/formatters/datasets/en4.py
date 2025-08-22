import logging
import os
from dask import array as da
import iris

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

def load_and_prepare_cube(fullpath, raw_var, glob_attrs, cmor_info):
    
    cubes = iris.load(fullpath, raw_var)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    if cube.units == 'K':
        cube.convert_units('degC')

    cube.coord('depth').units = 'm'   # Was metres but gets changed when running recipe
    cube = utils.fix_coords(cube)   
    utils.fix_var_metadata(cube, cmor_info)   
    utils.set_global_atts(cube, glob_attrs)
    
    return cube

def extract_surface_var(cube, srf_var, cmor_info):
    logger.info('Extracting surface level')

    depth0 = iris.Constraint(depth=cube.coord('depth').points[0])
    surface_cube = cube.extract(depth0)

    utils.fix_var_metadata(surface_cube, cmor_info) 

    return surface_cube

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORization main function."""

    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']
    filename = cfg['filename']
    fullpath = os.path.join(in_dir, filename)

    for var, var_info in cfg['variables'].items():
        logger.info(f'Loading {fullpath}')

        glob_attrs["mip"] = var_info["mip"]

        raw_var = var_info['raw_var']
        srf_var = var_info['srf_var']

        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        cmor_info_srf = cmor_table.get_variable(var_info['mip'], srf_var)

        cube = load_and_prepare_cube(fullpath, raw_var, glob_attrs, cmor_info)
        surface_cube = extract_surface_var(cube, srf_var, cmor_info_srf)
        logger.info(f'Saving for {var}')
        utils.save_variable(cube, var, out_dir, glob_attrs)

        logger.info(f'Saving for {srf_var}')
        utils.save_variable(surface_cube, srf_var, out_dir, glob_attrs)
