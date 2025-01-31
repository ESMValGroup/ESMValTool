import logging
import re
from pathlib import Path
from pprint import pformat

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)

def _add_longitude_coordinate(cube):
    """Add longitude coordinate to cube."""
    lon_coord = iris.coords.DimCoord(
        [180.0], bounds=[[0.0, 360.0]],
        var_name='lon', standard_name='longitude',
        long_name='longitude', units='degrees_east'
    )
    new_aux_coords = [(c, cube.coord_dims(c)) for c in cube.aux_coords]
    new_dim_coords = [(c, cube.coord_dims(c)) for c in cube.dim_coords]
    new_dim_coords.append((lon_coord, cube.ndim))
    new_metadata = cube.metadata
    cube = iris.cube.Cube(
        cube.core_data()[..., None],
        dim_coords_and_dims=new_dim_coords,
        aux_coords_and_dims=new_aux_coords
    )
    cube.metadata = new_metadata
    return cube

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
    logger.debug("Found input files:\n%s", pformat([str(f['filename']) for f in input_files_dicts]))
    return input_files_dicts

def _extract_variable(short_name, var, cfg, file_dict, out_dir):
    """Extract variable."""
    raw_name = var.get('raw_name', short_name)
    mip = var['mip']
    cmor_info = cfg['cmor_table'].get_variable(mip, short_name)
    filename = file_dict['filename']

    # Extract data
    cube = iris.load_cube(filename, utils.var_name_constraint(raw_name))

    # Fix units
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)
    
    # Add longitude coordinate if missing
    if 'longitude' not in [coord.standard_name for coord in cube.coords()]:
        logger.info("Adding 'longitude' coordinate to variable '%s'", short_name)
        cube = _add_longitude_coordinate(cube)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = mip
    attrs['version'] = f"{file_dict['version']}-{file_dict['resolution']}"
    if 'version_suffix' in var:
        attrs['version'] += f"-{var['version_suffix']}"
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube, short_name, out_dir, attrs, unlimited_dimensions=['time'])

def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization function."""
    input_files_dicts = _get_input_file_dicts(in_dir, cfg)
    
    for file_dict in input_files_dicts:
        for (short_name, var) in cfg['variables'].items():
            logger.info("CMORizing variable '%s' from file '%s'", short_name, str(file_dict['filename']))
            _extract_variable(short_name, var, cfg, file_dict, out_dir)