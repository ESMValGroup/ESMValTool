"""ESMValTool CMORizer for HadCRUT5 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://crudata.uea.ac.uk/cru/data/temperature

Last access
    20210113

Download and processing instructions
    Download the following files:
        infilling
            [Source]/HadCRUT.5.0.0.0.analysis.anomalies.ensemble_mean.nc
        no-infilling
            [Source]/HadCRUT.5.0.0.0.anomalies.ensemble_mean.nc
        climatology
            [Source]/absolute_v5.nc
"""

import copy
import logging
import os

import iris
import numpy as np
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, version, filename, cfg, in_dir,
                      out_dir):
    """Extract variable."""
    # load data
    filepath = os.path.join(in_dir, filename)
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # load climatology
    filepath_clim = os.path.join(in_dir, cfg['climatology']['filename'])
    raw_var = var.get('raw_clim', short_name)
    clim_cube = iris.load_cube(filepath_clim,
                               utils.var_name_constraint(raw_var))

    # fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    for cub in [cube, clim_cube]:
        if cub.units != cmor_info.units:
            cub.convert_units(cmor_info.units)

    # derive absolute temperatures
    clim_data = clim_cube.data
    clim_data = np.tile(clim_data, [cube.shape[0] // 12, 1, 1])
    if cube.shape[0] % 12 != 0:
        for i in range(cube.shape[0] % 12):
            clim_data = np.vstack([clim_data, clim_data[i:i + 1]])

    cube.data = cube.data + clim_data

    # fix time units
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    # Fix coordinates
    utils._fix_dim_coordnames(cube)
    cube_coord = cube.coord('longitude')
    if cube_coord.points[0] < 0. and \
            cube_coord.points[-1] < 181.:
        cube_coord.points = \
            cube_coord.points + 180.
        utils._fix_bounds(cube, cube_coord)
        cube.attributes['geospatial_lon_min'] = 0.
        cube.attributes['geospatial_lon_max'] = 360.
        nlon = len(cube_coord.points)
        utils._roll_cube_data(cube, nlon // 2, -1)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata and  update version information
    attrs = copy.deepcopy(cfg['attributes'])
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


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        for (version, filename) in cfg['filenames'].items():
            logger.info("CMORizing variable '%s' '%s'", short_name, version)
            _extract_variable(short_name, var, version, filename, cfg, in_dir,
                              out_dir)
