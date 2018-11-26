"""Convenience functions for writing netcdf files."""
import logging
import os
from datetime import datetime

import iris

logger = logging.getLogger(__name__)


def save_iris_cube(cube, path, cfg):
    """Save `iris.cube.Cube` and append ESMValTool information.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube to be saved.
    path : str
        Desired path.
    cfg : dict
        Diagnostic script configuration.

    """
    if not cfg['write_netcdf']:
        logger.warning("Could not write netcdf file '%s', 'write_netcdf' is "
                       "set to 'False' in user configuration file.", path)
        return
    attr = {
        'created_by':
        'ESMValTool version {}'.format(cfg['version']) +
        ', diagnostic {}'.format(cfg['script']),
        'creation_date':
        datetime.utcnow().isoformat(' ') + 'UTC'
    }
    cube.attributes.update(attr)
    iris.save(cube, path)
    logger.info("Wrote %s", path)


def save_scalar_data(data, filename, cfg, **cube_atts):
    """Save scalar data for multiple datasets.

    Create 1D cube with the auxiliary dimension `dataset` and save scalar data
    for every appearing dataset.

    Parameters
    ----------
    data : dict
        Scalar data (values) and corresponding datasets (keys).
    filename : str
        Desired name of the file (without extension).
    cfg : dict
        Diagnostic script configuration.
    cube_atts : optional keyword arguments
        Attributes for the cube (`var_name`, `standard_name`, `long_name`,
        `units`, etc.)

    """
    dataset_coord = iris.coords.AuxCoord(list(data), long_name='dataset')
    cube = iris.cube.Cube(
        list(data.values()),
        aux_coords_and_dims=[(dataset_coord, 0)],
        **cube_atts)
    path = os.path.join(cfg['work_dir'], filename + '.nc')
    save_iris_cube(cube, path, cfg)
