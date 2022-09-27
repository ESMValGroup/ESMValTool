"""ESMValTool CMORizer for Earth System Data Cube data.

<We will add some useful info here later>
"""
import logging
from pathlib import Path

import xarray as xr
from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def fix_cube(short_name, var, cube, cfg):
    """General fixes for all cubes"""
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    utils.fix_coords(cube)
    cube.convert_units(cmor_info.units)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)
    cube = monthly_statistics(cube, operator="mean")  # Fix frequency
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize the dataset."""
    # 1. find the input data
    logger.debug("in_dir: '%s'", in_dir)
    logger.debug("cfg: '%s'", cfg)
    logger.debug("cfg_user: '%s'", cfg_user)

    if start_date:
        logger.warning("start_date set to %s, but will be ignored", start_date)
    if end_date:
        logger.warning("end_date set to %s, but will be ignored", end_date)

    attributes = cfg['attributes']
    variables = cfg['variables']

    version = attributes['version']

    filename_pattern = cfg['filename'].format(grid=attributes['grid'],
                                              chunking=attributes['chunking'],
                                              version=version)

    def cmorize_cube(path):
        logger.info('Opening zarr in "%s"', path)
        try:
            # dataset = xr.open_zarr(path, consolidated=True)
            # TODO: Where to put consolidated, kwargs or backend_kwargs?
            dataset = xr.open_dataset(path, engine='zarr', consolidated=True)
        except KeyError as exception:
            logger.info('Could not open zarr dataset "%s": "KeyError: %s"',
                        path, exception)
            return

        for short_name, var in variables.items():
            all_attributes = {
                **attributes,
                **var
            }  # add the mip to the other attributes
            raw_name = var['raw']
            cube_xr = dataset[raw_name]

            cube_iris = cube_xr.to_iris()

            # 2. apply the necessary fixes
            cube = fix_cube(short_name, var, cube_iris, cfg)

            # 3. store the data with the correct filename

            utils.save_variable(cube=cube,
                                var=short_name,
                                outdir=out_dir,
                                attrs=all_attributes)

    matches = list(Path(in_dir, f'v{version}').glob(filename_pattern))
    logger.debug('Pattern %s matched: %s',
                 Path(in_dir, version, filename_pattern), matches)

    if len(matches) != 0:
        for match in matches:
            cmorize_cube(match)
    else:
        logger.info(
            "Dataset not found locally, attempting connection to the cloud.")
        if '*' in filename_pattern:
            logger.warning(
                "For cloud connection, \"%s\" shouldn't contain wildcards",
                filename_pattern)
        cmorize_cube(f'{attributes["source"]}/v{version}/{filename_pattern}')
