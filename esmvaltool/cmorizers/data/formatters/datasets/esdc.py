"""ESMValTool CMORizer for Earth System Data Cube data.

<We will add some useful info here later>
"""
import logging
from pathlib import Path

import xarray as xr
from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def fix_time_coordinate(cube):
    time = cube.coord(axis='T')
    time.convert_units('days since 1850-1-1 00:00:00.0')


def fix_cube(short_name, var, cube, cfg):
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    # fix_time_coordinate(cube)
    utils.fix_coords(cube)
    logger.info(cmor_info.units)
    cube.convert_units(cmor_info.units)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)
    cube = monthly_statistics(cube, operator="mean")  # Fix frequency
    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize the dataset."""

    # This is where you'll add the cmorization code
    # 1. find the input data
    logger.debug("in_dir: '%s'", in_dir)
    logger.debug("cfg: '%s'", cfg)
    logger.debug("cfg_user: '%s'", cfg_user)

    attributes = cfg['attributes']
    variables = cfg['variables']

    version = attributes['version']

    filename_pattern = cfg['filename'].format(grid=attributes['grid'],
                                              chunking=attributes['chunking'],
                                              version=version)

    def cmorize_cube(path):
        logger.info(f'Opening zarr in "{path}"')
        try:
            ds = xr.open_zarr(path, consolidated=True)
        except KeyError as e:
            logger.info(
                f'Could not open zarr dataset "{path}": "KeyError: {e}"')
            return

        for short_name, var in variables.items():
            all_attributes = {
                **attributes,
                **var
            }  # add the mip to the other attributes
            raw_name = var['raw']
            cube_xr = ds[raw_name]

            cube_iris = cube_xr.to_iris()

            # 2. apply the necessary fixes
            cube = fix_cube(short_name, var, cube_iris, cfg)

            # 3. store the data with the correct filename

            utils.save_variable(cube=cube,
                                var=short_name,
                                outdir=out_dir,
                                attrs=all_attributes)

    matches = [m for m in Path(in_dir, f'v{version}').glob(filename_pattern)]
    logger.debug(
        f'Pattern {Path(in_dir, version, filename_pattern)} matched: {matches}'
    )

    if len(matches) != 0:
        for match in matches:
            cmorize_cube(match)
    else:
        logger.info(
            "Dataset not found locally, attempting connection to the cloud.")
        if '*' in filename_pattern:
            logger.warning(
                f"For cloud connection, \"{filename_pattern}\" shouldn't contain wildcards"
            )
        cmorize_cube(f'{attributes["source"]}/v{version}/{filename_pattern}')
