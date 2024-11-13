"""ESMValTool CMORizer for Earth System Data Cube data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://data.rsc4earth.de/EarthSystemDataCube/

Last access
    20230126

Download and processing instructions
    It is not necessary to download the data, as the cmorizer script can access
    it directly from the cloud if it is not available locally.

    To download a dataset, the dataset folder can be explored on the source
    website, and downloaded using wget:
        ```wget -m -nH -np -R "index.html*" http://data.rsc4earth.de/EarthSystemDataCube/v3.0.1/```
"""  # noqa: E501
import logging
from copy import deepcopy
from pathlib import Path

import cf_units
import iris.std_names
import xarray as xr
from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _fix_cube(var, cube, cfg):
    """General fixes for all cubes."""
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], var['short_name'])

    # Set correct names
    cube.var_name = cmor_info.short_name
    if cmor_info.standard_name:
        cube.standard_name = cmor_info.standard_name
    cube.long_name = cmor_info.long_name

    # Set calendar to gregorian instead of proleptic gregorian
    old_unit = cube.coord('time').units
    if old_unit.calendar == 'proleptic_gregorian':
        logger.info("Converting time units to gregorian")
        cube.coord('time').units = cf_units.Unit(old_unit.origin,
                                                 calendar='gregorian')
    cube = utils.fix_coords(cube)
    cube.convert_units(cmor_info.units)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)
    # Conversion from 8-d to monthly frequency
    cube = monthly_statistics(cube, operator="mean")

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    return cube


def _open_zarr(path):
    """Open zarr dataset."""
    logger.info('Opening zarr in "%s"', path)
    try:
        zarr_dataset = xr.open_dataset(path, engine='zarr')
        return zarr_dataset
    except KeyError as exception:
        # Happens when the zarr folder is missing metadata, e.g. when
        # it is a zarr array instead of a zarr dataset.
        logger.error('Could not open zarr dataset "%s": "KeyError: %s"', path,
                     exception)
        raise exception


def _extract_variable(zarr_path, var, cfg, out_dir):
    """Open and cmorize cube."""
    attributes = deepcopy(cfg['attributes'])
    all_attributes = {
        **attributes,
        **var
    }  # add the mip to the other attributes
    raw_name = var['raw']
    zarr_dataset = _open_zarr(zarr_path)
    cube_xr = zarr_dataset[raw_name]

    # Invalid standard names must be removed before converting to iris
    standard_name = cube_xr.attrs.get('standard_name', None)
    if (standard_name is not None
            and standard_name not in iris.std_names.STD_NAMES):
        del cube_xr.attrs['standard_name']
        logger.info('Removed invalid standard name "%s".', standard_name)

    cube_iris = cube_xr.to_iris()
    cube = _fix_cube(var, cube_iris, cfg)

    utils.save_variable(cube=cube,
                        var=var['short_name'],
                        outdir=out_dir,
                        attrs=all_attributes,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize the dataset."""
    if start_date:
        logger.warning('start_date set to "%s", but will be ignored',
                       start_date)
    if end_date:
        logger.warning('end_date set to "%s", but will be ignored', end_date)

    attributes = cfg['attributes']
    variables = cfg['variables']
    version = attributes['version']
    filename_pattern = cfg['filename'].format(grid=attributes['grid'],
                                              chunking=attributes['chunking'],
                                              version=version)

    local_path = Path(in_dir)
    in_files = list(local_path.glob(filename_pattern))
    logger.debug('Pattern %s matched: %s', Path(local_path, filename_pattern),
                 in_files)

    if len(in_files) > 1:
        logger.warning(
            'Pattern has matched "%i" files, '
            'but only the first one will be used.', len(in_files))
        logger.warning('The following files will be ignored.: "%s"',
                       in_files[1:])
        zarr_path = in_files[0]
    elif len(in_files) == 0:
        logger.info(
            'No local matches for pattern "%s", '
            'attempting connection to the cloud.',
            Path(local_path, filename_pattern))
        if '*' in filename_pattern:
            logger.warning(
                'Detected a wildcard character in path (*), '
                'online connection to \"%s\" may not work', filename_pattern)
        zarr_path = f'{attributes["source"]}/v{version}/{filename_pattern}'

    for short_name, var in variables.items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        _extract_variable(zarr_path, var, cfg, out_dir)
