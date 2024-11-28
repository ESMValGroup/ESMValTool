"""ESMValTool CMORizer for RSS data.

Tier
    Tier 2: other freely available dataset.

Source
    ftp://ftp.remss.com/vapor/monthly_1deg/

Last access
    20201119

Download and processing instructions
    A registration is required for downloading the data, but no licence
    agreement necessary. Download vapor, here ncdf3 file is used.

Modification history
   20170419-A_gier_bettina: written.
   20201201-A_weigel_katja: portet to ESMValTool v2.
   20231214-A_weigel_katja: update to current ESMValTool version.
   20241122-A_weigel_katja: changed to python.

"""
import glob
import logging
import os
import xarray as xr
from copy import deepcopy

from dask import array as da
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

def _load_cube(in_files, var):
    
    dataset = xr.open_mfdataset(in_files, engine="netcdf4")
    da_tmp = dataset[var['raw']]
    cube = da_tmp.to_iris()

    cube.data = da.ma.masked_equal(cube.core_data(), -500)

    return cube


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    axis2def = {'T': 'time', 'X': 'longitude', 'Y': 'latitude'}
    axes = ['T', 'X', 'Y']
    for axis in axes:
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)

            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name
            coord.points = coord.core_points().astype('float64')

            if len(coord.points) > 1:
                coord.guess_bounds()
        else:
            raise logger.error("Bounds for coordinate %s "
                               "cannot be guessed", coord)

    return cube


def _extract_variable(in_files, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input files '%s'",
                var['short_name'], ', '.join(in_files))
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    attributes['raw'] = var['raw']

    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = _load_cube(in_files, var)

    # keep the following raw cube attributes
    attrs_to_keep = [
        "institution", "Institution",
        "institute_id", "VersionID",
        "experiment_id",
        "source", "Source",  # overrides empty string default
        "model_id", "ModelID",
        "contact", "Contact",
        "references",
        "tracking_id",
        "mip_specs",  # described by "mip" already
        "source_id", "SourceID",
        "product", "Product",
        "frequency", "Frequency",
        "creation_date",
        "project_id", "ProjectID",
        "table_id", "TableID",
        "title", "Title",
        "modeling_realm",
        "doi",
        "VersionID",  # described by "version" already
    ]

    attrs_to_keep_exist = [
        att for att in cube.attributes if att in attrs_to_keep
    ]
    for att in attrs_to_keep_exist:
        attributes[att] = cube.attributes[att]

    utils.set_global_atts(cube, attributes)

    # Set correct names
    cube.var_name = definition.short_name
    cube.long_name = definition.long_name

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Roll longitude
    cube.coord('longitude').points = cube.coord('longitude').points + 180.
    nlon = len(cube.coord('longitude').points)
    cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=-1)

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)
    
    cube.coord('latitude').attributes = None
    cube.coord('longitude').attributes = None

    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    utils.save_variable(cube, cube.var_name,
                        out_dir, attributes,
                        unlimited_dimensions=['time'])
    logger.info("Finished CMORizing %s", ', '.join(in_files))


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Run CMORizer for RSS."""
    cfg.pop('cmor_table')
    
    for short_name, var in cfg['variables'].items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        logger.info("short_name, var KW")
        logger.info(short_name, var)
        logger.info(var['file'])
        logger.info(var['file'].format())
        # Now get list of files
        filepattern = os.path.join(in_dir, var['file'].format())
        in_files = glob.glob(filepattern)
        if not in_files:
                logger.warning('%s does not exist', filepattern)
                continue
        _extract_variable(in_files, var, cfg, out_dir)
