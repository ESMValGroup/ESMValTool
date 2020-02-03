"""ESMValTool CMORizer for MLS data.

Tier
    Tier 3: restricted dataset.

Source
    https://disc.gsfc.nasa.gov/datasets/ML2RHI_004/summary?keywords=ML2RHI_004

Last access
    20200203

Download and processing instructions
    Select Data Access -> Subset/Get Data -> Get Data and followe the
    "Instructions for downloading". All *.he5 files need to be saved in the
    $RAWOBS/Tier3/MLS directory, where $RAWOBS refers to the RAWOBS directory
    defined in the user configuration file.
    A registration is required for downloading the data.

"""

import glob
import logging
import os
import shutil
from pprint import pformat

import iris
import iris.coord_categorisation
import netCDF4 as nc
import numpy as np
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(cmor_info, attrs, var_file, out_dir, cfg):
    """Extract variable."""
    grid_file = gdal.Open(var_file)
    array = grid_file.ReadAsArray()
    for missing_value in cfg['missing_values']:
        array = np.ma.masked_equal(array, missing_value)
    array = array.astype(np.float)
    np.ma.set_fill_value(array, 1e20)
    array = np.ma.expand_dims(array, 0)
    time = iris.coords.DimCoord([183.0],
                                bounds=[0.0, 366.0],
                                units=Unit('days since 2000-01-01 00:00:00'),
                                standard_name='time',
                                var_name='time',
                                long_name='time')
    lats = iris.coords.DimCoord(
        90.0 - np.arange(array.shape[1]) * cfg['delta_degrees'],
        standard_name='latitude',
        var_name='lat',
        long_name='latitude')
    lons = iris.coords.DimCoord(
        180.0 + np.arange(array.shape[2]) * cfg['delta_degrees'],
        standard_name='longitude',
        var_name='lon',
        long_name='longitude')
    cube = iris.cube.Cube(array,
                          dim_coords_and_dims=[(time, 0), (lats, 1),
                                               (lons, 2)],
                          units=Unit('t ha-1'))
    cube.convert_units('kg m-2')
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(cube,
                        cmor_info.short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        files = glob.glob(os.path.join(in_dir, var_info['filename']))

        # Extract cubes
        cubes = iris.cube.CubeList()
        for filename in files[:10]:
            dataset = nc.Dataset(filename, mode='r')
            rhi_group = dataset.groups['HDFEOS'].groups['SWATHS'].groups['RHI']
            print(rhi_group)

        # print(pformat(files))
        # glob_attrs['mip'] = var_info['mip']
        # cmor_info = cmor_table.get_variable(var_info['mip'], var)
        # logger.info("Found input file '%s' for variable '%s'", var_file, var)
        # _extract_variable(cmor_info, glob_attrs, var_file, out_dir, cfg)
