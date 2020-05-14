"""ESMValTool CMORizer for FLUXCOM GPP data.

Tier
    Tier 3: restricted dataset.

Source
    http://www.bgc-jena.mpg.de/geodb/BGI/Home

Last access
    20190727

Download and processing instructions
    From the website, select FLUXCOM as the data choice and click download.
    Two files will be displayed. One for Land Carbon Fluxes and one for
    Land Energy fluxes. The Land Carbon Flux file (RS + METEO) using
    CRUNCEP data file has several data files for different variables.
    The data for GPP generated using the
    Artificial Neural Network Method will be in files with name:
    GPP.ANN.CRUNCEPv6.monthly.*.nc
    A registration is required for downloading the data.
    Users in the UK with a CEDA-JASMIN account may request access to the jules
    workspace and access the data.
    Note : This data may require rechunking of the netcdf files.
    This constraint will not exist once iris is updated to
    version 2.3.0 Aug 2019
"""
import logging
import os
import re
import numpy as np
import iris
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepath(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    regex = re.compile(basename)

    all_files = [
        f for f in os.listdir(in_dir)
        if os.path.isfile(os.path.join(in_dir, f))
    ]
    for filename in all_files:
        if regex.match(filename):
            return os.path.join(in_dir, basename)
    raise OSError(
        f"Cannot find input file matching pattern  '{basename}' in '{in_dir}'")


def _extract_variable(cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    logger.info("Var is %s", var)
    cubes = iris.load(filepath)
    for cube in cubes:
        # convert data from gc/m2/day to kg/m2/s
        cube = cube / (1000 * 86400)
        cube.units = 'kg m-2 s-1'

        # The following two lines are needed for iris.util.guess_coord_axis
        cube.coord('lat').standard_name = 'latitude'
        cube.coord('lon').standard_name = 'longitude'
        utils.fix_var_metadata(cube, cmor_info)
        utils.convert_timeunits(cube, 1950)
        utils.fix_coords(cube)
        utils.set_global_atts(cube, attrs)
        utils.flip_dim_coord(cube, 'latitude')
        coord = cube.coord('latitude')
        coord.bounds = np.flip(coord.bounds, axis=1)
        logger.info("Saving file")
        utils.save_variable(cube,
                            var,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = _get_filepath(in_dir, cfg['filename'])
    logger.info("Found input file '%s'", filepath)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        logger.info(var_info['mip'])
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        _extract_variable(cmor_info, glob_attrs, filepath, out_dir)
