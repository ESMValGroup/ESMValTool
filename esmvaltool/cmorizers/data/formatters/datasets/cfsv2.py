"""
ESMValTool CMORizer for CFSv2 data.

Tier
    Tier 2: other freely-available dataset.

Source
    Research Data Archive (RDA):
    https://rda.ucar.edu/datasets/ds094.2/

Last access
    20230403

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/cfsv2.py
"""

import copy
import glob
import logging
import os

from datetime import datetime
from cf_units import Unit

import iris
import pygrib

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_coords(sel_grb):
    """Get correct coordinates for cube."""
    # Time
    point = datetime(year=sel_grb.year, month=sel_grb.month, day=sel_grb.day)
    time_units = Unit('days since 1950-01-01 00:00:00', calendar='standard')
    time_coord = iris.coords.DimCoord(time_units.date2num(point),
                                      var_name='time',
                                      standard_name='time',
                                      long_name='time',
                                      units=time_units)

    # Latitude
    lat_coord = iris.coords.DimCoord(sel_grb.distinctLatitudes,
                                     standard_name='latitude',
                                     long_name='latitude',
                                     var_name='lat',
                                     units='degrees_north')

    # Longitude
    lon_coord = iris.coords.DimCoord(sel_grb.distinctLongitudes,
                                     standard_name='longitude',
                                     long_name='longitude',
                                     var_name='lon',
                                     units='degrees_east')

    return time_coord, [(lat_coord, 0), (lon_coord, 1)]


def _load_cfsv2_grib2(infile, var):
    """Load data from GRIB file and return list of cubes."""
    cubelist = iris.cube.CubeList()
    # create list of files (needed in case 'infile' contains wildcards)
    listing = sorted(glob.glob(infile), key=os.path.basename)
    for fname in listing:
        logger.info("Reading file '%s'", fname)

        # reading grib with pygrib
        grbs = pygrib.open(fname)
        grbs.seek(0)
        sel_grb = grbs.select(shortName=var.get('raw'))
        if len(sel_grb) > 1:
            grb = [grb for grb in sel_grb
                   if grb.typeOfLevel == var.get('level')][0]
        else:
            grb = sel_grb[0]

        logger.info("Reading: %s", grb)

        raw_data, lats, lons = grb.data()

        # Get coordinates
        time_coord, coords = _get_coords(grb)

        # Build cube
        tmp_cube = iris.cube.Cube(raw_data, dim_coords_and_dims=coords)
        tmp_cube.add_aux_coord(time_coord)

        cubelist.append(tmp_cube)

    cube = cubelist.merge_cube()

    return cube


def _extract_variable(short_name, var, in_files, cfg, out_dir):
    """Extract variable."""
    # load data (returns a list of cubes)

    cube = _load_cfsv2_grib2(in_files, var)

    # Fix metadata
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)

#    # fix z-coordinate (if present)
#    for coord in cube.dim_coords:
#        coord_type = iris.util.guess_coord_axis(coord)
#        if coord_type == 'Z':
#            coord.standard_name = 'air_pressure'
#            coord.long_name = 'pressure'
#            coord.var_name = 'plev'
#            coord.attributes['positive'] = 'down'
#            if coord.units == "hPa":
#                coord.convert_units('Pa')
#            utils.flip_dim_coord(cube, coord.standard_name)

    utils.fix_dim_coordnames(cube)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        short_name = var['short_name']

        filename = os.path.join(in_dir, var['file'])

        logger.info("CMORizing variable '%s' from file '%s'", short_name,
                    filename)
        _extract_variable(short_name, var, filename, cfg, out_dir)
