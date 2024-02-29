"""ESMValTool CMORizer for ESACCI-PERMAFROST data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/permafrost/data

Last access
   20240227

Download and processing instructions
   Download the data from:
     active_layer_thickness/L4/area4/pp/v03.0
     ground_temperature/L4/area4/pp/v03.0
     permafrost_extent/L4/area4/pp/v03.0
   Put all files in a single directory.
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime
from dateutil import relativedelta

from netCDF4 import Dataset
from cdo import *
import os.path

import cf_units
import iris
import numpy as np
from dask import array as da
from esmvalcore.cmor.table import CMOR_TABLES
#from esmvalcore.preprocessor import regrid
#from esmvalcore.preprocessor.regrid_schemes import (
#    ESMPyAreaWeighted, ESMPyLinear, ESMPyNearest, UnstructuredNearest)
#from esmf_regrid.schemes import ESMFAreaWeighted, ESMFBilinear

from esmvaltool.cmorizers.data.utilities import (
    save_variable, set_global_atts)

#from iris import NameConstraint
from iris.cube import Cube

logger = logging.getLogger(__name__)
#cdo = Cdo()


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    axis2def = {'T': 'time', 'X': 'longitude', 'Y': 'latitude'}
    axes = ['T', 'X', 'Y']

    for axis in axes:
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)
            if axis == 'T':
                coord.convert_units('days since 1850-1-1 00:00:00.0')
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name
            coord.points = coord.core_points().astype('float64')
            if len(coord.points) > 1:
                if coord.bounds is not None:
                    coord.bounds = None
                coord.guess_bounds()

    return cube


def _regrid_infile(infile, outfile):
    """Regrid infile to 0.5 deg x 0.5 deg grid using cdo."""
    cdo = Cdo()
    # ESACCI-PERMAFROST v3.0 dimensions of raw input data
    xsize = 14762
    ysize = 10353
    totalsize = xsize * ysize

    # description of ESACCI-PERMAFROST v3.0 polar stereographic
    # grid for cdo
    esagrid = (f'gridtype  = projection\n'
               f'gridsize  = {totalsize}\n'
               f'xsize     = {xsize}\n'
               f'ysize     = {ysize}\n'
               f'xname     = x\n'
               f'xlongname = "x coordinate of projection"\n'
               f'xunits    = "m"\n'
               f'yname     = y\n'
               f'ylongname = "y coordinate of projection"\n'
               f'yunits    = "m"\n'
               f'xfirst    = -6111475.22239475\n'
               f'xinc      = 926.625433138333\n'
               f'yfirst    = 4114895.09469662\n'
               f'yinc      = -926.625433138333\n'
               f'grid_mapping = polar_stereographic\n'
               f'grid_mapping_name = polar_stereographic\n'
               f'straight_vertical_longitude_from_pole = 0.\n'
               f'false_easting = 0.\n'
               f'false_northing = 0.\n'
               f'latitude_of_projection_origin = 90.\n'
               f'standard_parallel = 71.\n'
               f'longitude_of_prime_meridian = 0.\n'
               f'semi_major_axis = 6378137.\n'
               f'inverse_flattening = 298.257223563\n'
               f'proj_params = "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0'
               f' +k=1" +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"\n')

    esagrid_file = "./esacci_grid.txt"

    # write grid description to ASCII file
    file = open(esagrid_file, "w")
    file.write(esagrid)
    file.close()

    # define dimensions of target grid (regular lat-lon grid)
    target_dimx = 720  # delta_lon = 0.5 deg
    target_dimy = 360  # delta_lat = 0.5 deg
    target_grid = f'r{target_dimx}x{target_dimy}'

    # check if suitable weights file already exists
    # (e.g. from previous call to _regrid_file)

#    weightsfile = "./esa_weights.nc"
    weightsfile = "/work/bd0854/b380103/esmvaltool_output/esa_weights.nc"
    weightsfile_ok = False

    if os.path.isfile(weightsfile):
        weights = Dataset(weightsfile, "r")
        # make sure dimensions of source and target grids match
        # expected values
        src = weights.variables['src_grid_dims']
        dst = weights.variables['dst_grid_dims']
        if (xsize == src[0] and ysize == src[1] and
            target_dimx == dst[0] and target_dimy == dst[1]):
            logger.info("Using matching weights file %s for regridding.",
                        weightsfile)
            weightsfile_ok = True
        weights.close()

    # if no suitable weights file, generate new weights for regridding

    if not weightsfile_ok:
        logger.info("Generating regridding weights. This will take"
                    " about 5-10 minutes (or more)...")
        cdo.genbil(f"{target_grid} -setgrid,{esagrid_file}",
                   input=infile, output=weightsfile, options="-f nc")

    # now regrid data to 0.5 deg x 0.5 deg
    cdo.remap(f"{target_grid},{weightsfile} -setgrid,{esagrid_file}",
              input=infile, output=outfile, options='-f nc')
    return


def _extract_variable(in_file, var, cfg, out_dir, year):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']
    attributes['raw'] = var['raw']
    cmor_table = CMOR_TABLES[attributes['project_id']]
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    # regrid input file using cdo
    # (using the preprocessor or ESMF regrid is too slow)
    
    regridded_file = f"./{year}_{var['short_name']}.nc"
    print(regridded_file)
    _regrid_infile(in_file, regridded_file)
    exit()

    # load input file
    cube = iris.load_cube(regridded_file)

    # --> drop attributes that differ among input files
    # global attributes to remove
    drop_attrs = [
        'source', 'date_created', 'history', 'tracking_id',
        'id', 'time_coverage_start', 'time_coverage_end', 'platform',
        'sensor', 'keywords'
    ]
    # variable attributes to remove
    drop_var_attrs = [
        'flag_meanings', 'flag_values', 'grid_mapping', 'actual_range',
        'ancillary_variables'
    ]
    for attr in drop_attrs:
        if attr in cube.attributes.keys():
            cube.attributes.pop(attr)
    for attr in drop_var_attrs:
        if attr in cube.attributes.keys():
            cube.attributes.pop(attr)

    set_global_atts(cube, attributes)

#    iris.util.unify_time_units(cube)
    cube.coord('time').points = cube.coord('time').core_points().astype(
        'float64')

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

#    # Fix units
#    # input variable for snc (sncf) reports 'percent' --> rename to '%'
#    # input variable for snw (swe) reports 'mm' --> rename to 'kg m-2'
#    cube.units = definition.units

    # Fix data type
    cube.data = cube.core_data().astype('float32')

    # Fix coordinates
#    cube = _fix_coordinates(cube, definition)
#    cube.coord('latitude').attributes = None
#    cube.coord('longitude').attributes = None

#    cube.coord('projection_y_coordinate').standard_name = "latitude"
#    cube.coord('projection_x_coordinate').standard_name = "longitude"
#    print(cube)

    dlon = 0.5
    dlat = 0.5
    mid_dlon, mid_dlat = dlon / 2, dlat / 2

    latdata = np.linspace(-90.0, 90.0, int(180.0 / dlat) + 1)
    londata = np.linspace(0.0, 360.0 - dlon, int(360.0 / dlon))

    lats = iris.coords.DimCoord(latdata,
                                standard_name='latitude',
                                units='degrees_north',
                                var_name='lat',
                                circular=False)

    lons = iris.coords.DimCoord(londata,
                                standard_name='longitude',
                                units='degrees_east',
                                var_name='lon',
                                circular=False)

    lats.guess_bounds()
    lons.guess_bounds()

    # Construct the resultant stock cube, with dummy data.
    shape = (latdata.size, londata.size)
    dummy = np.empty(shape, dtype=np.dtype('int8'))
    coords_spec = [(lats, 0), (lons, 1)]
    target_grid_cube = Cube(dummy, dim_coords_and_dims=coords_spec)

#    regridded_cube = cube.regrid(target_grid_cube, ESMPyLinear())
#    regridded_cube = cube.regrid(target_grid_cube, ESMPyAreaWeighted())
#    regridded_cube = cube.regrid(target_grid_cube, ESMPyNearest())
#    regridded_cube = cube.regrid(target_grid_cube, ESMFAreaWeighted())
#    regridded_cube = cube.regrid(target_grid_cube, ESMFBilinear())
#    regridded_cube = cube.regrid(target_grid_cube, iris.analysis.PointInCell())
#    regridded_cube = cube.regrid(target_grid_cube, iris.analysis.Nearest())
#    regridded_cube = cube.regrid(target_grid_cube, iris.analysis.UnstructuredNearest())
#    regridded_cube = cube.regrid(target_grid_cube, UnstructuredNearest())

#    # regridding from polar stereographic to 0.5x0.5
#    cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')
#    cube = regrid(cube, target_grid='0.5x0.5', scheme='nearest')
#    cube = regrid(cube, target_grid='0.5x0.5', scheme='linear')

#    cube.attributes.update({"geospatial_lon_resolution": "0.5",
#                            "geospatial_lat_resolution": "0.5",
#                            "spatial_resolution": "0.5"})

    # Save results
    logger.debug("Saving cube\n%s", regridded_cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    save_variable(regridded_cube, regridded_cube.var_name,
                  out_dir, attributes,
                  unlimited_dimensions=['time'])
    logger.info("Finished CMORizing %s", in_file)
    
    exit()


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORize ESACCI-PERMAFROST dataset."""
    glob_attrs = cfg['attributes']

    logger.info("Starting CMORization for tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing ESACCI-PERMAFROST version %s",
                glob_attrs['version'])

    if start_date is None:
        start_date = datetime(2003, 1, 1)
    if end_date is None:
        end_date = datetime(2019, 12, 31)

    loop_date = start_date
    while loop_date <= end_date:
        for short_name, var in cfg['variables'].items():
            if 'short_name' not in var:
                var['short_name'] = short_name
            filepattern = os.path.join(
                in_dir, var['file'].format(year=loop_date.year)
                )
            in_file = glob.glob(filepattern)[0]
            if not in_file:
                logger.info(f'{loop_date.year}: no data not found for '
                            f'variable {short_name}')
            else:
                _extract_variable(in_file, var, cfg, out_dir, loop_date.year)

        loop_date += relativedelta.relativedelta(years=1)
