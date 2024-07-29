"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
    https://www.copernicus.eu/en/access-data/copernicus-services-catalogue/arctic-ocean-sea-and-ice-surface-temperature
   

Last access
   2024/07

Download and processing instructions

Modification history

"""

import logging
import os
from warnings import catch_warnings, filterwarnings

import iris
import glob
from cf_units import Unit
from datetime import date
from iris.time import PartialDateTime
from esmvalcore.cmor._fixes.shared import  add_scalar_typesi_coord


from esmvaltool.cmorizers.data.utilities import (
    constant_metadata,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _fix_data(cube, var, version):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")

    with constant_metadata(cube):
        if var in ['siconc']:
            cube *= 100.  # Convert from fraction to percentage
            cube.units = '%'

            add_scalar_typesi_coord(cube, value='sea_ice')

    return cube

def collect_files(in_dir, var, cfg):
    """Compose input file list and download if missing."""
    file_list = []
    var_dict = cfg['variables'][var]
    in_dir = os.path.join(in_dir, var_dict['raw_var'])
    print("> Looking for files:", os.path.join( in_dir, var_dict['files']) )
    file_list = glob.glob( os.path.join( in_dir, var_dict['files']) )

    print("> Reading files:", file_list)
    assert len(file_list)>0, "Error: no input files found"

    return file_list

def _select_timeslice(cube, select):
    """Slice a cube along its time axis."""
    if select.any():
        coord = cube.coord('time')
        time_dims = cube.coord_dims(coord)
        if time_dims:
            time_dim = time_dims[0]
            slices = tuple(select if i == time_dim else slice(None)
                           for i in range(cube.ndim))
            cube_slice = cube[slices]
        else:
            cube_slice = cube
    else:
        cube_slice = None
    return cube_slice
    
def get_month_cubes( cube ):

    # get axis of time coordinate
    axis_time = None
    for i,c in enumerate( cube.coords() ):
        if c.name()=='time': axis_time = i ; break
    assert axis_time is not None

    time_coord = cube.coord('time')

    month_cubes = []

    start_day = 1
    start_year = time_coord.units.num2date(time_coord.points[0] ).year

    for i in range(1,13):
        start_month = i
        end_month = i
        if i<12: end_year = start_year ; next_month = start_month+1
        else: end_year = start_year + 1 ; next_month = 1

        # calculate days in month
        end_day = ( date(end_year, next_month, start_day) - date(start_year, start_month, start_day) ).days

        start_datetime = PartialDateTime(year=int(start_year), month=int(start_month), day=int(start_day))
        end_datetime = PartialDateTime(year=int(end_year), month=int(end_month), day=int(end_day))
        dates = time_coord.units.num2date(time_coord.points)

        # slice cube
        select = (dates >= start_datetime) & (dates < end_datetime)
        cube_slice = _select_timeslice(cube, select)

        month_cubes.append( cube_slice.collapsed('time', iris.analysis.MEAN) )

    cubes = iris.cube.CubeList( month_cubes )
    return cubes

def extract_variable(in_files, out_dir, attrs, raw_info, cmor_table):
    """Extract variables and create OBS dataset."""
    var = raw_info['var']
    var_info = cmor_table.get_variable(raw_info['mip'], var)
    rawvar = raw_info['raw_var']

    # read in all in_files
    # for each file, get monthly average
    cubes = []
    valid_min, valid_max = 99999,-99999
    for f in sorted( in_files ):
        cube_year = iris.load_cube(f)
        cube_year_months = get_month_cubes( cube_year )

        iris.util.unify_time_units(cube_year_months)
        iris.util.equalise_attributes(cube_year_months)
        cube_year_months_merged = cube_year_months.merge_cube()

        # ignore this attribute, which differs between all cubes
        time_coord  = cube_year_months_merged.coord('time')
        # record min and max values
        if time_coord.attributes['valid_min'] < valid_min:
            valid_min = time_coord.attributes['valid_min']
        if time_coord.attributes['valid_max'] > valid_max:
            valid_max = time_coord.attributes['valid_max']
        time_coord.attributes.pop('valid_max')
        time_coord.attributes.pop('valid_min')

        cubes.append( cube_year_months_merged )

    cubes = iris.cube.CubeList( cubes )
    iris.util.unify_time_units(cubes)
  
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    # add back attributes
    cube.coord('time').attributes['valid_min'] = valid_min
    cube.coord('time').attributes['valid_max'] = valid_max

    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    _fix_data(cube, var, attrs['version'])
    set_global_atts(cube, attrs)
    save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        in_files = collect_files(in_dir, var, cfg)
        logger.info("CMORizing var %s from input set %s", var, vals['raw_var'])
        raw_info = cfg['variables'][var]
        raw_info.update({
            'var': var,
        })

        glob_attrs['mip'] = vals['mip']
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
