"""ESMValTool CMORizer for ESACCI-CLOUD data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/L3U/AVHRR-PM/

Last access
    20230619

Download and processing instructions
    see downloading script
"""

import copy
import glob
import logging
import os

import cf_units
import iris
import numpy as np
from cf_units import Unit
from iris import NameConstraint
from calendar import monthrange
from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.data import utilities as utils
from esmvalcore.preprocessor import regrid

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, year, month, cfg, in_dir,
                      out_dir):
    """Extract variable."""

    def adjust_dim(dim):
        if dim == 0:
            return dim
        return dim + 1

    cube = iris.cube.CubeList()
    filename = f'{year}{month:02}*' + var['file']
    #filename = f'{year}{month:02}01' + var['file']
    filelist = glob.glob(os.path.join(in_dir, filename))
    print(filelist)
    num_days = monthrange(year, month)[1]
    #for filename in sorted(filelist):
    for iday in range(1, num_days+1):

        print(iday)
        ifile = glob.glob(os.path.join(in_dir, f'{year}{month:02}' + f'{iday:02}' + var['file']))
        #ifile = os.path.join(in_dir, f'{year}{month:02}' + f'{iday:02}' + var['file'])
        print(ifile)

        #if ifile in filelist:
        #if os.path.isfile(ifile[0]):
        if ifile:
            logger.info("CMORizing file %s", ifile)

            # load data
            raw_var = var.get('raw', short_name)
            daily_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_var))
            daily_cube.attributes.clear()

            # Fix coordinates
            daily_cube = utils.fix_coords(daily_cube)
            utils.convert_timeunits(daily_cube, 1950)

            if iday == 1:
                fill_cube = daily_cube
                print(fill_cube.coord('time').points)

        else:

            logger.info("Fill missing day %s in month %s and year %s", iday, month, year)

            #daily_cube = fill_cube

            #daily_cube = iris.cube.Cube(fill_cube.data)
            #daily_cube.coords = fill_cube.coords
            #daily_cube.metadata = fill_cube.metadata
            #daily_cube.data = np.ma.zeros(daily_cube.data)

            dim_coords = fill_cube.coords(dim_coords=True)
            aux_coords = fill_cube.coords(dim_coords=False)
            dim_coords_and_dims = [(coord, adjust_dim(fill_cube.coord_dims(coord)[0]))
                                   for coord in dim_coords]
            aux_coords_and_dims = [(coord, (adjust_dim(d)
                                            for d in fill_cube.coord_dims(coord)))
                                   for coord in aux_coords]
            
            daily_cube = iris.cube.Cube(
                    fill_cube.data,
                    fill_cube.standard_name,
                    fill_cube.long_name,
                    fill_cube.var_name,
                    fill_cube.units,
                    fill_cube.attributes,
                    fill_cube.cell_methods,
                    dim_coords_and_dims,
                    aux_coords_and_dims,
            )
            
            # Read dataset time unit and calendar from file
            dataset_time_unit = str(fill_cube.coord('time').units)
            dataset_time_calender = fill_cube.coord('time').units.calendar
            # Convert datetime
            newtime = datetime(year=year, month=month, day=iday)
            newtime = cf_units.date2num(newtime, dataset_time_unit,
                                        dataset_time_calender)
            daily_cube.coord('time').points = float(newtime)

            #coord_time = fill_cube.coord('time')
            #time = coord_time.points
            #print(time)
            #times = datetime(year=year, month=month, day=iday)
            #print(times)
            #print(time.strftime("%Y%m%d"))
            #print(daily_cube.coord('time').points + relativedelta.relativedelta(days=1))
            #daily_cube.coord('time').points = float(time.strftime("%Y%m%d"))
            #daily_cube.coord('time').points = float(time_step) + 1.
            #daily_cube.coord('time').points += relativedelta.relativedelta(days=1)
            #new_time = time + float(iday) - 1.
            #daily_cube.coord('time') = iris.coords.DimCoord(time_units.date2num(times),
            #                                                var_name='time',
            #                                                standard_name='time',
            #                                                long_name='time',
            #                                                units=time_units)

            #daily_cube = fill_cube.filled(np.nan)
            #daily_cube.data = 0. #np.nan
            #print(fill_cube.coord('time').points)

        ## Fix coordinates
        utils.fix_coords(daily_cube)

        print(daily_cube.coord('time'))
        print(daily_cube.dim_coords)
        cube.append(daily_cube)

    print(cube)
    print([c.dim_coords[0] for c in cube])
    cube = cube.concatenate_cube()

    # regridding from 0.05x0.05 to 0.5x0.5
    cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    # Fix units
    if short_name == 'clt':
        daily_cube.data = 100. * daily_cube.data
    else:
        if 'raw_units' in var:
            cube.units = var['raw_units']
        cube.convert_units(cmor_info.units)

    utils.convert_timeunits(cube, 1950)

    ## Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata and  update version information
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        for year in range(glob_attrs['start_year'],
                          glob_attrs['end_year'] + 1):
            for month in range(1,13):
                print(month)
                _extract_variable(short_name, var, year, month, cfg, in_dir,
                                     out_dir)
