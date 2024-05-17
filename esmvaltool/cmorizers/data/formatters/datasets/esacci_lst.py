"""ESMValTool CMORizer for ESACCI-LST.

Tested for the V3 of the data
Works for both day and night files
"""

import datetime
import logging
import iris
import cf_units as unit
import numpy as np

from ...utilities import fix_coords

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    # variable_list contains the variable list
    # variable_keys has the short 'code' as a key for the variables.
    # both these lists are in made in the same order

    # vals has the info from the yml file
    # var is set up in the yml file
    for var, vals in cfg['variables'].items():
        glob_attrs['mip'] = vals['mip']

        # loop over years and months
        # get years from start_year and end_year
        for year in range(vals['start_year'], vals['end_year'] + 1):
            for month in range(1, 13):
                logger.info(year)
                logger.info(month)
                try:
                    cubes = load_cubes(in_dir,
                                       vals['file'],
                                       year,
                                       month,
                                       vals['raw'],
                                   )
                except:
                    logger.info(f'Problem with Month {month} in {year}')
                    continue

                #  make time coords
                time_units = 'hours since 1970-01-01 00:00:00'
                time_point = unit.date2num(datetime.datetime(year, month, 1),
                                           time_units,
                                           unit.CALENDAR_STANDARD
                                           )

                time_coord = iris.coords.DimCoord(time_point,
                                                  standard_name='time',
                                                  long_name='time',
                                                  var_name='time',
                                                  units=time_units,
                                                  bounds=None,
                                                  attributes=None,
                                                  coord_system=None,
                                                  circular=False
                                                  )

                cubes.attributes = {}
                cubes.attributes['var'] = var

                try:
                    cubes.remove_coord('time')
                except:
                    logger.info('Coord fix issue %s' % cubes.long_name)

                cubes.add_dim_coord(time_coord, 0)

                if cubes.long_name == 'land surface temperature':
                    cubes.long_name = 'surface_temperature'
                    cubes.standard_name = 'surface_temperature'

                try:
                    cubes = fix_coords(cubes)
                except:
                    logger.info('skip fixing')
                    logger.info(cubes.long_name)

                # this is needed for V1 data, V3 data is ok
                try:
                    cubes.coords()[2].standard_name = 'longitude'
                except:
                    # No change needed
                    pass

                var_name = cubes.attributes['var']

                if cubes.var_name == 'lst':
                    cubes.var_name = 'ts'

                if 'Day' in var_name:
                    cubes.long_name += ' Day'
                    cubes.var_name += '_day'

                if 'Night' in var_name:
                    cubes.long_name += ' Night'
                    cubes.var_name += '_night'
                    
                # land cover class gives this error when loading in CMORised fiels
                # OverflowError: Python int too large to convert to C long error
                # Attempt to fix it
                if 'land cover' in cubes.long_name:
                    cubes.data.fill_value = 0
                    cubes.data = cubes.data.filled()
                    cubes.data = np.ma.masked_equal(cubes.data, 0)

                save_name = f'{out_dir}/OBS_ESACCI-LST_sat_3.00_Amon_{var_name}_{year}{month:02d}.nc'
                iris.save(cubes,
                          save_name
                          )


def load_cubes(in_dir, file, year, month, variable_list):
    """Load files into cubes based on variables wanted in variable_list."""
    logger.info(f'Loading {in_dir}/{file}{year}{month:02d}.nc')
    cube = iris.load_cube(f'{in_dir}/{file}{year}{month:02d}*.nc',
                          variable_list
                          )

    return cube
