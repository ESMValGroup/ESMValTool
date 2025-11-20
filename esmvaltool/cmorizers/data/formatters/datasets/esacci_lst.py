"""ESMValTool CMORizer for ESACCI-LST.

Modified from the CCI LST V3 CMORizer
for V5.11 Microwave sensors
Compatibility with V3 IR sensors kept
Should be exapandable to any single sensor CCI LST product, daily or monthly
"""

import datetime
import logging
import iris
import cf_units
import glob

from ...utilities import fix_coords
from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

# decimal hours
# ECT taken from https://space.oscar.wmo.int/satellites
# times with a comment are what is on the Oscar website (HHMM)
# times without a comment are the 12 hour differences
overpass_time = {'SSMI13': {
                    'ASC': 17.85,
                    'DES': 5.85, # 0551
                    },
                'SSMI17': {
                    'ASC': 18.58,
                    'DES': 6.58, # 0635
                    },
                'AMSR_E': {
                    'ASC': 2.30,
                    'DES': 14.30, # 1418
                    },
                'AMSR_2': {
                    'ASC': 1.50,
                    'DES': 13.50, # 1330
                    },
                'MODISA': {
                    'DAY': 14.30, # 1418
                    'NIGHT': 2.30,
                    },
                'MODIST': {
                    'DAY': 9.67, # 0940
                    'NIGHT': 23.67,
                }
            }

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    # vals has the info from the yml file
    # var is set up in the yml file
    logger.info(f"{cfg=}")
    logger.info(f"{glob_attrs=}")

    for var, vals in cfg['variables'].items():
        glob_attrs['mip'] = vals['mip']
        cmor_table = cfg["cmor_table"]

        all_files = glob.glob(f"{in_dir}/{vals['file_base']}-{vals['platform']}-{vals['spatial_resolution']}_{vals['temporal_resolution']}_*-*-*{cfg['attributes']['version']}.nc")
          
        # loop over years then files
        for year in range(vals['start_year'], vals['end_year'] + 1):
            for file in all_files:
                if f"-{year}" not in file:
                    continue

                logger.info(f"Loading {file}")
                cube = iris.load_cube(file, vals['raw'])
                
                cube.attributes = {}
                cube.attributes['var'] = var

                cmor_info = cmor_table.get_variable(vals["mip"], var)
                
                # this is needed for V1 data, V3 data is ok
                try:
                    cube.coords()[2].standard_name = 'longitude'
                    logger.info('V1 longitude issue fixed')
                except IndexError:
                    # No change needed
                    logger.info('No V1 longitude issue')

                # this gives time as hours since 19500101
                cube = fix_coords(cube)
                
                # this is need for when uncertainity variables added
                if cube.long_name == 'land surface temperature':
                    cube.long_name = 'surface_temperature'
                    cube.standard_name = 'surface_temperature'

                if cube.var_name == 'lst':
                    cube.var_name = 'ts'

                orbit_time = ''
                if 'DAY' in file: orbit_time = 'DAY'
                if 'NIGHT' in file: orbit_time = 'NIGHT'
                if 'ASC' in file: orbit_time = 'ASC'
                if 'DES' in file: orbit_time = 'DES'
                if orbit_time == '':
                    # can't work out which overpass time to use
                    # default to 0.0 for midnight
                    time_of_overpass = 0.0
                else:
                    time_of_overpass = overpass_time[vals['platform']][orbit_time]
                
                # make a new time coordinate with new times
                new_time_points = cube.coord('time').points + time_of_overpass
                new_time_coord = iris.coords.DimCoord(new_time_points,
                                                      standard_name='time',
                                                      long_name='time',
                                                      var_name='time',
                                                      units=cube.coord('time').units,
                                                      bounds=None,
                                                      attributes=None,
                                                      coord_system=None,
                                                      circular=False
                                                      )
                
                try:
                    cube.remove_coord('time')
                    cube.add_dim_coord(new_time_coord, 0)
                    logger.info('New time coord added')
                except iris.exceptions.CoordinateNotFoundError:
                    logger.info('Problem adding overpass time to time coord')
        
#  Leaving this here for when ESMValCore PR with new CMOR tables for uncertainity are avaible
#                 # Land cover class gives this error when
#                 # loading in CMORised files
#                 # OverflowError:
#                 # Python int too large to convert to C long error
#                 # This fixes it:
#                 if 'land cover' in cubes.long_name:
#                     cubes.data.fill_value = 0

#                     # Get rid of anything outide 0-255
#                     cubes.data[np.where(cubes.data < 0)] = 0
#                     cubes.data[np.where(cubes.data > 255)] = 0

#                     cubes.data = cubes.data.filled()

#                     # This is the line the ultimately solves things.
#                     # Will leave the other checks above in because they
#                     cubes.data = cubes.data * 1.0

                utils.fix_var_metadata(cube, cmor_info)

                # Fix global metadata
                glob_attrs['platform'] = vals['platform']
                glob_attrs['orbit_direction'] = orbit_time

                utils.set_global_atts(cube, glob_attrs)

                # This util funtion will give all files in the same year the same name
                # can not see a way to resolve this from the source code
                # esmvaltool/cmorizers/data/utilities.py
                # utils.save_variable(
                #     cube,
                #     var_name,
                #     out_dir,
                #     glob_attrs,
                #     unlimited_dimensions=["time"],
                #    )

                # Alternative saving function
                time_dt = cf_units.num2pydate(cube.coord('time').points[0],
                                            str(cube.coord('time').units),
                                            cube.coord('time').units.calendar
                                            )
                datestr = datetime.datetime.strftime(time_dt,'%Y%m%d%H%M00')

                
                save_name = f"{out_dir}/OBS_ESACCI-LST_sat_{cfg['attributes']['version']}_{vals['mip']}_" + \
                    f"{var}_{datestr}.nc"
                iris.save(cube, save_name)