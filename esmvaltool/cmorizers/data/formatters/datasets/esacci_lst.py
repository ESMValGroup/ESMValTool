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
import calendar
import numpy as np

from ...utilities import fix_coords
from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

# from CCI SNOW CMORISER
def _create_nan_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = np.full_like(nan_cube.data, np.nan, dtype=np.float32)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord('time').units)
    dataset_time_calender = nan_cube.coord('time').units.calendar

    # Convert datetime
    newtime = datetime.datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(newtime, dataset_time_unit,
                                dataset_time_calender)
    
    nan_cube.coord('time').points = np.float32(newtime)
      
    return nan_cube

def load_callback(cube, field, filename):

    cube.attributes = None

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    # vals has the info from the yml file
    # var is set up in the yml file
    for var, vals in cfg['variables'].items():
        glob_attrs['mip'] = vals['mip']
        # cmor_table = cfg["cmor_table"]

        variable_name = var.split('_')[0]

        logger.info(f"Starting CMORiser for {var}")

        if vals['ir_or_mw'] == 'mw':
            file1 = 'ASC'
            file2 = 'DES' 
        else:
            file1 = 'DAY'
            file2 = 'NIGHT'

        for ir_mw in [file1, file2]:
            start_date = datetime.datetime(vals['start_year'], 1, 1)
            end_date = datetime.datetime(vals['end_year'], 12, 31)
            current_date = start_date

            while current_date <= end_date:

                overpass_time = vals[f'{ir_mw.lower()}_overpass_time']

                month_cube_list = iris.cube.CubeList([])
                # put all cubes in here then concat
                for day in range(1, calendar.monthrange(current_date.year, current_date.month)[1]+1):  

                    datestr = datetime.datetime.strftime(current_date,
                                                        "%Y%m%d")
                    
                    file_pattern_1 = f"{in_dir}/{vals['file_base']}-{vals['platform']}-" + \
                    f"{vals['spatial_resolution']}_{vals['temporal_resolution']}" + \
                        f"_{ir_mw}-{datestr}000000-*{cfg['attributes']['version']}.nc"

                    logger.info(f"{file_pattern_1=}")
                    all_cubes = iris.cube.CubeList([])
                    for filename in [file_pattern_1]: #[file_pattern_1, file_pattern_2]:
                        logger.info(f"{filename}")
                        try:
                            cube = iris.load_cube(filename,
                                            vals['raw'],
                                            callback = load_callback)
                        except: # chek what the exception is
                            logger.info(f"{current_date=}  MISSING")
                            if len(month_cube_list) > 0:
                                cube = _create_nan_cube(month_cube_list[0], current_date.year,
                                                    current_date.month, current_date.day)
                            else:
                                # no first of the month to template from
                                logger.info("NO 1st of the month template")
                                continue
                            
                        cube = fix_coords(cube)

                        month_cube_list.append(cube)
                    
                    current_date = current_date + datetime.timedelta(days=1)

                try:
                    all_cubes =  month_cube_list.concatenate_cube()
                except ValueError:
                    # happens when no files in a month
                    logger.info("No files for this month")
                    continue

                # this is needed for V1 data, V3 data is ok
                # try:
                #     cube.coords()[2].standard_name = 'longitude'
                #     logger.info('V1 longitude issue fixed')
                # except IndexError:
                #     # No change needed
                #     logger.info('No V1 longitude issue')

                # this gives time as hours since 19500101
                all_cubes = fix_coords(all_cubes)

                # this is need for when uncertainity variables added
                if all_cubes.long_name == 'land surface temperature':
                    all_cubes.long_name = 'surface_temperature'
                    all_cubes.standard_name = 'surface_temperature'

                if all_cubes.var_name == 'lst':
                    all_cubes.var_name = 'ts'

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

                cmor_info = cfg["cmor_table"].get_variable(vals["mip"], vals['raw'])
                #utils.fix_var_metadata(cube, cmor_info)

                # Fix global metadata
                glob_attrs['platform'] = vals['platform']

                utils.set_global_atts(all_cubes, glob_attrs)

                old_version = glob_attrs["dataset_id"]
                glob_attrs["dataset_id"] = glob_attrs["dataset_id"] + f"_{glob_attrs['platform']}" + f"_{ir_mw}"

                all_cubes.attributes['overpass_time'] = overpass_time

                utils.save_variable(
                    all_cubes,
                    variable_name,
                    out_dir,
                    glob_attrs,
                    unlimited_dimensions=["time"],
                    zlib=True
                )

                # need to put this back
                glob_attrs["dataset_id"] = old_version
