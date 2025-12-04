"""ESMValTool CMORizer for ESACCI-LST.

Modified from the CCI LST V3 CMORizer
for V5.11 Microwave sensors
Compatibility with V3 IR sensors kept
Should be exapandable to any single sensor CCI LST product, daily or monthly
"""

import datetime
import logging
import glob
import iris
import cf_units
import calendar

from ...utilities import fix_coords
from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

# decimal hours
# ECT taken from https://space.oscar.wmo.int/satellites
# times with a comment are what is on the Oscar website (HHMM)
# times without a comment are the 12 hour differences
overpass_time = {'SSMI13': {
                    'ASC': 17.85,
                    'DES': 5.85,  # 0551
                    },
                'SSMI17': {
                    'ASC': 18.58,
                    'DES': 6.58,  # 0635
                    },
                'AMSR_E': {
                    'ASC': 2.30,
                    'DES': 14.30,  # 1418
                    },
                'AMSR_2': {
                    'ASC': 1.50,
                    'DES': 13.50,  # 1330
                    },
                'MODISA': {
                    'DAY': 14.30,  # 1418
                    'NIGHT': 2.30,
                    },
                'MODIST': {
                    'DAY': 9.67,  # 0940
                    'NIGHT': 23.67,
                }
}

def load_callback(cube, field, filename):

    cube.attributes = None

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    # vals has the info from the yml file
    # var is set up in the yml file
    for var, vals in cfg['variables'].items():
        glob_attrs['mip'] = vals['mip']
        cmor_table = cfg["cmor_table"]

        variable_name = var.split('_')[0]

        logger.info(f"Starting CMORiser for {var}")

        start_date = datetime.datetime(vals['start_year'], 1, 1)
        end_date = datetime.datetime(vals['end_year'], 12, 31)

        current_date = start_date

        if vals['ir_or_mw'] == 'mw':
            file1 = 'ASC'
            file2 = 'DES'
        else:
            file1 = 'DAY'
            file2 = 'NIGHT'


        this_month = 1
        while current_date <= end_date:

            month_cube_list = iris.cube.CubeList([])
            # put all cubes in here then concat
            for day in range(1, calendar.monthrange(current_date.year, current_date.month)[1]+1):
                

                datestr = datetime.datetime.strftime(current_date,
                                                    "%Y%m%d")
                
                ##### ADD A MONTH LOOP to make a file per month ######

                # ESACCI-LST-L3C-LST-SSMI17-0.125deg_1DAILY_DES-20171214000000-fv5.11.nc
                file_pattern_1 = f"{in_dir}/{vals['file_base']}-{vals['platform']}-" + \
                f"{vals['spatial_resolution']}_{vals['temporal_resolution']}" + \
                    f"_{file1}-{datestr}000000-*{cfg['attributes']['version']}.nc"
            
                file_pattern_2 = f"{in_dir}/{vals['file_base']}-{vals['platform']}-" + \
                f"{vals['spatial_resolution']}_{vals['temporal_resolution']}" + \
                    f"_{file2}-{datestr}000000-*{cfg['attributes']['version']}.nc"
                
                #all_files = glob.glob(file_pattern)

                all_cubes = iris.cube.CubeList([])
                for filename in [file_pattern_1, file_pattern_2]:
                    logger.info(f"{filename}")
                    cube = iris.load_cube(filename,
                                    vals['raw'],
                                    callback = load_callback)
                    logger.info(f"{cube.coord('time')=}")
                    cube = fix_coords(cube)
                    logger.info(f"{cube.coord('time')=}")
                    logger.info(f"{cube.coord('time').points=}")


                    for KEY in overpass_time.keys():
                        new_time = 0.0
                        if KEY in filename:
                            if 'ASC' in filename:
                                new_time = cube.coord('time').points + overpass_time[KEY]['ASC']/24.
                            if 'DES' in filename:
                                new_time = cube.coord('time').points + overpass_time[KEY]['DES']/24.

                            #new_time_points = cube.coord('time').points + new_time/24.
                            new_time_coord = iris.coords.DimCoord(new_time,
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

                            month_cube_list.append(cube)


                    


                    
                
                current_date = current_date + datetime.timedelta(days=1)

            logger.info(f"{month_cube_list=}")
            all_cubes =  month_cube_list.concatenate_cube()
            logger.info(f"{all_cubes=}")


        
    
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
            #glob_attrs['orbit_direction'] = orbit_time

            utils.set_global_atts(all_cubes, glob_attrs)

            # This util funtion will give all files in the same year the same name
            # can not see a way to resolve this from the source code
            # esmvaltool/cmorizers/data/utilities.py
            utils.save_variable(
                all_cubes,
                variable_name,
                out_dir,
                glob_attrs,
                unlimited_dimensions=["time"],
               )

            ######## work out how to get utils save to work with the multiple platforms
            # save_name = f"{out_dir}/OBS_ESACCI-LST_sat_{cfg['attributes']['version']}_{vals['mip']}_" + \
            #     f"{var}_{datestr}.nc"
            # iris.save(all_cubes, save_name)
           