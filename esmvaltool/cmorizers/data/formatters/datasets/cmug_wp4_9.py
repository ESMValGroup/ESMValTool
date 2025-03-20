"""ESMValTool CMORizer for ESACCI-LST-UNCERT data.
Tier 2 # Is this the right tier????
Source
   on Jasmin:
   /group_workspaces/jasmin2/esacci_lst/public
Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-LST-UNCERT
   BOTH DAY and NIGHT files are needed for each month
Currently set to work with only the MODIS AQUA L3 monthly data
Modification history
   20201222 Started by Robert King, based on the CMUG WP5.3 cmorizer with no uncertanties
"""

import datetime
import logging
from calendar import monthrange

import iris
import cf_units as Unit

#from . import utilities as utils
from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)

height_coord = iris.coords.AuxCoord(2, var_name='height', units='metres')

#def cmorization(in_dir, out_dir, cfg, _):
def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""

    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    # This loop makes it easier for the CMUG WP5.4 work
    # variable_list contains the variable list
    # variable_keys has the short 'code' as a key for the variables.
    # both these lists are in made in teh same order
    # variable_list = []
    # variable_keys = []

    # VAR_LIST = {'day':[],
    #             'night':[]
    #              }

    # VAR_KEYS = {'day':[],
    #             'night':[]
    #              }
    # vals has the info from the yml file
    # var is set up in the yml file
    print(cfg['variables'].items())
    for var, vals in cfg['variables'].items():
        print('#####################')
        print(var)
        print(vals)
    
        glob_attrs['mip'] = vals['mip']

        # loop over years and months
        # get years from start_year and end_year
        output = iris.cube.CubeList()

        cubes = load_cubes(in_dir,
                                  vals['file'],
                                  1,#year,
                                  1,#month,
                                  vals['raw'],#variable_list
                                  
                ) 
        print(cubes)

        for i in range(len(cubes.coord('time').points)):
            print(i)

            date = cubes.coord('time')[i]
            print(date)

            datetime_point = Unit.num2date(date.points[0],
                                           str(date.units),
                                           date.units.calendar
                                           )
            year = datetime_point.year
            month = datetime_point.month
            print(datetime_point)
            print(year, month)
            this_cube = cubes[i:i+1,:,:]

            cubes.attributes = {}
            cubes.attributes['var'] = var

            # use CMORizer utils

            # this_cube.coord('longitude').bounds = None
            # this_cube.coord('longitude').points = this_cube.coord('longitude').points - 180

            # this_cube = iris.util.reverse(this_cube, 'latitude')


            print(this_cube.coord('latitude').var_name)
            print(this_cube.coord('longitude').var_name)

            this_cube.coord('latitude').var_name = 'lat'
            this_cube.coord('longitude').var_name = 'lon'
            try:
                this_cube = fix_coords(this_cube)
                print('FIXED')
            except:
                print('NOT NOT NOT FIXED')
                logger.info('skip fixing')
                logger.info(this_cube.long_name)

            
            # BOunds seemto cause problems later!!
            this_cube.coord('latitude').bounds = None
            this_cube.coord('longitude').bounds = None
            
            if var == 'ta':
                this_cube.add_aux_coord(height_coord)
                save_var='tas'
            else:
                save_var = var

            
            iris.save(this_cube,
                      '%s/OBS_CMUG_WP4_9_sat_1.00_Amon_%s_%s.nc' % 
                      (out_dir, save_var, '%s%02d'%(year,month))
            )
            print('############### SAVED %s' % var)

# #             # Not sure how to use the util saver with cube lists
# #             # utils.save_variable(
# #             #     monthly_cubes,
# #             #     var,
# #             #     out_dir,
# #             #     glob_attrs,
# #             # )


# no need to pass variable into any more, use global variable_list
def load_cubes(in_dir, file, year, month, variable_list):
    """
    variable = land surface temperature
    platform = AQUA not used for now
               but in place for future expansion to all ESC CCI LST plaforms
    """
    logger.info('Loading %s/%s%s%s*.nc', in_dir, file, year, month)

    cube = iris.load_cube('%s/%s*.nc' % (in_dir, file),
                          variable_list)


    return cube

