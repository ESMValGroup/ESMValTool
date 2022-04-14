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
        for year in range(vals['start_year'], vals['end_year'] + 1):
            print('*******')
            print(year)

            this_years_cubes = iris.cube.CubeList()
            for month in range(1,13):
                logger.info(year)
                logger.info(month)
                print(vals['file'])
                cubes = load_cubes(in_dir,
                                  vals['file'],
                                  year,
                                  month,
                                  vals['raw'],#variable_list
                                  
                ) 
                                
                #  make time coords
                time_point = Unit.date2num(datetime.datetime(year,month,1),
                                           'hours since 1970-01-01 00:00:00',
                                           Unit.CALENDAR_STANDARD)

                time_coord = iris.coords.DimCoord(time_point, 
                                                  standard_name='time',
                                                  long_name='time', 
                                                  var_name='time', 
                                                  units='hours since 1970-01-01',
                                                  bounds=None,#bounds, 
                                                  attributes=None, 
                                                  coord_system=None, 
                                                  circular=False
                                              )


# #                 for i,cube in enumerate(cubes):
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

                 # use CMORizer utils

# #                 for cube in output:
                try:
                    cubes = fix_coords(cubes)
                    print('FIXED')
                except:
                    print('NOT NOT NOT FIXED')
                    logger.info('skip fixing')
                    logger.info(cubes.long_name)

                # Using  utils save :This seems to save files all with the same name!
                # so need to write this

                try:
                    cubes.coords()[2].standard_name = 'longitude'
                    print('YES LST LON CHNAGE')
                except:
                    print('NO LAT LON CHANGE')

                var_name = cubes.attributes['var']

                if cubes.var_name == 'lst':
                    cubes.var_name = 'ts'

                iris.save(cubes,
                          '%s/OBS_ESACCI_LST_UNCERTS_sat_1.00_Amon_%s_%s.nc' % 
                          (out_dir, var_name, '%s%02d'%(year,month))
                )
                print('############### SAVED %s' % var_name)

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
    cube = iris.load_cube('%s/%s%s%02d*.nc' % (in_dir, file,
                                          year, month),
                     variable_list)

    return cube

