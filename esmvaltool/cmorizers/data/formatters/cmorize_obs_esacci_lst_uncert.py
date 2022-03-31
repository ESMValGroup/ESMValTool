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
from utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)

def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""

    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    # This loop makes it easier for the CMUG WP5.4 work
    # variable_list contains the variable list
    # variable_keys has the short 'code' as a key for the variables.
    # both these lists are in made in teh same order
    variable_list = []
    variable_keys = []
    # vals has the info from the yml file
    # var is set up in the yml file

    for var, vals in cfg['variables'].items():
    
        glob_attrs['mip'] = vals['mip']
        
        for key in vals.keys():
            logger.info("%s %s", key, vals[key])

        variable_list.append(vals['raw'])
        variable_keys.append(var)

    # loop over years and months
    # get years from start_year and end_year
    for year in range(vals['start_year'], vals['end_year'] + 1):
        print('****************************')
        print(variable_list)
        this_years_cubes = iris.cube.CubeList()
        for month in range(1,13):
            logger.info(year)
            logger.info(month)
            print(vals['file'])
            cubes = load_cubes(in_dir,
                              vals['file'],
                              '',
                              year,
                              month,
                              variable_list
            ) 
            print(cubes)
            #day_cube_lst   = day_cube.extract('land surface temperature')
            #night_cube_lst = night_cube.extract('land surface temperature')

            #monthly_cubes = make_monthly_average(day_cube_lst, night_cube_lst,
            #                                    year, month)

            output = iris.cube.CubeList()
            # cube list should be in same order as the variable_keys list

            # make time coords
            time_point = Unit.date2num(datetime.datetime(year,month,1),
                                       'hours since 1970-01-01 00:00:00',
                                       Unit.CALENDAR_STANDARD)
            #day_point = time_point + 12
            #night_point = time_point
            #all_point = time_point + 18
          
            
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

            
            for i,cube in enumerate(cubes):
                #if cube.long_name == 'land surface temperature': continue
                cube.attributes = {}
                cube.attributes['var'] = variable_keys[i]
                
                try:
                    cube.remove_coord('time')
                except:
                    logger.info('Coord fix issue %s' % cube.long_name)

                cube.add_dim_coord(time_coord, 0)
                    
                output.append(cube)

            output = output.merge()
            output = output.concatenate()

            # use CMORizer utils
            ########################## make this a loop over cube list

            for cube in output:
                try:
                    cube = utils.fix_coords(cube)
                except:
                    logger.info('skip fixing')
                    logger.info(cube.long_name)
                    pass
            print('Get here $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        # Use utils save
        # This seems to save files all with the same name!
            logger.info(out_dir)
            for cube in output:
                #print(cube)
                var_name = cube.attributes['var']
                print('************************ %s' % cube.long_name)
                print(cube)
                if var_name == 'tsLSSysErrDay' or var_name == 'tsLSSysErrNight':
                    print('BOOO')
                else:
                    print('HELLOE!')
                    cube.coords()[2].standard_name = 'longitude'
                iris.save(cube,
                          '%s/OBS_ESACCI_LST_UNCERTS_sat_1.00_Amon_%s_%s.nc' % 
                          (out_dir, var_name, '%s%02d'%(year,month))
                          )

        # Not sure how to use the util saver with cube lists
        # utils.save_variable(
        #     monthly_cubes,
        #     var,
        #     out_dir,
        #     glob_attrs,
        # )

# no need to pass variable into any more, use global variable_list
def load_cubes(in_dir, file, file_night, year, month, variable_list):
    """
    variable = land surface temperature
    platform = AQUA not used for now
               but in place for future expansion to all ESC CCI LST plaforms
    """
    logger.info('Loading %s/%s%s%s*.nc', in_dir, file, year, month)
    cube = iris.load('%s/%s%s%02d*.nc' % (in_dir, file,
                                                   year, month),
                         variable_list)

    # logger.info('Loading %s/%s%s%s*.nc', in_dir, file_night, year, month)
    # night_cube = iris.load('%s/%s%s%02d*.nc' % (in_dir, file_night,
    #                                                  year, month),
    #                        variable_list)

    return cube#, night_cube


def make_monthly_average(day_cube, night_cube, year, month):
    """Make the average LST form the day time and night time files."""
    for cube in day_cube:
        cube.attributes.clear()
        try:
            cube.coords()[2].var_name = 'longitude'
            cube.coords()[2].standard_name = 'longitude'
            cube.coords()[2].long_name = 'longitude'
        except:
            pass

    for cube in night_cube:
        cube.attributes.clear()
        try:
            cube.coords()[2].var_name = 'longitude'
            cube.coords()[2].standard_name = 'longitude'
            cube.coords()[2].long_name = 'longitude'
        except:
            pass

    # This is an arbitary change to one cube's time point
    # to allow the concatenation of two cubes for the 'same' time
    co_time = night_cube[0].coord('time') # zeroth cube should be the first one because of list
    co_time.points = co_time.points + 100.0
 
    # make a monthly LST value
    mean_lst_cubes = iris.cube.CubeList([day_cube[0], night_cube[0]]).concatenate_cube()
    print(mean_lst_cubes)
    monthly_lst = mean_lst_cubes.collapsed('time', iris.analysis.MEAN)

    monthly_lst.attributes = {}
    monthly_lst.attributes = {'information':
                              'Mean of Day and Night Aqua MODIS monthly LST'
                              }
    monthly_lst.long_name = 'land surface temperature'
    monthly_lst.attributes['var']  = 'ts'
    monthly_lst.cell_methods = None

    return monthly_lst
