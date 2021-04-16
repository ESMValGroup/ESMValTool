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

from . import utilities as utils

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
        # not currently used, but referenced for future
        # platform = 'MODISA'

    # loop over years and months
    # get years from start_year and end_year
    for year in range(vals['start_year'], vals['end_year'] + 1):
        
        this_years_cubes = iris.cube.CubeList()
        for month in range(1,3):
            logger.info(year)
            logger.info(month)

            day_cube, night_cube = load_cubes(in_dir,
                                              vals['file_day'],
                                              vals['file_night'],
                                              year,
                                              month,
                                              variable_list
                                              ) 

            day_cube_lst   = day_cube.extract('land surface temperature')
            night_cube_lst = night_cube.extract('land surface temperature')

            monthly_cubes = make_monthly_average(day_cube_lst, night_cube_lst,
                                                year, month)

            output = iris.cube.CubeList()
            # cube list should be in same order as the variable_keys list

            time_point = (datetime.datetime(year, month, 1, 0, 0) -
                          datetime.datetime(1970, 1, 1, 0, 0, 0)).total_seconds()
            print('QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ')
            print(time_point)

            print(time_point/24/3600)
            print(datetime.datetime(1970, 1, 1, 0, 0, 0) + datetime.timedelta(hours=time_point/24/3600))
            # time_point is the 1200 point some how ?!
            day_point = time_point 
            night_point = time_point - 12
            all_point = time_point + 8
          
            # THINK ABOUT BOUNDS ON TIME COORDS!!!
            num_days = monthrange(year, month)[1]
            bounds = [time_point-12,
                      time_point + ((num_days - 1) * 24)+11]

            day_time_coord = iris.coords.DimCoord(day_point, 
                                                  standard_name='time',
                                                  long_name='time', 
                                                  var_name='time', 
                                                  units='hours since 1970-01-01',
                                                  bounds=None,#bounds, 
                                                  attributes=None, 
                                                  coord_system=None, 
                                                  circular=False
            )

            night_time_coord = iris.coords.DimCoord(night_point, 
                                                    standard_name='time',
                                                    long_name='time', 
                                                    var_name='time', 
                                                    units='hours since 1970-01-01',
                                                    bounds=None,#bounds, 
                                                    attributes=None, 
                                                    coord_system=None, 
                                                    circular=False
            )

            all_time_coord = iris.coords.DimCoord(all_point, 
                                                  standard_name='time',
                                                  long_name='time', 
                                                  var_name='time', 
                                                  units='hours since 1970-01-01',
                                                  bounds=bounds, 
                                                  attributes=None, 
                                                  coord_system=None, 
                                                  circular=False
            )

            
            # DAY TIME VALUES given a 1200 time
            # NIGHT TIME VALUES given a 0000 time
            # ALL TIME VALUES given a 1800 time
            # all on the 1st of the month
            monthly_cubes.remove_coord('time')
            monthly_cubes.add_aux_coord(all_time_coord)
            output.append(monthly_cubes)

            for i,cube in enumerate(day_cube):
                cube.attributes = {}
                #cube.attributes['time'] = 'day'
                cube.attributes['var'] = variable_keys[i]
                
                try:
                    cube.remove_coord('time')
                except:
                    print(cube.long_name)

                cube.add_dim_coord(day_time_coord, 0)
                    
                output.append(cube)

            for i,cube in enumerate(night_cube):
                cube.attributes = {}
                #cube.attributes['time'] = 'night'
                cube.attributes['var'] = variable_keys[i]
                
                try:
                    cube.remove_coord('time')
                except:
                    print(cube.long_name)

                cube.add_dim_coord(night_time_coord, 0)
                    
                output.append(cube)

            output = output.merge()
            output = output.concatenate()
            print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
            X = output.extract('uncertainty from large-scale systematic errors')
            print(X[0].coord('time'))
            print(X[1].coord('time'))

            Y = X.merge_cube()
            print(Y)
            


            # use CMORizer utils
            ########################## make this a loop over cube list
            for cube in output:
                try:
                    cube = utils.fix_coords(cube)
                except:
                    logger.info('skip fixing')
                    logger.info(cube.long_name)
                    pass

            print(output)
            print(0/0)
        # Use utils save
        # This seems to save files all with the same name!
        # Fixed by making yearly files

        # need to do this only for the LST variable
        # this_years_cubes = this_years_cubes.merge_cube()
        # this_years_cubes.long_name = 'Surface Temperature'
        # this_years_cubes.standard_name = 'surface_temperature'
            logger.info(out_dir)
            for cube in output:
                var_name = cube.attributes['var']
                iris.save(cube,
                          '%s/OBS_ESACCI-LST-UNCERTS_sat_1.00_Amon_%s_%s.nc' % 
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
def load_cubes(in_dir, file_day, file_night, year, month, variable_list):
    """
    variable = land surface temperature
    platform = AQUA not used for now
               but in place for future expansion to all ESC CCI LST plaforms
    """
    logger.info('Loading %s/%s%s%s*.nc', in_dir, file_day, year, month)
    day_cube = iris.load('%s/%s%s%02d*.nc' % (in_dir, file_day,
                                                   year, month),
                         variable_list)

    logger.info('Loading %s/%s%s%s*.nc', in_dir, file_night, year, month)
    night_cube = iris.load('%s/%s%s%02d*.nc' % (in_dir, file_night,
                                                     year, month),
                           variable_list)

    return day_cube, night_cube


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

    co_time = night_cube[0].coord('time') # zeroth cube should be the first one because of list
    co_time.points = co_time.points + 100.0
    # maybe the arbitary difference should go on day cubes to
    # take the timestamp to 12Z?
    # not really an issue when using monthly files

    # Make a time suitable time coordinate
    # make a monthly LST value, need to do this here to get coordinate values
    mean_lst_cubes = iris.cube.CubeList([day_cube[0], night_cube[0]]).concatenate_cube()
    monthly_lst = mean_lst_cubes.collapsed('time', iris.analysis.MEAN)

    # fix time coordinate bounds
    # monthly_co_time = monthly_lst.coord('time')

    # time_point = (datetime.datetime(year, month, 1, 0, 0) -
    #               datetime.datetime(1981, 1, 1, 0, 0, 0)).total_seconds()
    # monthly_co_time.points = time_point

    # num_days = monthrange(year, month)[1]
    # monthly_co_time.bounds = [time_point,
    #                           time_point + ((num_days - 1) * 24 * 3600)]
    # # should this be num_days or num_days-1 ### question for Valeriu or Axel
    # # or 23:59:59 ???
    monthly_lst.attributes = {}
    monthly_lst.attributes = {'information':
                               'Mean of Day and Night Aqua MODIS monthly LST'
                               }
    monthly_lst.long_name = 'surface temperature'

    #monthly_lst.attributes['time'] = 'all'
    monthly_lst.attributes['var']  = 'ts'

    return monthly_lst
