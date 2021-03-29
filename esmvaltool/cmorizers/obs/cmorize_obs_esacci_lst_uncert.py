"""ESMValTool CMORizer for ESACCI-LST data.
Tier 2 # Is this the right tier????
Source
   /group_workspaces/jasmin2/esacci_lst/public
Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-LST
   BOTH DAY and NIGHT files are needed for each month
Currently set to work with only the MODIS AQUA L3 monthly data
Modification history
   20201222 Started by Robert King, based on the CMUG WP5.3 cmorizer with no uncertanties
"""

import datetime
import logging
from calendar import monthrange

import iris

from . import utilities as utils

logger = logging.getLogger(__name__)

# Need a list of the uncertainity information variable names to load
variable_list = ['land surface temperature',
                 'uncertainty from locally correlated errors on atmospheric scales',
                 # 'uncertainty from locally correlated errors on surface scales',
                 # 'uncertainty from uncorrelated errors',
                 # 'uncertainty from large-scale systematic errors',
                 # 'land surface temperature total uncertainty',
                 # 'land cover class'
]


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    # cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    print('1111111111111111111111111111')
    print(cmor_table)
    print('222222222222222222222222222222')
    # run the cmorization

    # The loop of variables like the WP5.3 isnt needed as
    # variable_list contains the variable list

    # # vals has the info from the yml file
    # # var is set up in the yml file
    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in
        # the future for getting other info
        # like uncertainty information from the original files
        print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        print(var,vals)
        print(cmor_table.get_variable(vals['mip'], var))
        print("QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ")
        glob_attrs['mip'] = vals['mip']

        for key in vals.keys():
            logger.info("%s %s", key, vals[key])

        variable = vals['raw']
        # not currently used, but referenced for future
        # platform = 'MODISA'

    # loop over years and months
    # get years from start_year and end_year
    # note 2003 doesn't start until July so not included at this stage
    for year in range(vals['start_year'], vals['end_year'] + 1):

        this_years_cubes = iris.cube.CubeList()
        for month in range(1,13):
            logger.info(month)
            logger.info(in_dir)
            logger.info(vals['file_day'])
            logger.info(vals['file_night'])
            logger.info(year)


            day_cube, night_cube = load_cubes(in_dir,
                                              vals['file_day'],
                                              vals['file_night'],
                                              year,
                                              month
                                              ) # variable moved
            logger.info("KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK")
            print(day_cube, night_cube)
            monthly_cubes = make_monthly_average(day_cube, night_cube,
                                                year, month)

            # use CMORizer utils
            ########################## make this a loop over cube list
            for cube in monthly_cubes:
                try:
                    cube = utils.fix_coords(cube)
                except:
                    logger.info('skip fixing')
                    logger.info(cube)
                    pass

        # Use utils save
        # This seems to save files all with the same name!
        # Fixed by making yearly files

        # need to do this only for the LST variable
        # this_years_cubes = this_years_cubes.merge_cube()
        # this_years_cubes.long_name = 'Surface Temperature'
        # this_years_cubes.standard_name = 'surface_temperature'
            logger.info(out_dir)
            iris.save(monthly_cubes, '%s/test_%s_%02d.nc' % (out_dir,year,month))

        # Not sure how to use the util saver with cube lists
        # utils.save_variable(
        #     monthly_cubes,
        #     var,
        #     out_dir,
        #     glob_attrs,
        # )

# no need to pass variable into any more, use global variable_list
def load_cubes(in_dir, file_day, file_night, year, month):
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
    monthly_co_time = monthly_lst.coord('time')

    time_point = (datetime.datetime(year, month, 1, 0, 0) -
                  datetime.datetime(1981, 1, 1, 0, 0, 0)).total_seconds()
    monthly_co_time.points = time_point

    num_days = monthrange(year, month)[1]
    monthly_co_time.bounds = [time_point,
                              time_point + ((num_days - 1) * 24 * 3600)]
    # should this be num_days or num_days-1 ### question for Valeriu or Axel
    # or 23:59:59 ???


    # put all cubes here
    output = iris.cube.CubeList()

    
    monthly_lst.attributes = {'information':
                               'Mean of Day and Night Aqua MODIS monthly LST'
                               }
    monthly_lst.long_name = 'surface temperature'

    output.append(monthly_lst)


    for cube in day_cube:
        cube.long_name = 'day ' + cube.long_name
        output.append(cube)

    for cube in night_cube:
        cube.long_name = 'night ' + cube.long_name
        output.append(cube)

    logger.info(output)

    return output
