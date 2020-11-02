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
   20201015 Started by Robert King
   20201029 Day/Night averaging added along with CMOR utils
"""

import glob
import os

import datetime
from calendar import monthrange

import iris

import logging
from . import utilities as utils

logger = logging.getLogger(__name__) # from OC example, dont know what this does!

# this function started from OC example
# This function is REQUIRED and it calls the other functions I write
def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes'] # this is needed later, not sure how it works though
    
    # run the cmorization

    # vals has the info from the yml file
    # var is set up in the yml file
    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in future for getting other info
        # like uncertainty information from the original files
        
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        
        for KEY in vals.keys():logger.info("%s %s" % (KEY,vals[KEY]))
        
        variable = vals['raw']
        platform = 'MODISA' # not currently used

        # loop over years and months
        # get years from start_year and end_year
        # not 2003 doesn't start until July
        for YEAR in range(2004,2019): # Change this in final version
            this_years_cubes = iris.cube.CubeList()
            for month in range(12): # Change this in final version
                MONTH = month + 1
                logger.info(MONTH)
                day_cube, night_cube = load_cubes(in_dir,
                                                  vals['file_day'],
                                                  vals['file_night'],
                                                  YEAR,
                                                  MONTH,
                                                  variable,
                                                  platform
                )

                monthly_cube = make_monthly_average(day_cube, night_cube, YEAR, MONTH)
                # day_cube and night_cube not need anymore this loop
                
                # use CMORizer utils
                monthly_cube = utils.fix_coords(monthly_cube)
                
                this_years_cubes.append(monthly_cube)
                
                # is there anything else needed for CMOR?????
                
            # Use utils save
            # This seems to save files all with the same name!! Need to fix!!!
            this_years_cubes = this_years_cubes.merge_cube()
            utils.save_variable(
                this_years_cubes, # monthly_cube,
                var,
                out_dir,
                glob_attrs,
            )

            
def load_cubes(in_dir,file_day,file_night,YEAR,MONTH,variable,platform):
    """
    variable = land surface temperature
    platform = AQUA not used for now but in place for future expansion to all ESC CCI LST plaforms
    """
    logger.info('%s/%s%s%02d*.nc' % (in_dir,file_day,YEAR,MONTH))
    day_cube = iris.load_cube('%s/%s%s%02d*.nc' % (in_dir,file_day,YEAR,MONTH),
                              variable)
    night_cube = iris.load_cube('%s/%s%s%02d*.nc' % (in_dir,file_night,YEAR,MONTH),
                                variable)

    return day_cube, night_cube


def make_monthly_average(day_cube, night_cube, YEAR, MONTH):
    """
    Make the average LST form the day time and night time files
    uses techniques from Patrick at AVD email 9/10/20
    """

    #for cube in day_cube + night_cube:
    day_cube.attributes.clear()
    night_cube.attributes.clear()

    #for cube in night_cube:
    co_time = night_cube.coord('time')
    co_time.points = co_time.points + 100.0
    # maybe the arbitary difference should go on day cubes to take the timestamp to 12Z?
    # not really an issue when using monthly files
 
    result = iris.cube.CubeList([day_cube,night_cube]).concatenate_cube()

    # This corrects the lonitude coord name issue
    logger.info("Longitude coordinate correction being applied")
    result.coords()[2].var_name = 'longitude'
    result.coords()[2].standard_name = 'longitude'
    result.coords()[2].long_name = 'longitude'
    
    monthly_cube = result.collapsed('time',iris.analysis.MEAN)

    
    # fix time coordinate bounds
    monthly_co_time = monthly_cube.coord('time')
    
    time_point = (datetime.datetime(YEAR,MONTH,1,0,0) - datetime.datetime(1981,1,1,0,0,0)).total_seconds()
    monthly_co_time.points = time_point

    num_days = monthrange(YEAR, MONTH)[1]
    monthly_co_time.bounds = [time_point, time_point+((num_days-1)*24*60*60)]
    # should this be num_days or num_days-1 ### question for Valeriu or Axel
    # or 23:59:59 ???

    monthly_cube.attributes = {'information' :'Mean of Day and Night Aqua MODIS monthly LST'}


    return monthly_cube
