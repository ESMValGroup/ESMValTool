"""ESMValTool CMORizer for ESA LAI
Source
   
Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/CMUG_WP4_10
"""

import datetime
import logging
from calendar import monthrange
import glob

import iris
import iris.coord_categorisation as icc
import cf_units as Unit

#from . import utilities as utils
from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)

#def cmorization(in_dir, out_dir, cfg, _):
def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""

    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']

    # vals has the info from the yml file
    # var is set up in the yml file
    print(cfg['variables'].items())
    for var, vals in cfg['variables'].items():
        print('#####################')
        print(var)
        print(vals)
    
        glob_attrs['mip'] = vals['mip']

        var_name = var
        # loop over years and months
        # get years from start_year and end_year
        for year in range(vals['start_year'], vals['end_year'] + 1):
            print('*******')
            print(year)

            for month in range(1,13):
                print(month)

                output = iris.cube.CubeList()
                print(vals['file'])
                cubes = load_cubes(in_dir,
                                   vals['file'],
                                   vals['raw'],#variable_list,
                                   year,
                                   month
                               ) 

                
                cubes = cubes.concatenate_cube()

                cube_mean = cubes.collapsed('time', iris.analysis.MEAN)
                cube_mean.attributes['num_files'] = len(cubes.coord('time').points)

                cube_mean.units = 1

                print(cube_mean)
                print(cubes.coord('time'))
                
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
             


                try:
                    cube_mean.remove_coord('time')
                except:
                    logger.info('Coord fix issue %s' % cube_mean.long_name)

                cube_mean.add_aux_coord(time_coord)
                cube_mean = iris.util.new_axis(cube_mean, 'time')
                print(cube_mean)

                # use CMORizer utils
                  
                try:
                    cubemean = fix_coords(cubes)
                    print('FIXED')
                except:
                    print('NOT NOT NOT FIXED')
                    logger.info('skip fixing')
                    logger.info(cube_mean.long_name)

                iris.save(cube_mean,
                          '%s/OBS_CMUG_WP4_10_sat_1.00_Lmon_%s_%s.nc' % 
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
def load_cubes(in_dir, file, variable_list, year,month):
    """
    Descrip
    """
    filelist = glob.glob('%s/%s*%s%02d*.nc' % (in_dir, file, year,month))
    
    loaded_cubes = iris.cube.CubeList()
    print(filelist)

    for file in filelist:
        try:
            cubes = iris.load(file,
                              variable_list)
            for cube in cubes:
                cube.attributes = None
                loaded_cubes.append(cube)
        except:
            pass




    return loaded_cubes
