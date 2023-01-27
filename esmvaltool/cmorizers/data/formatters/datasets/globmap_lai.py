"""ESMValTool CMORizer for GLOBMAP LAI data.
Tier
   Tier 2: other freely-available dataset.
Source
   https://zenodo.org/record/4700264
Download and processing instructions
   Original data files are hdf4 files that are download in rar packages.
"""

import datetime
import logging
from pyhdf.SD import SD, SDC
import glob
import numpy as np
import iris

from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""

    glob_attrs = cfg['attributes']

    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in the future for getting
        # other info like uncertainty information from the original files
        glob_attrs['mip'] = vals['mip']
        variable = vals['raw']
        # loop over years
        for year in range(vals['start_year'], vals['end_year'] + 1):

            loaded_cubes = load_and_make_cubes(in_dir,
                                               vals['file'],
                                               year,
                                               variable
                                               )

            # might be better to redefine how coordinates are made?
            # or just do this anyway to make sure nothing subtle is missed
            fixed_cubes = fix_coords(loaded_cubes)

            # save this year's data
            save_variable(fixed_cubes,
                          var,
                          out_dir,
                          glob_attrs
                          )
            

def load_and_make_cubes(in_dir, filepath, year, variable):
    """Load GLOBMAP HDF4 files.

    Also make cubes of the data and return a single merged cube
    of every file. 
    Note frequency is variable, because of a resolution change in the 
    original data, and files with errors/issues are ignored.
    """

    logger.info(f'Loading {in_dir}/{filepath}{year}*.hdf')

    filelist = glob.glob(f'{in_dir}/{filepath}{year}*.hdf')

    # CubeList for storing data from each file this year.
    cubes_output = iris.cube.CubeList()

    for filename in filelist:
        hdf_data = SD(filename, SDC.READ)
        print(filename)
        try:
            data_from_file  = hdf_data.select(variable)
        except: #HDF4Error:
            print('ERROR!')
            hdf_data.end()
            continue
        data_array = data_from_file[:,:]

        lon_points = np.array([hdf_data.attributes()['PROJ_UL_XY'][0] + i*hdf_data.attributes()['PIXEL_SIZE']
                               for i in range(data_array.shape[1])]
                             )
        lat_points = np.array([hdf_data.attributes()['PROJ_UL_XY'][1] - i*hdf_data.attributes()['PIXEL_SIZE']
                               for i in range(data_array.shape[0])]
                             )

        lat_coord = iris.coords.DimCoord(lat_points,
                                         standard_name='latitude',
                                         units='degrees'
                                         )
    
        lon_coord = iris.coords.DimCoord(lon_points,
                                         standard_name='longitude',
                                         units='degrees'
                                        )

        # Date comes form filename, easiest place to get it
        time_str = filename.split('.')[-3].split('A')[1]
        year = int(time_str[:4])
        days = int(time_str[4:])

        # use days since Jan 1st 1950 - this covers the whole LAI dataset even though
        # this fits CMOR fix???
        datum = datetime.datetime(1950,1,1)
        this_date = datetime.datetime(year,1,1) + datetime.timedelta(days=days)
        num_days = (this_date - datum).days
    
        time_coord = iris.coords.AuxCoord(num_days, 
                                     standard_name = 'time',
                                     units = 'days since 1950-1-1'
                                     )

        cube = iris.cube.Cube(data_array,
                          long_name = 'leaf area index',
                          dim_coords_and_dims = ((lat_coord,0),
                                                 (lon_coord,1)
                                                 )
                          )

        cube.add_aux_coord(time_coord)

        data_from_file.endaccess()
        hdf_data.end()

        cubes_output.append(cube)
    
    print(cubes_output)
    cubes_output = cubes_output.merge_cube()
    print(cubes_output)

    return cubes_output

        # My original look at this:
        # Leave this as int16 and apply factor in analysis
        # smaller files and lets me check the factor is correct
        # need to decide what to do in this contect

        # make a cube
        # make lat lon DimCoords

#     cube = iris.load(f'{in_dir}/{filepath}{year}{month:02d}*.nc',
#                      variable,
#                      callback=attr_ancil_callback
#                      )

#     cube = cube.concatenate_cube()

#     return 


# def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
#     """Cmorization func call.
#     """

#     glob_attrs = cfg['attributes']

#     # Run the cmorization
#     # vals has the info from the yml file
#     # var is set up in the yml file
#     print(cfg)  # while testing
#     for var, vals in cfg['variables'].items():
#         # leave this loop in as might be useful in the future for getting
#         # other info like uncertainty information from the original files
        
#         # loop over years
#         for year in range(vals['start_year'], vals['end_year'] + 1):
#             for month in range(1, 3):  # just a couple of months while testing
#                 # do this montnthly due to data size
               

#                 # now sort coordinates
#                 fixed_cubes = fix_coords(loaded_cubes)

#                 # save this year's data
#                 save_variable(fixed_cubes,
#                               var,
#                               out_dir,
#                               glob_attrs
#                               )


# def attr_ancil_callback(cube, field, filename):
#     """Callback function for loading cubes.
#     Removes cube attributes and unwanted ancillary variables.
#     """
#     cube.attributes = None



