#  CMORiser for ESA CCI SOil Moisture V7.1 

"""ESMValTool CMORizer for ESACCI-Soil Moisture data.
Tier
   Tier 2: other freely-available dataset.
Source

Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-SM

This is for the daily COMBINED product.
"""

import datetime
import logging
from calendar import monthrange

import iris

from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization

    # vals has the info from the yml file
    # var is set up in the yml file
    print(cfg)

    for var, vals in cfg['variables'].items():
        print('##############')
        print(var)
        print('............................')
        print(vals)
        # leave this loop in as might be useful in
        # the future for getting other info
        # like uncertainty information from the original files

        glob_attrs['mip'] = vals['mip']

        for key in vals.keys():
            logger.info("%s %s", key, vals[key])

        variable = vals['raw']

        # loop over years and months
        # get years from start_year and end_year
        for year in range(vals['start_year'], vals['end_year'] + 1):
            print(year)
            print(in_dir)
            cube = load_cubes(in_dir,
                              vals['file'],
                              year,
                              variable
                              )


def attribute_callback(cube, field, filename):
    cube.attributes = None

    for coord in cube.ancillary_variables():
        cube.remove_ancillary_variable(coord)

    #cube.remove_ancillary_variable('Observation Timestamp')
    #cube.remove_ancillary_variable('soil_moisture_content standard_error')
    #cube.remove_ancillary_variable('soil_moisture_content status_flag')


def load_cubes(in_dir, filepath, year, variable):
    """Need to write description.
    """
    logger.info('Loading %s/%s%s*.nc', in_dir, filepath, year)
    cube = iris.load(f'{in_dir}/{filepath}{year}01*.nc',
                          variable,
                          callback = attribute_callback
                      )
    cube = cube.concatenate_cube()
    print(cube)


    return cube
