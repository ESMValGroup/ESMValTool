"""ESMValTool CMORizer for ESACCI-Soil Moisture data.

Tier
   Tier 2: other freely-available dataset.
Source
   ftp.geo.tuwien.ac.at
   An account can be obtained from https://esa-soilmoisture-cci.org/

Download and processing instructions
   Put all the daily files under a single directory with no subdirectories
   with years in ${RAWOBS}/Tier2/ESACCI-SM

This is for the daily COMBINED product.
"""
import logging
import iris
from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']

    # Run the cmorization
    # vals has the info from the yml file
    # var is set up in the yml file

    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in the future for getting
        # other info like uncertainty information from the original files
        glob_attrs['mip'] = vals['mip']

        for key in vals.keys():
            logger.info("%s %s", key, vals[key])

        variable = vals['raw']

        # loop over years
        for year in range(vals['start_year'], vals['end_year'] + 1):
            cubes = load_cubes(in_dir,
                               vals['file'],
                               year,
                               variable
                               )

            # now sort coordinates
            cubes = fix_coords(cubes)

            # save this year's data
            save_variable(cubes,
                          var,
                          out_dir,
                          glob_attrs
                          )


def attr_ancil_callback(cube, field, filename):
    """Callback function for loading cubes.

    Removes cube attributes and unwanted ancillary variables.
    """
    cube.attributes = None

    for coord in cube.ancillary_variables():
        cube.remove_ancillary_variable(coord)


def load_cubes(in_dir, filepath, year, variable):
    """Load ESA CCI SM netcdf files.

    Also remove attributes and ancillary variables to allow all of the
    year's data to concatenate into a single cube.
    """
    logger.info(f'Loading {in_dir}/{filepath}{year}*.nc')
    cube = iris.load(f'{in_dir}/{filepath}{year}*.nc',
                     variable,
                     callback=attr_ancil_callback
                     )
    cube = cube.concatenate_cube()

    return cube
