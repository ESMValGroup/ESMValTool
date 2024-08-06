"""ESMValTool CMORizer for ESACCI-LST data.

Tier
   Tier 2: other freely-available dataset.
Source
   CMSAF website www.cmsaf.eu
Download and processing instructions
   Order .tar files from CMSAF website
   Untar to get netcdf files, one for each hour of every day
   To fix non-CF complient naming of semiminor and semimajor axes in the
   netcdf files run this shell script in the directory:
   for file in ls *.nc;
   do
   ncrename -O -h -a grid_mapping@semiminoraxis,semi_minor_axis \
                  -a grid_mapping@semimajoraxis,semi_major_axis ${file};
   done

   This means that Iris will load the files sucessfully, otherwise it
   will raise ValueError: No ellipsoid specified
"""

import logging

##### NOTE esmvalcore custom ts 1hr .dat file added

import iris

from ...utilities import fix_coords, save_variable

logger = logging.getLogger(__name__)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call.
    """

    glob_attrs = cfg['attributes']

    # Run the cmorization
    # vals has the info from the yml file
    # var is set up in the yml file
    print(cfg)  # while testing
    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in the future for getting
        # other info like uncertainty information from the original files
        glob_attrs['mip'] = vals['mip']

        for key in vals.keys():
            logger.info("%s %s", key, vals[key])

        variable = vals['raw']
        # loop over years
        for year in range(vals['start_year'], vals['end_year'] + 1):
            for month in range(1, 13):  # just a couple of months while testing
                # do this montnthly due to data size
                loaded_cubes = load_cubes(in_dir,
                                          vals['file'],
                                          year,
                                          month,
                                          variable
                                          )

                # now sort coordinates
                # fix_coords breaks longitude!
                # fixed_cubes = fix_coords(loaded_cubes)

                # save this year's data
#               save_variable(fixed_cubes,
                save_variable(loaded_cubes,
                              var,
                              out_dir,
                              glob_attrs
                              )


def attr_ancil_callback(cube, field, filename):
    """Callback function for loading cubes.

    Removes cube attributes and unwanted ancillary variables.
    """
    cube.attributes = None


def load_cubes(in_dir, filepath, year, month, variable):
    """Load ESA CCI SM netcdf files.

    Also remove attributes and ancillary variables to allow all of the
    year's data to concatenate into a single cube.
    """
    logger.info(f'Loading {in_dir}/{filepath}{year}{month:02d}*.nc')
    cube = iris.load(f'{in_dir}/{filepath}{year}{month:02d}*.nc',
                     variable,
                     callback=attr_ancil_callback
                     )

    cube = cube.concatenate_cube()

    return cube
