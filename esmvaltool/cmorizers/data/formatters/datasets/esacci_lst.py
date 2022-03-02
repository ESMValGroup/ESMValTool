"""ESMValTool CMORizer for ESACCI-LST data.

Tier
   Tier 2: other freely-available dataset.

Source
   On CEDA-JASMIN
   /gws/nopw/j04/esacci_lst/public
   For access to this JASMIN group workspace please register at
   https://accounts.jasmin.ac.uk/services/group_workspaces/esacci_lst/

Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-LST
   BOTH DAY and NIGHT files are needed for each month

Currently set to work with only the MODIS AQUA L3 monthly data

Modification history
   20201015 Started by Robert King
   20201029 Day/Night averaging added along with CMOR utils
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
    for var, vals in cfg['variables'].items():
        # leave this loop in as might be useful in
        # the future for getting other info
        # like uncertainty information from the original files

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
            for month0 in range(12):  # Change this in final version
                month = month0 + 1
                logger.info(month)
                day_cube, night_cube = load_cubes(in_dir, vals['file_day'],
                                                  vals['file_night'], year,
                                                  month, variable)

                monthly_cube = make_monthly_average(day_cube, night_cube, year,
                                                    month)

                # use CMORizer utils
                monthly_cube = fix_coords(monthly_cube)

                this_years_cubes.append(monthly_cube)

            # Use utils save
            # This seems to save files all with the same name!
            # Fixed by making yearly files
            this_years_cubes = this_years_cubes.merge_cube()
            this_years_cubes.long_name = 'Surface Temperature'
            this_years_cubes.standard_name = 'surface_temperature'

            save_variable(
                this_years_cubes,
                var,
                out_dir,
                glob_attrs,
            )


def load_cubes(in_dir, file_day, file_night, year, month, variable):
    """Variable description.

    variable = land surface temperature
    platform = AQUA not used for now
               but in place for future expansion to all ESC CCI LST platforms
    """
    logger.info('Loading %s/%s%s%s*.nc', in_dir, file_day, year, month)
    day_cube = iris.load_cube(
        '%s/%s%s%02d*.nc' % (in_dir, file_day, year, month), variable)
    logger.info('Loading %s/%s%s%s*.nc', in_dir, file_night, year, month)
    night_cube = iris.load_cube(
        '%s/%s%s%02d*.nc' % (in_dir, file_night, year, month), variable)

    return day_cube, night_cube


def make_monthly_average(day_cube, night_cube, year, month):
    """Make the average LST form the day time and night time files."""
    day_cube.attributes.clear()
    night_cube.attributes.clear()

    co_time = night_cube.coord('time')
    co_time.points = co_time.points + 100.0
    # maybe the arbitrary difference should go on day cubes to
    # take the timestamp to 12Z?
    # not really an issue when using monthly files

    result = iris.cube.CubeList([day_cube, night_cube]).concatenate_cube()

    # This corrects the lonitude coord name issue
    # This should be fixed in the next version of the CCI data
    logger.info("Longitude coordinate correction being applied")
    result.coords()[2].var_name = 'longitude'
    result.coords()[2].standard_name = 'longitude'
    result.coords()[2].long_name = 'longitude'

    monthly_cube = result.collapsed('time', iris.analysis.MEAN)

    # fix time coordinate bounds
    monthly_co_time = monthly_cube.coord('time')

    time_point = (datetime.datetime(year, month, 1, 0, 0) -
                  datetime.datetime(1981, 1, 1, 0, 0, 0)).total_seconds()
    monthly_co_time.points = time_point

    num_days = monthrange(year, month)[1]
    monthly_co_time.bounds = [
        time_point, time_point + ((num_days - 1) * 24 * 3600)
    ]
    # should this be num_days or num_days-1 ### question for Valeriu or Axel
    # or 23:59:59 ???

    monthly_cube.attributes = {
        'information': 'Mean of Day and Night Aqua MODIS monthly LST'
    }

    return monthly_cube
