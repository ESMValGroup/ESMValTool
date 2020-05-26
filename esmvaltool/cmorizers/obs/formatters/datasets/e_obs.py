"""ESMValTool CMORizer for E-OBS data.

Tier
    Tier 2: other freely-available dataset.

Source
    http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles

Last access
    20200225

Download and processing instructions
    Download the ensemble mean files for:
        TG TN TX RR PP
"""

import logging
import os
import iris
import numpy as np
from cf_units import Unit

from esmvalcore.preprocessor import monthly_statistics
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def fix_coords_non_symetric_lon(cube):
    """Fix the time units and values to CMOR standards."""
    # first fix any completely missing coord var names
    utils._fix_dim_coordnames(cube)
    # fix individual coords
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == 'time':
            logger.info("Fixing time...")
            cube.coord('time').convert_units(
                Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
            utils._fix_bounds(cube, cube.coord('time'))

        # fix longitude
        if cube_coord.var_name == 'lon':
            logger.info("Fixing longitude...")
            if cube_coord.ndim == 1:
                if cube_coord.points[0] < 0. and \
                        cube_coord.points[-1] < 181.:
                    lon_coord = cube.coord('longitude').copy()
                    lons_below_0 = lon_coord.points[lon_coord.points < 0.] + \
                        360.
                    lons_above_0 = lon_coord.points[lon_coord.points >= 0.]
                    lons = np.hstack((lons_above_0, lons_below_0))
                    cube_coord.points = lons

                    utils._fix_bounds(cube, cube_coord)
                    cube.attributes['geospatial_lon_min'] = 0.
                    cube.attributes['geospatial_lon_max'] = 360.
                    utils._roll_cube_data(cube, len(lons_above_0), -1)

        # fix latitude
        if cube_coord.var_name == 'lat':
            logger.info("Fixing latitude...")
            utils._fix_bounds(cube, cube.coord('latitude'))

        # fix depth
        if cube_coord.var_name == 'lev':
            logger.info("Fixing depth...")
            utils._fix_bounds(cube, cube.coord('depth'))

        # fix air_pressure
        if cube_coord.var_name == 'air_pressure':
            logger.info("Fixing air pressure...")
            utils._fix_bounds(cube, cube.coord('air_pressure'))

    # remove CS
    cube.coord('latitude').coord_system = None
    cube.coord('longitude').coord_system = None

    return cube


def _extract_variable(short_name, var, res, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.units = var.get('raw_units', short_name)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    fix_coords_non_symetric_lon(cube)
    if 'height2m' in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    utils.fix_var_metadata(cube, cmor_info)
    attrs = cfg['attributes'].copy()
    attrs['version'] = 'v' + attrs['version'] + '-' + str(res)
    attrs.pop('resolution')
    attrs['mip'] = var['mip']
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    #####
    # also derive monthly data
    if 'add_mon' in var:
        if var['add_mon']:
            logger.info("Building monthly means")

            # Calc monthly
            cube = monthly_statistics(cube)
            cube.remove_coord('month_number')
            cube.remove_coord('year')

            # Fix metadata
            attrs['mip'] = 'Amon'

            # Fix coordinates
            fix_coords_non_symetric_lon(cube)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    ver = cfg['attributes']['version']
    for res in cfg['attributes']['resolution'].values():
        for (short_name, var) in cfg['variables'].items():
            logger.info("CMORizing variable '%s' on %s°x%s°",
                        short_name, res, res)
            raw_var = var.get('raw', short_name)
            filepath = raw_filepath.format(raw_name=raw_var, resolution=res,
                                           version=ver)
            _extract_variable(short_name, var, res, cfg, filepath, out_dir)
