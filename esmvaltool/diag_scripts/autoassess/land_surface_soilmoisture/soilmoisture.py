"""Run module for soil moisture metrics."""

import os
import logging
import numpy as np
import iris
from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean


logger = logging.getLogger(__name__)


def land_sm_top(run):
    """
    Calculate median absolute errors for soil mosture against CCI data.

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns:
        metrics - dictionary of metrics names and values

    """
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    seasons = ['djf', 'mam', 'jja', 'son']

    # Constants
    # density of water and ice
    rhow = 1000.
    rhoi = 917.
    # first soil layer depth
    dz1 = 0.1

    #   Work through each season
    metrics = dict()
    for season in seasons:
        fname = 'ecv_soil_moisture_{}.nc'.format(season)
        clim_file = os.path.join(run['climfiles_root'], fname)
        ecv_clim = iris.load_cube(clim_file)
        # correct invalid units
        if (ecv_clim.units == 'unknown' and
                'invalid_units' in ecv_clim.attributes):
            if ecv_clim.attributes['invalid_units'] == 'm^3m^-3':
                ecv_clim.units = 'm3 m-3'

        # m01s08i223
        # standard_name: mrsos
        smcl_run = get_supermean('moisture_content_of_soil_layer', season,
                                 supermean_data_dir)

        # m01s08i229i
        # standard_name: ???
        # TODO: uncomment when implemented
        # sthu_run = get_supermean(
        #     'mass_fraction_of_unfrozen_water_in_soil_moisture', season,
        #     supermean_data_dir)

        # m01s08i230
        # standard_name: ??? soil_frozen_water_content - mrfso
        # TODO: uncomment when implemented
        # sthf_run = get_supermean(
        #     'mass_fraction_of_frozen_water_in_soil_moisture', season,
        #     supermean_data_dir)

        # TODO: remove after correct implementation
        sthu_run = smcl_run
        sthf_run = smcl_run

        # extract top soil layer
        cubes = [smcl_run, sthu_run, sthf_run]
        for i, cube in enumerate(cubes):
            if cube.coord('depth').attributes['positive'] != 'down':
                logger.warning('Cube %s depth attribute is not down', cube)
            top_level = min(cube.coord('depth').points)
            topsoil = iris.Constraint(depth=top_level)
            cubes[i] = cube.extract(topsoil)
        smcl_run, sthu_run, sthf_run = cubes

        # Set all sea points to missing data np.nan
        smcl_run.data[smcl_run.data < 0] = np.nan
        sthu_run.data[sthu_run.data < 0] = np.nan
        sthf_run.data[sthf_run.data < 0] = np.nan

        # set soil moisture to missing data on ice points (i.e. no soil)
        sthu_plus_sthf = (dz1 * rhow * sthu_run) + (dz1 * rhoi * sthf_run)
        ice_pts = sthu_plus_sthf.data == 0
        sthu_plus_sthf.data[ice_pts] = np.nan

        # Calculate the volumetric soil moisture in m3/m3
        theta_s_run = smcl_run / sthu_plus_sthf
        vol_sm1_run = theta_s_run * sthu_run
        vol_sm1_run.units = "m3 m-3"
        vol_sm1_run.long_name = "Top layer Soil Moisture"

        # update the coordinate system ECV data with a WGS84 coord system
        # TODO: ask Heather why this is needed
        # TODO: who is Heather?
        # unify coord systems for regridder
        vol_sm1_run.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        vol_sm1_run.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('longitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)
        ecv_clim.coord('latitude').coord_system = \
            iris.coord_systems.GeogCS(semi_major_axis=6378137.0,
                                      inverse_flattening=298.257223563)

        # Interpolate to the grid of the climatology and form the difference
        vol_sm1_run = regrid(vol_sm1_run, ecv_clim, 'linear')
        # diff the cubes
        dff = vol_sm1_run - ecv_clim

        # Remove NaNs from data before aggregating statistics
        dff.data = np.ma.masked_invalid(dff.data)

        # save output
        iris.save(dff, os.path.join(run['dump_output'],
                                    'soilmoist_diff_{}.nc'.format(season)))
        name = 'soilmoisture MedAbsErr {}'.format(season)
        metrics[name] = float(np.ma.median(np.ma.abs(dff.data)))

    return metrics
