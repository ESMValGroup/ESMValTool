"""Module for snow metrics"""

import os
import numpy as np
import iris

from esmvaltool.preprocessor._regrid import regrid
from .supermeans import get_supermean


def land_swe_top(run):
    """
    Compute median-absolute difference of SWE against GlobSnow

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns:
        metrics - dictionary of metrics names and values

    """
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    snow_seasons = ['son', 'djf', 'mam']

    # Location of climatology
    # TODO (VPREDOI) Make sure this is a standard path
    # clim_swe_dir = os.path.join(run['clim_root'], 'GlobSnow')
    clim_swe_dir = os.path.join('/home/users/valeriu', 'GlobSnow')

    # Calculate rms errors for seasons with snow.
    metrics = dict()
    for season in snow_seasons:
        clim_file = os.path.join(clim_swe_dir, 'SWE_clm_{}.pp'.format(season))
        swe_clim = iris.load_cube(clim_file)

        # snowfall
        swe_run = get_supermean('snowfall_flux', season, supermean_data_dir)

        # Force same coord_system
        swe_run.coord('longitude').coord_system = swe_clim.coord(
            'longitude').coord_system
        swe_run.coord('latitude').coord_system = swe_clim.coord(
            'latitude').coord_system

        # Force the units for SWE to match the model
        swe_clim.units = swe_run.units

        # form the difference
        # active regridding here
        swe_run = regrid(swe_run, swe_clim, 'linear')
        dff = swe_run - swe_clim

        #  Calculate median absolute error of the difference
        name = "snow MedAbsErr {}".format(season)
        metrics[name] = float(np.ma.median(np.ma.abs(dff.data)))

    return metrics
