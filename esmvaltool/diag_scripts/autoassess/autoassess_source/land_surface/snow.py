"""Module for snow metrics"""

import os
import sys

import numpy as np

import iris

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
    SUPERMEAN_DATA_DIR = os.path.join(run['data_root'],
                                      run['runid'],
                                      run['_area'] + '_supermeans')
    print('SUPERMEAN', SUPERMEAN_DATA_DIR)
    # create: /group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/landsurface/recipe_autoassess_landsurface_20180821_161123/work/aa_strato/autoassess_strato_test_1/MPI-ESM-LR/land_surface_supermeans

    snow_seasons = ['son', 'djf', 'mam']

    # Location of climatology
    clim_swe_dir = os.path.join(run['clim_root'], 'GlobSnow')

    # Calculate rms errors for seasons with snow.
    metrics = dict()
    for season in snow_seasons:
        clim_file = os.path.join(clim_swe_dir,
                                 'SWE_clm_{}.pp'.format(season))
        swe_clim = iris.load_cube(clim_file)
        # Need to touch data to get it recognised as a Masked Array?
        swe_clim.data

        # Force the units for SWE to match the model
        swe_clim.units = "kg m-2"

        # m01s00i023
        swe_run = get_supermean('snowfall_amount', season, SUPERMEAN_DATA_DIR)

        # Interpolate to the grid of the climatology and form the difference
        dff = swe_run.regrid(swe_clim, iris.analysis.Linear()) - swe_clim

        #  Calculate median absolute error of the difference
        name = "snow MedAbsErr {}".format(season)
        metrics[name] = float(np.ma.median(np.ma.abs(dff.data)))

    return metrics
