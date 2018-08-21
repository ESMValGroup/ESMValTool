"""Module to compute surface radiation metrics"""

import os

import numpy as np

import iris

from .supermeans import get_supermean


def land_surf_rad(run):
    """
    Compute median absolute errors against CERES-EBAF data

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns:
        metrics - dictionary of metrics names and values

    """
    SUPERMEAN_DATA_DIR = os.path.join(run['data_root'],
                                      run['runid'],
                                      run['_area'] + '_supermeans')

    rad_seasons = ['ann', 'djf', 'mam', 'jja', 'son']
    rad_fld = ['SurfRadNSW', 'SurfRadNLW']

    # Location of climatology
    ebaf_dir = os.path.join(run['clim_root'], 'CERES', 'EBAF_CMIP5',
                            'supermeansnc')

    # Land mask: Use fractional mask for now.
    # Fraction of Land m01s03i395
    lnd = get_supermean('land_area_fraction', 'ann', SUPERMEAN_DATA_DIR)
    incld = np.where(lnd.data < 0.98)

    metrics = dict()
    for season in rad_seasons:

        # Set the seasonal file
        # TODO obs filename
        fname = 'CERES-EBAF_Ed2-6r_200003-201002.{}.nc'.format(season)
        clm_file = os.path.join(ebaf_dir, fname)


        for fld in rad_fld:
            if (fld == 'SurfRadNSW'):
                # Net (downward) SW from EBAF:
                cube_name = 'surface_upwelling_shortwave_flux_in_air'
                ebaf_rsus = iris.load_cube(clm_file, cube_name)
                cube_name = 'surface_downwelling_shortwave_flux_in_air'
                ebaf_rsds = iris.load_cube(clm_file, cube_name)

                # remove the redundant time axis by selecting 1st time
                ebaf_fld = ebaf_rsds[0]
                ebaf_fld.data -= ebaf_rsus[0].data

                # m01s01i201
                run_fld_rad = get_supermean('surface_net_downward_shortwave_flux',
                                            season,
                                            SUPERMEAN_DATA_DIR)

            elif (fld == 'SurfRadNLW'):
                # Net (downward) LW from EBAF:
                cube_name = 'surface_upwelling_longwave_flux_in_air'
                ebaf_rlus = iris.load_cube(clm_file, cube_name)
                cube_name = 'surface_downwelling_longwave_flux_in_air'
                ebaf_rlds = iris.load_cube(clm_file, cube_name)

                # remove the redundant time axis by selecting 1st time
                ebaf_fld = ebaf_rlds[0]
                ebaf_fld.data -= ebaf_rlus[0].data

                # m01s02i201
                run_fld_rad = get_supermean('surface_net_downward_longwave_flux',
                                            season,
                                            SUPERMEAN_DATA_DIR)

            else:
                raise Exception('Skipping unassigned case.')

            # Regrid both to land points and mask out where this is below
            # a threshold. Force the coordinate system on EBAF.
            ebaf_fld.coord('latitude').coord_system = \
                run_fld_rad.coord('latitude').coord_system
            ebaf_fld.coord('longitude').coord_system = \
                run_fld_rad.coord('longitude').coord_system
            # TODO would area weighted regridding be better?
            dff = run_fld_rad.regrid(lnd, iris.analysis.Linear()) - \
                ebaf_fld.regrid(lnd, iris.analysis.Linear())

            name = "{} MedAbsErr {}".format(fld, season)
            metrics[name] = float(np.median(np.abs(dff.data[incld])))

    return metrics
