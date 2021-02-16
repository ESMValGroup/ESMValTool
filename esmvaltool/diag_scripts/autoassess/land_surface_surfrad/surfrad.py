"""Module to compute surface radiation metrics."""

import os

import numpy as np

import iris

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared._supermeans import get_supermean


def land_surf_rad(run):
    """
    Compute median absolute errors against CERES-EBAF data.

    Arguments:
        run - dictionary containing model run metadata
              (see auto_assess/model_run.py for description)

    Returns:
        metrics - dictionary of metrics names and values.
    """
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    rad_seasons = ['ann', 'djf', 'mam', 'jja', 'son']
    rad_fld = ['SurfRadNSW', 'SurfRadNLW']

    # Land mask: Use fractional mask for now.
    # Fraction of Land m01s03i395
    # replaced with a constant sftlf mask; original was
    # lnd = get_supermean('land_area_fraction', 'ann', supermean_data_dir)
    cubes = iris.load(os.path.join(supermean_data_dir, 'cubeList.nc'))
    lnd = cubes.extract_cube(iris.Constraint(name='land_area_fraction'))

    metrics = dict()
    for season in rad_seasons:
        for fld in rad_fld:
            if fld == 'SurfRadNSW':
                ebaf_fld = get_supermean(
                    'Surface Net downward Shortwave Radiation', season,
                    run['clim_root'], obs_flag='CERES-EBAF')
                run_fld_rad = get_supermean(
                    'Surface Net downward Shortwave Radiation', season,
                    supermean_data_dir)

            elif fld == 'SurfRadNLW':
                ebaf_fld = get_supermean(
                    'Surface Net downward Longwave Radiation', season,
                    run['clim_root'], obs_flag='CERES-EBAF')
                run_fld_rad = get_supermean(
                    'Surface Net downward Longwave Radiation', season,
                    supermean_data_dir)

            else:
                raise Exception('Skipping unassigned case.')

            # Regrid both to land points and mask out where this is below
            # a threshold. Force the coordinate system on model.
            ebaf_fld.coord('latitude').coord_system = \
                run_fld_rad.coord('latitude').coord_system
            ebaf_fld.coord('longitude').coord_system = \
                run_fld_rad.coord('longitude').coord_system
            lnd.coord('latitude').coord_system = \
                run_fld_rad.coord('latitude').coord_system
            lnd.coord('longitude').coord_system = \
                run_fld_rad.coord('longitude').coord_system

            reg_run_fld = regrid(run_fld_rad, lnd, 'linear')
            reg_ebaf_fld = regrid(ebaf_fld, lnd, 'linear')

            # apply the mask
            reg_run_fld.data = np.ma.masked_array(
                reg_run_fld.data, mask=(lnd.data > 90.))
            reg_ebaf_fld.data = np.ma.masked_array(
                reg_ebaf_fld.data, mask=(lnd.data > 90.))

            # do a simple diff
            dff = reg_run_fld - reg_ebaf_fld

            name = "{} MedAbsErr {}".format(fld, season)
            metrics[name] = float(np.ma.median(np.abs(dff.data)))

    return metrics
