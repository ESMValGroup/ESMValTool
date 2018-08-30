"""Module to compute surface radiation metrics"""

import os

import numpy as np

import iris

from esmvaltool.preprocessor._regrid import regrid
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
    supermean_data_dir = os.path.join(run['data_root'], run['runid'],
                                      run['_area'] + '_supermeans')

    rad_seasons = ['ann', 'djf', 'mam', 'jja', 'son']
    rad_fld = ['SurfRadNSW', 'SurfRadNLW']

    # Land mask: Use fractional mask for now.
    # Fraction of Land m01s03i395
    # TODO VPREDOI this is replaced with a constant sftlf mask
    # lnd = get_supermean('land_area_fraction', 'ann', supermean_data_dir)
    name_constraint = iris.Constraint(name='land_area_fraction')
    cubes_path = os.path.join(supermean_data_dir, 'cubeList.nc')
    cubes = iris.load(cubes_path)
    lnd = cubes.extract_strict(name_constraint)

    metrics = dict()
    for season in rad_seasons:
        for fld in rad_fld:
            if fld == 'SurfRadNSW':
                # Net (downward) SW from EBAF:
                # rsus: 'surface_upwelling_shortwave_flux_in_air'
                # rsds: 'surface_downwelling_shortwave_flux_in_air'
                # rsns (derived, custom): rsds - rsus: ANNUAL
                # name original: 'surface_net_downward_shortwave_flux'
                # name custom CMOR: surface_net_downward_shortwave_radiation
                ebaf_fld = get_supermean(
                    'surface_net_downward_shortwave_radiation', season,
                    run['clim_root'])

                # m01s01i201
                run_fld_rad = get_supermean(
                    'surface_net_downward_shortwave_radiation', season,
                    supermean_data_dir)

            elif fld == 'SurfRadNLW':
                # Net (downward) LW from EBAF:
                # rlus: 'surface_upwelling_longwave_flux_in_air'
                # rlds: 'surface_downwelling_longwave_flux_in_air'
                # rlns (derived, custom): rlds - rlus: ANNUAL
                # name original: 'surface_net_downward_longwave_flux'
                # name custom CMOR: surface_net_downward_longwave_radiation
                ebaf_fld = get_supermean(
                    'surface_net_downward_longwave_radiation', season,
                    run['clim_root'])

                # m01s02i201
                run_fld_rad = get_supermean(
                    'surface_net_downward_longwave_radiation', season,
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

            # TODO would area weighted regridding be better?
            reg_run_fld = regrid(run_fld_rad, lnd, 'linear')
            reg_ebaf_fld = regrid(ebaf_fld, lnd, 'linear')

            # apply the mask
            reg_run_fld.data = np.ma.masked_array(
                reg_run_fld.data, mask=(lnd.data < 0.98))
            reg_ebaf_fld.data = np.ma.masked_array(
                reg_ebaf_fld.data, mask=(lnd.data < 0.98))

            # do a simple diff
            dff = reg_run_fld - reg_ebaf_fld

            name = "{} MedAbsErr {}".format(fld, season)
            metrics[name] = float(np.ma.median(np.abs(dff.data)))

    return metrics
