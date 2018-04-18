"""Module for radiation metrics"""

import os.path
from ..auto_assess_deprecated.supermeans import get_supermean
from ..utility.area_utils import area_average


def invertmask(mask):
    newmask = mask.copy()
    newmask.data = 1.0 - mask.data
    return newmask


def round_metrics(metrics_dict):
    for key, val in metrics_dict.items():
        metrics_dict[key] = round(float(val), 5)
    return metrics_dict


def global_rd(data_dir, run_id, area_name):
    """
    Compute global means of radiative quantities

    Arguments:
        run - dictionary containing model run metadata

    Returns:
        metrics - dictionary of metrics names and values

    """
    SUPERMEAN_DATA_DIR = os.path.join(data_dir, run_id, area_name + '_supermeans')

    glob_seasons = ['ann', 'djf', 'mam', 'jja', 'son']

    metrics = dict()
    for season in glob_seasons:

        # Surface net downward shortwave flux (STASH code m01s01i201)
        rsns = get_supermean('surface_net_downward_shortwave_flux', season, SUPERMEAN_DATA_DIR)
        rsns.convert_units('W m-2')
        rsnsg = area_average(rsns, weighted=True)
        metrics['Surface net downward shortwave flux Global' + season] = rsnsg.data

        # TOA Incoming shortwave flux (STASH code m01s01i207)
        rsdt = get_supermean('toa_incoming_shortwave_flux', season, SUPERMEAN_DATA_DIR)
        rsdt.convert_units('W m-2')
        rsdtg = area_average(rsdt, weighted=True)
        metrics['TOA incoming shortwave flux Global' + season] = rsdtg.data

        # TOA outgoing shortwave flux (STASH code m01s01i208)
        rsut = get_supermean('toa_outgoing_shortwave_flux', season, SUPERMEAN_DATA_DIR)
        rsut.convert_units('W m-2')
        rsutg = area_average(rsut, weighted=True)
        metrics['TOA outgoing shortwave flux Global' + season] = rsutg.data

        # TOA outgoing shortwave flux assuming clear sky(STASH code m01s01i209)
        rsutcs = get_supermean('toa_outgoing_shortwave_flux_assuming_clear_sky', season, SUPERMEAN_DATA_DIR)
        rsutcs.convert_units('W m-2')
        rsutcsg = area_average(rsutcs, weighted=True)
        metrics['TOA outgoing shortwave flux assuming clear sky Global' + season] = rsutcsg.data

        # surface_downwelling_shortwave_flux_in_air(STASH code m01s01i235)
        rsds = get_supermean('surface_downwelling_shortwave_flux_in_air', season, SUPERMEAN_DATA_DIR)
        rsds.convert_units('W m-2')
        rsdsg = area_average(rsds, weighted=True)
        metrics['Surface downwelling shortwave flux in air Global' + season] = rsdsg.data

        # surface_net_downward_longwave_flux(STASH code m01s02i201)
        rls = get_supermean('surface_net_downward_longwave_flux', season, SUPERMEAN_DATA_DIR)
        rls.convert_units('W m-2')
        rlsg = area_average(rls, weighted=True)
        metrics['Surface net downward longwave flux Global' + season] = rlsg.data

        # toa_outgoing_longwave_flux(STASH code m01s02i205)
        rlut = get_supermean('toa_outgoing_longwave_flux', season, SUPERMEAN_DATA_DIR)
        rlut.convert_units('W m-2')
        rlutg = area_average(rlut, weighted=True)
        metrics['TOA outgoing longwave flux Global' + season] = rlutg.data

        # toa_outgoing_longwave_flux_assuming_clear_sky(STASH code m01s02i206)
        rlutcs = get_supermean('toa_outgoing_longwave_flux_assuming_clear_sky', season, SUPERMEAN_DATA_DIR)
        rlutcs.convert_units('W m-2')
        rlutcsg = area_average(rlutcs, weighted=True)
        metrics['TOA outgoing longwave flux assuming clear sky Global' + season] = rlutcsg.data

        # surface_downwelling_longwave_flux_in_air(STASH code m01s02i207)
        # TODO cf standard name is 'surface_downwelling_longwave_flux_in_air'
        rlds = get_supermean('surface_downwelling_longwave_flux', season, SUPERMEAN_DATA_DIR)
        rlds.convert_units('W m-2')
        rldsg = area_average(rlds, weighted=True)
        metrics['Surface downwelling longwave flux in air Global ' + season] = rldsg.data

        # surface_upward_sensible_heat_flux(STASH code m01s03i217)
        hfss = get_supermean('surface_upward_sensible_heat_flux', season, SUPERMEAN_DATA_DIR)
        hfss.convert_units('W m-2')
        hfssg = area_average(hfss, weighted=True)
        metrics['Surface upward sensible heat flux Global' + season] = hfssg.data

        # surface_upward_latent_heat_flux(STASH code m01s03i234)
        hfls = get_supermean('surface_upward_sensible_heat_flux', season, SUPERMEAN_DATA_DIR)
        hfls.convert_units('W m-2')
        hflsg = area_average(hfls, weighted=True)
        metrics['Surface upward latent heat flux' + season] = hflsg.data

    return round_metrics(metrics)
