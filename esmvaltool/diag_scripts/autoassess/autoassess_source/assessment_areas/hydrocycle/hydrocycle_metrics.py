"""Module for global hydrocycle metrics"""

import os.path
import sys

import numpy as np

import iris
import iris.coord_categorisation
from ..auto_assess_deprecated.supermeans import get_supermean

from ..utility.area_utils import area_average
from ..auto_assess_deprecated.loaddata import load_run_ss

import hydrocycle_get_landfrac


def invertmask(mask):
    newmask = mask.copy()
    newmask.data = 1.0 - mask.data
    return newmask


def round_metrics(metrics_dict):
    for key, val in metrics_dict.items():
        metrics_dict[key] = round(float(val), 5)
    return metrics_dict


def global_hc(run):
    """
    Compute global means of hydrological quantities

    Arguments:
        run - dictionary containing model run metadata

    Returns:
        metrics - dictionary of metrics names and values

    """
    SUPERMEAN_DATA_DIR = os.path.join(run['data_root'],
                                      run['runid'],
                                      run['_area'] + '_supermeans')

    glob_seasons = ['ann', 'djf', 'mam', 'jja', 'son']

    # Latitude definitions for sub-regions
    ro_global = (-60., 90., False, True)

    metrics = dict()
    for season in glob_seasons:

        # Precipitation (STASH code m01s05i216)
        ppn = get_supermean('precipitation_flux', season, SUPERMEAN_DATA_DIR)
        ppn.convert_units('kg m-2 day-1')
        ppng = area_average(ppn, weighted=True)
        metrics['Precipitation Global ' + season] = ppng.data

        # Evaporation (STASH code m01s03i223)
        evap = get_supermean('surface_upward_water_flux', season, SUPERMEAN_DATA_DIR)
        evap.convert_units('kg m-2 day-1')
        evapg = area_average(evap, weighted=True)
        metrics['Evaporation Global ' + season] = evapg.data

        # P minus E
        pme = ppn - evap
        pmeg = area_average(pme, weighted=True)
        metrics['P minus E Global ' + season] = pmeg.data

        # Water vapour (STASH code m01s30i403 and m01s30i404)
        dmt = get_supermean('m01s30i403', season, SUPERMEAN_DATA_DIR)  # TOTAL COLUMN DRY MASS  RHO GRID
        wvti = get_supermean('atmosphere_mass_per_unit_area', season, SUPERMEAN_DATA_DIR)
        dmt.units = wvti.units
        wvt = wvti - dmt
        wvtg = area_average(wvt, weighted=True)
        metrics['Total water vapour Global ' + season] = wvtg.data

        # Cloud liquid water (STASH code m01s30i405)
        qlt = get_supermean('atmosphere_cloud_liquid_water_content', season, SUPERMEAN_DATA_DIR)
        qltg = area_average(qlt, weighted=True)
        metrics['Total cloud liquid water Global ' + season] = qltg.data

        # Cloud ice water (STASH code m01s30i406)
        qit = get_supermean('atmosphere_cloud_ice_content', season, SUPERMEAN_DATA_DIR)
        qitg = area_average(qit, weighted=True)
        metrics['Total cloud ice water Global ' + season] = qitg.data

        # Sensible heat (STASH code m01s03i217)
        shf = get_supermean('surface_upward_sensible_heat_flux', season, SUPERMEAN_DATA_DIR)
        shfg = area_average(shf, weighted=True)
        metrics['Sensible heat flux Global ' + season] = shfg.data

        # Latent heat (STASH code m01s03i234)
        lhf = get_supermean('surface_upward_latent_heat_flux', season, SUPERMEAN_DATA_DIR)
        lhfg = area_average(lhf, weighted=True)
        metrics['Latent heat flux Global ' + season] = lhfg.data

        # Get land and sea fractions
        try:
            # Get land fraction (STASH code m01s03i395) from supermeans
            landfrac = get_supermean('land_area_fraction', season, SUPERMEAN_DATA_DIR)
        except iris.exceptions.ConstraintMismatchError:
            print 'Cannot find land fraction in supermeans'
            landfrac = hydrocycle_get_landfrac.get_landfrac(ppn)
            # Remove bounding from land mask as it gets in the way
            landfrac.coord('longitude').bounds = None
            landfrac.coord('latitude').bounds = None
        seafrac = invertmask(landfrac)

        # Now for land and sea separately:
        #  create means weighted by land/sea fraction
        ppnlg = area_average(ppn, weighted=True, mask=landfrac)
        ppnsg = area_average(ppn, weighted=True, mask=seafrac)

        evaplg = area_average(evap, weighted=True, mask=landfrac)
        evapsg = area_average(evap, weighted=True, mask=seafrac)

        wvtlg = area_average(wvt, weighted=True, mask=landfrac)
        wvtsg = area_average(wvt, weighted=True, mask=seafrac)

        qltlg = area_average(qlt, weighted=True, mask=landfrac)
        qltsg = area_average(qlt, weighted=True, mask=seafrac)

        qitlg = area_average(qit, weighted=True, mask=landfrac)
        qitsg = area_average(qit, weighted=True, mask=seafrac)

        shflg = area_average(shf, weighted=True, mask=landfrac)
        shfsg = area_average(shf, weighted=True, mask=seafrac)

        lhflg = area_average(lhf, weighted=True, mask=landfrac)
        lhfsg = area_average(lhf, weighted=True, mask=seafrac)

        # Runoff (STASH code m01s08i234 and m01s08i235)
        sroi = get_supermean('surface_runoff_flux', season, SUPERMEAN_DATA_DIR)
        sroi.convert_units('kg m-2 day-1')

        ssroi = get_supermean('subsurface_runoff_flux', season, SUPERMEAN_DATA_DIR)
        ssroi.convert_units('kg m-2 day-1')

        # Convert run-offs to full gridbox means (output as per land fraction)
        sro = sroi.copy() * landfrac
        ssro = ssroi.copy() * landfrac

        # Calculate total runoff
        tro = sro + ssro

        # Runoff observations extend 60S to 90N
        kwargs = dict(weighted=True, mask=landfrac, latitude=ro_global)
        trolg = area_average(tro, **kwargs)
        srolg = area_average(sro, **kwargs)
        ssrolg = area_average(ssro, **kwargs)

        # Enter metrics into dictionary
        pmelg = ppnlg - evaplg
        pmesg = ppnsg - evapsg

        metrics['Precipitation Land ' + season] = ppnlg.data
        metrics['Precipitation Ocean ' + season] = ppnsg.data

        metrics['Evaporation Land ' + season] = evaplg.data
        metrics['Evaporation Ocean ' + season] = evapsg.data

        metrics['Total cloud liquid water Land ' + season] = qltlg.data
        metrics['Total cloud liquid water Ocean ' + season] = qltsg.data

        metrics['Total cloud ice water Land ' + season] = qitlg.data
        metrics['Total cloud ice water Ocean ' + season] = qitsg.data

        metrics['Latent heat flux Land ' + season] = lhflg.data
        metrics['Latent heat flux Ocean ' + season] = lhfsg.data

        metrics['Sensible heat flux Land ' + season] = shflg.data
        metrics['Sensible heat flux Ocean ' + season] = shfsg.data

        metrics['Total water vapour Land ' + season] = wvtlg.data
        metrics['Total water vapour Ocean ' + season] = wvtsg.data

        metrics['Total Runoff Land ' + season] = trolg.data

        metrics['Surface Runoff Land ' + season] = srolg.data

        metrics['Sub-surface Runoff Land ' + season] = ssrolg.data

        metrics['P minus E Land ' + season] = pmelg.data
        metrics['P minus E Ocean ' + season] = pmesg.data

    return round_metrics(metrics)


def regional_hc(run):
    """
    Compute regional means of hydrological quantities

    Arguments:
        run - dictionary containing model run metadata

    Returns:
        metrics - dictionary of metrics names and values

    """
    SUPERMEAN_DATA_DIR = os.path.join(run['data_root'],
                                      run['runid'],
                                      run['_area'] + '_supermeans')

    glob_seasons = ['ann', 'djf', 'mam', 'jja', 'son']

    # Latitude definitions for sub-regions
    tropics = (-30., 30., True, True)
    nhextrop = (30., 90., False, True)
    shextrop = (-90., -30., True, False)
    ro_shextrop = (-60., -30., False, False)

    metrics = dict()
    for season in glob_seasons:

        # Precipitation (STASH code m01s05i216)
        ppn = get_supermean('precipitation_flux', season, SUPERMEAN_DATA_DIR)
        ppn.convert_units('kg m-2 day-1')

        ptropavg = area_average(ppn, weighted=True, latitude=tropics)
        pnhavg = area_average(ppn, weighted=True, latitude=nhextrop)
        pshavg = area_average(ppn, weighted=True, latitude=shextrop)

        metrics["Mean Precip: Tropics " + season] = ptropavg.data
        metrics["Mean Precip: NH Extra-tropics " + season] = pnhavg.data
        metrics["Mean Precip: SH Extra-tropics " + season] = pshavg.data

        # Get land and sea fractions
        try:
            # Get land fraction (STASH code m01s03i395) from supermeans
            landfrac = get_supermean('land_area_fraction', season, SUPERMEAN_DATA_DIR)
        except iris.exceptions.ConstraintMismatchError:
            print 'Cannot find land fraction in supermeans'
            landfrac = hydrocycle_get_landfrac.get_landfrac(ppn)
            # Remove bounding from land mask as it gets in the way
            landfrac.coord('longitude').bounds = None
            landfrac.coord('latitude').bounds = None
        seafrac = invertmask(landfrac)

        # Evaporation (STASH code m01s03i223)
        evap = get_supermean('surface_upward_water_flux', season, SUPERMEAN_DATA_DIR)
        evap.convert_units('kg m-2 day-1')

        etropavg = area_average(evap, weighted=True, latitude=tropics)
        enhavg = area_average(evap, weighted=True, latitude=nhextrop)
        eshavg = area_average(evap, weighted=True, latitude=shextrop)

        metrics["Mean Evap: Tropics " + season] = etropavg.data
        metrics["Mean Evap: NH Extra-tropics " + season] = enhavg.data
        metrics["Mean Evap: SH Extra-tropics " + season] = eshavg.data

        # Latent heat (STASH code m01s03i234)
        lhf = get_supermean('surface_upward_latent_heat_flux', season, SUPERMEAN_DATA_DIR)

        lhtropavg = area_average(lhf, weighted=True, mask=seafrac, latitude=tropics)
        lhnhavg = area_average(lhf, weighted=True, mask=seafrac, latitude=nhextrop)
        lhshavg = area_average(lhf, weighted=True, mask=seafrac, latitude=shextrop)

        metrics["Mean Sfc Latent heat: Tropical Ocean " + season] = lhtropavg.data
        metrics["Mean Sfc Latent heat: NH Extra-tropical Ocean " + season] = lhnhavg.data
        metrics["Mean Sfc Latent heat: SH Extra-tropical Ocean " + season] = lhshavg.data

        # Sensible heat (STASH code m01s03i217)
        shf = get_supermean('surface_upward_sensible_heat_flux', season, SUPERMEAN_DATA_DIR)

        shtropavg = area_average(shf, weighted=True, mask=landfrac, latitude=tropics)
        shnhavg = area_average(shf, weighted=True, mask=landfrac, latitude=nhextrop)
        shshavg = area_average(shf, weighted=True, mask=landfrac, latitude=shextrop)

        metrics["Mean Sfc Sensible heat: Tropical Land " + season] = shtropavg.data
        metrics["Mean Sfc Sensible heat: NH Extra-tropical Land " + season] = shnhavg.data
        metrics["Mean Sfc Sensible heat: SH Extra-tropical Land " + season] = shshavg.data

        # Runoff (STASH code m01s08i234 and m01s08i235)
        sroi = get_supermean('surface_runoff_flux', season, SUPERMEAN_DATA_DIR)
        sroi.convert_units('kg m-2 day-1')

        ssroi = get_supermean('subsurface_runoff_flux', season, SUPERMEAN_DATA_DIR)
        ssroi.convert_units('kg m-2 day-1')

        # Convert run-offs to full gridbox means (output as per land fraction)
        sro = sroi.copy() * landfrac
        ssro = ssroi.copy() * landfrac

        # Calculate total runoff
        tro = sro + ssro

        trotropavg = area_average(tro, weighted=True, mask=landfrac, latitude=tropics)
        tronhavg = area_average(tro, weighted=True, mask=landfrac, latitude=nhextrop)
        troshavg = area_average(tro, weighted=True, mask=landfrac, latitude=ro_shextrop)

        metrics["Total runoff: Tropical Land " + season] = trotropavg.data
        metrics["Total runoff: NH Extra-tropical Land " + season] = tronhavg.data
        metrics["Total runoff: SH Extra-tropical Land " + season] = troshavg.data

        # Water vapour (STASH code m01s30i403 and m01s30i404)
        dmt = get_supermean('m01s30i403', season, SUPERMEAN_DATA_DIR)  # TOTAL COLUMN DRY MASS  RHO GRID
        wvti = get_supermean('atmosphere_mass_per_unit_area', season, SUPERMEAN_DATA_DIR)
        dmt.units = wvti.units
        wvt = wvti - dmt

        wvtropavg = area_average(wvt, weighted=True, latitude=tropics)
        wvnhavg = area_average(wvt, weighted=True, latitude=nhextrop)
        wvshavg = area_average(wvt, weighted=True, latitude=shextrop)

        metrics["Mean Water Vapour: Tropics " + season] = wvtropavg.data
        metrics["Mean Water Vapour: NH Extra-tropics " + season] = wvnhavg.data
        metrics["Mean Water Vapour: SH Extra-tropics " + season] = wvshavg.data

        # Cloud liquid water (STASH code m01s30i405)
        qlt = get_supermean('atmosphere_cloud_liquid_water_content', season, SUPERMEAN_DATA_DIR)

        qlttropavg = area_average(qlt, weighted=True, latitude=tropics)
        qltnhavg = area_average(qlt, weighted=True, latitude=nhextrop)
        qltshavg = area_average(qlt, weighted=True, latitude=shextrop)

        metrics["Mean Cloud liquid water: Tropics " + season] = qlttropavg.data
        metrics["Mean Cloud liquid water: NH Extra-tropics " + season] = qltnhavg.data
        metrics["Mean Cloud liquid water: SH Extra-tropics " + season] = qltshavg.data

        # Cloud ice water (STASH code m01s30i406)
        qit = get_supermean('atmosphere_cloud_ice_content', season, SUPERMEAN_DATA_DIR)

        qittropavg = area_average(qit, weighted=True, latitude=tropics)
        qitnhavg = area_average(qit, weighted=True, latitude=nhextrop)
        qitshavg = area_average(qit, weighted=True, latitude=shextrop)

        metrics["Mean Cloud ice water: Tropics " + season] = qittropavg.data
        metrics["Mean Cloud ice water: NH Extra-tropics " + season] = qitnhavg.data
        metrics["Mean Cloud ice water: SH Extra-tropics " + season] = qitshavg.data

    return round_metrics(metrics)


def regional_hc_sdev(run):

    """
    Compute regional interannual standard deviations of hydrological quantities
    Currently just for precipitation.

    Arguments:
        run - dictionary containing model run metadata

    Returns:
        metrics - dictionary of metrics names and values
    """

    metrics = dict()

    # Latitude definitions for sub-regions
    tropics = (-30., 30., True, True)
    nhextrop = (30., 90., False, True)
    shextrop = (-90., -30., True, False)

    # First use annual means
    ppn = load_run_ss(run, 'annual', 'precipitation_flux')  # m01s05i216
    ppn.convert_units('kg m-2 day-1')
    ann = ppn.copy()

    # Now use seasonal means
    ppn = load_run_ss(run, 'seasonal', 'precipitation_flux')  # m01s05i216
    ppn.convert_units('kg m-2 day-1')
    iris.coord_categorisation.add_season(ppn, 'time', name='clim_season')
    djf = ppn.extract(iris.Constraint(clim_season='djf'))
    mam = ppn.extract(iris.Constraint(clim_season='mam'))
    jja = ppn.extract(iris.Constraint(clim_season='jja'))
    son = ppn.extract(iris.Constraint(clim_season='son'))

    # Combine into one dictionary to loop over
    ppnseasons = dict(ann=ann, djf=djf, mam=mam, jja=jja, son=son)
    for (season, ppn) in ppnseasons.items():

        ppng = area_average(ppn, weighted=True)
        psd = ppng.collapsed('time', iris.analysis.STD_DEV)
        metrics["St dev Precip: Global " + season] = psd.data

        ptropavg = area_average(ppn, weighted=True, latitude=tropics)
        ptropsd = ptropavg.collapsed('time', iris.analysis.STD_DEV)
        metrics["St dev Precip: Tropics " + season] = ptropsd.data

        pnhavg = area_average(ppn, weighted=True, latitude=nhextrop)
        pnhsd = pnhavg.collapsed('time', iris.analysis.STD_DEV)
        metrics["St dev Precip: NH Extra-tropics " + season] = pnhsd.data

        pshavg = area_average(ppn, weighted=True, latitude=shextrop)
        pshsd = pshavg.collapsed('time', iris.analysis.STD_DEV)
        metrics["St dev Precip: SH Extra-tropics " + season] = pshsd.data

    return round_metrics(metrics)

