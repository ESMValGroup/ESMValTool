"""Module for OBS hydrocycle metrics"""

import os.path
import sys

import numpy as np

import iris
import iris.coord_categorisation

from utility.area_utils import area_average, invertmask
import hydrocycle_obs_filenames
import hydrocycle_get_landfrac


def global_obs_runoff_ann():
    """
    Compute global annual MEAN of runoff from Fekete obs, in mm/day
    It does NOT reproduce the latitudinally-binned values.

    NOTE that this is different from Pete Falloon's IDL code, which calculated
    a global TOTAL in Sverdrups.

    1 Sverdrup = 1.0e6 m3/s. 1 mm per day(year) = 0.001 m3 per m2 per day(year)
    So to convert from mm/year to Sverdrups, Pete divided by:
                (86400 seconds per day) x (360 days) x 0.001 x 1.0e6
    Pete's global total was 1.26164 Sv. A global area average is 9.639e-15 Sv
                which is equivalent to 0.8328 mm/day.

    Called from hydrocycle_obs (below) and values written to file
    along with other metrics.

    ANNUAL MEANS ONLY
    """
    # TODO 360d in docstring

    # Latitude definitions for sub-regions
    tr = (-30., 30., True, True)
    nh = (30., 90., False, True)
    ro_gl = (-60., 90., False, True)
    ro_sh = (-60., -30., False, False)

    metrics = dict()
    season = 'ann'

    # TODO local path
    obs_file = '/project/cma/clim/hydrocycle/cmp_ro.pp'
    obs = iris.load_cube(obs_file)

    obs_tot = area_average(obs, weighted=True, latitude=ro_gl)
    obstravg = area_average(obs, weighted=True, latitude=tr)
    obsnhavg = area_average(obs, weighted=True, latitude=nh)
    obsshavg = area_average(obs, weighted=True, latitude=ro_sh)

    # TODO 360d
    name = "Total Runoff Land " + season
    metrics[name] = float(obs_tot.data) / 360.

    name = "Total runoff: Tropical Land " + season
    metrics[name] = float(obstravg.data) / 360.

    name = "Total runoff: NH Extra-tropical Land " + season
    metrics[name] = float(obsnhavg.data) / 360.

    name = "Total runoff: SH Extra-tropical Land " + season
    metrics[name] = float(obsshavg.data) / 360.

    return metrics


def ppn_seasons(ppn):

    # TODO 28d months?
    # Set up constraint for 3 month period
    spans_three_months = lambda t: (t.bound[1] - t.bound[0]) > 3*28*24.0
    three_months_bound = iris.Constraint(time=spans_three_months)

    # Set up constraint for 12 month period
    spans_twelve_months = lambda t: (t.bound[1] - t.bound[0]) > 12*28*24.0
    twelve_months_bound = iris.Constraint(time=spans_twelve_months)

    # Extract time series of seasonal means
    iris.coord_categorisation.add_season(ppn, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(ppn, 'time', name='season_year')
    djf = ppn.extract(iris.Constraint(clim_season='djf'))
    mam = ppn.extract(iris.Constraint(clim_season='mam'))
    jja = ppn.extract(iris.Constraint(clim_season='jja'))
    son = ppn.extract(iris.Constraint(clim_season='son'))

    # Create time series of seasonal(annual) means, one for each season_year.
    # Use mdtol=0 (as in pp_avg) i.e. set mean to missing data if ANY
    #  contributing field is missing
    djfa1 = djf.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN, mdtol=0)
    mama1 = mam.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN, mdtol=0)
    jjaa1 = jja.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN, mdtol=0)
    sona1 = son.aggregated_by(['clim_season', 'season_year'],
                              iris.analysis.MEAN, mdtol=0)
    anna1 = ppn.aggregated_by(['season_year'], iris.analysis.MEAN, mdtol=0)

    # Remove all of the resultant 'times' which do not cover a three month
    # period (note: judged here as > 3*28 days). Similarly for 12 months
    # (since 'annual' will go from 1st Dec to 30th November)
    djfa = djfa1.extract(three_months_bound)
    mama = mama1.extract(three_months_bound)
    jjaa = jjaa1.extract(three_months_bound)
    sona = sona1.extract(three_months_bound)
    anna = anna1.extract(twelve_months_bound)

    return dict(ann=anna, djf=djfa, mam=mama, jja=jjaa, son=sona)


def ppn_stats(ppn):

    # First create timeseries of seasonal(annual) area means, then take st dev
    ppng = area_average(ppn, weighted=True)
    ppnsd = np.std(ppng.data, ddof=1)

    # For area mean of climatology, create climatology, THEN take area mean
    ppnavg = area_average(ppn, coords=['time'], mdtol=0)
    ppngavg = area_average(ppnavg, weighted=True)
    return ppnavg, ppngavg, ppnsd


def hydrocycle_obs():
    """
    Computes means (and some stdevs) of observed hydrological quantities

    NOTE that annual means run from 1st December to 30th November, as in
         the MetUM.
    NOTE that missing data tolerance is set to zero, as in pp_avg.pro and
         (I think) in ppmnvar (used for supermeans).

    Writes out a file of values.
    """

    with open('hydrocycle_obs.csv', 'w') as f:

        # Latitude definitions for sub-regions
        tr = (-30., 30., True, True)
        nh = (30., 90., False, True)
        sh = (-90., -30., True, False)

        # Read filenames from file
        filedict = hydrocycle_obs_filenames.obs_filenames()
        gpcp_file = filedict['gpcpm']
        cmap_file = filedict['cmapm']
        dasilva_file_lh = filedict['dasilvalh']
        dasilva_file_sh = filedict['dasilvash']
        nocs_file = filedict['nocss']

        # TODO STASH codes
        constr_noc_lh = iris.AttributeConstraint(STASH="m??s03i234")
        constr_noc_sh = iris.AttributeConstraint(STASH="m??s03i217")
        constr_noc_ev = iris.AttributeConstraint(STASH="m??s03i232")

        # Precipitation:  annual and seasonal means and stdev
        gpcpmon = iris.load_cube(gpcp_file)
        cmapmon = iris.load_cube(cmap_file)

        gppnseasons = ppn_seasons(gpcpmon)
        cppnseasons = ppn_seasons(cmapmon)

        # Now loop through the seasons calculating means and (some) standard
        # deviations, global and regional

        glob_seasons = ['ann', 'djf', 'mam', 'jja', 'son']

        for season in glob_seasons:

            # GPCP first
            gpcp = gppnseasons[season]
            (gpcpavg, gpcpgavg, gpcpsd) = ppn_stats(gpcp)

            # Same for CMAP
            cmap = cppnseasons[season]
            (cmapavg, cmapgavg, cmapsd) = ppn_stats(cmap)

            # Now write those values to the output csv file
            name = "Precipitation Global " + season
            value1 = float(gpcpgavg.data)
            value2 = float(cmapgavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            if season == 'ann':
                # Global, annual CLIMATOLOGY of evaporation taken as the same
                #  as precipitation
                name = "Evaporation Global ann"
                f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "St dev Precip: Global " + season
            value1 = float(gpcpsd)
            value2 = float(cmapsd)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Land only and sea only, using the climatologies we created above

            # First get land fraction on correct resolution

            landfracg = hydrocycle_get_landfrac.get_landfrac(gpcpavg)
            seafracg = invertmask(landfracg)

            landfracc = hydrocycle_get_landfrac.get_landfrac(cmapavg)
            seafracc = invertmask(landfracc)

            gpcplavg = area_average(gpcpavg, weighted=True, mask=landfracg)
            gpcpsavg = area_average(gpcpavg, weighted=True, mask=seafracg)

            cmaplavg = area_average(cmapavg, weighted=True, mask=landfracc)
            cmapsavg = area_average(cmapavg, weighted=True, mask=seafracc)

            name = "Precipitation Land " + season
            value1 = float(gpcplavg.data)
            value2 = float(cmaplavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "Precipitation Ocean " + season
            value1 = float(gpcpsavg.data)
            value2 = float(cmapsavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Tropics
            gpcptrseas = area_average(gpcp, weighted=True, latitude=tr)
            gpcptrsd = np.std(gpcptrseas.data, ddof=1)
            gpcptravga = area_average(gpcpavg, weighted=True, latitude=tr)

            cmaptrseas = area_average(cmap, weighted=True, latitude=tr)
            cmaptrsd = np.std(cmaptrseas.data, ddof=1)
            cmaptravga = area_average(cmapavg, weighted=True, latitude=tr)

            name = "Mean Precip: Tropics " + season
            value1 = float(gpcptravga.data)
            value2 = float(cmaptravga.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "St dev Precip: Tropics " + season
            value1 = float(gpcptrsd)
            value2 = float(cmaptrsd)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # N hem extratropics
            gpcpnhseas = area_average(gpcp, weighted=True, latitude=nh)
            gpcpnhsd = np.std(gpcpnhseas.data, ddof=1)
            gpcpnhavga = area_average(gpcpavg, weighted=True, latitude=nh)

            cmapnhseas = area_average(cmap, weighted=True, latitude=nh)
            cmapnhsd = np.std(cmapnhseas.data, ddof=1)
            cmapnhavga = area_average(cmapavg, weighted=True, latitude=nh)

            name = "Mean Precip: NH Extra-tropics " + season
            value1 = float(gpcpnhavga.data)
            value2 = float(cmapnhavga.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "St dev Precip: NH Extra-tropics " + season
            value1 = float(gpcpnhsd)
            value2 = float(cmapnhsd)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # S hem extratropics
            gpcpshseas = area_average(gpcp, weighted=True, latitude=sh)
            gpcpshsd = np.std(gpcpshseas.data, ddof=1)
            gpcpshavga = area_average(gpcpavg, weighted=True, latitude=sh)

            cmapshseas = area_average(cmap, weighted=True, latitude=sh)
            cmapshsd = np.std(cmapshseas.data, ddof=1)
            cmapshavga = area_average(cmapavg, weighted=True, latitude=sh)

            name = "Mean Precip: SH Extra-tropics " + season
            value1 = float(gpcpshavga.data)
            value2 = float(cmapshavga.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "St dev Precip: SH Extra-tropics " + season
            value1 = float(gpcpshsd)
            value2 = float(cmapshsd)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Other fields: annual and seasonal, global and regional means

            # DaSilva: surface LH fluxes over OCEAN only

            model_file = '{}{}.pp'.format(dasilva_file_lh, season)
            dasilval = iris.load_cube(model_file)
            dasilvalg = area_average(dasilval, weighted=True)
            dasilvaltravg = area_average(dasilval, weighted=True, latitude=tr)
            dasilvalnhavg = area_average(dasilval, weighted=True, latitude=nh)
            dasilvalshavg = area_average(dasilval, weighted=True, latitude=sh)

            # NOCS2.0: surface SH, LH fluxes and evaporation over OCEAN only

            # Latent heat (ocean)
            model_file = '{}{}.pp'.format(nocs_file, season)
            nocslh = iris.load_cube(model_file, constr_noc_lh)
            nocslhg = area_average(nocslh, weighted=True)
            nocslhtravg = area_average(nocslh, weighted=True, latitude=tr)
            nocslhnhavg = area_average(nocslh, weighted=True, latitude=nh)
            nocslhshavg = area_average(nocslh, weighted=True, latitude=sh)

            name = "Latent heat flux Ocean " + season
            value1 = float(dasilvalg.data)
            value2 = float(nocslhg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Sfc Latent heat: Tropical Ocean " + season
            value1 = float(dasilvaltravg.data)
            value2 = float(nocslhtravg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Sfc Latent heat: NH Extra-tropical Ocean " + season
            value1 = float(dasilvalnhavg.data)
            value2 = float(nocslhnhavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Sfc Latent heat: SH Extra-tropical Ocean " + season
            value1 = float(dasilvalshavg.data)
            value2 = float(nocslhshavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Sensible heat (ocean)
            # DaSilva
            model_file = '{}{}.pp'.format(dasilva_file_sh, season)
            dasilvas = iris.load_cube(model_file)
            dasilvasg = area_average(dasilvas, weighted=True)
            dasilvastravg = area_average(dasilvas, weighted=True, latitude=tr)
            dasilvasnhavg = area_average(dasilvas, weighted=True, latitude=nh)
            dasilvasshavg = area_average(dasilvas, weighted=True, latitude=sh)

            # NOCS2.0
            model_file = '{}{}.pp'.format(nocs_file, season)
            nocssh = iris.load_cube(model_file, constr_noc_sh)
            nocsshg = area_average(nocssh, weighted=True)
            nocsshtravg = area_average(nocssh, weighted=True, latitude=tr)
            nocsshnhavg = area_average(nocssh, weighted=True, latitude=nh)
            nocsshshavg = area_average(nocssh, weighted=True, latitude=sh)

            if season == 'ann':
                # Global, annual mean SH flux over ocean also from literature
                name = "Sensible heat flux Ocean ann"
                # Trenberth et al. (2009)
                value2 = max(float(nocsshg.data), float(dasilvasg.data), 12.0)
                # SOC (Josey et al 1999)
                valuel = min(float(nocsshg.data), float(dasilvasg.data), 7.0)
                f.write('{0},{1},{2}\n'.format(name, valuel, value2))
            else:
                name = "Sensible heat flux Ocean " + season
                value1 = float(nocsshg.data)
                value2 = float(dasilvasg.data)
                f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "Mean Sfc Sensible heat: Tropical Ocean " + season
            value1 = float(nocsshtravg.data)
            value2 = float(dasilvastravg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Sfc Sensible heat: NH Extra-tropical Ocean " + season
            value1 = float(nocsshnhavg.data)
            value2 = float(dasilvasnhavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Sfc Sensible heat: SH Extra-tropical Ocean " + season
            value1 = float(nocsshshavg.data)
            value2 = float(dasilvasshavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Evaporation (ocean)
            model_file = '{}{}.pp'.format(nocs_file, season)
            nocsev = iris.load_cube(model_file, constr_noc_ev)
            nocsev.units = 'kg m-2 s-1'
            nocsev.convert_units('kg m-2 day-1')
            nocsevg = area_average(nocsev, weighted=True)
            nocsevtravg = area_average(nocsev, weighted=True, latitude=tr)
            nocsevnhavg = area_average(nocsev, weighted=True, latitude=nh)
            nocsevshavg = area_average(nocsev, weighted=True, latitude=sh)

            if season == 'ann':
                # Global, annual mean evaporation over ocean also from
                #  literature
                name = "Evaporation Ocean ann"
                value2 = float(nocsevg.data)
                valuel = 2.97     # Yu 2007
                f.write('{0},{1},{2}\n'.format(name, valuel, value2))
            else:
                name = "Evaporation Ocean " + season
                value1 = float(nocsevg.data)
                value2 = float(nocsevg.data)
                f.write('{0},{1},{2}\n'.format(name, value1, value2))

            name = "Mean Evap: Tropical Ocean " + season
            value1 = float(nocsevtravg.data)
            value2 = float(nocsevtravg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Evap: NH Extra-tropical Ocean " + season
            value1 = float(nocsevnhavg.data)
            value2 = float(nocsevnhavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))
            name = "Mean Evap: SH Extra-tropical Ocean " + season
            value1 = float(nocsevshavg.data)
            value2 = float(nocsevshavg.data)
            f.write('{0},{1},{2}\n'.format(name, value1, value2))

            # Now write the P-E values where available
            if season == 'ann':
                # Set these to 0.0 explicitly so that obs have P-E balance
                value1 = 0.0
                value2 = 0.0
                name = "P minus E Global ann"
                f.write('{0},{1},{2}\n'.format(name, value1, value2))

                # Use range of values in Mueller et al. (2011) for land evap
                value1a = max(float(gpcplavg.data), float(cmaplavg.data))
                value1b = min(1.36000, 1.76000)
                value1 = value1a-value1b
                value2a = min(float(gpcplavg.data), float(cmaplavg.data))
                value2b = max(1.36000, 1.76000)
                value2 = value2a-value2b
                name = "P minus E Land ann"
                f.write('{0},{1},{2}\n'.format(name, value2, value1))

            if season == 'ann':
                # Use a fixed value of 2.97 from Yu 2007 as second value
                #  of ocean evap
                value1a = max(float(gpcpsavg.data), float(cmapsavg.data))
                value1b = min(float(nocsevg.data), 2.97)
                value1 = value1a-value1b
                value2a = min(float(gpcpsavg.data), float(cmapsavg.data))
                value2b = max(float(nocsevg.data), 2.97)
                value2 = value2a-value2b
            else:
                value1a = max(float(gpcpsavg.data), float(cmapsavg.data))
                value1b = float(nocsevg.data)
                value1 = value1a-value1b
                value2a = min(float(gpcpsavg.data), float(cmapsavg.data))
                value2b = float(nocsevg.data)
                value2 = value2a-value2b

            name = "P minus E Ocean " + season
            f.write('{0},{1},{2}\n'.format(name, value2, value1))

        # end of loop through seasons

        # Calculate annual global runoff metrics and write these to this
        #  file too

        runoff_metrics = global_obs_runoff_ann()
        for name, value in runoff_metrics.items():
            f.write('{0},{1},{1}\n'.format(name, value))

        # Now add other annual mean values from literature

        # Mueller et al. 2011 from range of obs
        name = "Evaporation Land ann"
        values = [1.36, 1.76]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Range of estimates quoted in Trenberth et al. (2009)
        name = "Latent heat flux Global ann"
        values = [80.0, 83.0]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Jimenez et al. (2011); excludes Antarctic, Greenland
        # Also Trenberth et al. (2009)
        name = "Latent heat flux Land ann"
        values = [38.0, 51.0]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Range of estimates quoted in Trenberth et al. (2009)
        name = "Sensible heat flux Global ann"
        values = [15.7, 18.9]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Jimenez et al. (2011), Trenberth et al. (2009) estimates
        name = "Sensible heat flux Land ann"
        values = [26.0, 47.0]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Trenberth (2011) (SSM/I)
        name = "Total water vapour Global ann"
        values = [24.2000, 24.2000]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Trenberth (2011) (SSM/I)
        name = "Total water vapour Land ann"
        values = [18.5000, 18.5000]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

        # Trenberth (2011) (SSM/I)
        name = "Total water vapour Ocean ann"
        values = [26.6000, 26.6000]
        f.write('{0},{1[0]},{1[1]}\n'.format(name, values))

if __name__ == '__main__':
    hydrocycle_obs()
