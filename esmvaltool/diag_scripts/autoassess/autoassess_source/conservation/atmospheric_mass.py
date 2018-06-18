"""module with routines to calculate atmospheric mass conservation"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

import iris
from .loaddata import load_run_ss
from .area_utils import area_average


def global_atmos_mass_conservation(run):
    """
    Purpose

    Compute and plot HadGEM3 dry and wet mass conservation (trend) metrics
    for AutoAssess.

    Inputs:
    Annual mean global model data for stash codes: 30403, 30404

    Notes:
        (1) The mass fixer should normally be applied in the run concerned, in
            which case dry mass should be conserved to numerical accuracy and
            the computed trend should not differ significantly from zero.
        (2) Wet mass is not formally conserved, but a significant trend could
            be indicative of unacceptable climate drift (i.e. a non-stationary
            state, for some reason), which is worth monitoring and could be
            related directly to surface P-E fluxes in principle.

    Program initially written in IDL by Tim Johns.
    23/12/2014 Porting of program by Jose Rodriguez.

    """
    metrics = dict()
    g = 9.80665  # standard gravity acceleration in m s^-2

    # calculate ATMOSPHERIC MASS METRIC
    # read atmospheric dry mass and mass data
    # VPREDOI
    # temporarily replacing these crazy variables with standard ones
    # so I can test the code flow
    # dry = load_run_ss(run, 'annual', 'm01s30i403')
    # TOTAL COLUMN DRY MASS  RHO GRID
    # wet = load_run_ss(run, 'annual', 'atmosphere_mass_per_unit_area')
    # m01s30i404
    dry = load_run_ss(run, 'monthly', 'eastward_wind', lbproc=192)
    wet = load_run_ss(run, 'monthly', 'eastward_wind', lbproc=192)

    # calculating global means
    # Get time series of global means:
    dryg = g * area_average(dry, weighted=True)
    wetg = g * area_average(wet, weighted=True)

    # COMPUTE metrics (first subtracting mean over the data period)
    #    Best fit linear trend and inter-annual s.d. for dry and wet mass
    #    [Note: One could also compute the trend significance | H0(trend=0.0)]

    # Compute metrics

    # The data used here are time series of dry and wet mass annual means.
    # The time mean is calculated from 1/12/year to 1/12/year+1 .
    # On the other hand, time coordinate units of cube are hours since
    # 1/1/1970, where the lower bound of coordinate is the initial
    # period of time meaning.
    # TODO gregorian calendar
    # Convert cube time units to year of annual mean, following Hadley Centre's
    # convention of naming that year as the lower bound of meaning:

    # Should be able to do this using netcdf dataetime?
    time_years = 1970 + (dryg.coord('time').bounds[:, 0] / 720. + 1.0) / 12
    ctime_years = ['{}'.format(int(year)) for year in time_years]

    deldry = dryg - dryg.collapsed('time', iris.analysis.MEAN)
    coeffs = np.polyfit(time_years, deldry.data, 1)
    slope_dry = coeffs[0]
    sd_dry = np.std(deldry.data)

    delwet = wetg - wetg.collapsed('time', iris.analysis.MEAN)
    coeffs = np.polyfit(time_years, delwet.data, 1)
    slope_wet = coeffs[0]
    sd_wet = np.std(delwet.data)

    # dry mass trend (Pa/sec)
    metric_dry_1 = slope_dry
    # log10(e-folding time, years)
    metric_dry_2 = np.log10(1.0e5 / np.abs(metric_dry_1))

    # wet mass trend (Pa/sec)
    metric_wet_1 = slope_wet
    # log10(e-folding time, years)
    metric_wet_2 = np.log10(1.0e5 / np.abs(metric_wet_1))

    # assign metric value:

    # Dry mass: 1/log10(e-folding time[years])
    metric_name = 'atmos. dry mass drift: 1/log!D10!N(e-folding time/years)'
    metrics[metric_name] = 1.0 / metric_dry_2

    # produce extra plots:

    xtickFormatter = FormatStrFormatter('%4d')
    ytickFormatter = FormatStrFormatter('%4.1e')
    titl1_temp = '{0}   Global {1} mass (deviation from time mean, Pa)'
    # TODO unit from cube
    # VPREDOI
    # see below at use of titl2_temp
    # titl2_temp = '$log_{{10}}$(e-fold time/year):{0:6.2f}s.d.:{1:9.2e} Pa'
    # TODO unit from cube?

    plt.figure(figsize=(8.27, 11.69))

    top_ax = plt.subplot(2, 1, 1)
    top_ax.xaxis.set_major_formatter(xtickFormatter)
    top_ax.yaxis.set_major_formatter(ytickFormatter)
    titl1 = titl1_temp.format(run['runid'], 'dry')
    # VPREDOI
    # Python3: not liking np array formatting
    # titl2 = titl2_temp.format(metric_dry_2, sd_dry)
    titl2 = 'dry titties'
    plt.plot(time_years, deldry.data, linewidth=2, color='black')
    plt.xticks(time_years)
    top_ax.set_xticklabels(
        ctime_years, ha='left', fontsize='x-small', rotation=315)
    plt.xlim([np.min(time_years), np.max(time_years)])
    plt.ylim([np.min(deldry.data) - 1.0e-4, np.max(deldry.data) + 1.0e-4])
    plt.title(titl1 + '\n' + titl2)
    plt.axhline(0.0, linestyle=':', color='black')

    bottom_ax = plt.subplot(2, 1, 2)
    bottom_ax.xaxis.set_major_formatter(xtickFormatter)
    bottom_ax.yaxis.set_major_formatter(ytickFormatter)
    titl1 = titl1_temp.format(run['runid'], 'wet')
    # VPREDOI
    # same array formatting issue for Python 3
    # titl2 = titl2_temp.format(metric_wet_2, sd_wet)
    titl2 = 'wet titties'
    plt.plot(time_years, delwet.data, linewidth=2, color='black')
    plt.xticks(time_years)
    bottom_ax.set_xticklabels(
        ctime_years, ha='left', fontsize='x-small', rotation=315)
    plt.xlim([np.min(time_years), np.max(time_years)])
    plt.ylim([np.min(delwet.data) - 1.0, np.max(delwet.data) + 1.0])
    plt.title(titl1 + '\n' + titl2)
    plt.axhline(0.0, linestyle=':', color='black')

    plt.savefig(run['runid'] + '_atmospheric_mass_conservation.png')

    return metrics
