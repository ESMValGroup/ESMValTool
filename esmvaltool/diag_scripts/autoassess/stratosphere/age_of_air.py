"""Stratospheric age-of-air assessment code."""
import datetime
import logging
import os
import warnings

import iris
import iris.analysis as iai
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.autoassess.loaddata import load_run_ss

from .strat_metrics_1 import weight_lat_ave

logger = logging.getLogger(__name__)


# Constant for number of seconds in a 360 day calendar year
# Wrong if gregorian calendar!
RSECS_PER_360DAY_YEAR = float(60 * 60 * 24 * 360)

# What is the source of the reference data???
# Diag 1
# SF6 based data
AGE_YRS = [
    0.25, 0.43, 0.80, 1.36, 1.67, 2.25, 2.46, 2.83, 3.06, 3.23, 3.59, 3.85,
    4.06, 4.03, 4.11, 4.07, 3.95
]
ZSF6_KM = [
    17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
    29.0, 30.0, 31.0, 32.0, 33.0
]
# CO2 based data
AGE_YRS2 = [
    0.37, 0.39, 0.59, 0.44, 0.54, 1.24, 1.62, 2.04, 2.61, 2.77, 2.82, 3.13,
    3.47, 3.50, 3.40, 3.53, 3.92, 3.73
]
ZCO2_KM = [
    16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0,
    28.0, 29.0, 30.0, 31.0, 32.0, 33.0
]

# Diag 2
# SF6 based data
AGE2_YRS = [
    -0.19, -0.01, 0.18, 0.36, 0.50, 0.77, 1.12, 1.52, 1.95, 2.50, 3.10, 3.66,
    4.02, 4.40, 4.68, 4.87, 4.95, 5.01, 5.10, 5.16, 5.23, 5.32, 5.40, 5.19,
    5.57, 5.90
]
Z2_KM = [
    10.16, 11.09, 12.06, 13.08, 14.05, 15.05, 16.07, 17.06, 18.05, 18.99,
    20.03, 21.04, 22.05, 23.03, 24.03, 25.04, 26.03, 27.01, 27.97, 29.01,
    30.03, 31.09, 31.98, 32.68, 34.13, 40.10
]
# CO2 based data
AGE2_YRS2 = [
    -0.36, 0.59, 0.78, 0.67, 0.79, 0.98, 1.25, 1.49, 1.88, 2.42, 3.03, 3.71,
    4.16, 4.48, 4.56, 4.71, 4.75, 4.85, 4.78, 4.79, 4.75, 4.82, 4.74, 4.93,
    4.78, 4.80
]
Z2_KM2 = [
    10.18, 11.06, 11.90, 13.02, 14.04, 15.04, 16.08, 17.07, 18.03, 19.06,
    20.09, 21.11, 22.04, 23.03, 23.98, 25.01, 26.00, 26.98, 27.91, 28.93,
    30.01, 31.10, 32.10, 33.01, 33.99, 35.07
]


def calculate_analysis_years(run):
    """Calculate years."""
    # 1) Discard first 10 years of run.
    analysis_start_year = int(run['start']) + 10
    analysis_end_year = int(run['start']) + int(run['nyear'])

    # 2) If at least 2 years remain, then age of air can be assessed.
    if int(run['nyear']) >= 12:

        # 3) If at least 5 years remain, then use the last 5 years.
        #    If less than 5 years remain, then use all remaining years.
        if int(run['nyear']) >= 15:
            analysis_start_year = analysis_end_year - 5

    else:
        raise ValueError()

    analysis_start_dt = datetime.datetime(analysis_start_year - 1, 12, 1)
    analysis_end_dt = datetime.datetime(analysis_end_year, 11, 1)

    return analysis_start_dt, analysis_end_dt


def age_of_air(run):
    """Calculate the age of air metrics."""
    # Create metrics dictionary with MDI incase age of air
    # diagnostics not available
    metrics = {
        'RMS error: tropical Age of Air': -10000.,
        'RMS error: NH midlatitude Age of Air': -10000.
    }

    try:
        # Set up to only run for 5 year period
        analysis_start_dt, analysis_end_dt = calculate_analysis_years(run)
        constraint = dict(
            from_dt=analysis_start_dt, to_dt=analysis_end_dt, lbproc=128)
        # Calculate age of air metrics if appropriate diagnostic available
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', '.*orography.*', UserWarning)
            agecube = load_run_ss(run, 'monthly', 'age_of_stratospheric_air',
                                  **constraint)  # m01s34i150
    except iris.exceptions.ConstraintMismatchError:
        logger.warning('Age of air fields absent.  Skipping this diagnostic.')
    except ValueError:
        logger.warning("Run length < 12 years: Can't assess age of air")
    else:
        # Create time/zonal means of age data
        agecube = agecube.collapsed(['longitude', 'time'], iris.analysis.MEAN)
        # Convert units of data from seconds to years
        agecube.data /= RSECS_PER_360DAY_YEAR
        # Create latitude bounds for area-weighted averaging
        agecube.coord('latitude').guess_bounds()
        # Convert units of height from metres to km
        agecube.coord('level_height').convert_units('km')

        # Calculate area-weighted means for tropics
        trop_cons = iris.Constraint(latitude=lambda lat: -10 <= lat <= 10)
        diag1 = weight_lat_ave(agecube.extract(trop_cons))
        diag1.var_name = 'tropical_age_of_air'

        # Calculate area-weighted means for mid-latitudes
        mlat_cons = iris.Constraint(latitude=lambda lat: 35 <= lat <= 45)
        diag2 = weight_lat_ave(agecube.extract(mlat_cons))
        diag2.var_name = 'midlat_age_of_air'

        # Write age of air data to CWD
        outfile = '{0}_age_of_air_{1}.nc'
        cubelist = iris.cube.CubeList([diag1, diag2])
        with iris.FUTURE.context(netcdf_no_unlimited=True):
            iris.save(cubelist, outfile.format(run['runid'], run.period))

        # Calculate metrics
        diag1sf6 = iai.Linear(diag1, [('level_height', ZSF6_KM)])
        diag1co2 = iai.Linear(diag1, [('level_height', ZCO2_KM)])
        diag2sf6 = iai.Linear(diag2, [('level_height', Z2_KM)])
        diag2co2 = iai.Linear(diag2, [('level_height', Z2_KM2)])
        metric1sf6 = np.sqrt(np.mean((diag1sf6.data - AGE_YRS)**2))
        metric1co2 = np.sqrt(np.mean((diag1co2.data - AGE_YRS2)**2))
        metric2sf6 = np.sqrt(np.mean((diag2sf6.data - AGE2_YRS)**2))
        metric2co2 = np.sqrt(np.mean((diag2co2.data - AGE2_YRS2)**2))
        trop_age = (metric1sf6 + metric1co2) / 2.
        midl_age = (metric2sf6 + metric2co2) / 2.

        # Add metrics to dictionary
        metrics['RMS error: tropical Age of Air'] = float(trop_age)
        metrics['RMS error: NH midlatitude Age of Air'] = float(midl_age)

    return metrics


def multi_age_plot(run):
    """
    Plot results.

    This function is plotting the results of the function age_of_air for each
    run against observations.
    """
    # Run age_of_air for each run.
    # Age_of_air returns metrics and writes results into an *.nc in the current
    # working directory.
    # To make this function independent of the previous call to age_of_air,
    # age_of_air is run again for each run in this function
    #
    # This behaviour is due to the convention that only metric_functions can
    # return metric values, multi_functions are supposed to
    # only produce plots (see __init__.py).

    ######################################

    # Set up constraints to deal with loading data
    trop_cons = iris.Constraint(
        cube_func=lambda c: c.var_name == 'tropical_age_of_air')
    midl_cons = iris.Constraint(
        cube_func=lambda c: c.var_name == 'midlat_age_of_air')

    # Set up generic input file name
    infile = '{0}_age_of_air_{1}.nc'

    # Create control filename
    cntlfile = infile.format(run['suite_id1'], run['period'])

    # Create experiment filename
    exptfile = infile.format(run['suite_id2'], run['period'])

    # If no control data then stop ...
    if not os.path.exists(cntlfile):
        logger.warning('Age of air for control absent. skipping ...')
        return

    # Create tropics plot
    fig = plt.figure()
    ax1 = plt.gca()
    # Plot OBS
    plt.plot(
        AGE_YRS,
        ZSF6_KM,
        linestyle='-',
        marker='s',
        color='black',
        label='SF6 obs')
    plt.plot(
        AGE_YRS2,
        ZCO2_KM,
        linestyle='-',
        marker='D',
        color='black',
        label='CO2 obs')
    # Plot control
    diag = iris.load_cube(cntlfile, trop_cons)
    levs = diag.coord('level_height').points
    plt.plot(diag.data, levs, label=run['suite_id1'])
    # Plot experiment
    if os.path.exists(exptfile):
        diag = iris.load_cube(exptfile, trop_cons)
        levs = diag.coord('level_height').points
        plt.plot(diag.data, levs, label=run['suite_id2'])
    ax1.set_title('Tropical mean age profile (10S-10N)')
    ax1.set_xlabel('Mean age (years)')
    ax1.set_ylabel('Height (km)')
    ax1.set_ylim(16, 34)
    ax1.legend(loc='upper left')
    fig.savefig('age_tropics.png')
    plt.close()

    # Create midlats plot
    fig = plt.figure()
    ax1 = plt.gca()
    # Plot OBS
    plt.plot(
        AGE2_YRS,
        Z2_KM,
        linestyle='-',
        marker='s',
        color='black',
        label='SF6 obs')
    plt.plot(
        AGE2_YRS2,
        Z2_KM2,
        linestyle='-',
        marker='D',
        color='black',
        label='CO2 obs')
    # Plot control
    diag = iris.load_cube(cntlfile, midl_cons)
    levs = diag.coord('level_height').points
    plt.plot(diag.data, levs, label=run['suite_id1'])
    # Plot experiment
    if os.path.exists(exptfile):
        diag = iris.load_cube(exptfile, midl_cons)
        levs = diag.coord('level_height').points
        plt.plot(diag.data, levs, label=run['suite_id2'])
    ax1.set_title('Midlatitude mean age profile (35N-45N)')
    ax1.set_xlabel('Mean age (years)')
    ax1.set_ylabel('Height (km)')
    ax1.set_ylim(16, 34)
    ax1.legend(loc='upper left')
    fig.savefig('age_midlatitudes.png')
    plt.close()
