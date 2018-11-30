"""
Compute Seasonal Jets.

This module is part of the ACSIS indicators computations. It is an
ESMValTool executable that can be called directly from the command line:
example:

acsisjet -m ERA-Interim_T3M_ua.nc -o FULL_PATH_TO_ACSIS_OUTPUT

Data availability:
------------------
This module uses data that must be supplied as comand line arguments
(see below for options). ACSIS data on CEDA-JASMIN can be found in:

group_workspaces
jasmin4
esmeval
acsis_data

or for OBS data in

group_workspaces
jasmin4
esmeval
obsdata

Specifically for ERA-Interim:
In obsdata - Tier3 - ERA-Interim -
OBS_ERA-Interim_reanaly_1_T3M_ua_197901-201412.nc

Usage:
------
Arguments:

-m --era-file     [OPTIONAL]   full path to the ERA file (netCDF);
-o --output       [OPTIONAL]   full path to output dir;
                               DEFAULT $HOME/$USER/ACSIS-JET;
-l --log-level    [OPTIONAL]   logging level;
                               OPTIONS: 'debug', 'info', 'warning', 'error';
                               DEFAULT: 'info';

Contact:
--------
Ed Hawkins, UREAD, e.hawkins@reading.ac.uk
Valeriu Predoi, UREAD, valeriu.predoi@ncas.ac.uk
First working version: November 2018
"""
import argparse
from datetime import datetime
import logging
import os
import numpy as np
import yaml
import iris
from esmvaltool.preprocessor import extract_region, \
    extract_levels, extract_season
from esmvaltool.utils.acsis._utils import _save_cubes, _set_logger

seasons = ['DJF', 'MAM', 'JJA', 'SON']

# set up logging
logger = logging.getLogger(__name__)

# print the header
HEADER = r"""
___________________________________________________________________________

 ESMValTool ACSIS Indicators Computation: ERA jet speeds and latitudes
___________________________________________________________________________

This module of the ACSIS indicators computation computes the seasonal values
for jet speeds and latitudes.

""" + __doc__


def get_args():
    """Define the `esmvaltool` command line."""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-m', '--era-file', type=str, nargs='+', help='ERA files')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default=os.path.join(os.environ['HOME'], 'ACSIS-JET'),
        help='Output directory [FULLPATH]')
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        choices=['debug', 'info', 'warning', 'error'])
    args = parser.parse_args()
    return args


def _make_dir_tree(base_dir):
    """Create output directory three for storing output."""
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    speed_dir = os.path.join(base_dir, 'jet_speeds')
    if not os.path.exists(speed_dir):
        os.makedirs(speed_dir)
    lats_dir = os.path.join(base_dir, 'jet_latitudes')
    if not os.path.exists(lats_dir):
        os.makedirs(lats_dir)
    return speed_dir, lats_dir


def _extract_u850(era_file):
    """
    Get jet speed and jet latitude.

    Extract eastern wind u at 850 hPa 15-75N lat.
    Extract mean 0-60W lon. Mean on LON. Extract season.
    Return each season's cube in a dictionary.
    """
    seasonal_data = {}
    # load ERA data
    era_cube = iris.load_cube(era_file)
    # extract 0-60W lon; 15-75N lat region
    era_cube = extract_region(era_cube, 0., 60., 15., 75.)
    # extract 850 hPa
    era_cube = extract_levels(era_cube, 85000., 'linear')
    # collapse-mean on lon
    era_cube = era_cube.collapsed(['longitude'], iris.analysis.MEAN)
    # extract seasons
    for season in seasons:
        seasonal_data[season] = extract_season(era_cube, season)
    return seasonal_data


def _get_jets(era_seasonal_data):
    """Take seasonal dictionary and get jets dicts (speeds and lats)."""
    jet_speeds = {}
    jet_lats = {}
    for season in seasons:
        era_cube = era_seasonal_data[season]
        # get the day when u850 is max ie max over t-axis(axis=0)
        jet_speed = np.amax(era_cube.data, axis=1)
        jet_speeds[season] = jet_speed
        t_ax, l_ax = era_cube.data.shape
        # get the jet latitudes
        jet_lat = np.zeros((t_ax, ))
        for t_idx in range(t_ax):
            for l_idx in range(l_ax):
                if era_cube.data[t_idx, l_idx] == jet_speed[t_idx]:
                    jet_lat[t_idx] = era_cube.coord('latitude').points[l_idx]
        jet_lats[season] = jet_lat
    return jet_speeds, jet_lats


def main():
    """Run the the meat of the jet code."""
    logger.info('Executing code...')
    args = get_args()

    # base dir needs timestamp
    now = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    base_output = '_'.join((args.output, now))

    # build output tree
    speed_dir, lats_dir = _make_dir_tree(base_output)

    # set logging
    log_level = args.log_level
    _set_logger(logging, base_output, 'acsisjets.log', log_level)
    logger.info(HEADER)
    logger.info("Running main function...")

    # get input files
    era_file = args.era_file

    # compute and save
    season_dict = _extract_u850(era_file)
    jet_speeds, jet_lats = _get_jets(season_dict)

    # save output
    # netCDF
    for season in seasons:
        _save_cubes('seasons', jet_speeds[season], 'u850',
                    speed_dir, '_'.join((season, 'u850_jet-speed.nc')))
        _save_cubes('seasons', jet_lats[season], 'latitude',
                    lats_dir, '_'.join((season, 'latitude_jet-latitude.nc')))

    # yaml files; for easy access
    # make them strings for output
    for season in seasons:
        jet_speeds[season] = str(list(jet_speeds[season]))
        jet_lats[season] = str(list(jet_lats[season]))
    jets_file = os.path.join(speed_dir, 'ERA-jets.yml')
    lats_file = os.path.join(lats_dir, 'ERA-jet-lats.yml')
    logger.info('Writing seasonal jet speeds in %s', jets_file)
    with open(jets_file, 'w') as jet_yamlfile:
        yaml.safe_dump(jet_speeds, jet_yamlfile)
    logger.info('Writing seasonal jet lats in %s', jets_file)
    with open(lats_file, 'w') as lats_yamlfile:
        yaml.safe_dump(jet_lats, lats_yamlfile)

    logger.info('Finished computing ACSIS jets!')


if __name__ == '__main__':
    main()
