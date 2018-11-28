"""
Compute MOC annual average and Greenland - Iceland airp DJF means.

This module is part of the ACSIS indicators computations. It is an
ESMValTool executable that can be called directly from the command line:
example:

acsismoc -m moc_transports.nc -s hadslp2r.nc -o FULL_PATH_TO_ACSIS_OUTPUT

Data availability:
------------------
This module uses data that must be supplied as comand line arguments
(see below for options). ACSIS data on CEDA-JASMIN can be found in:

group_workspaces
jasmin4
esmeval
acsis_data

with subdirs labelled by data creation date.

Usage:
------
Arguments:

-m --moc-file     [OPTIONAL]   full path to the moc_transport file (netCDF);
-s --hadslp-file  [OPTIONAL]   full path to the hadslp file (netCDF);
-o --output       [OPTIONAL]   full path to output dir;
                               DEFAULT $HOME/$USER/ACSIS;
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
import collections
import numpy as np
import iris
from esmvaltool.preprocessor import extract_region, extract_season

# decalre the global variables for analysis
# TODO Ed Hawkins - check these
MOC_VARIABLE = 'moc_mar_hc10'
AIRPRESS_VARIABLE = 'UM_0_fc8_vn405'
SEASON = 'DJF'

# set up logging
logger = logging.getLogger(__name__)

# print the header
HEADER = r"""
___________________________________________________________________________

 ESMValTool ACSIS Indicators Computation: annual MOC and Greenland/Iceland
___________________________________________________________________________

This module of the ACSIS indicators computation computes the annual mean
for MOC and Greenland - Iceland DJF mean difference for air pressure.

""" + __doc__


def get_args():
    """Define the `esmvaltool` command line."""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-m', '--moc-file', type=str, nargs='+', help='MOC files')
    parser.add_argument(
        '-s', '--hadslp-file', type=str, nargs='+', help='hadslp files')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default=os.path.join(os.environ['HOME'], 'ACSIS'),
        help='Output directory [FULLPATH]')
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        choices=['debug', 'info', 'warning', 'error'])
    args = parser.parse_args()
    return args


def _set_logger(logging, out_dir, log_file, log_level):
    # set logging for screen and file output
    root_logger = logging.getLogger()
    out_fmt = "%(asctime)s %(levelname)-8s %(name)s,%(lineno)s\t%(message)s"
    logging.basicConfig(
        filename=os.path.join(out_dir, log_file),
        filemode='a',
        format=out_fmt,
        datefmt='%H:%M:%S',
        level=logging.DEBUG)
    root_logger.setLevel(log_level.upper())
    logfmt = logging.Formatter(out_fmt)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logfmt)
    root_logger.addHandler(console_handler)


def _make_dir_tree(base_dir):
    """Create output directory three for storing output."""
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    moc_dir = os.path.join(base_dir, 'moc_annual')
    if not os.path.exists(moc_dir):
        os.makedirs(moc_dir)
    greenland_iceland_dir = os.path.join(base_dir, 'greenland_iceland')
    if not os.path.exists(greenland_iceland_dir):
        os.makedirs(greenland_iceland_dir)

    return moc_dir, greenland_iceland_dir


def _extract_yearly_data(cube):
    """Return yearly data."""
    # get time cells
    time_cells = [
        cube.coord('time').units.num2date(cell.point)
        for cell in cube.coord('time').cells()
    ]

    # extract years info
    all_years = []
    for date_obj in time_cells:
        real_date = datetime(date_obj.year, date_obj.month, date_obj.day, 0, 0)
        all_years.append(real_date.year)

    # compute means
    yearly_cube = {}
    measurements = {}
    for year in set(all_years):
        idx1 = all_years.index(year)
        idx2 = idx1 + len([s for s in all_years if s == year])
        measurements[str(year)] = idx2 - idx1
        yearly_cube[str(year)] = cube[idx1:idx2 - 1]

    return yearly_cube, measurements


def _yearly_mean(data_file, output_dir, variable, var_constraint):
    """Process needed data."""
    cube = iris.load(data_file, constraints=var_constraint)[0]
    datas, no_measurements = _extract_yearly_data(cube)
    data_dict = {}
    for y_r, ydat in datas.items():
        logger.info('For year %s : %s number of measurements: %s', y_r,
                    variable, no_measurements[y_r])
        moc_yearly_mean = np.ma.mean(ydat.data)
        logger.info('For year %s : %s yearly average: %s', y_r, variable,
                    moc_yearly_mean)
        data_dict[float(y_r)] = moc_yearly_mean

    # assemble output cube and save it
    data_dict = collections.OrderedDict(sorted(data_dict.items()))
    times = iris.coords.DimCoord(
        np.array([f for f in data_dict.keys()]),
        standard_name='time',
        units='years')
    cspec = [(times, 0)]
    stats_cube = iris.cube.Cube(
        np.array([f for f in data_dict.values()]),
        dim_coords_and_dims=cspec,
        long_name=variable)
    cube_name = '_'.join((variable, 'AnnualMean.nc'))
    logger.info('Saving %s', os.path.join(output_dir, cube_name))
    iris.save(stats_cube, os.path.join(output_dir, cube_name))


def _djf_greenland_iceland(data_file, output_dir, var_constraint, season):
    """Get the DJF mean for Greenland-Iceland."""
    season = SEASON
    cube = iris.load(data_file, constraints=var_constraint)[0]
    greenland_map = extract_region(cube, 25., 35., 30., 40.)
    iceland_map = extract_region(cube, 15., 25., 60., 70.)
    greenland = greenland_map.collapsed(['longitude', 'latitude'],
                                        iris.analysis.MEAN)
    iceland = iceland_map.collapsed(['longitude', 'latitude'],
                                    iris.analysis.MEAN)
    # get cubes of interest
    greenland_djf = extract_season(greenland, season)
    iceland_djf = extract_season(iceland, season)
    diff = greenland - iceland
    season_geo_diff = extract_season(diff, season)
    # save to disk
    logger.info('%s: Saving Greenland in %s', AIRPRESS_VARIABLE, output_dir)
    iris.save(
        greenland_map,
        os.path.join(output_dir,
                     ('').join(['Greenland_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Iceland in %s', AIRPRESS_VARIABLE, output_dir)
    iris.save(
        iceland_map,
        os.path.join(output_dir,
                     ('').join(['Iceland_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Greenland Mean in %s', AIRPRESS_VARIABLE,
                output_dir)
    iris.save(
        greenland,
        os.path.join(output_dir,
                     ('').join(['Greenland_Mean_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Iceland Mean in %s', AIRPRESS_VARIABLE, output_dir)
    iris.save(
        iceland,
        os.path.join(output_dir,
                     ('').join(['Iceland_Mean_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Greenland Mean %s in %s', AIRPRESS_VARIABLE,
                season, output_dir)
    iris.save(
        greenland_djf,
        os.path.join(output_dir, ('').join(
            ['Greenland_Mean_DJF_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Iceland Mean %s in %s', AIRPRESS_VARIABLE, season,
                output_dir)
    iris.save(
        iceland_djf,
        os.path.join(output_dir, ('').join(
            ['Iceland_Mean_DJF_', AIRPRESS_VARIABLE, '.nc'])))
    logger.info('%s: Saving Greenland - Iceland Mean %s in %s',
                AIRPRESS_VARIABLE, season, output_dir)
    iris.save(
        season_geo_diff,
        os.path.join(output_dir, ('').join(
            ['Greenland-Iceland_Mean_DJF_', AIRPRESS_VARIABLE, '.nc'])))
    return season_geo_diff


def main():
    """Run the the meat of the code."""
    args = get_args()

    # base dir needs timestamp
    now = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    base_output = '_'.join((args.output, now))

    # build output tree
    moc_dir, greice_dir = _make_dir_tree(base_output)

    # set logging
    log_level = args.log_level
    _set_logger(logging, base_output, 'acsismoc.log', log_level)
    logger.info(HEADER)
    logger.info("Running main function...")

    # get input files
    moc_file = args.moc_file
    hadslp_file = args.hadslp_file

    # compute means and write to file
    # moc
    moc_constraint = iris.Constraint(
        cube_func=(lambda c: c.var_name == MOC_VARIABLE))
    _yearly_mean(moc_file, moc_dir, MOC_VARIABLE, moc_constraint)

    # vn405
    vn_constraint = iris.Constraint(
        cube_func=(lambda c: c.var_name == AIRPRESS_VARIABLE))
    _yearly_mean(hadslp_file, moc_dir, AIRPRESS_VARIABLE, vn_constraint)

    # greenland-iceland
    gre_ic_djf = _djf_greenland_iceland(hadslp_file, greice_dir, vn_constraint,
                                        SEASON)
    logger.debug('Greenland - Iceland Mean DJF cube %s', gre_ic_djf)


if __name__ == '__main__':
    main()
