"""Compute annual average for moc_mar_hc10."""
import argparse
from datetime import datetime
import logging
import os
import numpy as np
import iris
from esmvaltool.preprocessor import extract_region, extract_season

# set up logging
logger = logging.getLogger(__name__)

# print the header
HEADER = r"""
______________________________________________________________________

            ESMValTool ACSIS Indicators Computation
______________________________________________________________________

""" + __doc__


def get_args():
    """Define the `esmvaltool` command line."""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-m',
        '--moc-file',
        type=str,
        nargs='+',
        help='MOC files')
    parser.add_argument(
        '-s',
        '--hadslp-file',
        type=str,
        nargs='+',
        help='hadslp files')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
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
        real_date = datetime(date_obj.year, date_obj.month,
                             date_obj.day, 0, 0)
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


def _yearly_mean(data_file, variable, var_constraint):
    """Process needed data."""
    cube = iris.load(data_file, constraints=var_constraint)[0]
    datas, no_measurements = _extract_yearly_data(cube)
    for y_r, ydat in datas.items():
        logger.info('For year %s : %s number of measurements: %s',
                    y_r, variable, no_measurements[y_r])
        moc_yearly_mean = np.ma.mean(ydat.data)
        logger.info('For year %s : %s yearly average: %s',
                    y_r, variable, moc_yearly_mean)


def _djf_greenland_iceland(data_file, var_constraint, season):
    """Get the DJF mean for Greenland-Iceland."""
    cube = iris.load(data_file, constraints=var_constraint)[0]
    greenland = extract_region(cube, 25., 35., 30., 40.)
    iceland = extract_region(cube, 15., 25., 60., 70.)
    greenland = greenland.collapsed(['longitude', 'latitude'],
                                    iris.analysis.MEAN)
    iceland = iceland.collapsed(['longitude', 'latitude'],
                                iris.analysis.MEAN)
    diff = greenland - iceland
    season_geo_diff = extract_season(diff, season)
    return season_geo_diff


def main():
    """Run the the meat of the code."""
    args = get_args()

    # set logging
    log_level = args.log_level
    _set_logger(logging, '.', 'moc.log',
                log_level)
    logger.info(HEADER)
    logger.info("Running main function...")

    # get command args
    moc_file = args.moc_file
    hadslp_file = args.hadslp_file

    # compute means and write to file

    # first moc
    moc_constraint = iris.Constraint(
        cube_func=(lambda c: c.var_name == 'moc_mar_hc10')
    )
    _yearly_mean(moc_file, 'moc_mar_hc10', moc_constraint)

    # then vn405
    vn_constraint = iris.Constraint(
        cube_func=(lambda c: c.var_name == 'UM_0_fc8_vn405')
    )
    _yearly_mean(hadslp_file, 'UM_0_fc8_vn405', vn_constraint)

    # now get greenland-iceland
    gre_ic_djf = _djf_greenland_iceland(hadslp_file, vn_constraint, 'DJF')
    iris.save(gre_ic_djf, os.path.join(args.output,
                                       'DJF_GreIce_UM_0_fc8_vn405.nc'))


if __name__ == '__main__':
    main()
