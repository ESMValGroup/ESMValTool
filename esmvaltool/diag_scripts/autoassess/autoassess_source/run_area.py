"""
(C) Crown Copyright 2017, the Met Office

Wrapper for existing assessment areas of AutoAssess.

Replicates some behaviours of the old version of AutoAssess
which are all deprecated.
All lot of this data is not required for data loading anymore,
but just to not break
the assessment area code.

"""

import sys
import os
import os.path
import argparse
import csv
import datetime
from pprint import pprint
import re
import tempfile

# use Agg backend for non-interactive runs; this is propagated to all funcs
# called from this module
import matplotlib
matplotlib.use('Agg')


def create_dir(path):
    """
    Make a dir
    """
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            pass
        else:
            raise


def create_run_object(args, area, suite_id):
    """
    Create run object containing all necessary information for Ass
    Areas of the previous version of AutoAssess.
    Use only information provided through command line options.

        # Private entries
        run['_area']           - assessment area name
        run['_start_date']     - start date
        run['_end_date']       - end date

        run['runid']           - name of the run (5-char UM runid,suite ID)
        run['data_root']       - dir with data to be assessed
        run['clim_root']       - dir for obs and climatologies
        run['ancil_root']      - directory for ancillary files
        run['nyear']           - length of the assessment period as full years
        run['start']           - start year

        run['from_monthly']    - Date ranges for different mean periods.
        run['to_monthly']        Climatologic years:
        run['from_daily']        all from_* start on 01/12/start_year-1
        run['to_daily']          The to_* date will always be:
        run['from_annual']       1st date + nyears - mean period length
        run['to_annual']         daily: 30/11/XX, monthly: 01/11/XX,
                                 seasonal: 01/09/XX
        run['from_seasonal']     annual: 01/12/XX-1
        run['to_seasonal']

        # only use by Stratosphere area
        run.id                 - suite ID
        run.title              - '# TITLE #'

    :param argparse.Namespace args: Command line arguments.
    :param str area: Name of assessment area.
    :param str suite_id: Model run suite ID.
    :returns: Run dictionary.
    :rtype: Dictionary with attributes
    """

    class run(dict):
        pass

    run = run()
    # added private entries required for replicating data retrieval with the
    # previous API (loaddata)
    run['_area'] = area
    run['_start_date'] = args.start_date
    run['_end_date'] = args.end_date

    run.title = '# TITLE #'  # stratosphere
    run['runid'] = suite_id
    run.id = run['runid']  # stratosphere
    run['data_root'] = args.data_dir
    run['clim_root'] = args.obs_dir
    run['ancil_root'] = args.ancil_dir

    start_year = int(args.start_date[0:4])
    end_year = int(args.end_date[0:4])
    run['start'] = start_year
    run['nyear'] = end_year - start_year
    run.period = '{:04d}_{:03d}'.format(run['start'],
                                        run['nyear'])  # stratosphere uses this

    year, month, day = map(int, args.start_date.split('/'))
    run['from_instantaneous'] = datetime.datetime(year, month, day)
    run['from_daily'] = datetime.datetime(year, month, day)
    run['from_monthly'] = datetime.datetime(year, month, day)
    run['from_seasonal'] = datetime.datetime(year, month, day)
    run['from_annual'] = datetime.datetime(year, month, day)

    year, month, day = map(int, args.end_date.split('/'))
    assert month == 12 and day == 1  # Climatological year
    run['to_instantaneous'] = datetime.datetime(year, 11, 30)
    run['to_daily'] = datetime.datetime(year, 11, 30)
    run['to_monthly'] = datetime.datetime(year, 11, 1)
    run['to_seasonal'] = datetime.datetime(year, 9, 1)
    run['to_annual'] = datetime.datetime(year - 1, 12, 1)
    return run


def parse_args(args):
    """
    Parse arguments in a function to facilitate testing. Contains all command
    line options.

    :param list args: Command line arguments from sys.argv.
    :returns: Checked command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--area', required=True, help='Name of assessment area.')
    parser.add_argument(
        '--suite-id1',
        required=True,
        help='The name of a Rose suite, something like: u-ab123')
    parser.add_argument(
        '--suite-id2',
        required=True,
        help='The name of a Rose suite, something like: u-ab123')
    parser.add_argument(
        '--start-date',
        required=True,
        help=('Start of assessment period. Must be the beginning '
              'of a climatologic year.'
              'Format: YYYY/12/01'))
    parser.add_argument(
        '--end-date',
        required=True,
        help=('End of assessment period. Beginning of the '
              'first climatologic year after the end of the '
              'assessment period.'
              'Format: YYYY/12/01'))
    parser.add_argument(
        '--obs-dir',
        required=True,
        help='Directory with Observations and Climatologies.')
    parser.add_argument(
        '--ancil-dir', required=True, help='Directory with Ancillary files.')
    parser.add_argument(
        '--out-dir', required=True, help='Write results into this directory.')
    parser.add_argument(
        '--data-dir', required=True, help='Directory tree with climate data.')
    parser.add_argument(
        '--tmp-dir',
        default='tmp',
        help='Place for temporary files. Default is $TMPDIR.')
    args = parser.parse_args(args)

    regex = '^[a-zA-Z_]+$'
    assert re.match(regex,
                    args.area), regex + ' does not match ' + str(args.area)

    regex = '^[a-z0-9-]+$'

    # at least ONE year
    year, month, day = map(int, args.start_date.split('/'))
    start = datetime.date(year, month, day)
    year, month, day = map(int, args.end_date.split('/'))
    end = datetime.date(year, month, day)
    assert end.year - start.year >= 1 and \
        end.month >= start.month and \
        end.day >= start.day, \
        'Assessment requires at least two years of data.'
    # As collapsing over a single year/season/month fails.
    # Cubes also loose the DimCoord 'time' if only one time point is left.

    # climatologic years
    regex = '^[0-9]{4}/12/01$'
    assert re.match(
        regex,
        args.start_date), regex + ' does not match ' + str(args.start_date)
    assert re.match(
        regex, args.end_date), regex + ' does not match ' + str(args.end_date)
    assert args.start_date < args.end_date, 'Start must be before end.'
    regex = '^(/[a-zA-Z0-9_-]*)+$'
    assert re.match(
        regex, args.obs_dir), regex + ' does not match ' + str(args.obs_dir)
    assert re.match(
        regex,
        args.ancil_dir), regex + ' does not match ' + str(args.ancil_dir)
    assert re.match(
        regex, args.out_dir), regex + ' does not match ' + str(args.out_dir)
    assert re.match(
        regex, args.data_dir), regex + ' does not match ' + str(args.data_dir)
    assert re.match(
        regex, args.tmp_dir), regex + ' does not match ' + str(args.tmp_dir)

    return args


def create_output_tree(out_dir, ref_suite_id, exp_suite_id, area):
    """
    Create directory tree for area output according to the following scheme:

        `out_dir`/`exp_suite_id`_vs_`ref_suite_id`/`area`

    If the leaf directory `area` exists raises OSError.

    :param str out_dir: Base directory for output.
    :param str suite_id1: Suite Id of reference model run.
    :param str suite_id2: Suite Id of test model run.
    :param str area: Name of asssessment area.
    :returns: Path to area output directory.
    :rtype: str
    :raises: OSError
    """
    assessment_name = exp_suite_id + '_vs_' + ref_suite_id
    # make sure out_dir exists in output folder
    _out_dir = os.path.join(out_dir, assessment_name)
    create_dir(out_dir)

    # create output folder for area
    area_out_dir = os.path.join(_out_dir, area)
    create_dir(area_out_dir)
    return area_out_dir


def create_tmp_dir(tmp_dir, ref_suite_id, exp_suite_id, area):
    """
    Create directory tree for temporary data according to the following scheme:

        `tmp_dir`/`exp_suite_id`_vs_`ref_suite_id`_random/`area`_random

    :param str tmp_dir: Base temporary directory.
    :param str suite_id1: Suite ID of reference model run.
    :param str suite_id2: Suite ID of test model run.
    :param str area: Name of asssessment area.
    :returns: Path to area temporary directory.
    :rtype: str
    """
    assessment_name = exp_suite_id + '_vs_' + ref_suite_id
    # create unique temporary folder in tmp dir
    _tmp_dir = tempfile.mkdtemp(prefix=assessment_name + '_', dir=tmp_dir)

    # create temporary folder for area
    area_tmp_dir = tempfile.mkdtemp(prefix=area + '_', dir=_tmp_dir)
    return area_tmp_dir


def run_area():
    """ """
    args = parse_args(sys.argv[1:])

    area = str(args.area).lower()

    # import area here to allow removal of areas
    if area == 'monsoon':
        import monsoon as area_package
    elif area == 'stratosphere':
        import stratosphere as area_package
    elif area == 'hydrocycle':
        import hydrocycle as area_package
    elif area == 'conservation':
        import conservation as area_package
    elif area == 'globaltrop':
        import globaltrop as area_package
    elif area == 'land_surface':
        import land_surface as area_package
    else:
        raise Exception('Unknown area: ' + str(area))

    area_tmp_dir = create_tmp_dir(args.tmp_dir, args.suite_id1, args.suite_id2,
                                  area)
    area_out_dir = create_output_tree(args.out_dir, args.suite_id1,
                                      args.suite_id2, area)

    # the areas write all output to the cwd
    os.chdir(area_out_dir)

    for suite_id in [args.suite_id1, args.suite_id2]:
        all_metrics = {}
        # run each metric function
        for metric_function in area_package.metrics_functions:
            print('# Call:', metric_function.__name__)
            run_object = create_run_object(args, area, suite_id)
            metrics = metric_function(run_object)
            print('# metrics: ', pprint(metrics))

            duplicate_metrics = set(all_metrics.keys()) & set(metrics.keys())
            if len(duplicate_metrics) != 0:
                raise AssertionError('Duplicate Metrics ' +
                                     str(duplicate_metrics))
            all_metrics.update(metrics)

        # write metrics to file
        create_dir(os.path.join(area_out_dir, suite_id))
        with open(os.path.join(area_out_dir, suite_id, 'metrics.csv'),
                  'w') as fh:
            writer = csv.writer(fh)
            for metric in all_metrics.items():
                writer.writerow(metric)

    # multimodel functions
    if hasattr(area_package, 'multi_functions'):
        ref_run_object = create_run_object(args, area, args.suite_id1)
        test_run_object = create_run_object(args, area, args.suite_id2)
        run_objects = [ref_run_object,
                       test_run_object]  # reference must be first
        for multi_function in area_package.multi_functions:
            multi_function(run_objects)
    else:
        print('# Area has no multi functions.')


if __name__ == '__main__':
    run_area()
