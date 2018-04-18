#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Command line tool to retrieve data from the Met Office Managed Archive Storage
System (a.k.a. MASS archive) using the MOOSE system.

Example usage:

./retriever.py --area Africa \\
               --start-date 1990/01/01 --end-date 1991/01/01 \\
               --suite-id u-ab663 \\
               --dest-dir=/scratch/username \\


This will look for query files inside the AutoAssess code tree. To run this
program independently of AutoAssess specify the directory containing query
files:


./retriever.py --area globaltrop_supermeans \\
               --suite-id u-ab642 \\
               --start-date 1990/12/01 --end-date 1991/12/01 \\
               --dest-dir pp_files \\
               --query-dir /abs/path/to/query_templates \\


- query template files
    For each area there must be template MOOSE query files that are sufficient
    to retrieve all required data. They must not contain a time range for
    retrieval.

    There is a naming convention for the filename:

    areaname.[daily|monthly|seasonal|annual].streamletter.txt

    Minimal example:

        monsoon.monthly.m.txt:
            begin
                stash=30201
            end


- destination directory tree structure
    For each query file template the following directory structure is created:

    DEST_DIR/suite_id/assessment_area/averaging_period/stream_letter

    for example:

    /scratch/tvoigt/u-ab663/africa/daily/a/*.pp
                                         e/*.pp

"""

import argparse
import glob
import logging
import multiprocessing
import os
import os.path
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from time import time
from autoassess import query_template_dir
from autoassess import create_dir

mp_logger = multiprocessing.log_to_stderr()
mp_logger.setLevel(logging.INFO)

# TODO query file should belong to each area


def add_time_range(query, start_date, end_date):
    """
    Add time range to MOOSE query string.

    :param str query: Content of a MOOSE query file that does not contain any
        time range specification.
    :param str start_date: Date string in the form: YYYY/MM/DD.
    :param str end_date: Date string in the form: YYYY/MM/DD.
    :returns: Content of a MOOSE query file with time range specification
        as global block.
    :rtype: string
    """
    regex = '^[0-9]{4}/(0[1-9]|1[0-2])/(0[1-9]|[1-2][0-9]|3[0-1])$'
    assert re.match(regex, start_date), 'Not a valid date: ' + start_date
    assert re.match(regex, end_date), 'Not a valid date: ' + end_date

    ############################################################################
    # T1 refers to the beginning of averaged periods, or time instants
    # T2 refers to the end of averaged periods, or in case of time instants to
    # the forecast data time (analysis time), which is the start of the
    # simulation for climate runs
    #
    # average periods:
    #
    #   START <= T1 <= END and END <= T2
    #                                            _____ _____ _____ _____
    #   one year seasons                        |_____|_____|_____|_____|
    #
    # instantaneous data (T2 is meaningless for selection):
    #
    #   START <= T1 <= END
    #
    #   one year one instance per season        |     |     |     |     |
    #
    ############################################################################

    year_block = (
    'begin_global\n'
    '  T1>={START_DATE 00:00}\n'
    '  T1<={END_DATE 00:00}\n'
    '  T2<={END_DATE 23:00}\n'
    'end_global\n'
    '\n')

    year_block = year_block.replace('START_DATE', str(start_date))
    year_block = year_block.replace('END_DATE', str(end_date))
    return year_block + query


def make_moose_command(suite_id, stream_letter, query_file, dest_dir):
    """
    Assemble MOOSE command from templates for the MOOSE URI and MOOSE command.

    :param str suite_id: Rose suite ID of the Model run, for example: u-ab663
    :param str stream_letter: Letter specifying the output stream of the UM.
    :param str query_file: Full path to query file for data retrieval.
    :param str dest_dir: Directory to store retrieved data.
    :returns: MOOSE command that can be run from the command line.
    :rtype: str
    """
    moose_uri_template = 'moose:/crum/SUITE_ID/apSTREAM_LETTER.pp -v'
    moo_select_command_template = 'moo select QUERY_FILE MOOSE_URI DEST_DIR'

    moose_uri = moose_uri_template.replace('SUITE_ID', suite_id)\
                                  .replace('STREAM_LETTER', stream_letter)

    command = moo_select_command_template.replace('QUERY_FILE', query_file)\
                                         .replace('MOOSE_URI', moose_uri)\
                                         .replace('DEST_DIR', dest_dir)
    return command


def make_target_path(suite_id, assessment_area, averaging_period, stream_letter,
                     dest_dir_base):
    """
    Create target path to store the data from individual retrieval queries.

    :param str suite_id: Rose suite ID of the Model run, for example: u-ab663
    :param str assessment_area: Name of the assessement area.
    :param str averaging_period: Averaging period of data. Either 'daily',
        'monthly', 'seasonal', or 'annual'.
    :param str stream_letter: Letter specifying the output stream of the UM.
    :param str dest_dir_base: Top directory into which to store retrieved data.
    :returns: Target path to store the data from individual retrieval queries.
    :rtype: str
    """
    valid_periods = ['daily', 'monthly', 'seasonal', 'annual']
    assert averaging_period in valid_periods, \
        'Invalid averaging period: ' + \
        str(averaging_period) + ' should be one of ' + str(valid_periods)

    dest_dir_template = 'SUITE_ID/ASSESSMENT_AREA/AVERAGING_PERIOD/STREAM_LETTER'

    dest_dir = dest_dir_template.replace('SUITE_ID', suite_id)\
                                .replace('ASSESSMENT_AREA', assessment_area)\
                                .replace('AVERAGING_PERIOD', averaging_period)\
                                .replace('STREAM_LETTER', stream_letter)
    return os.path.join(dest_dir_base, dest_dir)


def parse_args(args):
    """
    Parse arguments in a function to facilitate testing. Contains all command
    line options.

    :param list args: Command line arguments from sys.argv.
    :returns: Checked command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--area', required=True,
                        help='Name of assessment area.')
    parser.add_argument('--suite-id', required=True,
                        help='The name of a Rose suite, something like: u-ab123')
    parser.add_argument('--start-date', required=True,
                        help='Retrieve data at or after this date. Format: YYYY/MM/DD')
    parser.add_argument('--end-date', required=True,
                        help='Retrieve data before this date. Format: YYYY/MM/DD')
    parser.add_argument('--dest-dir', required=True,
                        help='Retrieve data into this directory.')
    parser.add_argument('--query-dir', default=query_template_dir(),
                        help='Provide alternative directory with query template files.')
    parser.add_argument('--tmp-dir', default=os.environ['TMPDIR'],
                        help='Place for temporary query files. Default is $TMPDIR.')
    parser.add_argument('--procs', type=int, default=2,
                        help='Number of MOOSE processes')
    args = parser.parse_args(args)

    regex = '^[a-zA-Z_]+$'
    assert re.match(regex, args.area), \
           regex + ' does not match ' + str(args.area)

    regex = '^[a-z0-9-]+$'
    assert re.match(regex, args.suite_id), \
           regex + ' does not match ' + str(args.suite_id)

    regex = '^[0-9]{4}/(0[1-9]|1[0-2])/(0[1-9]|[1-2][0-9]|3[0-1])$'
    assert re.match(regex, args.start_date), \
           regex + ' does not match ' + str(args.start_date)
    assert re.match(regex, args.end_date), \
           regex + ' does not match ' + str(args.end_date)
    assert args.start_date < args.end_date, \
           'The start date must be before the end date.'

    return args


def run_shell_cmd(cmd):
    """
    Run command in subprocess.

    :param str cmd: Full shell command, as one would run it in a shell.
    """
    t_start = time()
    print cmd
    subprocess.check_call(cmd, shell=True)
    print '### Command took {}s'.format(time() - t_start) + ' ' + cmd


def create_dated_query(query, template_name, start_date, end_date, tmp_dir):
    """
    Add date range to query and write to temporary file.

    :param str query: Query file template.
    :param str template_name: Name of the template file.
    :param str start_date: Start of the retrieval period; a date string in the
        form: YYYY/MM/DD.
    :param str end_date: End of the retrieval period; a date string in the
        form: YYYY/MM/DD.
    :param str tmp_dir: Full path to tempoarary directory to store the dated
        query file.
    :returns: Path to temporary query file with date range added.
    :rtype: str
    """
    query = add_time_range(query, start_date, end_date)
    query_file_path = os.path.join(tmp_dir, template_name)
    create_dir(tmp_dir)
    _, query_file_path = tempfile.mkstemp(prefix=template_name + '_',
                                          dir=tmp_dir)
    with open(query_file_path, 'w') as fh:
        fh.write(query)
    return query_file_path


def retrieve_from_mass():
    """Main function"""
    args = parse_args(sys.argv[1:])

    TMP_DIR = args.tmp_dir
    DEST_DIR = args.dest_dir

    # find query file templates for area
    query_template_paths = []
    for query_file in glob.glob(args.query_dir + '/*.txt'):
        query_filename = os.path.split(query_file)[1]
        if args.area.lower() == query_filename.split('.')[0]:
            query_template_paths.append(query_file)
    if not query_template_paths:
        raise Exception('No query templates found in: ' + str(args.query_dir))

    # assemble MOOSE command for each query file template
    moose_commands = []
    for template_path in query_template_paths:
        filename = os.path.split(template_path)[1]
        area_name, averaging_period, stream_letter, _txt = filename.split('.')

        # create target sub-directory structure
        dest_dir = make_target_path(args.suite_id, area_name,
                                    averaging_period, stream_letter,
                                    DEST_DIR)
        try:
            create_dir(dest_dir, safe=True)
        except OSError as e:
            print os.path.basename(__file__), e
            sys.exit(1)


        # make MOOSE commands
        with open(template_path, 'r') as fh:
            query = fh.read()
        template_name = os.path.split(template_path)[1]
        query_file_path = create_dated_query(query, template_name,
                                             args.start_date, args.end_date,
                                             TMP_DIR)
        command = make_moose_command(args.suite_id, stream_letter,
                                     query_file_path, dest_dir)
        moose_commands.append(command)

    # Run Moose commands in parallel
    t_start = time()
    pool = Pool(processes=args.procs)
    result_obj = pool.map_async(run_shell_cmd, moose_commands)
    _ = result_obj.get(timeout=sys.maxsize)
    # Without the timeout this blocking call ignores all signals, such as Ctrl-C.
    # They are only sent to the processes in the Queues of the Pool. Only the
    # result object of the map_async function has a get function.
    pool.close()
    pool.join()
    print '### Took {}s'.format(time() - t_start)

if __name__ == '__main__':
    retrieve_from_mass()
