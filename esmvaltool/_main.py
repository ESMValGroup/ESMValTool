"""ESMValTool - Earth System Model Evaluation Tool

http://www.esmvaltool.org

CORE DEVELOPMENT TEAM AND CONTACTS:
  Veronika Eyring (PI; DLR, Germany - veronika.eyring@dlr.de)
  Bjoern Broetz (DLR, Germany - bjoern.broetz@dlr.de)
  Niels Drost (NLESC, Netherlands - n.drost@esciencecenter.nl)
  Nikolay Koldunov (AWI, Germany - nikolay.koldunov@awi.de)
  Axel Lauer (DLR, Germany - axel.lauer@dlr.de)
  Benjamin Mueller (LMU, Germany - b.mueller@iggf.geo.uni-muenchen.de)
  Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
  Mattia Righi (DLR, Germany - mattia.righi@dlr.de)
  Javier Vegas-Regidor (BSC, Spain - javier.vegas@bsc.es)

For further help, please read the documentation at
http://esmvaltool.readthedocs.io. Have fun!
"""

# ESMValTool main script
#
# Authors:
# Bouwe Andela (NLESC, Netherlands - b.andela@esciencecenter.nl)
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

from __future__ import print_function

import argparse
import datetime
import errno
import glob
import logging
import os
import shutil
import sys
from multiprocessing import cpu_count

from . import __version__
from ._config import configure_logging, read_config_user_file
from ._namelist import read_namelist_file
from ._task import resource_usage_logger

# set up logging
logger = logging.getLogger(__name__)

HEADER = r"""
______________________________________________________________________
          _____ ____  __  ____     __    _ _____           _
         | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
         |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
         | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
         |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
______________________________________________________________________

""" + __doc__


def get_args():
    """Define the `esmvaltool` command line"""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__version__,
        help="return ESMValTool's version number and exit")
    parser.add_argument(
        '-n',
        '--namelist-file',
        help='Path to the namelist file',
        required=True)
    parser.add_argument(
        '-c',
        '--config-file',
        default=os.path.join(os.path.dirname(__file__), 'config-user.yml'),
        help='Config file')
    parser.add_argument(
        '-s',
        '--synda-download',
        action='store_true',
        help='Download input data using synda. This requires a working '
        'synda installation.')
    args = parser.parse_args()
    return args


def main(args):
    """Define the `esmvaltool` program"""
    namelist_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.namelist_file)))
    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    # Read user config file
    if not os.path.exists(config_file):
        print("ERROR: config file {} does not exist".format(config_file))

    namelist_name = os.path.splitext(os.path.basename(namelist_file))[0]
    cfg = read_config_user_file(config_file, namelist_name)

    # Create run dir
    if os.path.exists(cfg['run_dir']):
        print("ERROR: run_dir {} already exists, aborting to "
              "prevent data loss".format(cfg['output_dir']))
    os.makedirs(cfg['run_dir'])

    # configure logging
    log_files = configure_logging(
        output=cfg['run_dir'], console_log_level=cfg['log_level'])

    # log header
    logger.info(HEADER)

    logger.info("Using config file %s", config_file)
    logger.info("Writing program log files to:\n%s", "\n".join(log_files))

    cfg['synda_download'] = args.synda_download

    resource_log = os.path.join(cfg['run_dir'], 'resource_usage.txt')
    with resource_usage_logger(pid=os.getpid(), filename=resource_log):
        process_namelist(namelist_file=namelist_file, config_user=cfg)


def process_namelist(namelist_file, config_user):
    """Process namelist"""
    if not os.path.isfile(namelist_file):
        raise OSError(errno.ENOENT, "Specified namelist file does not exist",
                      namelist_file)

    timestamp1 = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d --  %H:%M:%S"

    logger.info(
        "Starting the Earth System Model Evaluation Tool v%s at time: %s ...",
        __version__, timestamp1.strftime(timestamp_format))

    logger.info(70 * "-")
    logger.info("NAMELIST   = %s", namelist_file)
    logger.info("RUNDIR     = %s", config_user['run_dir'])
    logger.info("WORKDIR    = %s", config_user["work_dir"])
    logger.info("PREPROCDIR = %s", config_user["preproc_dir"])
    logger.info("PLOTDIR    = %s", config_user["plot_dir"])
    logger.info(70 * "-")

    logger.info("Running tasks using at most %s processes",
                config_user['max_parallel_tasks'] or cpu_count())

    logger.info(
        "If your system hangs during execution, it may not have enough "
        "memory for keeping this number of tasks in memory. In that case, "
        "try reducing 'max_parallel_tasks' in your user configuration file.")

    # copy namelist to run_dir for future reference
    shutil.copy2(namelist_file, config_user['run_dir'])

    # parse namelist
    namelist = read_namelist_file(namelist_file, config_user)
    logger.debug("Namelist summary:\n%s", namelist)

    # run
    namelist.run()

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info(
        "Ending the Earth System Model Evaluation Tool v%s at time: %s",
        __version__, timestamp2.strftime(timestamp_format))
    logger.info("Time for running namelist was: %s", timestamp2 - timestamp1)

    # Remind the user about reference/acknowledgement file
    out_refs = glob.glob(
        os.path.join(config_user['output_dir'], '*', '*',
                     'references-acknowledgements.txt'))
    logger.info("For the required references/acknowledgements of these "
                "diagnostics see:\n%s", '\n'.join(out_refs))


def run():
    """Run the `esmvaltool` program, logging any exceptions."""
    args = get_args()
    try:
        main(args)
    except:  # noqa
        logger.exception(
            "Program terminated abnormally, see stack trace "
            "below for more information",
            exc_info=True)
        sys.exit(1)
    else:
        logger.info("Run was succesful")
