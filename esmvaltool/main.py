#! /usr/bin/env python
r"""
 ______________________________________________________________________
           _____ ____  __  ____     __    _ _____           _
          | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
          |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
          | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
          |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
 ______________________________________________________________________

 ESMValTool - Earth System Model Evaluation Tool
 http://www.esmvaltool.org/

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
 ______________________________________________________________________

 For further help, check the doc/-folder for pdfs
 and references therein. Have fun!
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
import logging.config
import os
import shutil
import sys
from multiprocessing import cpu_count

import yaml

# Hack to make this file executable
if __name__ == '__main__':  # noqa
    sys.path.insert(0,
                    os.path.dirname(
                        os.path.dirname(os.path.abspath(__file__))))  # noqa

from esmvaltool.interface_scripts.auxiliary import ncl_version_check
from esmvaltool.namelist import read_namelist_file
from esmvaltool.version import __version__

# set up logging
if __name__ == '__main__':
    logger = logging.getLogger('ESMValTool')
    logger.addHandler(logging.NullHandler())
else:
    logger = logging.getLogger(__name__)


def configure_logging(cfg_file=None, output=None, console_log_level=None):
    """Set up logging"""
    if cfg_file is None:
        cfg_file = os.path.join(
            os.path.dirname(__file__), 'config-logging.yml')

    if output is None:
        output = os.getcwd()

    cfg_file = os.path.abspath(cfg_file)
    with open(cfg_file) as file_handler:
        cfg = yaml.safe_load(file_handler)

    for handler in cfg['handlers'].values():
        if 'filename' in handler:
            if not os.path.isabs(handler['filename']):
                handler['filename'] = os.path.join(output, handler['filename'])
        if console_log_level is not None and 'stream' in handler:
            if handler['stream'] in ('ext://sys.stdout', 'ext://sys.stderr'):
                handler['level'] = console_log_level.upper()

    logging.config.dictConfig(cfg)


def read_config_file(config_file, namelist_name):
    """Read config file and store settings in a dictionary."""
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    # set defaults
    defaults = {
        'write_plots': True,
        'write_netcdf': True,
        'exit_on_warning': False,
        'max_data_filesize': 100,
        'output_file_type': 'ps',
        'output_dir': './output_dir',
        'save_intermediary_cubes': False,
        'max_parallel_tasks': 1,
        'run_diagnostic': True,
    }

    for key in defaults:
        if key not in cfg:
            logger.warning("No %s specification in config file, "
                           "defaulting to %s", key, defaults[key])
            cfg[key] = defaults[key]

    # expand ~ to /home/username in directory names and normalize paths
    cfg['output_dir'] = os.path.abspath(os.path.expanduser(cfg['output_dir']))

    for key in cfg['rootpath']:
        cfg['rootpath'][key] = os.path.abspath(
            os.path.expanduser(cfg['rootpath'][key]))

    # insert a directory date_time_namelist_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((namelist_name, now))
    cfg['output_dir'] = os.path.join(cfg['output_dir'], new_subdir)

    # create subdirectories
    cfg['preproc_dir'] = os.path.join(cfg['output_dir'], 'preproc')
    cfg['work_dir'] = os.path.join(cfg['output_dir'], 'work')
    cfg['plot_dir'] = os.path.join(cfg['output_dir'], 'plots')
    cfg['run_dir'] = os.path.join(cfg['output_dir'], 'tmp')

    return cfg


def main():
    """Define the `esmvaltool` program"""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
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

    namelist_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.namelist_file)))
    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    ####################################################
    # Set up logging before anything else              #
    ####################################################

    # configure logging

    if not os.path.exists(config_file):
        print("ERROR: config file {} does not exist".format(config_file))

    namelist_name = os.path.splitext(os.path.basename(namelist_file))[0]
    cfg = read_config_file(config_file, namelist_name)

    # Create run dir
    if os.path.exists(cfg['run_dir']):
        print("ERROR: run_dir {} already exists, aborting to "
              "prevent data loss".format(cfg['output_dir']))
    os.makedirs(cfg['run_dir'])

    configure_logging(
        output=cfg['run_dir'], console_log_level=cfg['log_level'])

    # log header
    logger.info(__doc__)

    logger.info("Using config file %s", config_file)

    # check NCL version
    ncl_version_check()

    cfg['synda_download'] = args.synda_download

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
    try:
        main()
    except:  # noqa
        logger.exception(
            "Program terminated abnormally, see stack trace "
            "below for more information",
            exc_info=True)
        sys.exit(1)
    else:
        logger.info("Run was succesful")


if __name__ == '__main__':
    run()
