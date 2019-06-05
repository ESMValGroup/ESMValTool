"""ESMValTool - Earth System Model Evaluation Tool.

http://www.esmvaltool.org

CORE DEVELOPMENT TEAM AND CONTACTS:
  Veronika Eyring (PI; DLR, Germany - veronika.eyring@dlr.de)
  Bouwe Andela (NLESC, Netherlands - b.andela@esciencecenter.nl)
  Bjoern Broetz (DLR, Germany - bjoern.broetz@dlr.de)
  Lee de Mora (PML, UK - ledm@pml.ac.uk)
  Niels Drost (NLESC, Netherlands - n.drost@esciencecenter.nl)
  Nikolay Koldunov (AWI, Germany - nikolay.koldunov@awi.de)
  Axel Lauer (DLR, Germany - axel.lauer@dlr.de)
  Benjamin Mueller (LMU, Germany - b.mueller@iggf.geo.uni-muenchen.de)
  Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
  Mattia Righi (DLR, Germany - mattia.righi@dlr.de)
  Manuel Schlund (DLR, Germany - manuel.schlund@dlr.de)
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
from ._config import configure_logging, read_config_user_file, DIAGNOSTICS_PATH
from ._recipe import TASKSEP, read_recipe_file
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
    """Define the `esmvaltool` command line."""
    # parse command line args
    parser = argparse.ArgumentParser(
        description=HEADER,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('recipe', help='Path or name of the yaml recipe file')
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__version__,
        help="return ESMValTool's version number and exit")
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
    parser.add_argument(
        '--max-datasets',
        type=int,
        help='Try to limit the number of datasets used to MAX_DATASETS.')
    parser.add_argument(
        '--max-years',
        type=int,
        help='Limit the number of years to MAX_YEARS.')
    parser.add_argument(
        '--skip-nonexistent',
        action='store_true',
        help="Skip datasets that cannot be found.")
    parser.add_argument(
        '--diagnostics',
        nargs='*',
        help="Only run the named diagnostics from the recipe.")
    args = parser.parse_args()
    return args


def main(args):
    """Define the `esmvaltool` program."""
    recipe = args.recipe
    if not os.path.exists(recipe):
        installed_recipe = os.path.join(
            DIAGNOSTICS_PATH, 'recipes', recipe)
        if os.path.exists(installed_recipe):
            recipe = installed_recipe
    recipe = os.path.abspath(os.path.expandvars(os.path.expanduser(recipe)))

    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    # Read user config file
    if not os.path.exists(config_file):
        print("ERROR: config file {} does not exist".format(config_file))

    recipe_name = os.path.splitext(os.path.basename(recipe))[0]
    cfg = read_config_user_file(config_file, recipe_name)

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

    cfg['skip-nonexistent'] = args.skip_nonexistent
    cfg['diagnostics'] = {
        pattern if TASKSEP in pattern else pattern + TASKSEP + '*'
        for pattern in args.diagnostics or ()
    }
    cfg['synda_download'] = args.synda_download
    for limit in ('max_datasets', 'max_years'):
        value = getattr(args, limit)
        if value is not None:
            if value < 1:
                raise ValueError("--{} should be larger than 0.".format(
                    limit.replace('_', '-')))
            cfg[limit] = value

    resource_log = os.path.join(cfg['run_dir'], 'resource_usage.txt')
    with resource_usage_logger(pid=os.getpid(), filename=resource_log):
        process_recipe(recipe_file=recipe, config_user=cfg)
    return cfg


def process_recipe(recipe_file, config_user):
    """Process recipe."""
    if not os.path.isfile(recipe_file):
        raise OSError(errno.ENOENT, "Specified recipe file does not exist",
                      recipe_file)

    timestamp1 = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d %H:%M:%S"

    logger.info(
        "Starting the Earth System Model Evaluation Tool v%s at time: %s UTC",
        __version__, timestamp1.strftime(timestamp_format))

    logger.info(70 * "-")
    logger.info("RECIPE   = %s", recipe_file)
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

    if config_user['compress_netcdf']:
        logger.warning(
            "You have enabled NetCDF compression. Accesing .nc files can be "
            "much slower than expected if your access pattern does not match "
            "their internal pattern. Make sure to specify the expected "
            "access pattern in the recipe as a parameter to the 'save' "
            "preprocessor function. If the problem persists, try disabling "
            "NetCDF compression.")

    # copy recipe to run_dir for future reference
    shutil.copy2(recipe_file, config_user['run_dir'])

    # parse recipe
    recipe = read_recipe_file(recipe_file, config_user)
    logger.debug("Recipe summary:\n%s", recipe)

    # run
    recipe.run()

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info(
        "Ending the Earth System Model Evaluation Tool v%s at time: %s UTC",
        __version__, timestamp2.strftime(timestamp_format))
    logger.info("Time for running the recipe was: %s", timestamp2 - timestamp1)

    # Remind the user about reference/acknowledgement file
    out_refs = glob.glob(
        os.path.join(config_user['output_dir'], '*', '*',
                     'references-acknowledgements.txt'))
    logger.info(
        "For the required references/acknowledgements of these "
        "diagnostics see:\n%s", '\n'.join(out_refs))


def run():
    """Run the `esmvaltool` program, logging any exceptions."""
    args = get_args()
    try:
        conf = main(args)
    except:  # noqa
        if not logger.handlers:
            # Add a logging handler if main failed to do so.
            logging.basicConfig()
        logger.exception(
            "Program terminated abnormally, see stack trace "
            "below for more information",
            exc_info=True)
        logger.info(
            "If you suspect this is a bug or need help, please open an issue "
            "on https://github.com/ESMValGroup/ESMValTool/issues and attach "
            "the run/recipe_*.yml and run/main_log_debug.txt files from the "
            "output directory.")
        sys.exit(1)
    else:
        if conf["remove_preproc_dir"]:
            logger.info("Removing preproc containing preprocessed data")
            logger.info("If this data is further needed, then")
            logger.info("set remove_preproc_dir to false in config")
            shutil.rmtree(conf["preproc_dir"])
        logger.info("Run was successful")
