#! /usr/bin/env python
r"""ESMValTool - Earth System Model Evaluation Tool

______________________________________________________
  _____ ____  __  ____     __    _ _____           _
 | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
 |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
 | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
 |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
______________________________________________________

 http://www.esmvaltool.org/
______________________________________________________________________

CORE DEVELOPMENT TEAM AND CONTACTS:
  Veronika Eyring (PI; DLR, Germany - veronika.eyring@dlr.de)
  Bjoern Broetz (DLR, Germany - bjoern.broetz@dlr.de)
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
import copy
import datetime
import errno
import logging
import logging.config
import os
import shutil
import sys

import iris
import yaml

# Hack to make this file executable
if __name__ == '__main__':  # noqa
    sys.path.insert(0,
                    os.path.dirname(
                        os.path.dirname(os.path.abspath(__file__))))  # noqa

from esmvaltool.interface_scripts.auxiliary import ncl_version_check
from esmvaltool.interface_scripts.data_interface import write_settings
from esmvaltool.interface_scripts.preprocess import run_executable
from esmvaltool.namelist import read_namelist_file

# Define ESMValTool version
__version__ = "2.0.0"

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
        'run_dir': './',
        'preproc_dir': './preproc/',
        'work_dir': './work/',
        'plot_dir': './plots/',
        'save_intermediary_cubes': False,
        'run_diagnostic': True,
        'user_tag': '',
    }

    for key in defaults:
        if key not in cfg:
            logger.warning("No %s specification in config file, "
                           "defaulting to %s", key, defaults[key])
            cfg[key] = defaults[key]

    # expand ~ to /home/username in directory names and normalize paths
    for key in cfg:
        if key.endswith('_dir'):
            cfg[key] = os.path.abspath(os.path.expanduser(cfg[key]))

    for key in cfg['rootpath']:
        cfg['rootpath'][key] = os.path.abspath(
            os.path.expanduser(cfg['rootpath'][key]))

    # insert a directory date_time_namelist_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((now, namelist_name))
    if cfg['user_tag']:
        new_subdir += '_' + cfg['user_tag']

    for key in 'preproc_dir', 'work_dir', 'plot_dir':
        path = cfg[key]
        if path.startswith(cfg['run_dir']):
            cfg[key] = os.path.join(
                os.path.dirname(path), new_subdir, os.path.basename(path))
        else:
            cfg[key] = os.path.join(path, new_subdir)

    cfg['run_dir'] = os.path.join(cfg['run_dir'], new_subdir)

    return cfg


def create_interface_data_dir(project_info, diagnostic_name):
    """Create a directory for storing temporary files needed by esmvaltool.

    ESMValTool transfers information from the `esmvaltool` program to the
    diagnostic programs that it runs by storing that information in files,
    this function creates a directory to store those files.

    Returns
    -------
    str
        The name of the created directory.

    """
    interface_data = os.path.join(
        project_info['GLOBAL']['work_dir'],
        'interface_data',
        diagnostic_name,
    )
    os.makedirs(interface_data)
    return interface_data


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
        print("ERROR: run_dir {} already exists, aborting to prevent data loss"
              .format(cfg['run_dir']))
    os.makedirs(cfg['run_dir'])

    configure_logging(
        output=cfg['run_dir'], console_log_level=cfg['log_level'])

    # log header
    logger.info(__doc__)

    logger.info("Using config file %s", config_file)

    # check NCL version
    ncl_version_check()

    iris.FUTURE.netcdf_no_unlimited = True
    process_namelist(namelist_file=namelist_file, global_config=cfg)


def process_namelist(namelist_file, global_config):
    """Process namelist"""
    if not os.path.isfile(namelist_file):
        raise OSError(errno.ENOENT, "Specified namelist file does not exist",
                      namelist_file)

    logger.info("Processing namelist %s", namelist_file)

    # copy namelist to run_dir for future reference
    shutil.copy2(namelist_file, global_config['run_dir'])

    # parse namelist
    namelist = read_namelist_file(namelist_file, global_config)
    # run (only preprocessors for now)
    logger.debug("Namelist %s", namelist)
    namelist.run()

    os.environ['0_ESMValTool_version'] = __version__

    script_root = os.path.abspath(os.path.dirname(__file__))
    esmvaltool_root = os.path.dirname(script_root)

    # project_info is a dictionary with all info from the namelist.
    project_info = {}
    project_info['GLOBAL'] = global_config
    project_info['MODELS'] = namelist._models
    project_info['DIAGNOSTICS'] = namelist._diagnostics

    # FIX-ME: outdated, keep until standard logging is fully implemented
    project_info['GLOBAL']['verbosity'] = 1

    # additional entries to 'project_info'
    project_info['RUNTIME'] = {}
    project_info['RUNTIME']['yml'] = namelist_file
    project_info['RUNTIME']['yml_name'] = os.path.basename(namelist_file)

    # set references/acknowledgement file
    refs_acknows_file = str.replace(project_info['RUNTIME']['yml_name'],
                                    "namelist_", "refs-acknows_")
    refs_acknows_file = refs_acknows_file.split(os.extsep)[0] + ".log"
    out_refs = os.path.join(project_info["GLOBAL"]['run_dir'],
                            refs_acknows_file)
    project_info['RUNTIME']['out_refs'] = out_refs

    # print summary
    logger.info(70 * "-")
    logger.info("NAMELIST   = %s", project_info['RUNTIME']['yml_name'])
    logger.info("RUNDIR     = %s", project_info["GLOBAL"]['run_dir'])
    logger.info("WORKDIR    = %s", project_info["GLOBAL"]["work_dir"])
    logger.info("PREPROCDIR = %s", project_info["GLOBAL"]["preproc_dir"])
    logger.info("PLOTDIR    = %s", project_info["GLOBAL"]["plot_dir"])
    logger.info("LOGFILE    = %s", project_info['RUNTIME']['out_refs'])
    logger.info(70 * "-")

    # master references-acknowledgements file (hard coded)
    in_refs = os.path.join(esmvaltool_root, 'doc',
                           'MASTER_authors-refs-acknow.txt')
    project_info['RUNTIME']['in_refs'] = in_refs
    logger.info("Using references from %s", in_refs)

    # create workdir
    if not os.path.isdir(project_info['GLOBAL']['work_dir']):
        logger.info('Creating work_dir %s', project_info['GLOBAL']['work_dir'])
        os.makedirs(project_info['GLOBAL']['work_dir'])

    # create rundir
    if not os.path.isdir(project_info['GLOBAL']['run_dir']):
        logger.info('Creating run_dir %s', project_info['GLOBAL']['run_dir'])
        os.makedirs(project_info['GLOBAL']['run_dir'])

    # create refs-acknows file in run_dir (empty if existing)
    with open(out_refs, "w"):
        pass

    # current working directory
    project_info['RUNTIME']['cwd'] = os.getcwd()

    # summary to std-out before starting the loop
    timestamp1 = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d --  %H:%M:%S"

    logger.info(
        "Starting the Earth System Model Evaluation Tool v%s at time: %s ...",
        __version__, timestamp1.strftime(timestamp_format))

    for diag_name, diag in project_info['DIAGNOSTICS'].items():

        project_info['RUNTIME']['currDiag'] = copy.deepcopy(diag)

        for name, variables in diag['variables'].items():
            var = variables[0]
            project_info['RUNTIME']['derived_var'] = var['short_name']
            project_info['RUNTIME']['derived_field_type'] = var['field']
            project_info['ALLMODELS'] = variables

            executable = diag['script'][-1]
            logger.info("Running diag_script: %s", executable)
            interface_data = create_interface_data_dir(
                project_info, diag_name)
            ext = 'ncl' if executable.lower().endswith('.ncl') else 'yml'
            cfg_file = os.path.join(interface_data, 'settings.' + ext)
            logger.info("with configuration file: %s", cfg_file)
            write_settings(diag['settings'], cfg_file)
            project_info['RUNTIME']['currDiag']['cfg_file'] = cfg_file
            project_info['RUNTIME']['interface_data'] = interface_data
            run_executable(
                executable,
                project_info,
                project_info['GLOBAL']['verbosity'],
                project_info['GLOBAL']['exit_on_warning'],
            )

    # delete environment variable
    del os.environ['0_ESMValTool_version']

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info(
        "Ending the Earth System Model Evaluation Tool v%s at time: %s",
        __version__, timestamp2.strftime(timestamp_format))
    logger.info("Time for running namelist was: %s", timestamp2 - timestamp1)

    # Remind the user about reference/acknowledgement file
    logger.info("For the required references/acknowledgements of these "
                "diagnostics see: %s", out_refs)


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
