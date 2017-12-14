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
import copy
import datetime
import errno
import logging
import logging.config
import os
import shutil
import sys

import yaml

# Hack to make this file executable
if __name__ == '__main__':  # noqa
    sys.path.insert(0,
                    os.path.dirname(
                        os.path.dirname(os.path.abspath(__file__))))  # noqa

from esmvaltool.interface_scripts.auxiliary import ncl_version_check
from esmvaltool.interface_scripts.data_interface import write_settings,\
    get_backward_compatible_ncl_interface
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
        'output_dir': './output_dir',
        'save_intermediary_cubes': False,
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

    cfg['synda_download'] = args.synda_download

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
    # run (only preprocessors for now) # TODO: also run diagnostics from here
    logger.debug("Namelist summary:\n%s", namelist)
    namelist.run()

    os.environ['0_ESMValTool_version'] = __version__

    script_root = os.path.abspath(os.path.dirname(__file__))
    esmvaltool_root = os.path.dirname(script_root)

    # project_info is a dictionary with all info from the namelist.
    project_info = {}
    project_info['GLOBAL'] = global_config
    project_info['MODELS'] = namelist.models
    project_info['DIAGNOSTICS'] = namelist.diagnostics

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

    # create directories
    for key in project_info['GLOBAL']:
        if key.endswith('_dir'):
            if not os.path.isdir(project_info['GLOBAL'][key]):
                logger.info('Creating dir %s', project_info['GLOBAL'][key])
                os.makedirs(project_info['GLOBAL'][key])

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

    def _add_to_order(all_scripts, order, script_name):
        """Add scripts to order depth first"""
        # skip scripts that have already been added earlier
        if script_name in order:
            return

        # if the script has ancestors, add those to order first
        if all_scripts[script_name]['ancestors']:
            for ancestor_name in all_scripts[script_name]['ancestors']:
                _add_to_order(all_scripts, order, ancestor_name)

        # add script to order and remove from todo list
        order.append(script_name)
        all_scripts.pop(script_name)

    def get_order(all_scripts, order):
        """Determine the order in which scripts should be run."""
        while all_scripts:
            # select a random script
            script_name, script_cfg = all_scripts.popitem()
            all_scripts[script_name] = script_cfg

            # add ancestors and script to order
            _add_to_order(all_scripts, order, script_name)

    for diag_name, diag in project_info['DIAGNOSTICS'].items():

        script_order = []
        get_order(dict(diag['scripts']), script_order)
        for script_name in script_order:
            script_cfg = diag['scripts'][script_name]
            if not script_cfg['script']:
                continue

            logger.info("Running diagnostic script %s of diagnostic %s",
                        script_name, diag_name)
            project_info['RUNTIME']['currDiag'] = copy.deepcopy(diag)

            tmp, models = diag['variables'].popitem()
            diag['variables'][tmp] = models
            project_info['ALLMODELS'] = models

            script = script_cfg['script']
            logger.info("Running diag_script: %s", script)
            interface_data = os.path.join(project_info['GLOBAL']['run_dir'],
                                          diag_name, script_name)
            os.makedirs(interface_data)
            ext = 'ncl' if script.lower().endswith('.ncl') else 'yml'
            cfg_file = os.path.join(interface_data, 'settings.' + ext)
            logger.info("with configuration file: %s", cfg_file)
            settings = script_cfg['settings']
            if ext == 'ncl':
                settings = {'diag_script_info': settings}
                compatibility = get_backward_compatible_ncl_interface(
                    diag['variables'], global_config, namelist_file,
                    os.path.basename(script))
                settings.update(compatibility)
            write_settings(settings, cfg_file)
            if ext == 'ncl':
                os.rename(cfg_file,
                          os.path.join(interface_data, 'ncl.interface'))
            project_info['RUNTIME']['currDiag']['script'] = script
            project_info['RUNTIME']['currDiag']['cfg_file'] = cfg_file
            project_info['RUNTIME']['interface_data'] = interface_data
            run_executable(
                script,
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
