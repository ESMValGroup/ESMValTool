#! /usr/bin/env python
r"""
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

ESMValTool - Earth System Model Evaluation Tool

For further help, check the doc/-folder for pdfs
and references therein. Have fun!
"""

# ESMValTool main script
#
# Authors:
# Bouwe Andela (NLESC, Netherlands - b.andela@esciencecenter.nl)
# Valeriu Predoi (URead, UK - valeriu.predoi@ncas.ac.uk)
# Mattia Righi (DLR, Germany - mattia.righi@dlr.de)

import argparse
import copy
import datetime
import errno
import logging
import logging.config
import os
import re
import sys
import uuid
import yaml

import interface_scripts.namelistchecks as pchk
import interface_scripts.preprocess as pp
from interface_scripts.auxiliary import ncl_version_check
from interface_scripts.yaml_parser import Parser as Ps

logger = logging.getLogger('ESMValTool')
logger.addHandler(logging.NullHandler())

# Define ESMValTool version
__version__ = "2.0.0"


def configure_logging(cfg_file=None, output=None, console_log_level=None):
    """Set up logging"""
    if cfg_file is None:
        cfg_file = os.path.join(os.path.dirname(__file__), 'config-logging.yml')

    if output is None:
        output = os.getcwd()

    cfg_file = os.path.abspath(cfg_file)
    with open(cfg_file) as file:
        cfg = yaml.safe_load(file)

    for handler in cfg['handlers'].values():
        if 'filename' in handler:
            if not os.path.isabs(handler['filename']):
                handler['filename'] = os.path.join(output, handler['filename'])
        if console_log_level is not None and 'stream' in handler:
            if handler['stream'] in ('ext://sys.stdout', 'ext://sys.stderr'):
                handler['level'] = console_log_level.upper()

    logging.config.dictConfig(cfg)


def get_log_level_and_path(config_file):
    """ Get the log level and path from config.ini or config.yml,
        as logging needs to be configured as early as possible.
    """

    # set defaults
    cfg = {
        'log_level': 'INFO',
        'log_path': None,
    }

    # update defaults if possible
    if isinstance(config_file, str) and os.path.exists(config_file):

        with open(config_file, 'r') as file:
            cfg_file = file.read()

        options = {
            'log_level': r'(?i)^\s*log_level\s*(?:=|:)\s*(debug|info|warning|error)\s*(?:#.*|)$',
            'log_path': r'(?i)^\s*run_dir\s*(?:=|:)\s*(.*)\s*(?:#.*|)$',
        }

        for key, regex in options.items():
            match = re.search(regex, cfg_file, re.MULTILINE)
            if match:
                cfg[key] = match.groups()[0]

    return cfg['log_level'], cfg['log_path']


def read_config_file(config_file):
    """ Read config file and store settings in a dictionary
        
    """
    glob = {}
    glob = yaml.load(file(config_file, 'r'))

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
        'cmip5_dirtype': None,
        'run_diagnostic': True
    }

    for key in defaults:
        if not key in glob:
            logger.warning("No %s specification in config file, defaulting to %s" % (key, defaults[key]))
            glob[key] = defaults[key]

    return(glob)


def create_interface_data_dir(project_info, executable):
    """ Create a temporary directory for storing files needed to run
        executable, returns the name of the created directory.
    """
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    interface_data = os.path.join(
        project_info['GLOBAL']['work_dir'],
        'interface_data',
        project_info['RUNTIME']['currDiag'].id,
        now + '_' + os.path.splitext(os.path.basename(executable))[0],
    )
    os.makedirs(interface_data)
    return interface_data


def main():
    """ Run the program"""

    # parse command line args
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--namelist-file',
                        help='Namelist file')
    parser.add_argument('-c', '--config-file',
                        default=os.path.join(os.path.dirname(__file__), 'config.ini'),
                        help='Config file')
    args = parser.parse_args()

    namelist_file = os.path.abspath(args.namelist_file)
    config_file = os.path.abspath(args.config_file)

    ####################################################
    # Set up logging before anything else              #
    ####################################################

    # get the required console log level and location to save the logs
    # from config file

    log_level, run_dir = get_log_level_and_path(config_file)

    # if run_dir exists, don't overwrite it
    previous_run_dir = None
    if not run_dir is None:
        run_dir = os.path.abspath(run_dir)
        if os.path.isdir(run_dir):
            previous_run_dir = run_dir + '_' + uuid.uuid4().hex
            os.rename(run_dir, previous_run_dir)
        os.makedirs(run_dir)

    # configure logging
    configure_logging(output=run_dir, console_log_level=log_level)

    # log header
    logger.info(__doc__)

    if run_dir is None:
        logger.warning("Failed to retrieve run_dir from config file")

    if previous_run_dir:
        logger.info('Renamed existing run directory %s to %s',
                    run_dir, previous_run_dir)

    logger.info("Using config file %s", config_file)

    # check NCL version
    ncl_version_check()

    process_namelist(namelist_file=namelist_file, config_file=config_file)


def process_namelist(namelist_file, config_file):
    """Process namelist"""

    if not os.path.isfile(namelist_file):
        raise OSError(errno.ENOENT, "Specified namelist file does not exist", namelist_file)

    if not os.path.isfile(config_file):
        raise OSError(errno.ENOENT, "Specified config file does not exist", config_file)

    logger.info("Processing namelist %s", namelist_file)

    os.environ['0_ESMValTool_version'] = __version__

    preprocess_id = None
    script_root = os.path.abspath(os.path.dirname(__file__))

    # get namelist file
    yml_path = namelist_file

    # parse input namelist into project_info-dictionary.
    Project = Ps()

    # parse config file into GLOBAL_DICT info_dictionary
    GLOBAL_DICT = read_config_file(config_file)

    # project_info is a dictionary with all info from the namelist.
    project_info_0 = Project.load_namelist(yml_path)

    # project_info is a dictionary with all info from the namelist.
    project_info = {}
    project_info['GLOBAL'] = GLOBAL_DICT
    project_info['MODELS'] = project_info_0.MODELS
    project_info['DIAGNOSTICS'] = project_info_0.DIAGNOSTICS

    # FIX-ME: outdated, keep until standard logging is fully implemented
    project_info['GLOBAL']['verbosity'] = 1
    verbosity = 1 

    # additional entries to 'project_info'
    project_info['RUNTIME'] = {}
    project_info['RUNTIME']['yml'] = yml_path
    project_info['RUNTIME']['yml_name'] = os.path.basename(yml_path)

    # set references/acknowledgement file
    refs_acknows_file = str.replace(project_info['RUNTIME']['yml_name'], "namelist_", "refs-acknows_")
    refs_acknows_file = refs_acknows_file.split(os.extsep)[0] + ".log"
    out_refs = os.path.join(project_info["GLOBAL"]['run_dir'], refs_acknows_file)
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

    # perform options integrity checks
    logger.info('Checking integrity of namelist')
    tchk1 = datetime.datetime.utcnow()
    pchk.models_checks(project_info['MODELS'])
    pchk.diags_checks(project_info['DIAGNOSTICS'])
    pchk.preprocess_checks(project_info_0.PREPROCESS)
    tchk2 = datetime.datetime.utcnow()
    dtchk = tchk2 - tchk1
    logger.info('Namelist check successful! Time: %s', dtchk)

    # this will have to be purget at some point in the future
    project_info['CONFIG'] = project_info_0.CONFIG

    # tell the environment about regridding
    project_info['RUNTIME']['regridtarget'] = []

    # master references-acknowledgements file (hard coded)
    in_refs = os.path.join(script_root, 'doc', 'MASTER_authors-refs-acknow.txt')
    project_info['RUNTIME']['in_refs'] = in_refs

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

    logger.info("Starting the Earth System Model Evaluation Tool v%s at time: %s ...",
                __version__, timestamp1.strftime(timestamp_format))

    # variables needed for target variable, according to variable_defs
    if not os.path.isabs(project_info['CONFIG']['var_def_scripts']):
        project_info['CONFIG']['var_def_scripts'] = \
            os.path.join(script_root, project_info['CONFIG']['var_def_scripts'])

    # loop over all diagnostics defined in project_info and create/prepare
    # netCDF files for each variable
    for c in project_info['DIAGNOSTICS']:

        # set current diagnostic
        currDiag = project_info['DIAGNOSTICS'][c]

        # Are the requested variables derived from other, more basic, variables?
        requested_vars = currDiag.variables

        # get all models
        project_info['ADDITIONAL_MODELS'] = currDiag.additional_models
        project_info['ALLMODELS'] = project_info['MODELS'] + project_info['ADDITIONAL_MODELS']

        # initialize empty lists to hold preprocess cubes and file paths for each model
        models_cubes = []
        models_fullpaths = []

        # prepare/reformat model data for each model
        for model in project_info['ALLMODELS']:
            model_name = model['name']
            project_name = model['project']
            logger.info("MODEL = %s (%s)", model_name, project_name)

            # start calling preprocess
            op = pp.Diag()

            # FIX-ME old packaging of variable objects (legacy from initial 
            # version), this needs to be changed once we have the new variable
            # definition codes in place
            variable_defs_base_vars = op.add_base_vars_fields(
                requested_vars, model,
                project_info['CONFIG']['var_def_scripts'])

            # if not all variable_defs_base_vars are available, try to fetch
            # the target variable directly (relevant for derived variables)
            base_vars = op.select_base_vars(variable_defs_base_vars, model,
                                            currDiag, project_info)

            # process base variables
            for base_var in base_vars:
                if project_info_0.CONFIG['var_only_case'] > 0:
                    if op.id_is_explicitly_excluded(base_var, model):
                        continue
                logger.info("VARIABLE = %s (%s)", base_var.name, base_var.field)

                # rewrite netcdf to expected input format.
                logger.info("Calling preprocessor")
                # REFORMAT: for backwards compatibility we can revert to ncl reformatting
                # by changing cmor_reformat_type = 'ncl'
                # for python cmor_check one, use cmor_reformat_type = 'py'
                # PREPROCESS ID: extracted from variable dictionary
                if hasattr(base_var, 'preproc_id'):
                    try:
                        preprocess_id = base_var.preproc_id
                        for preproc_dict in project_info_0.PREPROCESS:
                            if preproc_dict['id'] == preprocess_id:
                                logger.info('Preprocess id: %s', preprocess_id)
                                project_info['PREPROCESS'] = preproc_dict
                    except AttributeError:
                        logger.error('preprocess_id is not an attribute of variable object. Exiting...')
                        sys.exit(1)
                prep = pp.preprocess(project_info, base_var, model, currDiag, cmor_reformat_type = 'py')

                # add only if we need multimodel statistics
                if project_info['PREPROCESS']['multimodel_mean'] is True:
                    models_cubes.append(prep[0])
                    models_fullpaths.append(prep[1])

        # before we proceed more, we call multimodel operations
        if project_info['PREPROCESS']['multimodel_mean'] is True:
            pp.multimodel_mean(models_cubes, models_fullpaths)

        vardicts = currDiag.variables
        variables = []
        field_types = []
        # hack: needed by diags that
        # perform regridding on ref_model
        # and expect that to be a model attribute
        ref_models = []
        ###############
        for v in vardicts:
            variables.append(v['name'])
            field_types.append(v['field'])
            ref_models.append(v['ref_model'][0])

        project_info['RUNTIME']['currDiag'] = currDiag

        for derived_var, derived_field, refmodel in zip(variables, field_types, ref_models):

            # needed by external diag to perform regridding
            for model in project_info['ALLMODELS']:
                model['ref'] = refmodel
            logger.info("External diagnostic will use ref_model: %s", model['ref'])

            project_info['RUNTIME']['derived_var'] = derived_var
            project_info['RUNTIME']['derived_field_type'] = derived_field

            # iteration over diagnostics for each variable
            scrpts = copy.deepcopy(currDiag.scripts)
            for i in range(len(currDiag.scripts)):

                # because the NCL environment is dumb, it is important to
                # tell the environment which diagnostic script to use
                project_info['RUNTIME']['currDiag'].scripts = [scrpts[i]]

                # this is hardcoded, maybe make it an option
                executable = os.path.join(script_root, "interface_scripts",
                                          "derive_var.ncl")
                logger.info("Calling %s for '%s'", executable, derived_var)
                interface_data = create_interface_data_dir(project_info, executable)
                project_info['RUNTIME']['interface_data'] = interface_data
                pp.run_executable(executable, project_info, verbosity,
                                  project_info['GLOBAL']['exit_on_warning'])

                # run diagnostics...
                # ...unless run_diagnostic is specifically set to False
                if project_info['GLOBAL']['run_diagnostic'] is True:
                    executable = os.path.join(script_root, "diag_scripts",
                                              scrpts[i]['script'])
                    configfile = os.path.join(script_root, scrpts[i]['cfg_file'])
                    logger.info("Running diag_script: %s", executable)
                    logger.info("with configuration file: %s", configfile)
                    interface_data = create_interface_data_dir(project_info, executable)
                    project_info['RUNTIME']['interface_data'] = interface_data                    
                    pp.run_executable(executable,
                                      project_info,
                                      verbosity,
                                      project_info['GLOBAL']['exit_on_warning'],
                                      launcher_arguments=None)

    # delete environment variable
    del os.environ['0_ESMValTool_version']

    #End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info("Ending the Earth System Model Evaluation Tool v%s at time: %s",
                __version__, timestamp2.strftime(timestamp_format))
    logger.info("Time for running namelist was: %s", timestamp2 - timestamp1)

    # Remind the user about reference/acknowledgement file
    logger.info("For the required references/acknowledgements of these "
                "diagnostics see: %s", out_refs)


def run():
    """ Run main, logging any exceptions."""
    try:
        main()
    except:
        logger.exception("Program terminated abnormally, see stack trace "
                         "below for more information", exc_info=True)
        sys.exit(1)
    else:
        logger.info("Run was succesful")


if __name__ == '__main__':
    run()
