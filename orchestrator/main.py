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

# Completely rewritten wrapper to be able to deal with
# the new yaml parser and simplified interface_scripts
# toolbox. Author: Valeriu Predoi, University of Reading,
# Initial version: August 2017
# contact: valeriu.predoi@ncas.ac.uk

import argparse
import ConfigParser
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
        cfg_file = os.path.join(os.path.dirname(__file__), 'logging.yml')

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


class configFile:
    """
    Class to hold info from the configuration file (if present)
    """
    
    def s2b(self, s):
        """
        convert a string to boolean arg
        """
        if s == 'True':
            return True
        elif s == 'False':
            return False
        else:
            raise ValueError

    def get_par_file(self, params_file):
        """
        Gets the params file
        """
        # ---- Check the params_file exists
        if not os.path.isfile(params_file):
            logger.error("non existent configuration file %s", params_file)
        cp = ConfigParser.ConfigParser()
        cp.read(params_file)
        return cp

    # parse GLOBAL section
    # the beauty of ConfigParser is that
    # one can add sections, so here is
    # where you define them
    def GLOBAL(self, params_file):
        """
        Function to build dictionary containing
        the GLOBAL attributes. Takes ConfigParser object cp
        """
        cp = self.get_par_file(params_file)
        glob = {}
        if cp.has_option('GLOBAL','ini-version') :
            ini_version = int(cp.get('GLOBAL','ini-version'))
            glob['ini-version'] = ini_version
        if cp.has_option('GLOBAL','write_plots') :
            write_plots = self.s2b(cp.get('GLOBAL','write_plots'))
            glob['write_plots'] = write_plots
        else:
            logger.warning("no write_plots in config, set to True")
            glob['write_plots'] = True
        if cp.has_option('GLOBAL','write_netcdf') :
            write_netcdf = self.s2b(cp.get('GLOBAL','write_netcdf'))
            glob['write_netcdf'] = write_netcdf
        else:
            logger.warning("no write_netcdf in config, set to True")
            glob['write_netcdf'] = True
        if cp.has_option('GLOBAL','log_level') :
            glob['log_level'] = cp.get('GLOBAL','log_level').upper()
        else:
            glob['log_level'] = 'INFO'
        if cp.has_option('GLOBAL','verbosity') :
            #TODO: remove verbosity once it is no longer used
            verbosity = int(cp.get('GLOBAL','verbosity'))
            glob['verbosity'] = verbosity
            logger.warning("Option verbosity in %s is deprecated and will be "
                           "removed in the future, use log_level instead",
                           params_file)
        else:
            glob['verbosity'] = 1
        if cp.has_option('GLOBAL','exit_on_warning') :
            eow = self.s2b(cp.get('GLOBAL','exit_on_warning'))
            glob['exit_on_warning'] = eow
        else:
            logger.warning("no exit_on_warning in config, set to True")
            glob['exit_on_warning'] = False
        if cp.has_option('GLOBAL','output_file_type') :
            output_type = cp.get('GLOBAL','output_file_type')
            glob['output_file_type'] = output_type
        else:
            logger.warning("no output_file_type in config, set to ps")
            glob['output_file_type'] = 'ps'
        if cp.has_option('GLOBAL','preproc_dir') :
            preproc_dir = cp.get('GLOBAL','preproc_dir')
            glob['preproc_dir'] = preproc_dir
        else:
            logger.warning("no preproc_dir in config, set to ./preproc/")
            glob['preproc_dir'] = './preproc/'
        if cp.has_option('GLOBAL','work_dir') :
            work_dir = cp.get('GLOBAL','work_dir')
            glob['work_dir'] = work_dir
        else:
            logger.warning("no work_dir in config, set to ./work/")
            glob['work_dir'] = './work/'
        if cp.has_option('GLOBAL','plot_dir') :
            plot_dir = cp.get('GLOBAL','plot_dir')
            glob['plot_dir'] = plot_dir
        else:
            logger.warning("no plot_dir in config, set to ./plots/")
            glob['plot_dir'] = './plots/'
        if cp.has_option('GLOBAL','max_data_filesize') :
            mdf = int(cp.get('GLOBAL','max_data_filesize'))
            glob['max_data_filesize'] = mdf
        else:
            logger.warning("no max_data_filesize in config, set to 100")
            glob['max_data_filesize'] = 100
        if cp.has_option('GLOBAL','run_dir') :
            run_dir = cp.get('GLOBAL','run_dir')
            glob['run_dir'] = run_dir
        else:
            logger.warning("no run_dir in config, assuming .")
            glob['run_dir'] = '.'
        if cp.has_option('GLOBAL','save_intermediary_cubes') :
            save_intermediary_cubes = self.s2b(cp.get('GLOBAL','save_intermediary_cubes'))
            glob['save_intermediary_cubes'] = save_intermediary_cubes
        else:
            logger.warning("no save_intermediary_cubes in config, assuming False")
            glob['save_intermediary_cubes'] = False
        if cp.has_option('GLOBAL','rootpath_CMIP5'):
            mp = cp.get('GLOBAL','rootpath_CMIP5')
            glob['rootpath_CMIP5'] = mp
        else:
            logger.error("Model root path for CMIP5 not defined")
            sys.exit(1)
        if cp.has_option('GLOBAL','rootpath_OBS'):
            op = cp.get('GLOBAL','rootpath_OBS')
            glob['rootpath_OBS'] = op
        else:
            logger.error("Observations root path not defined")
            sys.exit(1)
        if cp.has_option('GLOBAL','drs_CMIP5'):
            drs = cp.get('GLOBAL','drs_CMIP5')
            glob['drs_CMIP5'] = drs
            # check if drs_CMIP5 has host_root fro BADC and DKRZ
            if drs == 'BADC' or drs == 'DKRZ':
                if cp.has_option('GLOBAL','host_root'):
                    hroot = cp.get('GLOBAL','host_root')
                    glob['host_root'] = hroot
                else:
                    logger.error("Using database DRS for CMIP5: database root path not defined")
                    sys.exit(1)
        else:
            logger.warning("No DRS specifier, assuming None and continuing")
        if cp.has_option('GLOBAL','run_diagnostic') :
            run_diagnostic = self.s2b(cp.get('GLOBAL','run_diagnostic'))
            glob['run_diagnostic'] = run_diagnostic
        else:
            logger.warning("no run_diagnostic in config, assuming True")
            glob['run_diagnostic'] = True
        return glob


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

    # Log header
    logger.info(__doc__)

    if run_dir is None:
        logger.warning("Failed to retrieve run_dir from config file")

    if previous_run_dir:
        logger.info('Renamed existing run directory %s to %s',
                    run_dir, previous_run_dir)

    # Check NCL version
    ncl_version_check()

    process_namelist(namelist_file=namelist_file, config_file=config_file)


def process_namelist(namelist_file, config_file):
    """Process namelist"""

    if not os.path.isfile(namelist_file):
        raise OSError(errno.ENOENT, "Specified namelist file does not exist", namelist_file)

    if not os.path.isfile(config_file):
        raise OSError(errno.ENOENT, "Specified config file does not exist", config_file)

    os.environ['0_ESMValTool_version'] = __version__
    # we may need this
    # os.environ["ESMValTool_force_calc"] = "True"

    preprocess_id = None

    # Get namelist file
    yml_path = namelist_file

    # Parse input namelist into project_info-dictionary.
    Project = Ps()
    
    # Parse config file into GLOBAL_DICT info_dictionary
    confFileClass = configFile()
    GLOBAL_DICT = confFileClass.GLOBAL(config_file)
    
    # Project_info is a dictionary with all info from the namelist.
    project_info_0 = Project.load_namelist(yml_path)
    verbosity = GLOBAL_DICT['verbosity']
    preproc_dir = GLOBAL_DICT['preproc_dir']
    exit_on_warning = GLOBAL_DICT.get('exit_on_warning', False)
    
    # Project_info is a dictionary with all info from the namelist.
    project_info = {}
    project_info['GLOBAL'] = GLOBAL_DICT
    project_info['MODELS'] = project_info_0.MODELS
    project_info['DIAGNOSTICS'] = project_info_0.DIAGNOSTICS
    
    # Additional entries to 'project_info'. The 'project_info' construct
    # is one way by which Python passes on information to the NCL-routines.
    project_info['RUNTIME'] = {}
    project_info['RUNTIME']['yml'] = yml_path
    project_info['RUNTIME']['yml_name'] = os.path.basename(yml_path)
    
    # Set references/acknowledgement file
    refs_acknows_file = str.replace(project_info['RUNTIME']['yml_name'], "namelist_", "refs-acknows_")
    refs_acknows_file = refs_acknows_file.split(os.extsep)[0] + ".log"
    out_refs = os.path.join(project_info["GLOBAL"]['run_dir'], refs_acknows_file)
    project_info['RUNTIME']['out_refs'] = out_refs
    
    # Print summary
    logger.info("NAMELIST   = %s", project_info['RUNTIME']['yml_name'])
    logger.info("RUNDIR     = %s", project_info["GLOBAL"]['run_dir'])
    logger.info("WORKDIR    = %s", project_info["GLOBAL"]["work_dir"])
    logger.info("PREPROCDIR = %s", project_info["GLOBAL"]["preproc_dir"])
    logger.info("PLOTDIR    = %s", project_info["GLOBAL"]["plot_dir"])
    logger.info("LOGFILE    = %s", project_info['RUNTIME']['out_refs'])
    logger.info(70 * "_")
    
    #    logger.info("REFORMATTING THE OBSERVATIONAL DATA...", vv, 1)
    
    # perform options integrity checks
    logger.info('Checking integrity of namelist')
    tchk1 = datetime.datetime.now()
    pchk.models_checks(project_info['MODELS'])
    pchk.diags_checks(project_info['DIAGNOSTICS'])
    pchk.preprocess_checks(project_info_0.PREPROCESS)
    tchk2 = datetime.datetime.now()
    dtchk = tchk2 - tchk1
    logger.info('Namelist check successful! Time: %s', dtchk)
    
    # this will have to be purget at some point in the future
    project_info['CONFIG'] = project_info_0.CONFIG
    
    # tell the environment about regridding
    project_info['RUNTIME']['regridtarget'] = []
    
    # Master references-acknowledgements file (hard coded)
    in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
    project_info['RUNTIME']['in_refs'] = in_refs
    
    # Create workdir
    if not os.path.isdir(project_info['GLOBAL']['work_dir']):
        logger.info('Creating work_dir %s', project_info['GLOBAL']['work_dir'])
        os.makedirs(project_info['GLOBAL']['work_dir'])

    # Create rundir
    if not os.path.isdir(project_info['GLOBAL']['run_dir']):
        logger.info('Creating run_dir %s', project_info['GLOBAL']['run_dir'])
        os.makedirs(project_info['GLOBAL']['run_dir'])

    # Create refs-acknows file in run_dir (empty if existing)
    with open(out_refs, "w"):
        pass

    # Current working directory
    project_info['RUNTIME']['cwd'] = os.getcwd()
    
    # Summary to std-out before starting the loop
    timestamp1 = datetime.datetime.now()
    timestamp_format = "%Y-%m-%d --  %H:%M:%S"
    
    logger.info("Starting the Earth System Model Evaluation Tool v%s at time: %s ...",
                __version__, timestamp1.strftime(timestamp_format))
    
    # Loop over all diagnostics defined in project_info and
    # create/prepare netCDF files for each variable
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
    
        # Prepare/reformat model data for each model
        for model in project_info['ALLMODELS']:
            #currProject = model['project']
            model_name = model['name']
            project_name = model['project']
            logger.info("MODEL = %s (%s)", model_name, project_name)
    
            # variables needed for target variable, according to variable_defs
            var_def_dir = project_info_0.CONFIG['var_def_scripts']
    
            # start calling preprocess
            op = pp.Diag()
    
            # old packaging of variable objects (legacy from intial version)
            # this needs to be changed once we have the new variable definition
            # codes in place
            variable_defs_base_vars = op.add_base_vars_fields(requested_vars, model, var_def_dir)
            # if not all variable_defs_base_vars are available, try to fetch
            # the target variable directly (relevant for derived variables)
            base_vars = op.select_base_vars(variable_defs_base_vars,
                                                  model,
                                                  currDiag,
                                                  project_info)
    
            # process base variables
            for base_var in base_vars:
                if project_info_0.CONFIG['var_only_case'] > 0:
                    if op.id_is_explicitly_excluded(base_var, model):
                        continue
                logger.info("VARIABLE = %s (%s)", base_var.name, base_var.field)
    
                # Rewrite netcdf to expected input format.
                logger.info("Calling preprocessing to check/reformat model data, and apply preprocessing steps")
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
                executable = "./interface_scripts/derive_var.ncl"
                logger.info("Calling %s for '%s'", executable, derived_var)
                pp.run_executable(executable, project_info, verbosity,
                                  exit_on_warning)
    
                # run diagnostics...
                # ...unless run_diagnostic is specifically set to False
                if project_info['GLOBAL']['run_diagnostic'] is True:
                    executable = "./diag_scripts/" + scrpts[i]['script']
                    configfile = scrpts[i]['cfg_file']
                    logger.info("Running diag_script: %s", executable)
                    logger.info("with configuration file: %s", configfile)
                    pp.run_executable(executable,
                                      project_info,
                                      verbosity,
                                      exit_on_warning,
                                      launcher_arguments=None)
    
    # delete environment variable
    del os.environ['0_ESMValTool_version']
    
    #End time timing
    timestamp2 = datetime.datetime.now()
    logger.info("Ending the Earth System Model Evaluation Tool v%s at time: %s",
                __version__, timestamp2.strftime(timestamp_format))
    logger.info("Time for running namelist was: %s", timestamp2 - timestamp1)
    
    # Remind the user about reference/acknowledgement file
    logger.info("For the required references/acknowledgements of these "
                "diagnostics see: %s", out_refs)


if __name__ == '__main__':

    try:
        main()
    except:
        logger.exception("Program terminated abnormally, see stack trace "
                         "below for more information", exc_info=True)
        sys.exit(1)
    else:
        logger.info("Run was succesful")
