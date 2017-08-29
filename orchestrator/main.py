#! /usr/bin/env python

"""
Completely rewritten wrapper to be able to deal with
the new yaml parser and simplified interface_scripts
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
import sys
sys.path.append("./interface_scripts")
import subprocess
import getopt
import logging
import logging.config
import yaml
from auxiliary import info, error, print_header, ncl_version_check
import datetime
import os
import pdb
from yaml_parser import Parser as Ps
import preprocess as pp
import copy
import ConfigParser
import namelistchecks as pchk
import uuid

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
                handler['level'] = console_log_level

    logging.config.dictConfig(cfg)

# Configure logging
#TODO: specify a good output directory for log files
#TODO: (maybe) allow changing the console log level from the command line
configure_logging()

# Get a logger
logger = logging.getLogger('ESMValTool')

logger.debug("Example debug message")
logger.info("Example info message")
logger.warn("Example warning message")
logger.error("Example error message")

# Define ESMValTool version
version = "2.0.0"
os.environ['0_ESMValTool_version'] = version

# Check NCL version
ncl_version_check()

# print usage
def usage():
    msg = """ 
              python main.py [OPTIONS]
              ESMValTool - Earth System Model Evaluation Tool.
              You can run with these command-line options:

              -h --help             [HELP]     Print help message and exit.
              -n --namelist-file    [REQUIRED] Namelist file.
              -c --config-file      [REQUIRED] Configuration file. 

              For further help, check the doc/-folder for pdfs 
              and references therein. Have fun!

    """
    print(msg)

# class to hold information for configuration
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
        GLOB = {}
        if cp.has_option('GLOBAL','ini-version') :
            ini_version = int(cp.get('GLOBAL','ini-version'))
            GLOB['ini-version'] = ini_version
        if cp.has_option('GLOBAL','write_plots') :
            write_plots = self.s2b(cp.get('GLOBAL','write_plots'))
            GLOB['write_plots'] = write_plots
        else:
            logger.warn("no write_plots in config, set to True")
            GLOB['write_plots'] = True
        if cp.has_option('GLOBAL','write_netcdf') :
            write_netcdf = self.s2b(cp.get('GLOBAL','write_netcdf'))
            GLOB['write_netcdf'] = write_netcdf
        else:
            logger.warn("no write_netcdf in config, set to True")
            GLOB['write_netcdf'] = True
        if cp.has_option('GLOBAL','verbosity') :
            verbosity = int(cp.get('GLOBAL','verbosity'))
            GLOB['verbosity'] = verbosity
        else:
            logger.warn("no verbosity in config, set to 1")
            GLOB['verbosity'] = 1
        if cp.has_option('GLOBAL','exit_on_warning') :
            eow = self.s2b(cp.get('GLOBAL','exit_on_warning'))
            GLOB['exit_on_warning'] = eow
        else:
            logger.warn("no exit_on_warning in config, set to True")
            GLOB['exit_on_warning'] = False
        if cp.has_option('GLOBAL','output_file_type') :
            output_type = cp.get('GLOBAL','output_file_type')
            GLOB['output_file_type'] = output_type
        else:
            logger.warn("no output_file_type in config, set to ps")
            GLOB['output_file_type'] = 'ps'
        if cp.has_option('GLOBAL','preproc_dir') :
            preproc_dir = cp.get('GLOBAL','preproc_dir')
            GLOB['preproc_dir'] = preproc_dir
        else:
            logger.warn("no preproc_dir in config, set to ./preproc/")
            GLOB['preproc_dir'] = './preproc/'
        if cp.has_option('GLOBAL','work_dir') :
            work_dir = cp.get('GLOBAL','work_dir')
            GLOB['work_dir'] = work_dir
        else:
            logger.warn("no work_dir in config, set to ./work/")
            GLOB['work_dir'] = './work/'
        if cp.has_option('GLOBAL','plot_dir') :
            plot_dir = cp.get('GLOBAL','plot_dir')
            GLOB['plot_dir'] = plot_dir
        else:
            logger.warn("no plot_dir in config, set to ./plots/")
            GLOB['plot_dir'] = './plots/'
        if cp.has_option('GLOBAL','max_data_filesize') :
            mdf = int(cp.get('GLOBAL','max_data_filesize'))
            GLOB['max_data_filesize'] = mdf
        else:
            logger.warn("no max_data_filesize in config, set to 100")
            GLOB['max_data_filesize'] = 100
        if cp.has_option('GLOBAL','run_dir') :
            run_dir = cp.get('GLOBAL','run_dir')
            GLOB['run_dir'] = run_dir
        else:
            logger.warn("no run_dir in config, assuming .")
            GLOB['run_dir'] = '.'
        if cp.has_option('GLOBAL','save_intermediary_cubes') :
            save_intermediary_cubes = self.s2b(cp.get('GLOBAL','save_intermediary_cubes'))
            GLOB['save_intermediary_cubes'] = save_intermediary_cubes
        else:
            logger.warn("no save_intermediary_cubes in config, assuming False")
            GLOB['save_intermediary_cubes'] = False
        if cp.has_option('GLOBAL','model_rootpath'): ## Use this for all classes except the ones for obs_rootpath
            mp = cp.get('GLOBAL','model_rootpath')
            GLOB['model_rootpath'] = mp
        else:
            logger.error("Model root path not defined")
            sys.exit(1)
        if cp.has_option('GLOBAL','obs_rootpath'):  ## Use this for OBS, obs4mips, ana4mips classes
            op = cp.get('GLOBAL','obs_rootpath')
            GLOB['obs_rootpath'] = op
        else:
            logger.error("Observations root path not defined")
            sys.exit(1)
        if cp.has_option('GLOBAL','rawobs_rootpath'):  ## For reformat_obs only, to be used later
            rop = cp.get('GLOBAL','rawobs_rootpath')
            GLOB['rawobs_rootpath'] = rop
        if cp.has_option('GLOBAL','cmip5_dirtype') :
            ddd = cp.get('GLOBAL','cmip5_dirtype')
            GLOB['cmip5_dirtype'] = ddd
            permitted_values = ['default', 'badc', 'dkrz', 'ethz', 'smhi', 'None']
            if ddd not in permitted_values:
                logger.error("Unrecognized option for cmip5_dirtype in config: %s", ddd)
                sys.exit(1)
        else:
            logger.warn("no cmip5_dirtype in config, "
                        "assuming None (unstructured data directory)")
            GLOB['cmip5_dirtype'] = 'None'
        return GLOB


print_header()

# start parsing command line args
# options initialize and descriptor
namelist_file    = None
config_file      = None
preprocess_id    = None

shortop = "hp:n:c:"
# ---- Long form.
longop = [
   "help",
   "namelist-file=",
   "config-file="]

# get command line arguments
try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

command_string = 'python main.py '

# parse options
for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(0)
    elif o in ("-n", "--namelist-file"):
        namelist_file = a
        command_string = command_string + ' -n ' + a
    elif o in ("-c", "--config-file"):
        config_file = a
        command_string = command_string + ' -c ' + a
    else:
        logger.error("Unknown option: %s", o)
        usage()
        sys.exit(1)

# condition options
if not namelist_file:
    logger.error("No namelist file specified.")
    logger.error("Use --namelist-file to specify it.")
    sys.exit(1)
if not config_file:
    logger.error("No configuration file specified.")
    logger.error("Use --config-file to specify it.")
    sys.exit(1)

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

# if run_dir exists, don't overwrite it
if os.path.isdir(project_info['GLOBAL']['run_dir']):
    suf = uuid.uuid4().hex
    newdir = project_info['GLOBAL']['run_dir'] + '_' + suf
    mvd = 'mv ' + project_info['GLOBAL']['run_dir'] + ' ' + newdir
    proc = subprocess.Popen(mvd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    logger.info('Renamed existing run directory to %s', newdir)

# tell the environment about regridding
project_info['RUNTIME']['regridtarget'] = []

# Master references-acknowledgements file (hard coded)
in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
project_info['RUNTIME']['in_refs'] = in_refs

# Open refs-acknows file in run_dir (delete if existing)
if not os.path.isdir(project_info['GLOBAL']['run_dir']):
    mkd = 'mkdir -p ' + project_info['GLOBAL']['run_dir']
    proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    logger.info('Created run directory %s', project_info['GLOBAL']['run_dir'])

if (os.path.isfile(out_refs)):
    os.remove(out_refs)
f = open(out_refs, "w")
f.close()


# Current working directory
project_info['RUNTIME']['cwd'] = os.getcwd()

# Summary to std-out before starting the loop
timestamp1 = datetime.datetime.now()
timestamp_format = "%Y-%m-%d --  %H:%M:%S"

logger.info("Starting the Earth System Model Evaluation Tool v%s at time: %s ...",
            version, timestamp1.strftime(timestamp_format))

# Loop over all diagnostics defined in project_info and
# create/prepare netCDF files for each variable
for c in project_info['DIAGNOSTICS']:
    currDiag = project_info['DIAGNOSTICS'][c]

    # Are the requested variables derived from other, more basic, variables?
    requested_vars = currDiag.variables

    # Prepare/reformat model data for each model
    for model in project_info['MODELS']:
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
                    logger.info('preprocess_id is not an attribute of variable object. Exiting...')
                    sys.exit(1)
            pp.preprocess(project_info, base_var, model, currDiag, cmor_reformat_type = 'py')

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

        # needed by external diag to perform refridding
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

            # run diagnostics
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
del(os.environ['0_ESMValTool_version'])

#End time timing
timestamp2 = datetime.datetime.now()
logger.info("Ending the Earth System Model Evaluation Tool v%s at time: %s",
            version, timestamp2.strftime(timestamp_format))
logger.info("Time for running namelist was: %s", timestamp2 - timestamp1)

# Remind the user about reference/acknowledgement file
logger.info("For the required references/acknowledgements of these diagnostics see: ")
logger.info(out_refs)
