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
from auxiliary import info, error, print_header, ncl_version_check
import datetime
import os
import pdb
from yaml_parser import Parser as Ps
import preprocess as pp
import copy
import ConfigParser

# Define ESMValTool version
version = "2.0.0"
os.environ['0_ESMValTool_version'] = version

# Check NCL version
ncl_version_check()

# print usage
def usage():
    msg = """ 
              -----------------------------------------------------------
              ------------------- Starting ESMValTool -------------------
              -----------------------------------------------------------
              python main.py [OPTIONS]
              ESMValTool - Earth System Model Evaluation Tool.
              You can run with these command-line options:

              -h --help             [HELP]     Print help message and exit.
              -n --namelist-file    [REQUIRED] Namelist file.
              -c --config-file      [REQUIRED] Configuration file. 

              For further help, check the doc/-folder for pdfs 
              and references therein. Have fun!
              -------------------------------------------------------------

    """
    print >> sys.stderr, msg

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
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> non existent configuration file: "
            sys.exit(1)
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
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no write_plots in config "
            GLOB['write_plots'] = False
        if cp.has_option('GLOBAL','write_netcdf') :
            write_netcdf = self.s2b(cp.get('GLOBAL','write_netcdf'))
            GLOB['write_netcdf'] = write_netcdf
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no write_netcdf in config "
            GLOB['write_netcdf'] = False
        if cp.has_option('GLOBAL','verbosity') :
            verbosity = int(cp.get('GLOBAL','verbosity'))
            GLOB['verbosity'] = verbosity
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no verbosity in config "
            GLOB['verbosity'] = 1
        if cp.has_option('GLOBAL','exit_on_warning') :
            eow = self.s2b(cp.get('GLOBAL','exit_on_warning'))
            GLOB['exit_on_warning'] = eow
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no exit_on_warning in config "
            GLOB['exit_on_warning'] = False
        if cp.has_option('GLOBAL','output_file_type') :
            output_type = cp.get('GLOBAL','output_file_type')
            GLOB['output_file_type'] = output_type
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no output_file_type in config "
            GLOB['output_file_type'] = 'ps'
        if cp.has_option('GLOBAL','climo_dir') :
            climo_dir = cp.get('GLOBAL','climo_dir')
            GLOB['climo_dir'] = climo_dir
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no climo_dir in config "
            GLOB['climo_dir'] = '.'
        if cp.has_option('GLOBAL','wrk_dir') :
            work_dir = cp.get('GLOBAL','wrk_dir')
            GLOB['wrk_dir'] = work_dir
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no wrk_dir in config "
            GLOB['wrk_dir'] = '.'
        if cp.has_option('GLOBAL','plot_dir') :
            plot_dir = cp.get('GLOBAL','plot_dir')
            GLOB['plot_dir'] = plot_dir
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no plot_dir in config "
            GLOB['plot_dir'] = '.'
        if cp.has_option('GLOBAL','max_data_filesize') :
            mdf = int(cp.get('GLOBAL','max_data_filesize'))
            GLOB['max_data_filesize'] = mdf
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no max_data_filesize in config "
            GLOB['max_data_filesize'] = 100
        if cp.has_option('GLOBAL','force_processing') :
            fop = self.s2b(cp.get('GLOBAL','force_processing'))
            GLOB['force_processing'] = fop
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no force_processing in config "
            GLOB['force_processing'] = True
        if cp.has_option('GLOBAL','preprocess_id') :
            preprocess_id = cp.get('GLOBAL','preprocess_id')
            GLOB['preprocess_id'] = preprocess_id
        else:
            print >> sys.stderr,"PY  WARNING:  >>> main.py >>> no preprocess_id in config "
            GLOB['preprocess_id'] = 'Default'
        return GLOB

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
        print >> sys.stderr, "Unknown option:", o
        usage()
        sys.exit(1)

# condition options
if not namelist_file:
    print >> sys.stderr, "PY  ERROR:  >>> main.py >>> No namelist file specified."
    print >> sys.stderr, "PY  ERROR:  >>> main.py >>> Use --namelist-file to specify it."
    sys.exit(1)
if not config_file:
    print >> sys.stderr, "PY  ERROR:  >>> main.py >>> No configuration file specified."
    print >> sys.stderr, "PY  ERROR:  >>> main.py >>> Use --config-file to specify it."
    sys.exit(1)

# Get namelis file
yml_path = namelist_file

# Parse input namelist into project_info-dictionary.
Project = Ps()

# Parse config file into GLOBAL_DICT info_dictionary
confFileClass = configFile()
GLOBAL_DICT = confFileClass.GLOBAL(config_file)

# Project_info is a dictionary with all info from the namelist.
project_info_0 = Project.load_namelist(yml_path)
verbosity = GLOBAL_DICT['verbosity']
climo_dir = GLOBAL_DICT['climo_dir']
exit_on_warning = GLOBAL_DICT.get('exit_on_warning', False)

# Project_info is a dictionary with all info from the namelist.
project_info = {}
project_info['GLOBAL'] = GLOBAL_DICT
project_info['MODELS'] = project_info_0.MODELS
project_info['DIAGNOSTICS'] = project_info_0.DIAGNOSTICS
# this will have to be purget at some point in the future
project_info['CONFIG'] = project_info_0.CONFIG

# identify which Preprocess to perfom
preprocess_id = 'Preprocess_' + GLOBAL_DICT['preprocess_id']
info('>>> main.py >>> Preprocess type: ' + preprocess_id, verbosity, required_verbosity=1)
project_info['PREPROCESS'] = project_info_0.PREPROCESS[preprocess_id]

# Additional entries to 'project_info'. The 'project_info' construct
# is one way by which Python passes on information to the NCL-routines.
project_info['RUNTIME'] = {}

# tell the environment about regridding
project_info['RUNTIME']['regridtarget'] = []

# Input xml path/file
project_info['RUNTIME']['xml'] = yml_path
input_xml_file = os.path.basename(yml_path)
project_info['RUNTIME']['xml_name'] = input_xml_file

# Master references-acknowledgements file (hard coded)
in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
project_info['RUNTIME']['in_refs'] = in_refs

# Create refs-acknows file in workdir (delete if existing)
wrk_dir = project_info['GLOBAL']['wrk_dir']
if not os.path.isdir(wrk_dir):
    mkd = 'mkdir -p ' + wrk_dir
    proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    info(' >>> main.py >>> Created work directory ' + wrk_dir, verbosity, required_verbosity=1)

# Prepare writing of references/acknowledgementes to file
refs_acknows_file = str.replace(input_xml_file, "namelist_", "refs-acknows_")
refs_acknows_file = refs_acknows_file.split(os.extsep)[0] + ".log"

out_refs = os.path.join(wrk_dir, refs_acknows_file)
if (os.path.isfile(out_refs)):
    os.remove(out_refs)
f = open(out_refs, "w")
f.close()
project_info['RUNTIME']['out_refs'] = out_refs

# Current working directory
project_info['RUNTIME']['cwd'] = os.getcwd()

# Summary to std-out before starting the loop
timestamp1 = datetime.datetime.now()
timestamp_format = "%Y-%m-%d --  %H:%M:%S"

print_header(project_info)
info(" >>> main.py >>> Starting the Earth System Model Evaluation Tool v" + version + " at time: "
     + timestamp1.strftime(timestamp_format) + "...", verbosity, 1)

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
        info(" >>> main.py >>> ", verbosity, 1)
        info(" >>> main.py >>> MODEL = " + model_name + " (" + project_name + ")", verbosity, 1)

        # variables needed for target variable, according to variable_defs
        var_def_dir = project_info_0.CONFIG['var_def_scripts']

        # start calling preprocess
        op = pp.Diag()
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
            info(" >>> main.py >>> VARIABLE = " + base_var.name + " (" + base_var.field + ")",
                 verbosity, 1)

            # Rewrite netcdf to expected input format.
            info(" >>> main.py >>> Calling preprocessing to check/reformat model data, and apply preprocessing steps",
                 verbosity, 2)
            # for backwards compatibility we can revert to ncl reformatting
            # by changing cmor_reformat_type = 'ncl'
            pp.preprocess(project_info, base_var, model, currDiag, cmor_reformat_type = 'py')

    vardicts = currDiag.variables
    variables = []
    field_types = []
    for v in vardicts:
        variables.append(v['name'])
        field_types.append(v['field'])

    project_info['RUNTIME']['currDiag'] = currDiag
    for derived_var, derived_field in zip(variables, field_types):

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
            info(" >>> main.py >>> ", verbosity, required_verbosity=1)
            info(" >>> main.py >>> Calling " + executable + " for '" + derived_var + "'",
                 verbosity, required_verbosity=1)
            pp.run_executable(executable, project_info, verbosity,
                              exit_on_warning)

            # run diagnostics
            executable = "./diag_scripts/" + scrpts[i]['script']
            configfile = scrpts[i]['cfg_file']
            info(" >>> main.py >>> ", verbosity, required_verbosity=1)
            info(" >>> main.py >>> Running diag_script: " + executable, verbosity, required_verbosity=1)
            info(" >>> main.py >>> with configuration file: " + configfile, verbosity,
                 required_verbosity=1)
            pp.run_executable(executable,
                              project_info,
                              verbosity,
                              exit_on_warning,
                              launcher_arguments=None)
            

# delete environment variable
del(os.environ['0_ESMValTool_version'])

#End time timing
timestamp2 = datetime.datetime.now()
info(" >>> main.py >>> ", verbosity, 1)
info(" >>> main.py >>> Ending the Earth System Model Evaluation Tool v" + version + " at time: "
     + timestamp2.strftime(timestamp_format), verbosity, 1)
info(" >>> main.py >>> Time for running namelist was: " + str(timestamp2 - timestamp1), verbosity, 1)

# Remind the user about reference/acknowledgement file
info(" >>> main.py >>> ", verbosity, 1)
info(" >>> main.py >>> For the required references/acknowledgements of these diagnostics see: ",
     verbosity, 1)
info(" >>> main.py >>> " + out_refs, verbosity, 1)
