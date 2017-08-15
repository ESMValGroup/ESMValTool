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
from auxiliary import info, error, print_header, ncl_version_check
from optparse import OptionParser
import datetime
import os
import pdb
from yaml_parser import Parser as Ps
import preprocess as pp
import copy

# Define ESMValTool version
version = "2.0.0"
os.environ['0_ESMValTool_version'] = version

# Check NCL version
ncl_version_check()

# Check command arguments.
usage = "%prog nml/namelist-file.yml"
description = """ESMValTool - Earth System Model Evaluation Tool.
For further help, check the doc/-folder for pdfs and references therein."""

parser = OptionParser(usage=usage, description=description)
parser.add_option("-d", "--dummy",
                  action="store_true", dest="dummy", default=False,
                  help="dummy: does nothing")
options, args = parser.parse_args()
if len(args) == 0:
    parser.print_help()
    sys.exit(0)

# Get command arguments.
yml_path = args[0]

# Parse input namelist into project_info-dictionary.
Project = Ps()

# Project_info is a dictionary with all info from the namelist.
project_info_0 = Project.load_namelist(yml_path)
verbosity = project_info_0.GLOBAL['verbosity']
climo_dir = project_info_0.GLOBAL['climo_dir']
exit_on_warning = project_info_0.GLOBAL.get('exit_on_warning', False)

# Project_info is a dictionary with all info from the namelist.
project_info = {}
project_info['GLOBAL'] = project_info_0.GLOBAL
project_info['MODELS'] = project_info_0.MODELS
project_info['DIAGNOSTICS'] = project_info_0.DIAGNOSTICS
project_info['PREPROCESS'] = project_info_0.PREPROCESS
project_info['CONFIG'] = project_info_0.CONFIG

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
            pp.preprocess(project_info, base_var, model, currDiag)

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
