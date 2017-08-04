#! /usr/bin/env python

# main python script
#
# 2008-06-27  CAF
# 2008-08-27  CAF  Added checks for existing files, fixed bug when variable
#                  name is changed U->ua
# 2008-11-06  CAF  added ensemble id, output dir for plots, and project name
# 2008-11-18  CAF  added climo and plot directories.
# 2008-11-20  CAF  added support for multiple required variables for
#                  derived fields.
# 2008-12-01  CAF  added namelist option on command line, switched
#                  system call to ncl to popen
#                  skip variable if it's missing
# 2009-01-27  CAF  changed climo dir structure
# 2009-03-04  CAF  added logic to handle different file numbers for
#                  required variables
# 2009-05-22  CAF  changed method of field_number creation, var_att
#                  files can now have *2*z
#                  *2*s *3* formatted Requires line
# 2009-05-27  CAF  seach for filenames w/wo date range
# 2009-06-12  CAF  Don't do climatology for I and D
# 2010-12-02  HS   added environment variables to set maximum data sizes
# 2010-12-14  HS   added termination call when error occurs with reformat
# 2012-06-08  HS   set environment variable (OPTIONS_FILE) if
#                  namelist.global_vars.write_plot_vars is defined
# 2015-06-30  a_laue_ax: added ESMValTool version

import sys
sys.path.append("./interface_scripts")
import preprocess_0 as preproc

# path to ESMVal python toolbox library
sys.path.append("./diag_scripts/lib/python")
sys.path.append("./diag_scripts")

from auxiliary import info, error, print_header, ncl_version_check
## from climate import climate
from optparse import OptionParser
import datetime
import projects
import os
import pdb
from parser_0 import Parser as Ps
import subprocess
import re


# Define ESMValTool version
version = "X.X.0"
os.environ['0_ESMValTool_version'] = version

# Check NCL version
ncl_version_check()

# Check command arguments.
usage = "%prog nml/namelist-file.ini"
description = """ESMValTool - Earth System Model Evaluation Tool.
For further help, check the doc/-folder for pdfs and references therein."""

parser = OptionParser(usage=usage, description=description)
options, args = parser.parse_args()
if len(args) == 0:
    parser.print_help()
    sys.exit(0)

# Get command arguments.
yml_path = args[0]

# Parse input namelist into project_info-dictionary.
Project = Ps()

# Project_info is a dictionary with all info from the namelist.
project_info = Project.load_namelist(yml_path)
verbosity = project_info.GLOBAL['verbosity']
climo_dir = project_info.GLOBAL['preproc_dir']
exit_on_warning = project_info.GLOBAL.get('exit_on_warning', False)

# Additional entries to 'project_info'. The 'project_info' construct
# is one way by which Python passes on information to the NCL-routines.
project_info.CONFIG['RUNTIME'] = {}

# Input yml path/file
project_info.CONFIG['RUNTIME']['yml'] = yml_path
input_yml_file = os.path.basename(yml_path)
project_info.CONFIG['RUNTIME']['yml_name'] = input_yml_file

# Master references-acknowledgements file (hard coded)
in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
project_info.CONFIG['RUNTIME']['in_refs'] = in_refs

# Create refs-acknows file in workdir (delete if existing)
wrk_dir = project_info.GLOBAL['wrk_dir']
if not os.path.isdir(wrk_dir):
    os.mkdir(wrk_dir)

# Prepare writing of references/acknowledgementes to file
refs_acknows_file = str.replace(input_yml_file, "namelist_", "refs-acknows_")
refs_acknows_file = refs_acknows_file.split(os.extsep)[0] + ".log"

out_refs = os.path.join(wrk_dir, refs_acknows_file)
if (os.path.isfile(out_refs)):
    os.remove(out_refs)
f = open(out_refs, "w")
f.close()
project_info.CONFIG['RUNTIME']['out_refs'] = out_refs

# Current working directory
project_info.CONFIG['RUNTIME']['cwd'] = os.getcwd()

# Summary to std-out before starting the loop
timestamp1 = datetime.datetime.now()
timestamp_format = "%Y-%m-%d --  %H:%M:%S"

print_header(project_info,'testing new parser')
info("Starting the Earth System Model Evaluation Tool v" + version + " at time: "
     + timestamp1.strftime(timestamp_format) + "...", verbosity, 1)

# Loop over all diagnostics defined in project_info
for c in project_info.DIAGNOSTICS:
    currDiag = project_info.DIAGNOSTICS[c]
    # Are the requested variables derived from other, more basic, variables?
    requested_var = currDiag['variable']
    # Prepare/reformat model data for each model
    for model in project_info.MODELS:
        model_name = model['name']
        project_name = model['project']
        info("", verbosity, 1)
        info("MODEL = " + model_name + " (" + project_name + ")", verbosity, 1)
        info("VARIABLE NAME = " + requested_var['name'], verbosity, 1)
        info("VARIABLE FIELD = " + requested_var['field'], verbosity, 1)
        # variables needed for target variable, according to variable_defs
        """
        V Predoi: NOTE1
        This bit about reformatting variables should be done here
        with an option in the ini file [GLOBAL] reformat=True/False
        """
        #STARTNOTE1
        #variable_defs_base_vars = reformat_var.add_base_vars_fields(currDiag, requested_vars, model)
        # if not all variable_defs_base_vars are available, try to fetch
        # the target variable directly (relevant for derived variables)
        #base_vars = reformat_var.select_base_vars(currDiag,variable_defs_base_vars,
        #                                      model,
        #                                      currProject,
        #                                      project_info)
        #ENDNOTE1

        # process base variables
        #for base_var in currDiag['variable']:
        """
        V Predoi: NOTE2
        Do we still need this bit?
        """
            #STARTNOTE2
            #if currDiag.id_is_explicitly_excluded(base_var, model):
            #    continue
            #info("VARIABLE = " + base_var.var + " (" + base_var.fld + ")",
            #     verbosity, 1)
            #ENDNOTE2

            # Rewrite netcdf to expected input format.
        info("Calling cmor_reformat.py to check/reformat model data",
             verbosity, 2)
        infiles = preproc.get_cmip_cf_infile(project_info, currDiag, requested_var['name'])
        info('NETCDF FILES LIST:' + str(infiles), verbosity, 1)
        info('PREPROCESSING ID: ' + currDiag['preprocess']['id'],  verbosity, 1)
        if currDiag['preprocess'] != 'None':
            #preproc.cmor_reformat(project_info, project_info['PREPROCESS'][ppidx]['id'])
            print('should do preprocessing here....')   

    project_info.CONFIG['RUNTIME']['variable_def_dir'] = './variable_defs'
    project_info.CONFIG['RUNTIME']['currDiag'] = currDiag
    project_info.CONFIG['RUNTIME']['derived_var'] = requested_var['name']
    project_info.CONFIG['RUNTIME']['derived_field_type'] = requested_var['field']
    with open('tempvars.ncl','w') as file:
        file.write('variable_def_dir=' + '"./variable_defs"' + '\n')
        file.write('derived_var=' + '"' + requested_var['name'] + '"' + '\n')
        file.write('variables=' + '"' + requested_var['name'] + '"' + '\n')
    executable = "./interface_scripts/derive_var.ncl"
    info("", verbosity, required_verbosity=1)
    info("Calling " + executable + " for '" + requested_var['name'] + "'",
         verbosity, required_verbosity=1)
    # just use launchers here
    import launchers_0 as lu
    lunch = lu.ncl_launcher()
    lunch.execute(executable, project_info.CONFIG, verbosity,exit_on_warning) 

    # diag running
    ppidx = project_info['DIAGNOSTICS'].index(currDiag)
    executable = currDiag['scripts'][ppidx]['script']
    configfile = currDiag['scripts'][ppidx]['cfg_file']
    info("", verbosity, required_verbosity=1)
    info("Running diag_script: " + executable, verbosity, required_verbosity=1)
    info("with configuration file: " + configfile, verbosity,
         required_verbosity=1)
    nclexecute(executable, project_info, verbosity,exit_on_warning)

# delete environment variable
del(os.environ['0_ESMValTool_version'])

#End time timing
timestamp2 = datetime.datetime.now()
info("", verbosity, 1)
info("Ending the Earth System Model Evaluation Tool v" + version + " at time: "
     + timestamp2.strftime(timestamp_format), verbosity, 1)
info("Time for running namelist was: " + str(timestamp2 - timestamp1), verbosity, 1)

# Remind the user about reference/acknowledgement file
info("", verbosity, 1)
info("For the required references/acknowledgements of these diagnostics see: ",
     verbosity, 1)
info(out_refs, verbosity, 1)
