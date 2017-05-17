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
# 2015-06-26  SR   Added optional <ESGF> section quality check and print
# 2015-06-30  a_laue_ax: added ESMValTool version

import sys
sys.path.append("./interface_scripts")

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
import reformat
import xml.sax
import xml_parsers

# Define ESMValTool version
version = "1.1.0"
os.environ['0_ESMValTool_version'] = version

# Check NCL version
ncl_version_check()

# Check command arguments.
usage = "%prog nml/namelist-file.xml"
description = """ESMValTool - Earth System Model Evaluation Tool.
For further help, check the doc/-folder for pdfs and references therein."""

parser = OptionParser(usage=usage, description=description)
parser.add_option("-r", "--reformat",
                  action="store_true", dest="reformat", default=False,
                  help="run reformat scripts for the observations according to namelist")
options, args = parser.parse_args()
if len(args) == 0:
    parser.print_help()
    sys.exit(0)

# Get command arguments.
input_xml_full_path = args[0]

# Parse input namelist into project_info-dictionary.
Project = xml_parsers.namelistHandler()
parser = xml.sax.make_parser()
parser.setContentHandler(Project)
parser.parse(input_xml_full_path)

# Project_info is a dictionary with all info from the namelist.
project_info = Project.project_info

if options.reformat:
	if 'REFORMAT' not in project_info.keys():
		error('No REFORMAT tag specified in {0}'.format(input_xml_full_path))
	if len(project_info['REFORMAT']) == 0:
		info('No reformat script specified',1,1)
        print_header({}, options.reformat)
	for k,v in project_info['REFORMAT'].iteritems():
		if not os.path.exists(v):
			error('Path {0} does not exist'.format(v))
		projects.run_executable(v,
	                                project_info,
	                                1,
	                                False,
	                                write_di=False)
	sys.exit(0)

verbosity = project_info['GLOBAL']['verbosity']
climo_dir = project_info['GLOBAL']['climo_dir']
exit_on_warning = project_info['GLOBAL'].get('exit_on_warning', False)

# Additional entries to 'project_info'. The 'project_info' construct
# is one way by which Python passes on information to the NCL-routines.
project_info['RUNTIME'] = {}

# Input xml path/file
project_info['RUNTIME']['xml'] = input_xml_full_path
input_xml_file = os.path.basename(input_xml_full_path)
project_info['RUNTIME']['xml_name'] = input_xml_file

# Master references-acknowledgements file (hard coded)
in_refs = os.path.join(os.getcwd(), 'doc/MASTER_authors-refs-acknow.txt')
project_info['RUNTIME']['in_refs'] = in_refs

# Create refs-acknows file in workdir (delete if existing)
wrk_dir = project_info['GLOBAL']['wrk_dir']
if not os.path.isdir(wrk_dir):
    os.mkdir(wrk_dir)

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

print_header(project_info, options.reformat)
info("Starting the Earth System Model Evaluation Tool v" + version + " at time: "
     + timestamp1.strftime(timestamp_format) + "...", verbosity, 1)

# Load ESGF config info (if specified in namelist)
print project_info['ESGF']
if 'ESGF' in project_info and 'config_file' in project_info['ESGF']:
    esgf_config_file = project_info['ESGF']['config_file']
    if os.path.isfile(esgf_config_file):
        info("Loading ESGF config file", verbosity, 2)
        esgf_config_handler = xml_parsers.ESGFConfigHandler()
        parser = xml.sax.make_parser()
        parser.setContentHandler(esgf_config_handler)
        parser.parse(esgf_config_file)
        esgf_config = esgf_config_handler.current_tag.config
        esgf_config.config_file_name = esgf_config_file
        info("ESGF config file loaded", verbosity, 2)

        # Summary of ESGF config info (if exists) to std-out
        info('', verbosity, 3)
        info('Summary of ESGF config information:',\
            verbosity, 3)
        msg = str(esgf_config)
        for msg_line in msg.split('\n'):
            info(msg_line, verbosity, 3)
        info("", verbosity, 3)

        # Perform quality check on ESGF config
        projects.ESGF.quality_check(esgf_config)
        info("ESGF config file passed quality check.", verbosity, 2)

        # Store ESGF config in project_info
        project_info['ESGF']['config'] = esgf_config

        # Also, add fullpath of namelist file to
        # ESGF section of project info, so it can be
        # included in report/console output
        project_info['ESGF']['namelist_fullpath']\
            = input_xml_full_path

    else:
        msg = "Cannot find ESGF config file '%s'" % esgf_config_file
        raise IOError(msg)

# Loop over all diagnostics defined in project_info and
# create/prepare netCDF files for each variable
for currDiag in project_info['DIAGNOSTICS']:

    # Are the requested variables derived from other, more basic, variables?
    requested_vars = currDiag.get_variables_list()

    # Update currDiag-specific models
    project_info['MODELS'] = projects.remove_diag_specific_models(
        project_info['MODELS'])
    diag_specific_models = currDiag.get_diag_models()
    projects.add_model(project_info, diag_specific_models)

    # Prepare/reformat model data for each model
    for model in project_info['MODELS']:
        currProject = getattr(vars()['projects'], model.split_entries()[0])()
        model_name = currProject.get_model_name(model)
        project_name = currProject.get_project_name(model)
        info("", verbosity, 1)
        info("MODEL = " + model_name + " (" + project_name + ")", verbosity, 1)

        # variables needed for target variable, according to variable_defs
        variable_defs_base_vars = currDiag.add_base_vars_fields(requested_vars, model)
        # if not all variable_defs_base_vars are available, try to fetch
        # the target variable directly (relevant for derived variables)
        base_vars = currDiag.select_base_vars(variable_defs_base_vars,
                                              model,
                                              currProject,
                                              project_info)

        # process base variables
        for base_var in base_vars:
            if currDiag.id_is_explicitly_excluded(base_var, model):
                continue
            info("VARIABLE = " + base_var.var + " (" + base_var.fld + ")",
                 verbosity, 1)

            # Rewrite netcdf to expected input format.
            info("Calling cmor_reformat.py to check/reformat model data",
                 verbosity, 2)
            reformat.cmor_reformat(currProject, project_info, base_var, model)

    variables = currDiag.get_variables()
    field_types = currDiag.get_field_types()

    project_info['RUNTIME']['currDiag'] = currDiag
    for derived_var, derived_field in zip(variables, field_types):
        project_info['RUNTIME']['derived_var'] = derived_var
        project_info['RUNTIME']['derived_field_type'] = derived_field

        executable = "./interface_scripts/derive_var.ncl"
        info("", verbosity, required_verbosity=1)
        info("Calling " + executable + " for '" + derived_var + "'",
             verbosity, required_verbosity=1)
        projects.run_executable(executable, project_info, verbosity,
                                exit_on_warning)
    project_info['RUNTIME']['derived_var'] = "Undefined"

    executable = "./diag_scripts/" + currDiag.get_diag_script()
    configfile = currDiag.get_diag_script_cfg()
    info("", verbosity, required_verbosity=1)
    info("Running diag_script: " + executable, verbosity, required_verbosity=1)
    info("with configuration file: " + configfile, verbosity,
         required_verbosity=1)

    projects.run_executable(executable,
                            project_info,
                            verbosity,
                            exit_on_warning,
                            launcher_arguments=currDiag.get_launcher_arguments())

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
