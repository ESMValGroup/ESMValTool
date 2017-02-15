#! /usr/bin/env python

# Parse namelist and list + check that expected input files
# exists.
#
# This scripts requires the "projects"-module so make sure
# it is executed from the correct path (ESMValTool root)

import glob
import os
import pdb
import re
import xml.sax

from argparse import ArgumentParser

import sys
sys.path.append("./interface_scripts")
import projects
import reformat
import xml_parsers


def get_first_year(infile):
    regex = re.compile('.*_([0-9][0-9][0-9][0-9])[0-9][0-9][0-9]*-[0-9][0-9][0-9][0-9][0-9][0-9][0-9]*.nc')
    return min(get_year(infile, regex))


def get_last_year(infile):
    regex = re.compile('.*_[0-9][0-9][0-9][0-9][0-9][0-9][0-9]*-([0-9][0-9][0-9][0-9])[0-9][0-9][0-9]*.nc')
    return max(get_year(infile, regex))


def get_year(infile, regex):
    files = glob.glob(infile)
    years = [regex.search(fil).group(1) for fil in files]
    return years


# Check command arguments.
description = """drs-to-esmval
Recursively searches the given <data_root_path> for matching entries
and rewrites these to ESMValTool <model>-tags."""

parser = ArgumentParser(description=description)
parser.add_argument("-n", "--nml",
                    dest="namelist",
                    type=str,
                    required=True,
                    help="Namelist to check")
parser.add_argument("-v", "--validate",
                    action="store_true",
                    default=False,
                    help="Validate that nml files exists on the FS")
parser.add_argument("-l", "--list",
                    action="store_true",
                    default=True,
                    help="List path to nml-expected files")
args = parser.parse_args()


# Parse input namelist into project_info-dictionary.
Project = xml_parsers.namelistHandler()
parser = xml.sax.make_parser()
parser.setContentHandler(Project)
parser.parse(args.namelist)

# Project_info is a dictionary with all info from the namelist.
project_info = Project.project_info
project_info['RUNTIME'] = {}
verbosity = project_info['GLOBAL']['verbosity'] = 0

if args.validate:
    prefix_not_ok = "MISSING but required: "
    prefix_ok = "required and present: "
else:
    prefix_not_ok = ""
    prefix_ok = ""

sys.stdout.write('\n')

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
        project_basename = currProject.get_project_basename()
        project_info['RUNTIME']['model'] = model_name
        project_info['RUNTIME']['project'] = project_name
        project_info['RUNTIME']['project_basename'] = project_basename

        # variables needed for target variable, according to variable_defs
        variable_defs_base_vars = currDiag.add_base_vars_fields(requested_vars, model)

        # if not all variable_defs_base_vars are available, try to fetch
        # the target variable directly (relevant for derived variables)
        for base_var in variable_defs_base_vars:
            # Check if this variable should be excluded for this model-id
            if currDiag.id_is_explicitly_excluded(base_var, model):
                continue

            # first try: use base variables provided by variable_defs script
            os.environ['__ESMValTool_base_var'] = base_var.var
            infile = reformat.infile(currProject,
                                     project_info,
                                     base_var,
                                     model)

            nml_syear = currProject.get_model_subsection(model, 'start_year')
            nml_eyear = currProject.get_model_subsection(model, 'end_year')

            if args.validate:
                year_string = ""
                if len(glob.glob(infile)) != 0:
                    fs_syear = get_first_year(infile)
                    fs_eyear = get_last_year(infile)

                    if nml_syear < fs_syear:
                        year_string = "requested initial year "\
                                      + str(nml_syear)\
                                      + " < "\
                                      + str(fs_syear)

                    if nml_eyear > fs_eyear:
                        if len(year_string) != 0:
                            year_string = year_string + ", "

                        year_string = year_string\
                                      + "requested final year "\
                                      + str(nml_eyear)\
                                      + " > "\
                                      + str(fs_eyear)

                    if len(year_string) == 0:
                        sys.stdout.write(prefix_ok + infile + '\n')
                    else:
                        sys.stdout.write(prefix_not_ok + infile)
                        if len(year_string) != 0:
                            sys.stdout.write(" (" + year_string + " )")
                        sys.stdout.write('\n')

                else:
                    sys.stdout.write(prefix_not_ok + infile + '\n')
            else:
                sys.stdout.write(infile + '\n')
