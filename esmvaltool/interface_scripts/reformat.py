#!/usr/bin/python

# Converts netcdf file CF 1.1 complaint netcdf file

# 2008/08/22 C. Fischer

# Changes
# 2008/08/27 CAF  Added check to see if files exist, if they do then
#                 don't process
# 2008/11/06 CAF  Changed shell call from csh script to ncl script,
#                 added ensemble id
# 2008-11-18 CAF  Added climo directory
# 2008-11-18 CAF  replace system call with popen
# 2009-01-27 CAF  added end_year and start_year to enviroment variables
#                 in order to create output filename.
#                 also changed climo dir structure
# 2009-05-28 CAF  fixed problem with current working directory
# 2009-06-11 CAF  added date range check to input file
# 2010-12-14 HS/SWW changed input directory and input files & output
# 2011-02-01 ZL   re-introduced one of input filename for ccsm
# 2012-05-11 HS   added a call to test_timerange.ncl,
#                 check time range for T*M*
# 2014-11-27 BvU  added optional horizontal and vertical grid files

# Known issues:

# What it does:
#  -  Sets environment variables then calls a csh script to convert the file
#

from auxiliary import info
import exceptions
import os
import pdb
import projects


def infile(currProject, project_info, variable, model):
    model_name = currProject.get_model_name(model)
    project_name = currProject.get_project_name(model)
    project_basename = currProject.get_project_basename()
    project_info['RUNTIME']['model'] = model_name
    project_info['RUNTIME']['project'] = project_name
    project_info['RUNTIME']['project_basename'] = project_basename

    # Build input and output file names
    indir, infile = currProject.get_cf_infile(project_info,
                                              model,
                                              variable.fld,
                                              variable.var,
                                              variable.mip,
                                              variable.exp)
    return(os.path.join(indir, infile))


def cmor_reformat(currProject, project_info, variable, model):
    model_name = currProject.get_model_name(model)
    project_name = currProject.get_project_name(model)
    project_basename = currProject.get_project_basename()
    project_info['RUNTIME']['model'] = model_name
    project_info['RUNTIME']['project'] = project_name
    project_info['RUNTIME']['project_basename'] = project_basename
    verbosity = project_info["GLOBAL"]["verbosity"]
    exit_on_warning = project_info['GLOBAL'].get('exit_on_warning', False)

    # Variable put in environment to be used for the (optional)
    # wildcard syntax in the model path, ".../${VARIABLE}/..."
    # in the namelist
    os.environ['__ESMValTool_base_var'] = variable.var

    # Build input and output file names
    indir, infile = currProject.get_cf_infile(project_info,
                                              model,
                                              variable.fld,
                                              variable.var,
                                              variable.mip,
                                              variable.exp)

    fullpath = currProject.get_cf_fullpath(project_info,
                                           model,
                                           variable.fld,
                                           variable.var,
                                           variable.mip,
                                           variable.exp)
#    print "indir = %s" % indir
#    print "infile = %s" % infile
#    print "fullpath = %s" % fullpath

    if (not os.path.isdir(os.path.dirname(fullpath))):
        os.makedirs(os.path.dirname(fullpath))

    # Area file name for ocean grids
    areafile_path = currProject.get_cf_areafile(project_info, model)

    # Land-mask file name for land variables
    lmaskfile_path = currProject.get_cf_lmaskfile(project_info, model)
    omaskfile_path = currProject.get_cf_omaskfile(project_info, model)

    # Porosity file name for land variables
    porofile_path = currProject.get_cf_porofile(project_info, model)

    # Additional grid file names for ocean grids, if available (ECEARTH)
    hgridfile_path = False
    zgridfile_path = False
    lsmfile_path = False
    if hasattr(currProject, "get_cf_hgridfile"):
        hgridfile_path = currProject.get_cf_hgridfile(project_info, model)
    if hasattr(currProject, "get_cf_zgridfile"):
        zgridfile_path = currProject.get_cf_zgridfile(project_info, model)
    if hasattr(currProject, "get_cf_lsmfile"):
        lsmfile_path = \
            currProject.get_cf_lsmfile(project_info, model, variable.fld)

    # General fx file name entry
    fx_file_path = False
    if hasattr(currProject, "get_cf_fx_file"):
        fx_file_path = currProject.get_cf_fx_file(project_info, model)

    project, name, ensemble, start_year, end_year, dir\
        = currProject.get_cf_sections(model)
    info("project is " + project, verbosity, required_verbosity=4)
    info("ensemble is " + ensemble, verbosity, required_verbosity=4)
    info("dir is " + dir, verbosity, required_verbosity=4)

    # Check if the current project has a specific reformat routine,
    # otherwise use default
    if (os.path.isdir("reformat_scripts/" + project)):
        which_reformat = project
    else:
        which_reformat = 'default'

    reformat_script = os.path.join("reformat_scripts",
                                   which_reformat,
                                   "reformat_" + which_reformat + "_main.ncl")

    # Set enviroment variables
    project_info['TEMPORARY'] = {}
    project_info['TEMPORARY']['indir_path'] = indir
    project_info['TEMPORARY']['outfile_fullpath'] = fullpath
    project_info['TEMPORARY']['infile_path'] = os.path.join(indir, infile)
    project_info['TEMPORARY']['areafile_path'] = areafile_path
    project_info['TEMPORARY']['lmaskfile_path'] = lmaskfile_path
    project_info['TEMPORARY']['omaskfile_path'] = omaskfile_path
    project_info['TEMPORARY']['porofile_path'] = porofile_path
    project_info['TEMPORARY']['start_year'] = start_year
    project_info['TEMPORARY']['end_year'] = end_year
    project_info['TEMPORARY']['ensemble'] = ensemble
    project_info['TEMPORARY']['variable'] = variable.var
    project_info['TEMPORARY']['field'] = variable.fld

    # FX file path
    if fx_file_path:
        project_info['TEMPORARY']['fx_file_path'] = fx_file_path

    # Special cases
    if 'realm' in currProject.get_model_sections(model):
        project_info['TEMPORARY']['realm'] = \
            currProject.get_model_sections(model)["realm"]
    if 'shift_year' in currProject.get_model_sections(model):
        project_info['TEMPORARY']['shift_year'] = \
            currProject.get_model_sections(model)["shift_year"]
    if 'case_name' in currProject.get_model_sections(model):
        project_info['TEMPORARY']['case_name'] = \
            currProject.get_model_sections(model)["case_name"]

    if hgridfile_path and zgridfile_path:
        project_info['TEMPORARY']['hgridfile_path'] = hgridfile_path
        project_info['TEMPORARY']['zgridfile_path'] = zgridfile_path
    if lsmfile_path:
        project_info['TEMPORARY']['lsmfile_path'] = lsmfile_path

    # Execute the ncl reformat script
    if ((not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath']))
            or project_info['GLOBAL']['force_processing']):

        info("  Calling " + reformat_script + " to check/reformat model data",
             verbosity,
             required_verbosity=1)

        projects.run_executable(reformat_script, project_info, verbosity,
                                exit_on_warning)
    if 'NO_REFORMAT' in reformat_script:
        pass
    else:
        if (not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath'])):
            raise exceptions.IOError(2, "Expected reformatted file isn't available: ",
                                     project_info['TEMPORARY']['outfile_fullpath'])
    del(project_info['TEMPORARY'])
