"""
New module to enable ESMValToolKit to deal with
the new yaml parser and simplified interface_scripts 
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
import os
import pdb
import subprocess
import re
import sys
import data_interface as dint
import auxiliary
from auxiliary import info, error, print_header, ncl_version_check
import exceptions
import launchers

#######################################################
### This script contains basic functionalities
### It contains all variables to call other
### more specialized modules for prerocessing purposes
########################################################

def write_data_interface(executable, project_info):
    """ @brief Write Python data structures to target script format interface
        @param executable String pointing to the script/binary to execute
        @param project_info Current namelist in dictionary format

        Data structures in Python are rewritten to the interface folder in
        a format appropriate for the target script/binary
    """
    suffix = os.path.splitext(executable)[1][1:]
    currInterface = vars(dint)[suffix.title() + '_data_interface'](project_info)
    currInterface.write_data_to_interface()

def run_executable(string_to_execute,
                   project_info,
                   verbosity,
                   exit_on_warning,
                   launcher_arguments=None,write_di=True):
    """ @brief Executes script/binary
        @param executable String pointing to the script/binary to execute
        @param project_info Current namelist in dictionary format
        @param verbosity The requested verbosity level
        @param exit_on_warning Boolean defining whether the wrapper should
                               crash on warnings

        Check the type of script/binary from the executable string suffix and
        execute the script/binary properly.
    """

    if write_di:
        write_data_interface(string_to_execute, project_info)

    suffix = os.path.splitext(string_to_execute)[1][1:]
    currLauncher = vars(launchers)[suffix + '_launcher']()
    if launcher_arguments is not None:
        currLauncher.arguments = launcher_arguments
    currLauncher.execute(string_to_execute,
                         project_info,
                         verbosity,
                         exit_on_warning)

def get_figure_file_names(project_info, model):
    """ @brief Returns names for plots
        @param project_info Current namelist in dictionary format
        @param some model from namelist
    """
    return "_".join([model['project'], model['name'], model['mip'], model['exp'], model['ensemble'],
                     str(model['start'])]) + "-" + str(model['end'])

def get_cf_fullpath(project_info, model, field, variable):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in

            This function specifies the full output path (directory + file) to
            the outupt file to use in the reformat routines and in climate.ncl
    """
    outdir = get_cf_outpath(project_info, model)
    outfile = get_cf_outfile(model, field, variable)
    return os.path.join(outdir, outfile)

def get_cf_outpath(project_info, model):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (output path)

            Standard path: dir_output/climo_dir/projectname/projectname_expname_ens_field_var_yrstart-yrend.nc
    """
    outdir1 = project_info['GLOBAL']['climo_dir']
    outdir2 = model['project']
    # let's check if directories exists, if not create them
    if not os.path.isdir(outdir1):
        mkd = 'mkdir -p ' + project_info['GLOBAL']['climo_dir']
        proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
    if not os.path.isdir(os.path.join(outdir1, outdir2)):
        mkd = 'mkdir -p ' + os.path.join(outdir1, outdir2)
        proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
    return os.path.join(outdir1, outdir2)

def get_cf_outfile(model, field, variable):
    """ @brief Returns the path and output file used in reformat
            @param variable Current variable
            @param climo_dir Where processed (reformatted) input files reside
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in
                            project_info)
            @return A string (output path)

            This function specifies the output file to use in reformat
            and in climate.ncl
    """
    outfile = '_'.join([model['project'],
                        model['name'],
                        model['mip'],
                        model['exp'],
                        model['ensemble'],
                        field,
                        variable,
                        str(model['start'])]) + '-' + str(model['end']) + '.nc'

    return outfile

def get_dict_key(model):
    """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the yaml namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
    """

    dict_key = "_".join([model['project'],
                         model['name'],
                         model['mip'],
                         model['exp'],
                         model['ensemble'],
                         str(model['start']),
                         str(model['end'])])
    return dict_key

def get_cmip_cf_infile(project_info, currentDiag, model):
    """@brief Path to input netCDF files for models
       @param project_info all info dictionary
       @param currDiag(dict) current diagnostic
       @param model(dict) what model we looking at
       This function supports multiple models and variables
       and operates on a single diagnostic
       Returns a dictionary {var1:[model_data], var2:[model_data], ...}
    """
    # FIXME introduce functionality to glob multiple nc files !!!
    verbosity = project_info['GLOBAL']['verbosity']
    # get variables
    data_files = {}
    variables = currentDiag.variables
    for var in variables:
        infiles = []
        # specify path completely
        rootdir = model['path']
        # use standard convention for CMIP file naming
        infile_id = '*' + '_'.join([var['name'], model['mip'],
                                    model['name'],
                                    model['exp'],
                                    model['ensemble']]) + '_*.nc*'
        # look for files: get paths
        srch = 'ls ' + rootdir + '/' + infile_id
        proc = subprocess.Popen(srch, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        fpaths = out.split('\n')
        # last element of ls will always be an empty string
        for fpath in fpaths[0:-1]:
            if os.path.exists(fpath.strip()):
                infiles.append(fpath.strip())
                info(" >>> preprocess.py >>> Using file " + fpath.strip(), verbosity, required_verbosity=1)
            else:
                fi = rootdir + '/' + infile_id
                info(" >>> preprocess.py >>> Could not find file type " + fi, verbosity, required_verbosity=1)
        data_files[var['name']] = infiles
        if len(infiles) > 1:
            # get pp.glob to glob the files into one single file
            import preprocessing_tools as pt
            info(" >>> preprocess.py >>> Found multiple netCDF files for current diagnostic, attempting to glob them", verbosity, required_verbosity=1)
            standard_name = rootdir + '/' + '_'.join([var['name'], model['mip'],
                                                     model['exp'],
                                                     model['name'],
                                                     model['ensemble']]) + '_GLOB.nc'
            globs = pt.glob(infiles, standard_name)
            if pt.glob(data_files[var['name']], standard_name) == 1:
                data_files[var['name']] = [standard_name]
            else:
                data_files[var['name']] = infiles
                info(" >>> preprocess.py >>> Could not glob files, keeping a list of files", verbosity, required_verbosity=1)
    return data_files

def get_cf_areafile(project_info, model):
    """ @brief Returns the path to the areacello file
        This function looks for the areafile of the ocean grid
    """
    areadir = model["path"]
    areafile = 'areacello_fx_' + model["name"] + "_" + model["exp"] + \
               "_r0i0p0.nc"

    return os.path.join(areadir, areafile)

def cmor_reformat(project_info, variable, model, currentDiag):
    model_name = model['name']
    project_name = model['project']
    #project_basename = currProject.get_project_basename()
    project_info['RUNTIME']['model'] = model_name
    project_info['RUNTIME']['project'] = project_name
    project_info['RUNTIME']['project_basename'] = 'OBS'
    verbosity = project_info["GLOBAL"]["verbosity"]
    exit_on_warning = project_info['GLOBAL'].get('exit_on_warning', False)

    # Variable put in environment to be used for the (optional)
    # wildcard syntax in the model path, ".../${VARIABLE}/..."
    # in the namelist
    os.environ['__ESMValTool_base_var'] = variable.name

    # Build input and output file names
    # FIXME - when globbing fails, we are looking only at a single (first) file currently
    infiles = get_cmip_cf_infile(project_info, currentDiag, model)[variable.name][0]
    outfilename = get_cf_outfile(model, variable.field, variable.name)
    info(' >>> preprocess.py >>> Reformatted file name: ' + outfilename, verbosity, required_verbosity=1)
    fullpath = get_cf_fullpath(project_info, model, variable.field, variable.name)
    info(' >>> preprocess.py >>> Reformatted target: ' + fullpath, verbosity, required_verbosity=1)

    # Area file name for ocean grids
    areafile_path = get_cf_areafile(project_info, model)

    # Land-mask file name for land variables
    # start FIXME get these funcs from projects
    #lmaskfile_path = get_cf_lmaskfile(project_info, model)
    #omaskfile_path = get_cf_omaskfile(project_info, model)
    # Porosity file name for land variables
    #porofile_path = get_cf_porofile(project_info, model)
    # end FIXME

    # Additional grid file names for ocean grids, if available (ECEARTH)
    # start FIXME these should be part of namelist
    #hgridfile_path = False
    #zgridfile_path = False
    #lsmfile_path = False
    #if hasattr(currProject, "get_cf_hgridfile"):
    #    hgridfile_path = currProject.get_cf_hgridfile(project_info, model)
    #if hasattr(currProject, "get_cf_zgridfile"):
    #    zgridfile_path = currProject.get_cf_zgridfile(project_info, model)
    #if hasattr(currProject, "get_cf_lsmfile"):
    #    lsmfile_path = \
    #        currProject.get_cf_lsmfile(project_info, model, variable.fld)
    # General fx file name entry
    #fx_file_path = False
    #if hasattr(currProject, "get_cf_fx_file"):
    #    fx_file_path = currProject.get_cf_fx_file(project_info, model)
    # end FIXME

    info(" >>> preprocess.py >>> Project is " + model['project'], verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Ensemble is " + model['ensemble'], verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Full IN path is " + infiles, verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Full OUT path is " + fullpath, verbosity, required_verbosity=1)

    # Check if the current project has a specific reformat routine,
    # otherwise use default
    if (os.path.isdir("reformat_scripts/" + model['project'])):
        which_reformat = model['project']
    else:
        which_reformat = 'default'

    reformat_script = os.path.join("reformat_scripts",
                                   which_reformat,
                                   "reformat_" + which_reformat + "_main.ncl")

    # Set enviroment variables
    project_info['TEMPORARY'] = {}
    # hardcoded to keep things tight; could be an option to namelist
    project_info['TEMPORARY']['indir_path'] = '.'
    project_info['TEMPORARY']['outfile_fullpath'] = fullpath
    project_info['TEMPORARY']['infile_path'] = infiles
    #project_info['TEMPORARY']['areafile_path'] = areafile_path
    #project_info['TEMPORARY']['lmaskfile_path'] = lmaskfile_path
    #project_info['TEMPORARY']['omaskfile_path'] = omaskfile_path
    #project_info['TEMPORARY']['porofile_path'] = porofile_path
    project_info['TEMPORARY']['start_year'] = model['start']
    project_info['TEMPORARY']['end_year'] = model['end']
    project_info['TEMPORARY']['ensemble'] = model['ensemble']
    project_info['TEMPORARY']['variable'] = variable.name
    project_info['TEMPORARY']['field'] = variable.field

    # FX file path
    # start FIXME add these in namelist
    #if fx_file_path:
    #    project_info['TEMPORARY']['fx_file_path'] = fx_file_path
    # end FIXME

    # Special cases
    # start FIXME - address these with correct variables
    #if 'realm' in currProject.get_model_sections(model):
    #    project_info['TEMPORARY']['realm'] = \
    #        currProject.get_model_sections(model)["realm"]
    #if 'shift_year' in currProject.get_model_sections(model):
    #    project_info['TEMPORARY']['shift_year'] = \
    #        currProject.get_model_sections(model)["shift_year"]
    #if 'case_name' in currProject.get_model_sections(model):
    #    project_info['TEMPORARY']['case_name'] = \
    #        currProject.get_model_sections(model)["case_name"]
    #if hgridfile_path and zgridfile_path:
    #    project_info['TEMPORARY']['hgridfile_path'] = hgridfile_path
    #    project_info['TEMPORARY']['zgridfile_path'] = zgridfile_path
    #if lsmfile_path:
    #    project_info['TEMPORARY']['lsmfile_path'] = lsmfile_path
    # end FIXME

    # Execute the ncl reformat script
    if ((not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath']))
            or project_info['GLOBAL']['force_processing']):

        info("  >>> preprocess.py >>>  Calling " + reformat_script + " to check/reformat model data",
             verbosity,
             required_verbosity=1)

        run_executable(reformat_script, project_info, verbosity,
                       exit_on_warning)
    if 'NO_REFORMAT' in reformat_script:
        pass
    else:
        if (not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath'])):
            raise exceptions.IOError(2, "Expected reformatted file isn't available: ",
                                     project_info['TEMPORARY']['outfile_fullpath'])
    del(project_info['TEMPORARY'])


###################################################################################
### a few definititory functions needed
###################################################################################
# variable class
class Var:
    """
    changed from original class diagdef.Var
    """
    def __init__(self, merged_dict, var0, fld0):
        # Write the variable attributes in 'merged_dict'
        # to class object attributes
        for name, value in merged_dict.iteritems():
            setattr(self, name, value)

        # Special cases, actually not sure what they do
        if (var0 == "none"):
            self.var0 = merged_dict['name']
        else:
            self.var0 = var0
        if (fld0 == "none"):
            self.fld0 = merged_dict['field']
        else:
            self.fld0 = fld0

# diagnostics class
class Diag:
    """
    changed from original diagdef.Diag
    """
    def add_base_vars_fields(self, variables, model, variable_def_dir):
        dep_vars = []
        model_der = False
        if hasattr(model, "attributes"):
            if "skip_derive_var" in model.attributes.keys():
                model_der = model.attributes["skip_derive_var"]

        for variable in variables:

            f = open(os.path.join(variable_def_dir, variable['name']
                                  + ".ncl"), 'r')
            for line in f:
                tokens = line.split()

                if "Requires:" in tokens:
                    # If 'none', return orig. field
                    if (tokens[2] == "none" or model_der == "True"):
                        dep_vars.append(Var(variable,
                                            "none",
                                            "none"))
                    else:
                        sub_tokens = tokens[2].split(",")
                        for sub in sub_tokens:
                            element = sub.split(":")

                            # Assume first digit is dimension in field type
                            offset = self.find_dimension_entry(element[1]) - 1

                            e_var = element[0]
                            e_fld = element[1]
                            e_fld = (variable['field'][0:offset + 1]
                                     + e_fld[1 + offset]
                                     + variable['field'][2 + offset]
                                     + e_fld[3 + offset:])

                            del keys
                            del vars
                            dep_var = copy.deepcopy(variable)
                            dep_var['name'] = e_var
                            dep_var['field'] = e_fld

                            dep_vars.append(Var(dep_var,
                                                variable.var,
                                                variable.fld))

        return dep_vars

    def select_base_vars(self, variables, model, currentDiag, project_info):

        verbosity = project_info['GLOBAL']['verbosity']

        for base_var in variables:
            # FIXME we should not need this anymore??
            # Check if this variable should be excluded for this model-id
            if self.id_is_explicitly_excluded(base_var, model):
                continue

            # first try: use base variables provided by variable_defs script
            os.environ['__ESMValTool_base_var'] = base_var.name
            infiles = get_cmip_cf_infile(project_info, currentDiag, model)[base_var.name]

            if len(infiles) == 0:
                info("  >>> preprocess.py >>> No input files found for " + base_var.name +
                     " (" + base_var.field + ")", verbosity, 1)

                base_var.var = base_var.var0
                base_var.fld = base_var.fld0

                # try again with input variable = base variable (non derived)
                infile = get_cmip_cf_infile(project_info, currentDiag, model)[base_var.name]

                if len(infile) == 0:
                    raise exceptions.IOError(2, "No input files found in ",
                                             infile)
                else:
                    info("  >>> preprocess.py >>> Using " + base_var.name + " (" + base_var.field +
                         ")", verbosity, 1)
                    base_vars = [base_var]
                    break  # discard other base vars

            else:
                base_vars = variables
        return base_vars

    def id_is_explicitly_excluded(self, var, model):
        """ Checks if variable attributes expclitly excludes
            use of this model
            Changed from original diagdef.Diag
        """
        exclude = False
        if hasattr(var, "exclude") and "id" in model.attributes:
            if model.attributes["id"] == var.exclude:
                    exclude = True

        if var.only != "None":
            exclude = True
            if "id" in model.attributes:
                if model.attributes["id"] == var.only:
                    exclude = False

        return exclude
