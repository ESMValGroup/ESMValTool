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
from regrid import regrid as rg
import iris

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
                     str(model['start_year'])]) + "-" + str(model['end_year'])

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

def get_regridded_cf_fullpath(project_info, model, field, variable, regrid_target_name):
    """ @brief Returns the path (only) to the output file used in reformat AFTER regridding
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in

            This function specifies the full output path (directory + file) to
            the outupt file to use in the reformat routines and in climate.ncl
    """
    outdir = get_cf_outpath(project_info, model) + '/REGRIDDED_ON_' + regrid_target_name
    if not os.path.isdir(outdir):
        mkd = 'mkdir -p ' + outdir
        proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
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
                        str(model['start_year'])]) + '-' + str(model['end_year']) + '.nc'

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
                         str(model['start_year']),
                         str(model['end_year'])])
    return dict_key

def get_cmip_cf_infile(project_info, currentDiag, model, currentVarName):
    """@brief Path to input netCDF files for models
       @param project_info all info dictionary
       @param currDiag(dict) current diagnostic
       @param model(dict) what model we looking at
       This function supports multiple models and variables
       and operates on a single diagnostic
       Returns a dictionary {var1:[model_data], var2:[model_data], ...}
    """
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
                if currentVarName == var['name']:
                    info(" >>> preprocess.py >>> Using file " + fpath.strip(), verbosity, required_verbosity=1)
            else:
                fi = rootdir + '/' + infile_id
                info(" >>> preprocess.py >>> Could not find file type " + fi, verbosity, required_verbosity=1)
        data_files[var['name']] = infiles
        if len(infiles) == 0:
            info(" >>> preprocess.py >>> Could not find any data files for " + model['name'], verbosity, required_verbosity=1)
        if len(infiles) > 1:
            # get pp.glob to glob the files into one single file
            import preprocessing_tools as pt
            info(" >>> preprocess.py >>> Found multiple netCDF files for current diagnostic, attempting to glob them; variable: " + var['name'], 
                 verbosity, required_verbosity=1)
            standard_name = rootdir + '/' + '_'.join([var['name'], model['mip'],
                                                     model['exp'],
                                                     model['name'],
                                                     model['ensemble']]) + '_GLOB.nc'
            # glob only if the GLOB.nc file doesnt exist
            if os.path.exists(standard_name) is False:
                globs = pt.glob(infiles, standard_name, verbosity)
                info(" >>> preprocess.py >>> Globbing files now...",
                     verbosity, required_verbosity=1)
            else:
                info(" >>> preprocess.py >>> Found GLOB file: " + standard_name,
                     verbosity, required_verbosity=1)
                globs = 1
            if globs == 1:
                data_files[var['name']] = [standard_name]
            else:
                data_files[var['name']] = infiles
                info(" >>> preprocess.py >>> Could not glob files, keeping a list of files", verbosity, required_verbosity=1)
    # VP-FIXME-question it is easy to implement a 'read-from-file' (cache file) functionality here
    # if people consider this to be useful (e.g. could use cmip5datafinder as a lookup data module)
    return data_files

def get_obs_cf_infile(project_info, currentDiag, obs_model, currentVarName):
    """@brief Function that returns the observation file for regridding
       Returns a full path dictionary keyed on variable name
       namelist field type: {name: ERA-Interim,  project: OBS,  type: reanaly,  
       version: 1,  start: 2000,  end: 2002,  path: /obspath/Tier3/ERA-Interim/}
       path type: test_data/OBS/Tier3/ERA-Interim/OBS_ERA-Interim_reanaly_1_T3M_ta_200001-200212.nc
    """
    verbosity = project_info['GLOBAL']['verbosity']
    # get variables
    data_files = {}
    variables = currentDiag.variables
    for var in variables:
        infiles = []
        # specify path completely
        rootdir = obs_model['path']
        # use standard convention for CMIP file naming
        infile_id = '*' + '_'.join([obs_model['project'], obs_model['name'],
                                    obs_model['type'], str(obs_model['version']),
                                    var['field'], var['name']]) + '_*.nc*'
        # look for files: get paths
        srch = 'ls ' + rootdir + '/' + infile_id
        proc = subprocess.Popen(srch, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        fpaths = out.split('\n')
        # last element of ls will always be an empty string
        for fpath in fpaths[0:-1]:
            if os.path.exists(fpath.strip()):
                infiles.append(fpath.strip())
                if currentVarName == var['name']:
                    info(" >>> preprocess.py >>> Using OBS file for regridding " + fpath.strip(), verbosity, required_verbosity=1)
            else:
                fi = rootdir + '/' + infile_id
                info(" >>> preprocess.py >>> Could not find OBS file type " + fi, verbosity, required_verbosity=1)
        data_files[var['name']] = infiles
        if len(infiles) == 0:
            info(" >>> preprocess.py >>> Could not find any OBS data files for " + obs_model['name'], verbosity, required_verbosity=1)
        if len(infiles) > 1:
            # get pp.glob to glob the files into one single file
            import preprocessing_tools as pt
            info(" >>> preprocess.py >>> Found multiple OBS netCDF files for current diagnostic, attempting to glob them", verbosity, required_verbosity=1)
            standard_name = rootdir + '/' + '_'.join([obs_model['project'], obs_model['name'],
                                                     obs_model['type'], obs_model['version'],
                                                     var['name'], var['field']]) + '_GLOB.nc'
            globs = pt.glob(infiles, standard_name, verbosity)
            if globs == 1:
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

# a couple functions needed by cmor reformatting (the new python one)
def get_attr_from_field_coord(ncfield, coord_name, attr):
    if coord_name is not None:
        attrs = ncfield.cf_group[coord_name].cf_attrs()
        attr_val = [value for (key, value) in attrs if key == attr]
        if attr_val:
            return attr_val[0]
    return None

# Use this callback to fix anything Iris tries to break!
# noinspection PyUnusedLocal
def merge_callback(raw_cube, field, filename):
    # Remove attributes that cause issues with merging and concatenation
    for attr in ['creation_date', 'tracking_id', 'history']:
        if attr in raw_cube.attributes:
            del raw_cube.attributes[attr]
    for coord in raw_cube.coords():
        # Iris chooses to change longitude and latitude units to degrees
        #  regardless of value in file, so reinstating file value
        if coord.standard_name in ['longitude', 'latitude']:
            units = get_attr_from_field_coord(field,
                                              coord.var_name,
                                              'units')
            if units is not None:
                coord.units = units

####################################################################################################################################

def preprocess(project_info, variable, model, currentDiag, cmor_reformat_type):

###############################################################################################################
#################### The Big Mean PREPROCESS Machine ##########################################################
###############################################################################################################

    # Starting to look at preprocessing needs
    verbosity = project_info["GLOBAL"]["verbosity"]
    info(' >>> preprocess.py >>> Looking at (namelist) PREPROCESS section now', verbosity, required_verbosity=1)
    # initialize the environment variables dictionary
    project_info['TEMPORARY'] = {}
    # key in the prperocess
    prp = project_info['PREPROCESS']

    #############################################################################
    ### PRERQUISITES: GET THE PARAMETERS From config file
    #############################################################################
    # VP-FIXME-question : these can be adjusted by simply changing the yaml config
    for k in prp.keys():
        if k == 'select_level':
            if prp[k] is not 'None':
                select_level = prp[k]
        if k == 'regrid':
            regrid = prp[k]
        if k == 'target_grid':
            if prp[k] is not 'None':
                target_grid = prp[k]
        if k == 'regrid_scheme':
            if prp[k] is not 'None':
                regrid_scheme = prp[k]
        if k == 'mask_fillvalues':
            if prp[k] is not False:
                mask_fillvalues = True
        if k == 'multimodel_mean':
            if prp[k] is not False:
                regrid_scheme = True
        if k == 'gridfile':
            if prp[k] is not 'None':
                areafile_path = prp[k]
                if areafile_path is not None:
                    project_info['TEMPORARY']['areafile_path'] = areafile_path
        if k == 'mask_landocean':
            # names and paths are taken from old projects.py
            # VP-FIXME-question : change?
            if prp[k] is not 'None': 
                lmaskdir = os.path.join(model["path"],
                                       model["exp"],
                                       'fx',
                                       'sftlf',
                                       model["name"],
                                       'r0i0p0')
                lmaskfile = 'sftlf_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                lmaskfile_path = os.path.join(lmaskdir, lmaskfile)
                if lmaskfile_path is not None:
                    project_info['TEMPORARY']['lmaskfile_path'] = lmaskfile_path
                omaskdir = model['path']
                omaskfile = 'sftof_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                omaskfile_path = os.path.join(omaskdir, omaskfile)
                if omaskfile_path is not None:
                    project_info['TEMPORARY']['omaskfile_path'] = omaskfile_path
        if k == 'porofile':
            if prp[k] is not 'None':
                pormaskdir = os.path.join(model["path"],
                                          model["exp"],
                                          'fx',
                                          'mrsofc',
                                          model["name"],
                                          'r0i0p0')
                pormaskfile = 'mrsofc_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                porofile_path = os.path.join(maskdir, maskfile)
                if porofile_path is not None:
                    project_info['TEMPORARY']['porofile_path'] = porofile_path

    # VP-FIXME-question : do we need these still?
    # if so, easily initialized in namelist file
    hgridfile_path = False
    zgridfile_path = False
    lsmfile_path = False
    fx_file_path = False
    indata_root = 'probably/set/this/to/a/common/fx/stuff/directory'
    for k in prp.keys():
        if k == 'get_cf_hgridfile':
            if prp[k] is not 'None':
                filepath = fx_files['nemo_hgrid_file'].get_fullpath()
                hgridfile_path = os.path.join(indata_root, filepath)
                if hgridfile_path is not None:
                    project_info['TEMPORARY']['hgridfile_path'] = hgridfile_path
        if k == 'get_cf_zgridfile':
            if prp[k] is not 'None':
                filepath = fx_files['nemo_zgrid_file'].get_fullpath()
                zgridfile_path = os.path.join(indata_root, filepath)
                if zgridfile_path is not None:
                    project_info['TEMPORARY']['zgridfile_path'] = zgridfile_path
        if k == 'get_cf_lsmfile':
            if prp[k] is not 'None':
                zgrid_path = fx_files['nemo_zgrid_file'].get_fullpath()
                if variable.field == "T3M":
                    filepath = fx_files['nemo_lsm3d_file'].get_fullpath()
                else:
                    filepath = fx_files['nemo_lsm_file'].get_fullpath()
                lsmfile_path = os.path.join(indata_root, filepath)
                if lsmfile_path is not None:
                    project_info['TEMPORARY']['lsmfile_path'] = lsmfile_path
        if k == 'get_fx_files':
            if prp[k] is not 'None':
                curr_fx = project_info["AUXILIARIES"]["FX_files"].fx_files
                new_fx = {}
                for key in curr_fx:
                    new_fx[key] = os.path.join(indata_root, curr_fx[key].get_fullpath())
                fx_file_path = os.path.abspath(new_fx)
                if fx_file_path:
                    project_info['TEMPORARY']['fx_file_path'] = fx_file_path

    # Specialy McSpecial cases
    if 'realm' in model.keys():
        project_info['TEMPORARY']['realm'] = model["realm"]
    if 'shift_year' in model.keys():
        project_info['TEMPORARY']['shift_year'] = model["shift_year"]
    if 'case_name' in model.keys():
        project_info['TEMPORARY']['case_name'] = model["case_name"]


    ######################################################################
    ### ENVIRONMENT VARIABLES
    # Initialize all needed files and variables (Leon: EEEvvryboooodyyy!!!) 
    #######################################################################
    model_name = model['name']
    project_name = model['project']
    project_info['RUNTIME']['model'] = model_name
    project_info['RUNTIME']['project'] = project_name
    project_info['RUNTIME']['project_basename'] = 'OBS'
    exit_on_warning = project_info['GLOBAL'].get('exit_on_warning', False)

    # Variable put in environment
    os.environ['__ESMValTool_base_var'] = variable.name

    # Build input and output file names
    infileslist = get_cmip_cf_infile(project_info, currentDiag, model, variable.name)[variable.name]
    # VP-FIXME-question : the code can glob multiple files, but how to handle the case when globbing fails; currently
    # the diagnostic is run on only the first file
    if len(infileslist) == 1:
        infiles = infileslist[0]
    else:
        info(" >>> preprocess.py >>> Found multiple netCDF files for current diagnostic", verbosity, required_verbosity=1)
        info(" >>> preprocess.py >>> netCDF globbing has failed", verbosity, required_verbosity=1)
        info(" >>> preprocess.py >>> Running diagnostic ONLY on the first file", verbosity, required_verbosity=1)
        infiles = infileslist[0]
    outfilename = get_cf_outfile(model, variable.field, variable.name)
    info(' >>> preprocess.py >>> Reformatted file name: ' + outfilename, verbosity, required_verbosity=1)
    # get full outpaths - original cmorized files that are preserved all through the process 
    fullpath = get_cf_fullpath(project_info, model, variable.field, variable.name)
    info(' >>> preprocess.py >>> Reformatted target: ' + fullpath, verbosity, required_verbosity=1)

    # indir is hardcoded to keep things tight; could be an option to namelist
    project_info['TEMPORARY']['indir_path'] = '.'
    project_info['TEMPORARY']['outfile_fullpath'] = fullpath
    project_info['TEMPORARY']['infile_path'] = infiles
    project_info['TEMPORARY']['start_year'] = model['start_year']
    project_info['TEMPORARY']['end_year'] = model['end_year']
    project_info['TEMPORARY']['ensemble'] = model['ensemble']
    project_info['TEMPORARY']['variable'] = variable.name
    project_info['TEMPORARY']['field'] = variable.field
    info(' >>> preprocess.py >>> Gathering runtime variables:', verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Project is " + model['project'], verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Model is " + model['name'], verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Ensemble is " + model['ensemble'], verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Full IN path is " + infiles, verbosity, required_verbosity=1)
    info(" >>> preprocess.py >>> Full OUT path is " + fullpath, verbosity, required_verbosity=1)

    ################## START CHANGING STUFF ###########################################

    ################## 0. CMOR_REFORMAT (NCL version) ##############
    # Legacy code that will be purged in the future
    # Check if the current project has a specific reformat routine,
    # otherwise use default
    # the cmor reformat is applied only once per variable
    if cmor_reformat_type == 'ncl':
        if (os.path.isdir("reformat_scripts/" + model['project'])):
            which_reformat = model['project']
        else:
            which_reformat = 'default'
    
        reformat_script = os.path.join("reformat_scripts",
                                       which_reformat,
                                       "reformat_" + which_reformat + "_main.ncl")
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

    ################## 0. CMOR_REFORMAT (NCL version) ##############
    # New code: cmor_check.py (by Javier Vegas)
    if cmor_reformat_type == 'py':
        # needed imports
        from cmor_check import CMORCheck as CC
        from cmor_check import CMORCheckError as CCE
        import warnings
        from variable_info import CMIP5Info

        variables_info = CMIP5Info()

        var_name = variable.name
        table = model['mip'] 


        try:
            # Load cubes for requested variable in given files
            # remember naming conbentions 
            # IN: infiles
            # OUT: project_info['TEMPORARY']['outfile_fullpath']
            files = infiles
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF measure variable',
                                        UserWarning)
                warnings.filterwarnings('ignore',
                                        'Missing CF-netCDF boundary variable',
                                        UserWarning)
                def cube_var_name(raw_cube):
                    return raw_cube.var_name == var_name
                var_cons = iris.Constraint(cube_func=cube_var_name)
                # force single cube; this function defaults a list of cubes
                reft_cube = iris.load(files, var_cons, callback=merge_callback)[0]

            # Concatenate data to single cube, i.e. merge time series
            # cube = cubes.concatenate_cube()
            # Create checker for loaded cube
            var_info = variables_info.get_variable(table, var_name)
            checker = CC(reft_cube, var_info, automatic_fixes=True)
            # Run checks
            checker.check_metadata()
            checker.check_data()
            # save reformatted cube
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

        except (iris.exceptions.ConstraintMismatchError,
                iris.exceptions.ConcatenateError,
                CCE) as ex:
            print(ex)

    #################### 1. MASK AND TIME-AREA OPERATIONA   ####################################################

    #################### 2. REGRID ####################################################
    if regrid is True and target_grid is not None:
        # we will regrid according to whatever regridding scheme and reference grids are needed
        # and create new regridded files from the original cmorized/masked ones
        # but preserving thos files (optionally)
        info("  >>> preprocess.py >>>  Calling regrid to regrid model data onto " + target_grid + " grid",
                 verbosity,
                 required_verbosity=1)
        # first let's see if the regrid source is a simple netCDF file (cube)
        if os.path.isfile(project_info['TEMPORARY']['outfile_fullpath']):
            info(' >>> preprocess.py >>> Preparing to regrid ' + project_info['TEMPORARY']['outfile_fullpath'], verbosity, required_verbosity=1)
            src_cube = iris.load_cube(project_info['TEMPORARY']['outfile_fullpath'])
        else:
            src_cube = infiles
        # try return a cube for regrid target
        # target_grid - could simply be a netCDF file or string model
        # descriptor eg 'ref_model'; currently netCDF and ref_model labels are implemented
        try:
            tgt_grid_cube = iris.load_cube(target_grid)
            if regrid_scheme:
                rgc = rg(src_cube, tgt_regrid_cube, regrid_scheme)
            else:
                info(' >>> preprocess.py >>> No regrid scheme specified, assuming linear', verbosity, required_verbosity=1)
                rgc = rg(src_cube, tgt_regrid_cube, 'linear')
            # save-append to outfile fullpath list to be further processed
            iris.save(rgc, project_info['TEMPORARY']['outfile_fullpath'])
        except (IOError, iris.exceptions.IrisError) as exc:
            info(" >>> preprocessing_tools.py >>> Target " + target_grid + " is not a netCDF file", "", verbosity)
            pass
            # continue to see what regridding we need
            # ref_model regrid string descriptor
            project_info['RUNTIME']['regridtarget'] = []
            if target_grid == 'ref_model':
                additional_models_dicts = currentDiag.additional_models
                # identify the current variable
                for var in currentDiag.variables:
                    if var['name'] == variable.name:
                        ref_model_list = var['ref_model']
                        for ref_model in ref_model_list:
                            for obs_model in additional_models_dicts:
                                if obs_model['name'] == ref_model:
                                    # add to environment variable
                                    project_info['RUNTIME']['regridtarget'].append(ref_model)
                                    info(' >>> preprocess.py >>> Regridding on ref_model ' + ref_model, verbosity, required_verbosity=1)
                                    tgt_nc_grid = get_obs_cf_infile(project_info, currentDiag, obs_model, variable.name)[variable.name][0]
                                    tgt_grid_cube = iris.load_cube(tgt_nc_grid)
                                    if regrid_scheme:
                                        rgc = rg(src_cube, tgt_grid_cube, regrid_scheme)
                                    else:
                                        info(' >>> preprocess.py >>> No regrid scheme specified, assuming linear', verbosity, required_verbosity=1)
                                        rgc = rg(src_cube, tgt_grid_cube, 'linear')
                                    # save specifically named regridded file
                                    newlyRegriddedCube = iris.save(rgc, get_regridded_cf_fullpath(project_info, model, variable.field, variable.name, ref_model))
                                    newlyRegriddedFilePath = get_regridded_cf_fullpath(project_info, model, variable.field, variable.name, ref_model)
                                    newlyRegriddedCube = iris.save(rgc, newlyRegriddedFilePath)
                                    info(' >>> preprocess.py >>> Running diagnostic on ' + newlyRegriddedFilePath, verbosity, required_verbosity=1)

    ############ FINISH all PREPROCESSING and delete environment
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
            infiles = get_cmip_cf_infile(project_info, currentDiag, model, base_var.name)[base_var.name]

            if len(infiles) == 0:
                info("  >>> preprocess.py >>> No input files found for " + base_var.name +
                     " (" + base_var.field + ")", verbosity, 1)

                base_var.var = base_var.var0
                base_var.fld = base_var.fld0

                # try again with input variable = base variable (non derived)
                infile = get_cmip_cf_infile(project_info, currentDiag, model, base_var.name)[base_var.name]

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
