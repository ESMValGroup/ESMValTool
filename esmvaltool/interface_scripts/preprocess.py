"""
New module to enable ESMValToolKit to deal with
the new yaml parser and simplified interface_scripts
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
import copy
import exceptions
import logging
import os
import subprocess

import iris
import iris.exceptions
import numpy as np

import data_interface as dint
import launchers
import preprocessing_tools as pt
from data_finder import get_input_filelist, get_output_file
from fixes.fix import Fix
from regrid import regrid, vertical_schemes, vinterp

logger = logging.getLogger(__name__)

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
    curr_interface = vars(dint)[suffix.title() + '_data_interface'](project_info)
    curr_interface.write_data_to_interface()


def run_executable(string_to_execute,
                   project_info,
                   verbosity,
                   exit_on_warning,
                   launcher_arguments=None, write_di=True):
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
    interface_data = project_info['RUNTIME']['interface_data']
    curr_launcher = vars(launchers)[suffix + '_launcher'](interface_data=interface_data)
    if launcher_arguments is not None:
        curr_launcher.arguments = launcher_arguments
    curr_launcher.execute(string_to_execute,
                          project_info,
                          verbosity,
                          exit_on_warning)

def get_figure_file_names(project_info, model):
    """ @brief Returns names for plots
        @param project_info Current namelist in dictionary format
        @param some model from namelist
    """
    #return "_".join([model['project'], model['name'], model['mip'], model['exp'], model['ensemble'],
    #                 str(model['start_year'])]) + "-" + str(model['end_year'])
    return "_".join([model['project'], model['name'], str(model['start_year'])]) + "-" + str(model['end_year'])

def get_cf_fullpath(project_info, model, variable):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @param field The field (see tutorial.pdf for available fields)
            @param variable The variable (defaults to the variable in

            This function specifies the full output path (directory + file) to
            the outupt file to use in the reformat routines and in climate.ncl
    """
    fullpath = get_output_file(project_info, model, variable)
    return fullpath


def get_cf_outpath(project_info, model):
    """ @brief Returns the path (only) to the output file used in reformat
            @param project_info Current namelist in dictionary format
            @param model One of the <model>-tags in the XML namelist file
            @return A string (output path)

            Standard path: dir_output/preproc_dir/projectname/projectname_expname_ens_field_var_yrstart-yrend.nc
    """
    outdir1 = project_info['GLOBAL']['preproc_dir']
    outdir2 = model['project']
    # let's check if directories exist, if not create them
    if not os.path.isdir(outdir1):
        mkd = 'mkdir -p ' + project_info['GLOBAL']['preproc_dir']
        proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
    if not os.path.isdir(os.path.join(outdir1, outdir2)):
        mkd = 'mkdir -p ' + os.path.join(outdir1, outdir2)
        proc = subprocess.Popen(mkd, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
    return os.path.join(outdir1, outdir2)

def get_dict_key(model):
    """ @brief Returns a unique key based on the model entries provided.
            @param model One of the <model>-tags in the yaml namelist file
            @return A string

            This function creates and returns a key used in a number of NCL
            scripts to refer to a specific dataset, see e.g., the variable
            'cn' in 'interface_scripts/read_data.ncl'
    """

    # allow for different projects

    # CMIP5
    if model['project'] == 'CMIP5':
        dict_key = "_".join([model['project'],
                             model['name'],
                             model['mip'],
                             model['exp'],
                             model['ensemble'],
                             str(model['start_year']),
                             str(model['end_year'])])
    else:
        dict_key = "_".join([model['project'],
                             model['name'],
                             str(model['start_year']),
                             str(model['end_year'])])

    return dict_key

def get_cf_infile(project_info, current_diag, model, current_var_dict):
    """@brief Path to input netCDF files for models
       @param project_info all info dictionary
       @param currDiag(dict) current diagnostic
       @param model(dict) what model we looking at
       This function supports multiple models and variables
       and operates on a single diagnostic
       Returns a dictionary {var1:[model_data], var2:[model_data], ...}
    """
    # assert verbosity level
    verbosity = project_info['GLOBAL']['verbosity']

    # look for files keyed by project name
    proj_name = model['project']

    # get variables; data files is a dictionary keyed on variable
    data_files = {}
    variables = current_diag.variables

    for var in variables:

        full_paths = get_input_filelist(project_info, model, var)
        standard_path = get_cf_outpath(project_info, model) # Need this call to create subdir in preproc (FIX-ME: get_cf_outpath and get_cf_fullpath shall be merged in a single function get_cf_outfile)
        standard_name = get_cf_fullpath(project_info, model, var)
        standard_name = standard_name.replace('.nc', '_GLOB.nc')

        if len(full_paths) == 0:
            logger.info("Could not find any data files for %s", model['name'])

        if len(full_paths) == 1:
            data_files[var['name']] = full_paths

        if len(full_paths) > 1:

            # get pp.glob to glob the files into one single file
            logger.info("Found multiple netCDF files for current diagnostic, attempting to glob them; variable: %s", var['name'])

            # glob only if the GLOB.nc file doesnt exist
            if os.path.exists(standard_name) is False:
                globs = pt.glob(full_paths, standard_name, var['name'], verbosity)
                logger.info("Globbing files now...")
            else:
                logger.info("Found GLOB file: %s", standard_name)
                globs = 1
            if globs == 1:
                data_files[var['name']] = [standard_name]
            else:
                data_files[var['name']] = full_paths
                logger.info("Could not glob files, keeping a list of files")

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

def preprocess(project_info, variable, model, current_diag, cmor_reformat_type):

###############################################################################################################
#################### The Big Mean PREPROCESS Machine ##########################################################
###############################################################################################################

    # Starting to look at preprocessing needs
    verbosity = project_info["GLOBAL"]["verbosity"]
    logger.info("Looking at (namelist) PREPROCESS section now")
    # initialize the environment variables dictionary
    project_info['TEMPORARY'] = {}
    # key in the prperocess
    prp = project_info['PREPROCESS']

    #############################################################################
    ### PRERQUISITES: GET THE PARAMETERS From config file
    #############################################################################

    # initialize variables
    save_intermediary_cubes = False
    mask_fillvalues = False
    multimodel_mean = False
    mask_land = False
    mask_ocean = False
    mask_poro = False

    # parse dictionary
    for k in prp.keys():
        if k == 'select_level':
            select_level = prp[k]
        if k == 'target_grid':
            target_grid = prp[k]
        if k == 'regrid_scheme':
            regrid_scheme = prp[k]
        if k == 'mask_fillvalues':
            mask_fillvalues = prp[k]
        if k == 'multimodel_mean':
            multimodel_mean = prp[k]
        if k == 'gridfile':

                areafile_path = prp[k]
                if areafile_path is not None:
                    project_info['TEMPORARY']['areafile_path'] = areafile_path


        # land (keeps only land regions)
        if k == 'mask_landocean':

            if prp[k] == 'land':
                mask_land = True
                lmaskdir = os.path.join(model["path"],
                                       model["exp"],
                                       'fx',
                                       'sftlf',
                                       model["name"],
                                       'r0i0p0')
                lmaskfile = 'sftlf_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                lmaskfile_path = os.path.join(lmaskdir, lmaskfile)

        # ocean (keeps only ocean regions)

            if prp[k] == 'ocean':
                mask_ocean = True
                if lmaskfile_path is not None:
                    project_info['TEMPORARY']['lmaskfile_path'] = lmaskfile_path
                omaskdir = model['path']
                omaskfile = 'sftof_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                omaskfile_path = os.path.join(omaskdir, omaskfile)
                if omaskfile_path is not None:
                    project_info['TEMPORARY']['omaskfile_path'] = omaskfile_path

        # poro (keeps only poro regions)
        if k == 'mask_poro':

            if prp[k] is not False:
                mask_poro = True
                pormaskdir = os.path.join(model["path"],
                                          model["exp"],
                                          'fx',
                                          'mrsofc',
                                          model["name"],
                                          'r0i0p0')
                pormaskfile = 'mrsofc_fx_' + model["name"] + "_" + model["exp"] + "_r0i0p0.nc"
                porofile_path = os.path.join(pormaskdir, pormaskfile)
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

    #############################'
    # need this since we are still using the old variable derivation through Var() object
    vari = {}
    vari['name'] = variable.name
    vari['field'] = variable.field
    ##############################

    # Build input and output file names
    infileslist = get_cf_infile(project_info, current_diag, model, vari)[variable.name]

    # VP-FIXME-question : the code can glob multiple files, but how to handle the case when globbing fails; currently
    # the diagnostic is run on only the first file
    if len(infileslist) == 1:
        infiles = infileslist[0]
    else:
        logger.info("Found multiple netCDF files for current diagnostic")
        logger.info("netCDF globbing has failed")
        logger.info("Running diagnostic ONLY on the first file")
        infiles = infileslist[0]

    outfilename = get_output_file(project_info, model, vari).split('/')[-1]
    logger.info("Reformatted file name: %s", outfilename)
    # get full outpaths - original cmorized files that are preserved all through the process
    fullpath = get_output_file(project_info, model, vari)
    logger.info("Reformatted target: %s", fullpath)

    # indir is hardcoded to keep things tight; could be an option to namelist
    project_info['TEMPORARY']['indir_path'] = project_info['GLOBAL']['run_dir']
    project_info['TEMPORARY']['outfile_fullpath'] = fullpath
    project_info['TEMPORARY']['infile_path'] = infiles
    project_info['TEMPORARY']['start_year'] = model['start_year']
    project_info['TEMPORARY']['end_year'] = model['end_year']
    if project_name == 'CMIP5':
        project_info['TEMPORARY']['ensemble'] = model['ensemble']
    project_info['TEMPORARY']['variable'] = variable.name
    project_info['TEMPORARY']['field'] = variable.field
    logger.info("Gathering runtime variables:")
    logger.info("Project is %s", model['project'])
    logger.info("Model is %s", model['name'])
    if project_name == 'CMIP5':
        logger.info("Ensemble is %s", model['ensemble'])
    logger.info("Full IN path is %s", infiles)
    logger.info("Full OUT path is %s", fullpath)

    ################## START CHANGING STUFF ###########################################

    # check if we want to save at each step
    save_intermediary_cubes = project_info['GLOBAL']['save_intermediary_cubes']
    # and initialize the latest_saver iterative variable
    latest_saver = project_info['TEMPORARY']['outfile_fullpath']

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
        if not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath']):

            logger.info(" Calling %s to check/reformat model data", reformat_script)

            run_executable(reformat_script, project_info, verbosity,
                           exit_on_warning)
        if 'NO_REFORMAT' in reformat_script:
            pass
        else:
            if (not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath'])):
                raise exceptions.IOError(2, "Expected reformatted file isn't available: ",
                                         project_info['TEMPORARY']['outfile_fullpath'])

        # load cube now, if available
        reft_cube = iris.load_cube(project_info['TEMPORARY']['outfile_fullpath'])

    ################## 0. CMOR_REFORMAT (PY version) ##############
    # New code: cmor_check.py (by Javier Vegas)
    elif cmor_reformat_type == 'py' and project_name == 'CMIP5':
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
            # remember naming conventions
            # IN: infiles
            # OUT: project_info['TEMPORARY']['outfile_fullpath']
            fixes = Fix.get_fixes(project_name, model_name, var_name)

            def apply_file_fixes(file_to_fix):
                for next_fix in fixes:
                    file_to_fix = next_fix.fix_file(file_to_fix)
                return file_to_fix

            # infiles is always a fullpath = single string
            files = apply_file_fixes(infiles)

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
                reft_cube_0 = iris.load(files, var_cons, callback=merge_callback)[0]

            for fix in fixes:
                reft_cube_0 = fix.fix_metadata(reft_cube_0)

            # Check metadata before any preprocessing start
            var_info = variables_info.get_variable(table, var_name)
            checker = CC(reft_cube_0, var_info, automatic_fixes=True)
            checker.check_metadata()

            # apply time gating so we minimize cube size
            yr1 = int(model['start_year'])
            yr2 = int(model['end_year'])
            reft_cube = pt.time_slice(reft_cube_0, yr1,1,1, yr2,12,31)

            for fix in fixes:
                reft_cube = fix.fix_data(reft_cube)

            # Check data after time (and maybe lat-lon slicing)
            checker = CC(reft_cube, var_info, automatic_fixes=True)
            checker.check_data()

            # save reformatted cube
            # saving the cube allows for more steps to be performed
            # on the (un)regridded files
            if save_intermediary_cubes is True:
                latest_saver = latest_saver.strip('.nc') + '_cmor.nc'
                iris.save(reft_cube, latest_saver)

            # save to default path (this file will change with
            # every iteration of preprocess step)
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

        except (iris.exceptions.ConstraintMismatchError,
                iris.exceptions.ConcatenateError,
                CCE) as ex:
            logger.warning("%s", ex)

    else:
        reft_cube = iris.load_cube(infiles)
        iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

    #################### 1. LAND/OCEAN/PORO MASK VIA sftlf/sftof/mrsofc FILE ####################################################

    # Land Mask
    if mask_land is True:
        if os.path.isfile(lmaskfile_path):
            logger.info(" Using mask file %s", lmaskfile_path)
            l_mask = iris.load_cube(lmaskfile_path)

            reft_cube = pt.fx_mask(reft_cube, l_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(" CMIP5 - checking cube after applying land mask...")
                checker = CC(reft_cube, var_info, automatic_fixes=True)
                checker.check_data()
            #######################

            if save_intermediary_cubes is True:
                latest_saver = latest_saver.strip('.nc') + '_land-mask.nc'
                iris.save(reft_cube, latest_saver)

            # save to default path (this file will change with
            # every iteration of preprocess step)
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

        else:
            logger.info("Could not find LAND mask file %s", lmaskfile)


    # Ocean Mask
    if mask_ocean is True:
        if os.path.isfile(omaskfile_path):
            logger.info("Using OCEAN mask file %s", omaskfile_path)
            o_mask = iris.load_cube(omaskfile_path)

            reft_cube = pt.fx_mask(reft_cube, o_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(" CMIP5 - checking cube after applying ocean mask...")
                checker = CC(reft_cube, var_info, automatic_fixes=True)
                checker.check_data()
            #######################

            if save_intermediary_cubes is True:
                latest_saver = latest_saver.strip('.nc') + '_ocean-mask.nc'
                iris.save(reft_cube, latest_saver)

            # save to default path (this file will change with
            # every iteration of preprocess step)
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

        else:
            logger.info("Could not find OCEAN mask file %s", omaskfile)

    # Poro Mask
    if mask_poro is True:
        if os.path.isfile(pormaskfile_path):
            logger.info(" Using PORO mask file %s", pormaskfile_path)
            por_mask = iris.load_cube(pormaskfile_path)

            reft_cube = pt.fx_mask(reft_cube, por_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(" CMIP5 - checking cube after applying poro mask...")
                checker = CC(reft_cube, var_info, automatic_fixes=True)
                checker.check_data()
            #######################

            if save_intermediary_cubes is True:
                latest_saver = latest_saver.strip('.nc') + '_poro-mask.nc'
                iris.save(reft_cube, latest_saver)

            # save to default path (this file will change with
            # every iteration of preprocess step)
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

        else:
            logger.info("Could not find OCEAN mask file %s", omaskfile)

    #################### 2. LEVEL/TIME/AREA OPS #################################################

    # SELECT LEVEL: PREPROCESS['select_level'] = DICT
    # DICT has keys: 'levels' and 'scheme'
    # DICT['levels']: could be a single element or a [list] with multiple values
    # DICT['scheme']: supports linear and nearest
    if select_level != 'None':

        # check if dictionary
        if isinstance(select_level, dict) is False:
            logger.warning("In namelist - select_level must be a dictionary with keys levels and scheme - no select level! ")
            nsl = 0
        else:
            nsl = 2

        # try get the parameters
        try:
            levels = select_level['levels']
            scheme = select_level['scheme']
        except (KeyError):
            logger.warning("select_level keys must be levels: and scheme: - no select level! ")
            nsl = 0

        # check scheme value
        if scheme not in vertical_schemes:
            logger.warning("Select level scheme should be one of the allowed ones %s - no select level! ", vertical_schemes)
            nsl = 0

        # check levels value
        if isinstance(levels, basestring):
            try:
                vlevels = levels.split(',')
                psl = 2
            except (AttributeError, "'int' object has no attribute 'split'"):
                logger.warning("Vertical levels must be either int, string or list %r", levels)
        elif isinstance(levels, int) or isinstance(levels, float):
            vlevels = [levels]
            psl = 2
        elif isinstance(levels, list):
            vlevels = levels
            psl = 2
        else:
            logger.warning("Vertical levels must be int, float, string or list (of ints or floats) %s - no select level!", levels)
            psl = 0

        if psl > 0 and nsl > 0:

            logger.info("Calling regrid to select vertical level %s Pa with scheme %s", vlevels, scheme)

            # check cube has 'air_pressure' coordinate
            ap = reft_cube.coord('air_pressure')
            if ap is None:
                logger.warning("Trying to select level but cube has no air_pressure coordinate ")
            else:
                ap_vals = ap.points
                logger.info("Cube air pressure values " + str(ap_vals))


                # warn if selected levels are outside data bounds
                # these cases will automatically make cmor.check() crash anyway
                if min(vlevels) < min(ap_vals):
                    logger.warning("Selected pressure level below lowest data point, expect large extrapolation errors! ")
                if max(vlevels) > max(ap_vals):
                    logger.warning("Selected pressure level above highest data point, expect large extrapolation errors! ")

                # call vinterp(interpolate)
                reft_cube = vinterp(reft_cube, vlevels, scheme)

                # check cube
                #########################
                if project_name == 'CMIP5':
                    logger.info(" CMIP5 - checking cube after selecting level(s)...")
                    checker = CC(reft_cube, var_info, automatic_fixes=True)
                    checker.check_data()
                #########################

                # save intermediary
                if save_intermediary_cubes is True:
                    vnames = '_'.join(['SL' + str(v) for v in vlevels])
                    latest_saver = latest_saver.strip('.nc') + '_' + vnames + '.nc'
                    iris.save(reft_cube, latest_saver)

                # save to default path (this file will change with
                # every iteration of preprocess step)
                iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

    #################### 3. REGRID ####################################################
    if target_grid != 'None':

        # we will regrid according to whatever regridding scheme and reference grids are needed
        # and create new regridded files from the original cmorized/masked ones
        # but preserving thos original files (optionally)
        logger.info("Calling regrid to regrid model data onto %s grid", target_grid)

        # get the floating cube
        src_cube = reft_cube

        logger.info("Source cube to be regridded --->\n%s", src_cube)

        # try return a cube for regrid target
        # target_grid = could simply be a netCDF file or string model
        # descriptor eg 'ref_model'; currently netCDF and ref_model labels are implemented
        try:
            tgt_grid_cube = iris.load_cube(target_grid)
            logger.info("Target regrid cube summary --->\n%s", tgt_grid_cube)
            if regrid_scheme:
                rgc = regrid(src_cube, tgt_regrid_cube, regrid_scheme)
            else:
                logger.info("No regrid scheme specified, assuming linear")
                rgc = regrid(src_cube, tgt_regrid_cube, 'linear')

            logger.info("Regridded cube summary --->\n%s", rgc)
            reft_cube = rgc
            print(reft_cube)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(" CMIP5 - checking cube after regridding on input netCDF file...")
                checker = CC(reft_cube, var_info, automatic_fixes=True)
                checker.check_data()
            #######################

            # save-append to outfile fullpath list to be further processed
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])
        except (IOError, iris.exceptions.IrisError) as exc:
            logger.info("Target %s is not a netCDF file", target_grid)
            pass
            # continue to see what regridding we need
            # ref_model regrid string descriptor
            if target_grid == 'ref_model':
                additional_models_dicts = current_diag.additional_models
                # identify the current variable
                for var in current_diag.variables:
                    if var['name'] == variable.name:
                        ref_model_list = var['ref_model']

                        # check if the ref_model list is populated
                        if len(ref_model_list) > 0:
                            # always regrid only on the first ref_model
                            ref_model = ref_model_list[0]
                            for obs_model in additional_models_dicts:
                                if obs_model['name'] == ref_model:

                                    logger.info("Regridding on ref_model %s", ref_model)
                                    tgt_nc_grid = get_cf_infile(project_info, current_diag, obs_model, vari)[variable.name][0]
                                    tgt_grid_cube = iris.load_cube(tgt_nc_grid)

                                    logger.info("Target regrid cube summary --->\n%s", tgt_grid_cube)

                                    if regrid_scheme:
                                        rgc = regrid(src_cube, tgt_grid_cube, regrid_scheme)
                                    else:
                                        logger.info("No regrid scheme specified, assuming linear")
                                        rgc = regrid(src_cube, tgt_grid_cube, 'linear')
                                    logger.info("Regridded cube summary --->\n%s", rgc)

                                    # check cube
                                    ###################
                                    reft_cube = rgc
                                    if project_name == 'CMIP5':
                                        logger.info(" CMIP5 - checking cube after regridding on REF model...")
                                        checker = CC(reft_cube, var_info, automatic_fixes=True)
                                        checker.check_data()
                                    ####################

                                    # save-append to outfile fullpath list to be further processed
                                    iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

                                    if save_intermediary_cubes is True:
                                        latest_saver = latest_saver.strip('.nc') + '_regridded_on_' + ref_model + '.nc'
                                        iris.save(reft_cube, latest_saver)


                        # otherwise don't do anything
                        else:
                            logger.info("No regridding model specified in variables[ref_model]. Skipping regridding.")

            else:
                # assume and check it is of XxY form
                if isinstance(target_grid.split('x')[0], basestring) and isinstance(target_grid.split('x')[1], basestring):
                    logger.info("Target regrid is XxY: %s", target_grid)
                    if regrid_scheme:
                        rgc = regrid(src_cube, target_grid, regrid_scheme)
                    else:
                        logger.info("No regrid scheme specified, assuming linear")
                        rgc = regrid(src_cube, target_grid, 'linear')

                    logger.info("Regridded cube summary --->\n%s", rgc)

                    # check cube
                    ######################
                    reft_cube = rgc
                    if project_name == 'CMIP5':
                        logger.info(" CMIP5 - checking cube after regridding on MxN cells...")
                        checker = CC(reft_cube, var_info, automatic_fixes=True)
                        checker.check_data()
                    #######################

                    # save-append to outfile fullpath list to be further processed
                    iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

                    if save_intermediary_cubes is True:
                        latest_saver = latest_saver.strip('.nc') + '_regrid_' + target_grid + '.nc'
                        iris.save(reft_cube, latest_saver)


    ############ 4. MASK FILL VALUES ##############################################
    if mask_fillvalues != 'None':

        # Brief explanation of functionality:
        # Mask fill values
        # -----------------
        # Suppose now we want to apply a mask that discards all the grid
        # points that have less than 20 data points in each 5 year block
        # c = mask_cube_counts(mycube, val_thr=1.0, counts_thr=20, time_window=5)[2]
        # Here: - mycube is the cube I am working on;
        #       - val_thr is a value threshold
        #       - counts_thr is a counts threshold for discarding
        #       - time_window is the window we count and apply counts_thr

        # check if dictionary
        if isinstance(mask_fillvalues, dict) is False:
            logger.warning("In namelist - mask_fillvalues must be a dictionary with keys min_value, threshold_percent and time_window - no mask_fillvalues! ")

        else:

            # try get the parameters
            try:
                val_thr = mask_fillvalues['min_value']
                percentage = mask_fillvalues['threshold_percent']
                time_window = int(mask_fillvalues['time_window'])

                logger.info(" >>> preprocess.py >>> Creating fillvalues mask...")
                # basic checks
                if percentage > 1.0:
                    logger.warning(" >>> preprocess.py >>> Percentage of missing values should be < 1.0")
                nr_time_points = len(reft_cube.coord('time').points)
                if time_window > nr_time_points:
                    logger.warning(" >>> preprocess.py >>> Time window (in time units) larger than total time span")

                # round to lower integer always
                max_counts_per_time_window = int(nr_time_points / time_window)
                count_thr = int(max_counts_per_time_window*percentage)

                # apply the mask
                reft_cube = pt.mask_cube_counts(reft_cube, val_thr, count_thr, time_window)[2]

                # save if needed
                if save_intermediary_cubes is True:
                    latest_saver = latest_saver.strip('.nc') + '_mfv.nc'
                    iris.save(reft_cube, latest_saver)

                # save-append to outfile fullpath list to be further processed
                iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

            except (KeyError):
                logger.warning("select_level keys must be levels: and scheme: - no select level! ")

    ############ 5. MULTIMODEL MEAN and anything else that needs ALL models #######
    # These steps are done outside the main Preprocess loop since they need ALL the models

    ############ FINISH all PREPROCESSING and delete environment
    del(project_info['TEMPORARY'])

    # return the latest cube, and latest path
    return reft_cube, latest_saver

# functions that perform collective models analysis
# apply multimodel means and/or anything else needed
# the order in which functions are called matters !!!
def multimodel_mean(cube_collection, path_collection):
    # time average
    means_list = [mycube.collapsed('time', iris.analysis.MEAN) for mycube in cube_collection]

    # global mean
    means = [np.mean(m.data) for m in means_list]
    logger.info("Multimodel global means: %s", means)

    # seasonal mean
    # need to fix this !
    #smeans_cubes = [pt.seasonal_mean(mycube) for mycube in cube_collection]
    #print(smeans_cubes)
    #smeans = [np.mean(c.data) for c in smeans_cubes]
    #logger.info(" >>> preprocess.py >>> Multimodel seasonal global means: %s" % str(smeans))


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

    def select_base_vars(self, variables, model, current_diag, project_info):

        for base_var in variables:
            # FIXME we should not need this anymore??
            # Check if this variable should be excluded for this model-id
            if self.id_is_explicitly_excluded(base_var, model):
                continue

            # first try: use base variables provided by variable_defs script
            os.environ['__ESMValTool_base_var'] = base_var.name

            # need this since we are still using the old variable derivation through Var() object
            vari = {}
            vari['name'] = base_var.name
            vari['field'] = base_var.field
            ##############################

            infiles = get_cf_infile(project_info, current_diag, model, vari)[base_var.name]

            if len(infiles) == 0:
                logger.info("No input files found for %s (%s)", base_var.name, base_var.field)

                base_var.var = base_var.var0
                base_var.fld = base_var.fld0

                # try again with input variable = base variable (non derived)
                infile = get_cf_infile(project_info, current_diag, model, vari)[base_var.name]

                if len(infile) == 0:
                    raise exceptions.IOError(2, "No input files found in ",
                                             infile)
                else:
                    logger.info("Using %s (%s)", base_var.name, base_var.field)
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
