"""
New module to enable ESMValToolKit to deal with
the new yaml parser and simplified interface_scripts
toolbox. Author: Valeriu Predoi, University of Reading,
Initial version: August 2017
contact: valeriu.predoi@ncas.ac.uk
"""
from __future__ import print_function

import copy
import logging
import os

import iris
import iris.exceptions
import numpy as np

from ..preprocessor.mask import fx_mask, mask_cube_counts
from ..preprocessor.regrid import regrid, vertical_schemes, vinterp
from ..preprocessor.time_area import time_slice
from .data_finder import get_input_filelist, get_output_file
from .fixes.fix import Fix
from .launchers import run_executable
from .preprocessing_tools import glob, merge_callback

logger = logging.getLogger(__name__)

#######################################################
# This script contains basic functionalities
# It contains all variables to call other
# more specialized modules for prerocessing purposes
########################################################


def get_cf_infile(project_info, current_diag, model, current_var_dict):
    """@brief Path to input netCDF files for models
       @param project_info all info dictionary
       @param currDiag(dict) current diagnostic
       @param model(dict) what model we looking at
       This function supports multiple models and variables
       and operates on a single diagnostic
       Returns a dictionary {var1:[model_data], var2:[model_data], ...}
    """
    # get variables; data files is a dictionary keyed on variable
    data_files = {}
    variables = current_diag.variables

    for var in variables:
        full_paths = get_input_filelist(project_info, model, var)
        if len(full_paths) == 0:
            logger.info("Could not find any data files for %s", model['name'])
        else:
            data_files[var['name']] = full_paths

    return data_files


def get_cf_areafile(project_info, model):
    """ @brief Returns the path to the areacello file
        This function looks for the areafile of the ocean grid
    """
    areadir = model["path"]
    areafile = 'areacello_fx_' + model["name"] + "_" + model["exp"] + \
               "_r0i0p0.nc"

    return os.path.join(areadir, areafile)


def preprocess(project_info, variable, model, current_diag,
               cmor_reformat_type):

    #################################################################
    # The Big Mean PREPROCESS Machine                               #
    #################################################################

    # Starting to look at preprocessing needs
    verbosity = project_info["GLOBAL"]["verbosity"]
    logger.info("Looking at (namelist) PREPROCESS section now")
    # initialize the environment variables dictionary
    project_info['TEMPORARY'] = {}
    # key in the prperocess
    prp = project_info['PREPROCESS']

    #################################################################
    # PRERQUISITES: GET THE PARAMETERS From config file             #
    #################################################################

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
                lmaskdir = os.path.join(
                    model["path"],
                    model["exp"],
                    'fx',
                    'sftlf',
                    model["name"],
                    'r0i0p0',
                )
                lmaskfile = 'sftlf_fx_{}_{}_r0i0p0.nc'.format(
                    model["name"], model["exp"])
                lmaskfile_path = os.path.join(lmaskdir, lmaskfile)

        # ocean (keeps only ocean regions)

            if prp[k] == 'ocean':
                mask_ocean = True
                if lmaskfile_path is not None:
                    project_info['TEMPORARY'][
                        'lmaskfile_path'] = lmaskfile_path
                omaskdir = model['path']
                omaskfile = 'sftof_fx_{}_{}_r0i0p0.nc'.format(
                    model["name"], model["exp"])
                omaskfile_path = os.path.join(omaskdir, omaskfile)
                if omaskfile_path is not None:
                    project_info['TEMPORARY'][
                        'omaskfile_path'] = omaskfile_path

        # poro (keeps only poro regions)
        if k == 'mask_poro':

            if prp[k] is not False:
                mask_poro = True
                pormaskdir = os.path.join(
                    model["path"],
                    model["exp"],
                    'fx',
                    'mrsofc',
                    model["name"],
                    'r0i0p0',
                )
                pormaskfile = 'mrsofc_fx_{}_{}_r0i0p0.nc'.format(
                    model["name"], model["exp"])
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
                    project_info['TEMPORARY'][
                        'hgridfile_path'] = hgridfile_path
        if k == 'get_cf_zgridfile':
            if prp[k] is not 'None':
                filepath = fx_files['nemo_zgrid_file'].get_fullpath()
                zgridfile_path = os.path.join(indata_root, filepath)
                if zgridfile_path is not None:
                    project_info['TEMPORARY'][
                        'zgridfile_path'] = zgridfile_path
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
                    new_fx[key] = os.path.join(indata_root,
                                               curr_fx[key].get_fullpath())
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

    #################################################################
    # ENVIRONMENT VARIABLES

    # Initialize all needed files and variables (Leon: EEEvvryboooodyyy!!!)

    model_name = model['name']
    project_name = model['project']
    project_info['RUNTIME']['model'] = model_name
    project_info['RUNTIME']['project'] = project_name
    project_info['RUNTIME']['project_basename'] = 'OBS'
    exit_on_warning = project_info['GLOBAL'].get('exit_on_warning', False)

    # Variable put in environment
    os.environ['__ESMValTool_base_var'] = variable.name

    #############################
    # need this since we are still using the old variable derivation through
    # Var() object
    vari = {}
    vari['name'] = variable.name
    vari['field'] = variable.field
    ##############################

    # Build input and output file names
    infiles = get_cf_infile(project_info, current_diag, model,
                            vari)[variable.name]

    outfilename = get_output_file(project_info, model, vari).split('/')[-1]
    logger.info("Reformatted file name: %s", outfilename)
    # get full outpaths - original cmorized files that are preserved
    # all through the process
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
    logger.info("Full IN path is %s", str(infiles))
    logger.info("Full OUT path is %s", fullpath)

    preproc_model_path = os.path.join(project_info['GLOBAL']['preproc_dir'],
                                      model['project'])
    if not os.path.isdir(preproc_model_path):
        os.makedirs(preproc_model_path)

    #################################################################
    # START CHANGING STUFF

    # check if we want to save at each step
    save_intermediary_cubes = project_info['GLOBAL']['save_intermediary_cubes']
    # and initialize the latest_saver iterative variable
    latest_saver = project_info['TEMPORARY']['outfile_fullpath']

    #################################################################
    # 0. CMOR_REFORMAT (NCL version)

    # Legacy code that will be purged in the future
    # Check if the current project has a specific reformat routine,
    # otherwise use default
    # the cmor reformat is applied only once per variable
    if cmor_reformat_type == 'ncl':
        if os.path.isdir("reformat_scripts/" + model['project']):
            which_reformat = model['project']
        else:
            which_reformat = 'default'

        reformat_script = os.path.join(
            "reformat_scripts", which_reformat,
            "reformat_" + which_reformat + "_main.ncl")
        if not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath']):

            logger.info(" Calling %s to check/reformat model data",
                        reformat_script)

            run_executable(reformat_script, project_info, verbosity,
                           exit_on_warning)
        if 'NO_REFORMAT' in reformat_script:
            pass
        else:
            if (not os.path.isfile(
                    project_info['TEMPORARY']['outfile_fullpath'])):
                raise IOError(2, "Expected reformatted file isn't available: ",
                              project_info['TEMPORARY']['outfile_fullpath'])

        # load cube now, if available
        reft_cube = iris.load_cube(
            project_info['TEMPORARY']['outfile_fullpath'])

    #################################################################
    # 0. CMOR_REFORMAT (PY version)

    # New code: cmor_check.py (by Javier Vegas)
    elif cmor_reformat_type == 'py':
        # needed imports
        from .cmor_check import CMORCheck
        from .cmor_check import CMORCheckError
        import warnings
        from .variable_info import CMIP5Info

        variables_info = CMIP5Info()

        var_name = variable.name
        if project_name == 'CMIP5':
            table = model['mip']
        else:
            table = project_info['MODELS'][0]['mip'] 

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

            # infiles are fullpaths
            files = [apply_file_fixes(infile) for infile in infiles]
            cfilelist = iris.cube.CubeList()

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
                for cfile in files:
                    reft_cube_i = iris.load(
                        cfile, var_cons, callback=merge_callback)[0]
                    cfilelist.append(reft_cube_i)

            # concatenate if needed
            if len(cfilelist) > 1:
                try:
                    iris.util.unify_time_units(cfilelist)
                    reft_cube_0 = cfilelist.concatenate()[0]
                except iris.exceptions.ConcatenateError as exc:
                    error_message = "Problem trying to concatenate cubes: %s"
                    logger.error(error_message, exc)
                    logger.debug('Cubes to concatenate:')
                    for cube in cfilelist:
                        logger.debug(str(cube))
                    return None, None

            else:
                reft_cube_0 = cfilelist[0]

            # fix again
            for fix in fixes:
                reft_cube_0 = fix.fix_metadata(reft_cube_0)

            # Check metadata before any preprocessing start
            var_info = variables_info.get_variable(table, var_name)
            checker = CMORCheck(reft_cube_0, var_info, automatic_fixes=True)
            checker.check_metadata()

            # apply time gating so we minimize cube size
            yr1 = int(model['start_year'])
            yr2 = int(model['end_year'])
            reft_cube = time_slice(reft_cube_0, yr1, 1, 1, yr2, 12, 31)

            for fix in fixes:
                reft_cube = fix.fix_data(reft_cube)

            # Check data after time (and maybe lat-lon slicing)
            checker = CMORCheck(reft_cube, var_info, automatic_fixes=True)
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

        except (iris.exceptions.ConcatenateError, CMORCheckError) as ex:
            logger.error("%s", ex)
            return None, None

    else:
        # check if we need to concatenate
        if len(infiles) > 1:
            reft_cube = glob(infiles, variable.name)
        else:
            reft_cube = iris.load_cube(infiles)
        iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])

    #################################################################
    # 1. LAND/OCEAN/PORO MASK VIA sftlf/sftof/mrsofc FILE

    # Land Mask
    if mask_land is True:
        if os.path.isfile(lmaskfile_path):
            logger.info(" Using mask file %s", lmaskfile_path)
            l_mask = iris.load_cube(lmaskfile_path)

            reft_cube = fx_mask(reft_cube, l_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(
                    " CMIP5 - checking cube after applying land mask...")
                checker = CMORCheck(reft_cube, var_info, automatic_fixes=True)
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

            reft_cube = fx_mask(reft_cube, o_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(
                    " CMIP5 - checking cube after applying ocean mask...")
                checker = CMORCheck(reft_cube, var_info, automatic_fixes=True)
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

            reft_cube = fx_mask(reft_cube, por_mask)

            # check cube
            ######################
            if project_name == 'CMIP5':
                logger.info(
                    " CMIP5 - checking cube after applying poro mask...")
                checker = CMORCheck(reft_cube, var_info, automatic_fixes=True)
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

    #################################################################
    # 2. LEVEL/TIME/AREA OPS

    # SELECT LEVEL: PREPROCESS['select_level'] = DICT
    # DICT has keys: 'levels' and 'scheme'
    # DICT['levels']: could be a single element or a [list]
    #                 with multiple values
    # DICT['scheme']: supports linear and nearest
    if select_level != 'None':

        # check if dictionary
        if isinstance(select_level, dict) is False:
            logger.warning("In namelist - select_level must be a dictionary "
                           "with keys levels and scheme - no select level!")
            nsl = 0
        else:
            nsl = 2

        # try get the parameters
        try:
            levels = select_level['levels']
            scheme = select_level['scheme']
        except KeyError:
            logger.warning("select_level keys must be levels: and scheme: "
                           "- no select level!")
            nsl = 0

        # check scheme value
        if scheme not in vertical_schemes:
            logger.warning(
                "Select level scheme should be one of the allowed ones %s "
                "- no select level!", vertical_schemes)
            nsl = 0

        # check levels value
        if isinstance(levels, str):
            try:
                vlevels = levels.split(',')
                psl = 2
            except (AttributeError, "'int' object has no attribute 'split'"):
                logger.warning("Vertical levels must be either int, string "
                               "or list %r", levels)
        elif isinstance(levels, int) or isinstance(levels, float):
            vlevels = [levels]
            psl = 2
        elif isinstance(levels, list):
            vlevels = levels
            psl = 2
        else:
            logger.warning("Vertical levels must be int, float, string or "
                           "list (of ints or floats) %s - no select level!",
                           levels)
            psl = 0

        if psl > 0 and nsl > 0:

            logger.info("Calling regrid to select vertical level %s Pa "
                        "with scheme %s", vlevels, scheme)

            # check cube has 'air_pressure' coordinate
            ap = reft_cube.coord('air_pressure')
            if ap is None:
                logger.warning("Trying to select level but cube has no "
                               "air_pressure coordinate")
            else:
                ap_vals = ap.points
                logger.info("Cube air pressure values " + str(ap_vals))

                # warn if selected levels are outside data bounds
                # these cases will automatically make cmor.check() crash anyway
                if min(vlevels) < min(ap_vals):
                    logger.warning("Selected pressure level below lowest "
                                   "data point, expect large extrapolation "
                                   "errors!")
                if max(vlevels) > max(ap_vals):
                    logger.warning("Selected pressure level above highest "
                                   "data point, expect large extrapolation "
                                   "errors!")

                # call vinterp(interpolate)
                reft_cube = vinterp(reft_cube, vlevels, scheme)

                # check cube
                #########################
                if project_name == 'CMIP5':
                    logger.info(
                        " CMIP5 - checking cube after selecting level(s)...")
                    checker = CMORCheck(reft_cube, var_info,
                                        automatic_fixes=True)
                    checker.check_data()
                #########################

                # save intermediary
                if save_intermediary_cubes is True:
                    vnames = '_'.join(['SL' + str(v) for v in vlevels])
                    latest_saver = latest_saver.strip(
                        '.nc') + '_' + vnames + '.nc'
                    iris.save(reft_cube, latest_saver)

                # save to default path (this file will change with
                # every iteration of preprocess step)
                iris.save(reft_cube,
                          project_info['TEMPORARY']['outfile_fullpath'])

    #################################################################
    # 3. REGRID

    if target_grid != 'None':

        # we will regrid according to whatever regridding scheme and
        # reference grids are needed and create new regridded files
        # from the original cmorized/masked ones but preserving those
        # original files (optionally)
        logger.info("Calling regrid to regrid model data onto %s grid",
                    target_grid)

        # get the floating cube
        src_cube = reft_cube

        logger.info("Source cube to be regridded --->\n%s", src_cube)

        # try return a cube for regrid target
        # target_grid = could simply be a netCDF file or string model
        # descriptor eg 'ref_model'; currently netCDF and ref_model
        # labels are implemented
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
                logger.info("CMIP5 - checking cube after regridding on "
                            "input netCDF file...")
                checker = CMORCheck(reft_cube, var_info, automatic_fixes=True)
                checker.check_data()
            #######################

            # save-append to outfile fullpath list to be further processed
            iris.save(reft_cube, project_info['TEMPORARY']['outfile_fullpath'])
        except (IOError, iris.exceptions.IrisError) as exc:
            logger.info("Target %s is not a netCDF file", target_grid)
            # continue to see what regridding we need
            # ref_model regrid string descriptor
            if target_grid == 'ref_model':
                try:
                    additional_models_dicts = current_diag.additional_models
                    if additional_models_dicts is None:
                        additional_models_dicts = project_info['ALLMODELS']
                except (AttributeError, "'Diagnostic' object has no attribute "
                                        "'additional_models'"):
                    logger.info("Regridding on one of the MODELS, "
                                "no ADDITIONAL MODELS specified")
                    additional_models_dicts = project_info['ALLMODELS']
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

                                    logger.info("Regridding on ref_model %s",
                                                ref_model)
                                    tgt_nc_grid = get_cf_infile(
                                        project_info, current_diag, obs_model,
                                        vari)[variable.name][0]

                                    # check if we need to concatenate
                                    if len([tgt_nc_grid]) > 1:
                                        tgt_grid_cube = glob(
                                            tgt_nc_grid, variable.name)
                                    else:
                                        tgt_grid_cube = iris.load_cube(
                                            tgt_nc_grid)

                                    logger.info(
                                        "Target regrid cube summary --->\n%s",
                                        tgt_grid_cube)

                                    if regrid_scheme:
                                        rgc = regrid(src_cube, tgt_grid_cube,
                                                     regrid_scheme)
                                    else:
                                        logger.info(
                                            "No regrid scheme specified, "
                                            "assuming linear")
                                        rgc = regrid(src_cube, tgt_grid_cube,
                                                     'linear')
                                    logger.info("Regridded cube summary"
                                                " --->\n%s", rgc)

                                    # check cube
                                    ###################
                                    reft_cube = rgc
                                    if project_name == 'CMIP5':
                                        logger.info(
                                            "CMIP5 - checking cube after "
                                            "regridding on REF model...")
                                        checker = CMORCheck(
                                            reft_cube,
                                            var_info,
                                            automatic_fixes=True)
                                        checker.check_data()
                                    ####################

                                    # save-append to outfile fullpath list
                                    # to be further processed
                                    iris.save(reft_cube, project_info[
                                        'TEMPORARY']['outfile_fullpath'])

                                    if save_intermediary_cubes is True:
                                        latest_saver = '_'.join([
                                            latest_saver.strip('.nc'),
                                            'regridded_on',
                                            ref_model + '.nc',
                                        ])
                                        iris.save(reft_cube, latest_saver)

                        # otherwise don't do anything
                        else:
                            logger.info(
                                "No regridding model specified in "
                                "variables[ref_model]. Skipping regridding.")

            else:
                # assume and check it is of XxY form
                if isinstance(target_grid.split('x')[0], str) and isinstance(
                        target_grid.split('x')[1], str):
                    logger.info("Target regrid is XxY: %s", target_grid)
                    if regrid_scheme:
                        rgc = regrid(src_cube, target_grid, regrid_scheme)
                    else:
                        logger.info(
                            "No regrid scheme specified, assuming linear")
                        rgc = regrid(src_cube, target_grid, 'linear')

                    logger.info("Regridded cube summary --->\n%s", rgc)

                    # check cube
                    ######################
                    reft_cube = rgc
                    if project_name == 'CMIP5':
                        logger.info("CMIP5 - checking cube after regridding "
                                    "on MxN cells...")
                        checker = CMORCheck(reft_cube, var_info,
                                            automatic_fixes=True)
                        checker.check_data()
                    #######################

                    # save-append to outfile fullpath list to be
                    # further processed
                    iris.save(reft_cube,
                              project_info['TEMPORARY']['outfile_fullpath'])

                    if save_intermediary_cubes is True:
                        latest_saver = latest_saver.strip(
                            '.nc') + '_regrid_' + target_grid + '.nc'
                        iris.save(reft_cube, latest_saver)

    #################################################################
    # 4. MASK FILL VALUES

    if mask_fillvalues != 'None':

        # Brief explanation of functionality:
        # Mask fill values
        # -----------------
        # Suppose now we want to apply a mask that discards all the grid
        # points that have less than 20 data points in each 5 year block
        # c = mask_cube_counts(
        #     mycube, val_thr=1.0, counts_thr=20, time_window=5)[2]
        # Here: - mycube is the cube I am working on;
        #       - val_thr is a value threshold
        #       - counts_thr is a counts threshold for discarding
        #       - time_window is the window we count and apply counts_thr

        # check if dictionary
        if isinstance(mask_fillvalues, dict) is False:
            logger.warning(
                "In namelist - mask_fillvalues must be a dictionary with "
                "keys min_value, threshold_percent and time_window - no "
                "mask_fillvalues!")

        else:

            # try get the parameters
            try:
                val_thr = mask_fillvalues['min_value']
                percentage = mask_fillvalues['threshold_percent']
                time_window = int(mask_fillvalues['time_window'])

                logger.info("Creating fillvalues mask...")
                # basic checks
                if percentage > 1.0:
                    logger.warning("Percentage of missing values should "
                                   "be < 1.0")
                nr_time_points = len(reft_cube.coord('time').points)
                if time_window > nr_time_points:
                    logger.warning("Time window (in time units) larger "
                                   "than total time span")

                # round to lower integer always
                max_counts_per_time_window = int(nr_time_points / time_window)
                count_thr = int(max_counts_per_time_window * percentage)

                # apply the mask
                reft_cube = mask_cube_counts(reft_cube, val_thr, count_thr,
                                             time_window)[2]

                # save if needed
                if save_intermediary_cubes is True:
                    latest_saver = latest_saver.strip('.nc') + '_mfv.nc'
                    iris.save(reft_cube, latest_saver)

                # save-append to outfile fullpath list to be further processed
                iris.save(reft_cube,
                          project_info['TEMPORARY']['outfile_fullpath'])

            except KeyError:
                logger.warning("select_level keys must be levels: and "
                               "scheme: - no select level!")

    #################################################################
    # 5. MULTIMODEL MEAN and anything else that needs ALL models

    # These steps are done outside the main Preprocess loop since
    # they need ALL the models

    #################################################################
    # FINISH all PREPROCESSING and delete environment

    del (project_info['TEMPORARY'])

    # return the latest cube, and latest path
    return reft_cube, latest_saver


# functions that perform collective models analysis
# apply multimodel means and/or anything else needed
# the order in which functions are called matters !!!
def multimodel_mean(cube_collection, path_collection):
    # time average
    means_list = [
        mycube.collapsed('time', iris.analysis.MEAN)
        for mycube in cube_collection if mycube is not None
    ]

    # global mean
    means = [np.mean(m.data) for m in means_list]
    logger.info("Multimodel global means: %s", means)

    # seasonal mean
    # need to fix this !
    # smeans_cubes = [seasonal_mean(mycube) for mycube in cube_collection]
    # print(smeans_cubes)
    # smeans = [np.mean(c.data) for c in smeans_cubes]
    # logger.info("Multimodel seasonal global means: %s", smeans)


#################################################################
# a few definititory functions needed
#################################################################
# variable class
class Var:
    """
    changed from original class diagdef.Var
    """

    def __init__(self, merged_dict, var0, fld0):
        # Write the variable attributes in 'merged_dict'
        # to class object attributes
        for name, value in merged_dict.items():
            setattr(self, name, value)

        # Special cases, actually not sure what they do
        if var0 == "none":
            self.var0 = merged_dict['name']
        else:
            self.var0 = var0
        if fld0 == "none":
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

            f = open(
                os.path.join(variable_def_dir, variable['name'] + ".ncl"), 'r')
            for line in f:
                tokens = line.split()

                if "Requires:" in tokens:
                    # If 'none', return orig. field
                    if tokens[2] == "none" or model_der == "True":
                        dep_vars.append(Var(variable, "none", "none"))
                    else:
                        sub_tokens = tokens[2].split(",")
                        for sub in sub_tokens:
                            element = sub.split(":")

                            # Assume first digit is dimension in field type
                            offset = self.find_dimension_entry(element[1]) - 1

                            e_var = element[0]
                            e_fld = element[1]
                            e_fld = (variable['field'][0:offset + 1] +
                                     e_fld[1 + offset] +
                                     variable['field'][2 + offset] +
                                     e_fld[3 + offset:])

                            del keys
                            del vars
                            dep_var = copy.deepcopy(variable)
                            dep_var['name'] = e_var
                            dep_var['field'] = e_fld

                            dep_vars.append(
                                Var(dep_var, variable.var, variable.fld))

        return dep_vars

    def select_base_vars(self, variables, model, current_diag, project_info):

        for base_var in variables:
            # FIXME we should not need this anymore??
            # Check if this variable should be excluded for this model-id
            if self.id_is_explicitly_excluded(base_var, model):
                continue

            # first try: use base variables provided by variable_defs script
            os.environ['__ESMValTool_base_var'] = base_var.name

            # need this since we are still using the old variable derivation
            # through Var() object
            vari = {}
            vari['name'] = base_var.name
            vari['field'] = base_var.field
            ##############################

            infiles = get_cf_infile(project_info, current_diag, model,
                                    vari)[base_var.name]

            if len(infiles) == 0:
                logger.info("No input files found for %s (%s)", base_var.name,
                            base_var.field)

                base_var.var = base_var.var0
                base_var.fld = base_var.fld0

                # try again with input variable = base variable (non derived)
                infile = get_cf_infile(project_info, current_diag, model,
                                       vari)[base_var.name]

                if len(infile) == 0:
                    raise IOError(2, "No input files found in ", infile)
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
