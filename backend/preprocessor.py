import logging
import projects
from auxiliary import info, error, print_header, ncl_version_check

# class PreProcessor(object):

def preprocess(project_info):
    logging.info("preprocessing")

    for diagnostic in project_info['DIAGNOSTICS']:
        #fetch and check preprocess options

        variables = []
        for variable in variables:


            #we changed the loop order. Is that ok?
            result_cubes = []
            for model in project_info['MODELS']:

                is_derived_var = False
                if is_derived_var:
                    cubes = iris.load(files) # all vars needed for derive

                    cubes = [extract_time(cube) for cube in cubes]

                    cube = derive_var(cubes)
                else:
                    cubes = load_cube() #just the var
                    cube = extract_time(cube)

                #Q: isn't this too late to do the fix?
                #use fix files
                cube = fix_format(cube)

                #should this be second
                #will raise exception if format is incorrect
                check_format(cube)

                #target grid should have right amount of vertical levels
                cube = regrid(cube, target_grid_cube, scheme)

                cube = mask(cube)

                #do we need this at all?
                cube = extract_region(cube)

                cube = calculate_time_statistics(cube)

                iris.save(cube, filename)

                if calculate_multi_model_statistics:
                    result_cubes.append(cube)

            #are there any multi-variable mms's?
            if (calculate_multi_model_statistics):
                mms_cube = calculate_multi_model_statistics(result_cubes)

                iris.save(mms_cube, filename)












    verbosity = project_info['GLOBAL']['verbosity']

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



                # # Rewrite netcdf to expected input format.
                # info("Calling cmor_reformat.py to check/reformat model data",
                #      verbosity, 2)
                # reformat.cmor_reformat(currProject, project_info, base_var, model)
                #
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

