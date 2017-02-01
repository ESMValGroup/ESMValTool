#
#  Wrapper to call ncl script that calculates climatology
#
#  2008-12-01  CAF

from auxiliary import info
import os
import pdb
import projects


def climate(currProject, currDiag, project_info, model, field, base_variable, variable_attributes):
    """ @brief Wrapper to call ncl script that calculates climatology
        @param currProject An instance of the current project
        @param project_info Dictionary with all info from the namelist
        @param currDiag Current diagnostic
    """
    infilename = currProject.get_cf_fullpath(project_info, model, field, base_variable, variable_attributes)
    cfield_m, cfield_s, cfield_a = currDiag.get_climate_field_type(field)

    file_monthly_clim = currProject.get_cf_fullpath(project_info, model, cfield_m, base_variable, variable_attributes)
    file_season_clim = currProject.get_cf_fullpath(project_info, model, cfield_s, base_variable, variable_attributes)
    file_annual_clim = currProject.get_cf_fullpath(project_info, model, cfield_a, base_variable, variable_attributes)

    if ((not os.path.isfile(file_monthly_clim)
         or not os.path.isfile(file_season_clim)
         or not os.path.isfile(file_annual_clim))
        or project_info['GLOBAL']['force_processing']):

        verbosity = project_info['GLOBAL']['verbosity']
        project_info['TEMPORARY'] = {}
        project_info['TEMPORARY']['infilename'] = infilename
        project_info['TEMPORARY']['mfile'] = file_monthly_clim
        project_info['TEMPORARY']['sfile'] = file_season_clim
        project_info['TEMPORARY']['afile'] = file_annual_clim
        project_info['TEMPORARY']['base_variable'] = base_variable

        info("  Calling climate.py", verbosity, required_verbosity=1)
        info("   INFILE = " + infilename,  1, verbosity)
        info("   OUTFILES = " + file_monthly_clim,  1, verbosity)
        info("              " + file_season_clim,  1, verbosity)
        info("              " + file_annual_clim,  1, verbosity)

        exit_on_warning = project_info['GLOBAL']['exit_on_warning']

        executable = "./interface_scripts/climate.ncl"
        projects.run_executable(executable, project_info, verbosity, exit_on_warning)

        del(project_info['TEMPORARY'])
