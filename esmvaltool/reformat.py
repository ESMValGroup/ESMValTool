"""Port reformat.py from v1 to v2 (life hack)."""
import logging
import exceptions
import os

logger = logging.getLogger(__name__)


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
    currLauncher = vars(launchers)[suffix + '_launcher']()
    if launcher_arguments is not None:
        currLauncher.arguments = launcher_arguments
    currLauncher.execute(string_to_execute,
                         project_info,
                         verbosity,
                         exit_on_warning)


def cmor_reformat(variable, data_file):
    """Run the oldskool v1 reformat scripts."""
    # define the oldschool dictionary
    project_info = {}

    # project_info['RUNTIME']['project_basename'] = project_basename
    project_info['RUNTIME']['model'] = variable['dataset']
    project_info['RUNTIME']['project'] = variable['project']

    # Build input and output file names
    indir = os.path.basedir(data_file)
    infile = os.path.basefile(data_file)
    fullpath = data_file
    project = variable['project']

    # Area file name for ocean grids
    # areafile_path = currProject.get_cf_areafile(project_info, model)

    # Land-mask file name for land variables
    # lmaskfile_path = currProject.get_cf_lmaskfile(project_info, model)
    # omaskfile_path = currProject.get_cf_omaskfile(project_info, model)

    # Porosity file name for land variables
    # porofile_path = currProject.get_cf_porofile(project_info, model)

    # Additional grid file names for ocean grids, if available (ECEARTH)
    # hgridfile_path = False
    # zgridfile_path = False
    # lsmfile_path = False
    # if hasattr(currProject, "get_cf_hgridfile"):
    #     hgridfile_path = currProject.get_cf_hgridfile(project_info, model)
    # if hasattr(currProject, "get_cf_zgridfile"):
    #     zgridfile_path = currProject.get_cf_zgridfile(project_info, model)
    # if hasattr(currProject, "get_cf_lsmfile"):
    #     lsmfile_path = \
    #         currProject.get_cf_lsmfile(project_info, model, variable.fld)

    # General fx file name entry
    # fx_file_path = False
    # if hasattr(currProject, "get_cf_fx_file"):
    #     fx_file_path = currProject.get_cf_fx_file(project_info, model)

    # Check if the current project has a specific reformat routine,
    # otherwise use default
    if (os.path.isdir("reformat_scripts/" + project)):
        # path was found; however on a case insensitive filesystem
        # like e.g. MacOS, this might be a problem, as
        # directories OBS and obs are considered to be the same
        # This causes errors in the namelist processing
        # a second check is therefore performed here
        if project in os.listdir('reformat_scripts'):
            which_reformat = project
        else:
            which_reformat = 'default'
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
    # project_info['TEMPORARY']['areafile_path'] = areafile_path
    # project_info['TEMPORARY']['lmaskfile_path'] = lmaskfile_path
    # project_info['TEMPORARY']['omaskfile_path'] = omaskfile_path
    # project_info['TEMPORARY']['porofile_path'] = porofile_path
    project_info['TEMPORARY']['start_year'] = variable['start_year']
    project_info['TEMPORARY']['end_year'] = variable['end_year']
    project_info['TEMPORARY']['ensemble'] = variable['ensemble']
    project_info['TEMPORARY']['variable'] = variable['short_name']
    project_info['TEMPORARY']['field'] = variable['field']

    # FX file path
    # if fx_file_path:
    #     project_info['TEMPORARY']['fx_file_path'] = fx_file_path

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

        logger.info("Calling %s to check/reformat model data", reformat_script)

        run_executable(reformat_script, project_info, verbosity,
                       exit_on_warning)
    if 'NO_REFORMAT' in reformat_script:
        pass
    else:
        if (not os.path.isfile(project_info['TEMPORARY']['outfile_fullpath'])):
            raise exceptions.IOError(2, "Expected reformatted file isn't available: ",
                                     project_info['TEMPORARY']['outfile_fullpath'])
    del(project_info['TEMPORARY'])
