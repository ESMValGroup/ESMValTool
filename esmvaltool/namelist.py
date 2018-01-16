"""Namelist parser"""
from __future__ import print_function

import copy
import fnmatch
import inspect
import logging
import os

import yaml
import yamale

from .interface_scripts.data_finder import (
    get_input_filelist, get_input_filename, get_output_file,
    get_start_end_year)
from .interface_scripts.data_interface import (get_legacy_ncl_env,
                                               write_legacy_ncl_interface)
from .preprocessor import (DEFAULT_ORDER, MULTI_MODEL_FUNCTIONS,
                           PREPROCESSOR_FUNCTIONS, PreprocessingTask)
from .preprocessor._download import synda_search
from .preprocessor._io import concatenate_callback
from .preprocessor._reformat import CMOR_TABLES
from .task import DiagnosticTask, get_independent_tasks, run_tasks

logger = logging.getLogger(__name__)

TASKSEP = os.sep


class NamelistError(Exception):
    """Namelist contains an error."""


def read_namelist_file(filename, config_user, initialize_tasks=True):
    """Read a namelist from file."""
    raw_namelist = check_namelist(filename)
    return Namelist(
        raw_namelist, config_user, initialize_tasks, namelist_file=filename)


def check_namelist_with_schema(filename):
    """Check if the namelist content matches schema."""
    schema_file = os.path.join(
        os.path.dirname(__file__), 'namelist_schema.yml')
    logger.debug("Checking namelist against schema %s", schema_file)
    namelist = yamale.make_data(filename)
    schema = yamale.make_schema(schema_file)
    yamale.validate(schema, namelist)


def check_namelist(filename):
    """Check a namelist file and return it in raw form."""
    # Note that many checks can only be performed after the automatically
    # computed entries have been filled in by creating a Namelist object.
    check_namelist_with_schema(filename)
    with open(filename, 'r') as file:
        raw_namelist = yaml.safe_load(file)

    # TODO: add more checks?
    check_preprocessors(raw_namelist['preprocessors'])
    check_diagnostics(raw_namelist['diagnostics'])
    return raw_namelist


def check_preprocessors(preprocessors):
    """Check preprocessors in namelist"""
    preprocessor_functions = set(DEFAULT_ORDER)
    for name, settings in preprocessors.items():
        invalid_functions = set(settings) - preprocessor_functions
        if invalid_functions:
            raise NamelistError(
                "Unknown function(s) {} in preprocessor {}, choose from: {}"
                .format(invalid_functions, name, ', '.join(DEFAULT_ORDER)))


def check_diagnostics(diagnostics):
    """Check diagnostics in namelist"""
    for name, diagnostic in diagnostics.items():
        if 'scripts' not in diagnostic:
            raise NamelistError("Missing scripts section in diagnostic {}"
                                .format(name))
        if diagnostic['scripts'] is None:
            continue
        for script_name, script in diagnostic['scripts'].items():
            if not script.get('script'):
                raise NamelistError(
                    "No script defined for script {} in diagnostic {}".format(
                        script_name, name))


def check_preprocessor_settings(settings):
    """Check preprocessor settings."""
    # The inspect functions getargspec and getcallargs are deprecated
    # in Python 3, but their replacements are not available in Python 2.
    # TODO: Use the new Python 3 inspect API
    for step in settings:
        if step not in DEFAULT_ORDER:
            raise NamelistError(
                "Unknown preprocessor function '{}', choose from: {}".format(
                    step, ', '.join(DEFAULT_ORDER)))
        function = PREPROCESSOR_FUNCTIONS[step]
        argspec = inspect.getargspec(function)
        args = argspec.args[1:]
        # Check for invalid arguments
        invalid_args = set(settings[step]) - set(args)
        if invalid_args:
            raise NamelistError(
                "Invalid argument(s) {} encountered for preprocessor "
                "function {}".format(invalid_args, step))
        # Check for missing arguments
        defaults = argspec.defaults
        end = None if defaults is None else -len(defaults)
        missing_args = set(args[:end]) - set(settings[step])
        if missing_args:
            raise NamelistError(
                "Missing required argument(s) {} for preprocessor "
                "function {}".format(missing_args, step))
        # Final sanity check in case the above fails to catch a mistake
        try:
            inspect.getcallargs(function, None, **settings[step])
        except TypeError:
            logger.error("Wrong preprocessor function arguments in "
                         "function '%s'", step)
            raise


def check_duplicate_models(models):
    """Check for duplicate models."""
    checked_models_ = []
    for model in models:
        if model in checked_models_:
            raise NamelistError(
                "Duplicate model {} in models section".format(model))
        checked_models_.append(model)


def check_variables(variables):
    """Check variables as derived from namelist"""
    generic = {
        'short_name', 'standard_name', 'field', 'model', 'project',
        'start_year', 'end_year', 'preprocessor', 'diagnostic'
    }
    project_specific = {
        # project: keys
        'CMIP5': {'ensemble', 'mip', 'exp'},
        'CCMVal1': {'ensemble', 'exp'},
        'CCMVal2': {'ensemble', 'exp'},
        'GFDL': {'ensemble', 'realm', 'shift'},
        'OBS': {'type', 'version', 'tier'},
    }
    for variable in variables:
        required = set(generic)
        project = variable.get('project')
        if project in project_specific:
            required.update(project_specific[project])
        missing = required - set(variable)
        if missing:
            raise NamelistError(
                "Missing keys {} from variable {} in diagnostic {}".format(
                    missing,
                    variable.get('short_name'), variable.get('diagnostic')))


def check_data_availability(input_files, start_year, end_year):
    """Check if the required input data is available"""
    if not input_files:
        raise NamelistError("No input files found")

    required_years = set(range(start_year, end_year + 1))
    available_years = set()
    for filename in input_files:
        start, end = get_start_end_year(filename)
        available_years.update(range(start, end + 1))

    missing_years = required_years - available_years
    if missing_years:
        raise NamelistError(
            "No input data available for years {} in files {}".format(
                ", ".join(str(year) for year in missing_years), input_files))


def _get_value(key, variables):
    """Get a value for key by looking at the other variables."""
    values = {variable[key] for variable in variables if key in variable}

    if len(values) > 1:
        raise NamelistError("Ambigous values {} for property {}".format(
            values, key))
    return values.pop()


def _add_cmor_info(variable, keys):
    """Add information from CMOR tables to variable."""
    if variable['project'] in CMOR_TABLES:
        variable_info = CMOR_TABLES[variable['project']].get_variable(
            variable['mip'], variable['short_name'])
        for key in keys:
            if key not in variable and hasattr(variable_info, key):
                value = getattr(variable_info, key)
                if value is not None:
                    variable[key] = value


def _update_target_grid(variable, all_variables, settings, config_user):
    """Replace the target grid model name with a filename if needed.

    This only works if the file is found by get_input_filelist.
    """
    if not settings.get('regrid') or 'target_grid' not in settings['regrid']:
        return

    target_grid = settings['regrid']['target_grid']

    if variable['model'] == target_grid:
        # No need to regrid model onto itself.
        del settings['regrid']
        return

    # Replace the target grid model name with a filename.
    for target_variable in all_variables:
        if target_variable['model'] == target_grid:
            files = get_input_filelist(
                variable=target_variable,
                rootpath=config_user['rootpath'],
                drs=config_user['drs'])
            target_file = files[0]
            settings['regrid']['target_grid'] = target_file
            return


def _get_default_settings(variable, config_user):
    """Get default preprocessor settings."""
    settings = {}

    # Set up downloading using synda if requested.
    if config_user['synda_download']:
        local_dir = os.path.dirname(
            get_input_filename(
                variable=variable,
                rootpath=config_user['rootpath'],
                drs=config_user['drs']))
        settings['download'] = {
            'dest_folder': local_dir,
        }

    # Configure loading
    settings['load'] = {
        'constraints': variable['standard_name'],
        'callback': concatenate_callback,
        'filename': variable['filename'],
    }

    # Configure fixes
    fix = {
        'project': variable['project'],
        'model': variable['model'],
        'short_name': variable['short_name'],
    }
    settings['fix_file'] = dict(fix)
    # Only supply mip if the CMOR check fixes are implemented.
    if variable.get('cmor_table'):
        fix['cmor_table'] = variable['cmor_table']
        fix['mip'] = variable['mip']
    settings['fix_metadata'] = dict(fix)
    settings['fix_data'] = dict(fix)

    # Configure time extraction
    settings['extract_time'] = {
        'yr1': variable['start_year'],
        'yr2': variable['end_year'] + 1,
        'mo1': 1,
        'mo2': 1,
        'd1': 1,
        'd2': 1,
    }

    # Configure final CMOR data check
    if variable.get('cmor_table'):
        settings['cmor_check_data'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
        }

    # Configure saving cubes to file
    settings['save'] = {}

    return settings


def _get_input_files(variable, config_user):
    """Get the input files for a single model"""
    # Find input files locally.
    input_files = get_input_filelist(
        variable=variable,
        rootpath=config_user['rootpath'],
        drs=config_user['drs'])

    # Set up downloading using synda if requested.
    # Do not download if files are already available locally.
    if config_user['synda_download'] and not input_files:
        input_files = synda_search(variable)

    logger.info("Using input files:\n%s", '\n'.join(input_files))
    check_data_availability(input_files, variable['start_year'],
                            variable['end_year'])

    return input_files


def _apply_preprocessor_settings(settings, profile_settings):
    """Apply settings from preprocessor profile."""
    for step, args in profile_settings.items():
        # Remove disabled preprocessor functions
        if args is False:
            if step in settings:
                del settings[step]
            continue
        # Enable/update functions without keywords
        if step not in settings:
            settings[step] = {}
        if isinstance(args, dict):
            settings[step].update(args)


def _get_preprocessor_settings(variables, preprocessors, config_user):
    """Get preprocessor settings for for a set of models."""
    all_settings = {}
    for variable in variables:
        # Start out with default settings
        settings = _get_default_settings(variable, config_user)
        # Get preprocessor settings from profile and apply those
        name = variable.get('preprocessor')
        if name not in preprocessors:
            raise NamelistError("Unknown preprocessor {} in variable {}"
                                .format(name, variable['short_name']))
        profile_settings = preprocessors[variable['preprocessor']]
        _apply_preprocessor_settings(settings, profile_settings)
        # if the target grid is a model name, replace it with a file name
        _update_target_grid(
            variable=variable,
            all_variables=variables,
            settings=settings,
            config_user=config_user)
        if 'derive' in settings:
            # create two pairs of settings?
            pass
        check_preprocessor_settings(settings)
        all_settings[variable['filename']] = settings

    _check_multi_model_settings(all_settings)
    return all_settings


def _check_multi_model_settings(all_settings):
    """Check that multi model settings are identical for all models."""
    multi_model_steps = (step for step in MULTI_MODEL_FUNCTIONS
                         if any(step in settings for settings in all_settings))
    for step in multi_model_steps:
        result = None
        for settings in all_settings:
            if result is None:
                result = settings[step]
            elif result != settings[step]:
                raise NamelistError(
                    "Unable to combine differing multi-model settings "
                    "{} and {}".format(result, settings[step]))


def _get_preprocessor_task(variables, preprocessors, config_user):
    """Create preprocessor tasks for a set of models."""
    # Configure preprocessor
    all_settings = _get_preprocessor_settings(
        variables=variables,
        preprocessors=preprocessors,
        config_user=config_user)

    # Input data, used by tasks without ancestors
    input_files = {
        v['filename']: _get_input_files(v, config_user)
        for v in variables
    }
    metadata = {'variables': variables}
    output_dir = os.path.dirname(next(v['filename'] for v in variables))
    task = PreprocessingTask(
        settings=all_settings,
        output_dir=output_dir,
        input_files=input_files,
        metadata=metadata,
        debug=config_user['save_intermediary_cubes'])

    return task


class Namelist(object):
    """Namelist object"""

    def __init__(self,
                 raw_namelist,
                 config_user,
                 initialize_tasks=True,
                 namelist_file=None):
        """Parse a namelist file into an object."""
        self._cfg = config_user
        self._namelist_file = namelist_file  # TODO: remove this dependency
        self._preprocessors = raw_namelist['preprocessors']
        self.models = raw_namelist['models']
        self.diagnostics = self._initialize_diagnostics(
            raw_namelist['diagnostics'])
        self.tasks = self.initialize_tasks() if initialize_tasks else None

    def _initialize_diagnostics(self, raw_diagnostics):
        """Define diagnostics in namelist"""
        logger.debug("Retrieving diagnostics from namelist")

        diagnostics = {}

        for name, raw_diagnostic in raw_diagnostics.items():
            diagnostic = {}
            diagnostic['name'] = name
            models = self._initialize_models(
                name, raw_diagnostic.get('additional_models'))
            diagnostic['models'] = models
            diagnostic['variable_collection'] = \
                self._initialize_variable_collection(
                    name, raw_diagnostic.get('variables'), models)
            diagnostic['scripts'] = self._initialize_scripts(
                name, raw_diagnostic.get('scripts'))
            diagnostics[name] = diagnostic

        return diagnostics

    def _initialize_models(self, diagnostic_name, raw_additional_models):
        """Define models in diagnostic"""
        logger.debug("Setting models for diagnostic %s", diagnostic_name)

        models = list(self.models)

        if raw_additional_models:
            models += raw_additional_models

        check_duplicate_models(models)
        return models

    def _initialize_variables(self, raw_variable, models):
        """Define variables for all models."""
        variables = []
        cmor_keys = ['standard_name', 'long_name', 'units']
        for model in models:
            variable = copy.deepcopy(raw_variable)
            variable.update(model)
            if ('cmor_table' not in variable
                    and variable.get('project') in CMOR_TABLES):
                variable['cmor_table'] = variable['project']
            _add_cmor_info(variable, cmor_keys)
            variables.append(variable)

        for variable in variables:
            for key in ['mip', 'cmor_table'] + cmor_keys:
                if key not in variable:
                    value = _get_value(key, variables)
                    if value is not None:
                        variable[key] = value

        check_variables(variables)

        for variable in variables:
            variable['filename'] = get_output_file(variable,
                                                   self._cfg['preproc_dir'])

        return variables

    def _initialize_variable_collection(self, diagnostic_name, raw_variables,
                                        models):
        """Define variables in diagnostic"""
        if not raw_variables:
            return {}

        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        variable_collection = {}

        for variable_name, raw_variable in raw_variables.items():
            if 'short_name' not in raw_variable:
                raw_variable['short_name'] = variable_name
            raw_variable['diagnostic'] = diagnostic_name
            raw_variable['preprocessor'] = str(raw_variable['preprocessor'])
            variable_collection[variable_name] = \
                self._initialize_variables(raw_variable, models)

        return variable_collection

    def _initialize_scripts(self, diagnostic_name, raw_scripts):
        """Define script in diagnostic"""
        if not raw_scripts:
            return {}

        logger.debug("Setting script for diagnostic %s", diagnostic_name)

        scripts = {}

        for script_name, raw_settings in raw_scripts.items():
            raw_script = raw_settings.pop('script')
            ancestors = []
            for id_glob in raw_settings.pop('ancestors', []):
                if TASKSEP not in id_glob:
                    id_glob = diagnostic_name + TASKSEP + id_glob
                ancestors.append(id_glob)
            settings = copy.deepcopy(raw_settings)
            settings['exit_on_ncl_warning'] = self._cfg['exit_on_warning']

            output_dir = os.path.join(
                self._cfg['output_dir'],
                diagnostic_name,
                script_name,
            )
            scripts[script_name] = {
                'script': raw_script,
                'output_dir': output_dir,
                'settings': settings,
                'ancestors': ancestors,
            }

        return scripts

    def _resolve_diagnostic_ancestors(self, diagnostic_tasks, later):
        """Resolve diagnostic ancestors"""
        for diagnostic_name, diagnostic in self.diagnostics.items():
            for script_name, script_cfg in diagnostic['scripts'].items():
                task_id = diagnostic_name + TASKSEP + script_name
                if diagnostic_tasks[task_id].ancestors is later:
                    logger.debug("Linking tasks for diagnostic %s script %s",
                                 diagnostic_name, script_name)
                    ancestors = []
                    for id_glob in script_cfg['ancestors']:
                        ancestor_ids = fnmatch.filter(diagnostic_tasks,
                                                      id_glob)
                        if not ancestor_ids:
                            raise NamelistError(
                                "Could not find any ancestors matching {}"
                                .format(id_glob))
                        logger.debug("Pattern %s matches %s", id_glob,
                                     ancestor_ids)
                        ancestors.extend(diagnostic_tasks[ancestor_id]
                                         for ancestor_id in ancestor_ids)
                    diagnostic_tasks[task_id].ancestors = ancestors

    def initialize_tasks(self):
        """Define tasks in namelist"""
        logger.info("Creating tasks from namelist")
        tasks = []

        diagnostic_tasks = {}
        later = object()
        for diagnostic_name, diagnostic in self.diagnostics.items():
            logger.info("Creating tasks for diagnostic %s", diagnostic_name)

            # Create preprocessor tasks
            preproc_tasks = []
            for variable_name in diagnostic['variable_collection']:
                logger.info("Creating preprocessor tasks for variable %s",
                            variable_name)
                task = _get_preprocessor_task(
                    variables=diagnostic['variable_collection'][variable_name],
                    preprocessors=self._preprocessors,
                    config_user=self._cfg)
                preproc_tasks.append(task)
            tasks.extend(preproc_tasks)

            # Create diagnostic tasks
            for script_name, script_cfg in diagnostic['scripts'].items():
                task_id = diagnostic_name + TASKSEP + script_name
                logger.info("Creating diagnostic task %s", task_id)
                # Ancestors will be resolved after creating all tasks,
                # unless there aren't any, in which case we set them to the
                # preprocessor tasks.
                ancestors = later if script_cfg['ancestors'] else preproc_tasks
                task = DiagnosticTask(
                    script=script_cfg['script'],
                    output_dir=script_cfg['output_dir'],
                    settings=script_cfg['settings'],
                    ancestors=ancestors,
                )
                diagnostic_tasks[task_id] = task
                # TODO: remove code below once new interface implemented
                os.makedirs(task.output_dir)
                write_legacy_ncl_interface(
                    variables=diagnostic['variable_collection'],
                    settings=task.settings,
                    config_user=self._cfg,
                    output_dir=task.output_dir,
                    namelist_file=self._namelist_file,
                    script=task.script)
                task.settings['env'] = get_legacy_ncl_env(
                    config_user=self._cfg,
                    output_dir=task.output_dir,
                    namelist_basename=os.path.basename(self._namelist_file))

        # Resolve diagnostic ancestors marked as 'later'
        self._resolve_diagnostic_ancestors(diagnostic_tasks, later)

        # TODO: check that no loops are created (will throw RecursionError)

        # Only add diagnostic tasks if enabled
        if self._cfg['run_diagnostic']:
            tasks.extend(diagnostic_tasks.values())

        # Return smallest possible set of tasks
        tasks = get_independent_tasks(tasks)

        return tasks

    def __str__(self):
        """Get human readable summary."""
        return '\n\n'.join(str(task) for task in self.tasks)

    def run(self):
        """Run all tasks in the namelist."""
        run_tasks(
            self.tasks, max_parallel_tasks=self._cfg['max_parallel_tasks'])
