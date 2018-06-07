"""Namelist parser"""
import copy
import fnmatch
import inspect
import logging
import os
import subprocess

import yamale
import yaml

from . import __version__, preprocessor
from ._data_finder import (get_input_filelist, get_input_filename,
                           get_output_file, get_start_end_year,
                           get_statistic_output_file)
from ._task import DiagnosticTask, get_independent_tasks, run_tasks, which
from .cmor.table import CMOR_TABLES
from .preprocessor._derive import get_required
from .preprocessor._download import synda_search
from .preprocessor._io import MODEL_KEYS, concatenate_callback
from .preprocessor._regrid import get_cmor_levels, get_reference_levels

logger = logging.getLogger(__name__)

TASKSEP = os.sep


class NamelistError(Exception):
    """Namelist contains an error."""


def read_namelist_file(filename, config_user, initialize_tasks=True):
    """Read a namelist from file."""
    raw_namelist = check_namelist(filename)
    return Namelist(
        raw_namelist, config_user, initialize_tasks, namelist_file=filename)


def check_ncl_version():
    """Check the NCL version"""
    ncl = which('ncl')
    if not ncl:
        raise NamelistError("Namelist contains NCL scripts, but cannot find "
                            "an NCL installation.")
    try:
        cmd = [ncl, '-V']
        version = subprocess.check_output(cmd, universal_newlines=True)
    except subprocess.CalledProcessError:
        logger.error("Failed to execute '%s'", ' '.join(' '.join(cmd)))
        raise NamelistError("Namelist contains NCL scripts, but your NCL "
                            "installation appears to be broken.")

    version = version.strip()
    logger.info("Found NCL version %s", version)

    major, minor = (int(i) for i in version.split('.')[:2])
    if major < 6 or (major == 6 and minor < 4):
        raise NamelistError("NCL version 6.4 or higher is required to run "
                            "a namelist containing NCL scripts.")


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
    preprocessor_functions = set(preprocessor.DEFAULT_ORDER)
    for name, settings in preprocessors.items():
        invalid_functions = set(settings) - preprocessor_functions
        if invalid_functions:
            raise NamelistError(
                "Unknown function(s) {} in preprocessor {}, choose from: {}"
                .format(invalid_functions, name, ', '.join(
                    preprocessor.DEFAULT_ORDER)))


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
        if step not in preprocessor.DEFAULT_ORDER:
            raise NamelistError(
                "Unknown preprocessor function '{}', choose from: {}".format(
                    step, ', '.join(preprocessor.DEFAULT_ORDER)))
        function = getattr(preprocessor, step)
        argspec = inspect.getargspec(function)
        args = argspec.args[1:]
        # Check for invalid arguments
        invalid_args = set(settings[step]) - set(args)
        if invalid_args:
            raise NamelistError(
                "Invalid argument(s): {} encountered for preprocessor "
                "function {}".format(', '.join(invalid_args), step))
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


def check_variable(variable, required_keys):
    """Check variables as derived from namelist"""
    required = set(required_keys)
    missing = required - set(variable)
    if missing:
        raise NamelistError(
            "Missing keys {} from variable {} in diagnostic {}".format(
                missing, variable.get('short_name'),
                variable.get('diagnostic')))


def check_data_availability(input_files, variable):
    """Check if the required input data is available"""
    if not input_files:
        raise NamelistError("No input files found for variable {}"
                            .format(variable))

    required_years = set(
        range(variable['start_year'], variable['end_year'] + 1))
    available_years = set()
    for filename in input_files:
        start, end = get_start_end_year(filename)
        available_years.update(range(start, end + 1))

    missing_years = required_years - available_years
    if missing_years:
        raise NamelistError(
            "No input data available for years {} in files {}".format(
                ", ".join(str(year) for year in missing_years), input_files))


def _get_value(key, models):
    """Get a value for key by looking at the other models."""
    values = {model[key] for model in models if key in model}

    if len(values) == 1:
        return values.pop()

    if len(values) > 1:
        raise NamelistError("Ambigous values {} for property {}".format(
            values, key))


def _update_from_others(variable, keys, models):
    """Get values for keys by copying from the other models."""
    for key in keys:
        if key not in variable:
            value = _get_value(key, models)
            if value is not None:
                variable[key] = value


def _update_cmor_table(table, mip, short_name):
    """Try to add an ESMValTool custom CMOR table file."""
    cmor_table = CMOR_TABLES[table]
    var_info = cmor_table.get_variable(mip, short_name)

    if var_info is None and hasattr(cmor_table, 'add_custom_table_file'):
        table_file = os.path.join(
            os.path.dirname(__file__), 'cmor', 'tables', 'custom',
            'CMOR_' + short_name + '.dat')
        if os.path.exists(table_file):
            logger.debug("Loading custom CMOR table from %s", table_file)
            cmor_table.add_custom_table_file(table_file, mip)
            var_info = cmor_table.get_variable(mip, short_name)

    if var_info is None:
        raise NamelistError(
            "Unable to load CMOR table '{}' for variable '{}' with mip '{}'"
            .format(table, short_name, mip))


def _add_cmor_info(variable, override=False):
    """Add information from CMOR tables to variable."""
    logger.debug("If not present: adding keys from CMOR table to %s", variable)

    if 'cmor_table' not in variable or 'mip' not in variable:
        logger.debug("Skipping because cmor_table or mip not specified")
        return

    if variable['cmor_table'] not in CMOR_TABLES:
        logger.warning("Unknown CMOR table %s", variable['cmor_table'])

    # Copy the following keys from CMOR table
    cmor_keys = ['standard_name', 'long_name', 'units']
    table_entry = CMOR_TABLES[variable['cmor_table']].get_variable(
        variable['mip'], variable['short_name'])

    for key in cmor_keys:
        if key not in variable or override and hasattr(table_entry, key):
            value = getattr(table_entry, key)
            if value is not None:
                variable[key] = value

    # Check that keys are available
    check_variable(variable, required_keys=cmor_keys)


def _update_target_levels(variable, variables, settings, config_user):
    """Replace the target levels model name with a filename if needed."""
    if (not settings.get('extract_levels')
            or 'levels' not in settings['extract_levels']):
        return

    levels = settings['extract_levels']['levels']

    # If levels is a model name, replace it by a dict with a 'model' entry
    if any(levels == v['model'] for v in variables):
        settings['extract_levels']['levels'] = {'model': levels}
        levels = settings['extract_levels']['levels']

    if not isinstance(levels, dict):
        return

    if 'cmor_table' in levels and 'coordinate' in levels:
        settings['extract_levels']['levels'] = get_cmor_levels(
            levels['cmor_table'], levels['coordinate'])
    elif 'model' in levels:
        if variable['model'] == levels['model']:
            del settings['extract_levels']
        else:
            filename = _model_to_file(levels['model'], variables, config_user)
            coordinate = levels.get('coordinate', 'air_pressure')
            settings['extract_levels']['levels'] = get_reference_levels(
                filename, coordinate)


def _update_target_grid(variable, variables, settings, config_user):
    """Replace the target grid model name with a filename if needed."""
    if not settings.get('regrid') or 'target_grid' not in settings['regrid']:
        return

    grid = settings['regrid']['target_grid']

    if variable['model'] == grid:
        del settings['regrid']
    elif any(grid == v['model'] for v in variables):
        settings['regrid']['target_grid'] = _model_to_file(
            grid, variables, config_user)


def _model_to_file(model, variables, config_user):
    """Find the first file belonging to model."""
    for variable in variables:
        if variable['model'] == model:
            files = get_input_filelist(
                variable=variable,
                rootpath=config_user['rootpath'],
                drs=config_user['drs'])
            return files[0]

    raise NamelistError(
        "Unable to find matching file for model {}".format(model))


def _get_default_settings(variable, config_user, derive=False):
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
    settings['load_cubes'] = {
        'callback': concatenate_callback,
        'filename': variable['filename'],
        'metadata': variable,
    }
    if not derive:
        settings['load_cubes']['constraints'] = variable['standard_name']
    # Configure merge
    settings['concatenate'] = {}

    # Configure fixes
    fix = {
        'project': variable['project'],
        'model': variable['model'],
        'short_name': variable['short_name'],
    }
    # File fixes
    fix_dir = variable['filename'] + '_fixed'
    if not derive:
        settings['fix_file'] = dict(fix)
        settings['fix_file']['output_dir'] = fix_dir
    # Cube fixes
    # Only supply mip if the CMOR check fixes are implemented.
    if variable.get('cmor_table'):
        fix['cmor_table'] = variable['cmor_table']
        fix['mip'] = variable['mip']
    settings['fix_data'] = dict(fix)
    settings['fix_metadata'] = dict(fix)

    # Configure time extraction
    settings['extract_time'] = {
        'yr1': variable['start_year'],
        'yr2': variable['end_year'] + 1,
        'mo1': 1,
        'mo2': 1,
        'd1': 1,
        'd2': 1,
    }

    # Configure CMOR metadata check
    if variable.get('cmor_table'):
        settings['cmor_check_metadata'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
        }
    # Configure final CMOR data check
    if variable.get('cmor_table'):
        settings['cmor_check_data'] = {
            'cmor_table': variable['cmor_table'],
            'mip': variable['mip'],
            'short_name': variable['short_name'],
        }

    # Clean up fixed files
    if not config_user['save_intermediary_cubes']:
        settings['cleanup'] = {
            'remove': [fix_dir],
        }

    # Configure saving cubes to file
    settings['save_cubes'] = {}

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

    logger.info("Using input files for variable %s of model %s:\n%s",
                variable['short_name'], variable['model'],
                '\n'.join(input_files))
    check_data_availability(input_files, variable)

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


def _update_multi_model_statistics(variables, settings, preproc_dir):
    """Configure multi model statistics."""
    if settings.get('multi_model_statistics', False):
        if settings['multi_model_statistics'] is True:
            settings['multi_model_statistics'] = {}
        stat_settings = settings['multi_model_statistics']

        variable = variables[0]

        # Define output files
        stat_settings['filenames'] = {}
        for statistic in stat_settings['statistics']:
            stat_settings['filenames'][statistic] = get_statistic_output_file(
                variable, statistic, preproc_dir)

        # Define models to exclude
        exclude_models = set(stat_settings.get('exclude', {}))
        for key in 'reference_model', 'alternative_model':
            if key in exclude_models and key in variable:
                exclude_models.remove(key)
                exclude_models.add(variable[key])
        exclude_files = {
            v['filename']
            for v in variables if v['model'] in exclude_models
        }
        logger.debug('Multimodel excludes files %s', exclude_files)
        stat_settings['exclude'] = {'_filename': exclude_files}


def _get_preprocessor_settings(variables, profile, config_user):
    """Get preprocessor settings for a set of models."""
    all_settings = {}
    profile = copy.deepcopy(profile)
    _update_multi_model_statistics(variables, profile,
                                   config_user['preproc_dir'])

    for variable in variables:
        derive = 'derive' in profile
        settings = _get_default_settings(variable, config_user, derive=derive)
        _apply_preprocessor_settings(settings, profile)
        # if the target grid is a model name, replace it with a file name
        # TODO: call _update_target_grid only once per variable?
        _update_target_levels(
            variable=variable,
            variables=variables,
            settings=settings,
            config_user=config_user)
        _update_target_grid(
            variable=variable,
            variables=variables,
            settings=settings,
            config_user=config_user)
        check_preprocessor_settings(settings)
        all_settings[variable['filename']] = settings

    _check_multi_model_settings(all_settings)
    return all_settings


def _check_multi_model_settings(all_settings):
    """Check that multi model settings are identical for all models."""
    multi_model_steps = (step for step in preprocessor.MULTI_MODEL_FUNCTIONS
                         if any(step in settings
                                for settings in all_settings.values()))
    for step in multi_model_steps:
        result = None
        for filename, settings in all_settings.items():
            if result is None:
                result = settings[step]
            elif result != settings[step]:
                raise NamelistError(
                    "Unable to combine differing multi-model settings "
                    "{} and {} for output file {}".format(
                        result, settings[step], filename))


def _get_single_preprocessor_task(variables,
                                  profile,
                                  config_user,
                                  ancestors=None):
    """Create preprocessor tasks for a set of models."""
    # Configure preprocessor
    all_settings = _get_preprocessor_settings(
        variables=variables, profile=profile, config_user=config_user)

    # Input files, used by tasks without ancestors
    input_files = None
    if not ancestors:
        input_files = [
            filename for variable in variables
            for filename in _get_input_files(variable, config_user)
        ]

    output_dir = os.path.dirname(variables[0]['filename'])

    task = preprocessor.PreprocessingTask(
        settings=all_settings,
        output_dir=output_dir,
        ancestors=ancestors,
        input_files=input_files,
        debug=config_user['save_intermediary_cubes'])

    return task


def _get_preprocessor_task(variables,
                           profiles,
                           config_user,
                           write_ncl_interface=False):
    """Create preprocessor task(s) for a set of models."""
    # First set up the preprocessor profile
    variable = variables[0]
    preproc_name = variable.get('preprocessor')
    if preproc_name not in profiles:
        raise NamelistError(
            "Unknown preprocessor {} in variable {} of diagnostic {}".format(
                preproc_name, variable['short_name'], variable['diagnostic']))
    profile = copy.deepcopy(profiles[variable['preprocessor']])
    logger.info("Creating preprocessor '%s' task for variable '%s'",
                variable['preprocessor'], variable['short_name'])

    # Create preprocessor task(s)
    derive_tasks = []
    if variable.get('derive'):
        # Create tasks to prepare the input data for the derive step
        derive_profile, profile = preprocessor.split_settings(
            profile, 'derive')
        profile['derive'] = {'short_name': variable['short_name']}

        derive_input = {}
        for variable in variables:
            _update_cmor_table(
                table=variable['cmor_table'],
                mip=variable['mip'],
                short_name=variable['short_name'])
            _add_cmor_info(variable)
            if not variable.get('force_derivation') and get_input_filelist(
                    variable=variable,
                    rootpath=config_user['rootpath'],
                    drs=config_user['drs']):
                # No need to derive, just process normally up to derive step
                short_name = variable['short_name']
                if short_name not in derive_input:
                    derive_input[short_name] = []
                derive_input[short_name].append(variable)
            else:
                # Process input data needed to derive variable
                for short_name, field in get_required(variable['short_name'],
                                                      variable['field']):
                    if short_name not in derive_input:
                        derive_input[short_name] = []
                    variable = copy.deepcopy(variable)
                    variable['short_name'] = short_name
                    variable['field'] = field
                    variable['filename'] = get_output_file(
                        variable, config_user['preproc_dir'])
                    _add_cmor_info(variable, override=True)
                    derive_input[short_name].append(variable)

        for derive_variables in derive_input.values():
            task = _get_single_preprocessor_task(derive_variables,
                                                 derive_profile, config_user)
            derive_tasks.append(task)

    # Add CMOR info
    for variable in variables:
        _add_cmor_info(variable)

    # Create (final) preprocessor task
    profile['extract_metadata'] = {'write_ncl': write_ncl_interface}
    task = _get_single_preprocessor_task(
        variables, profile, config_user, ancestors=derive_tasks)

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
        self._namelist_file = os.path.basename(namelist_file)
        self._preprocessors = raw_namelist['preprocessors']
        if raw_namelist.get('models'):
            self.models = raw_namelist['models']
        else:
            self.models = []
        self._support_ncl = self._need_ncl(raw_namelist['diagnostics'])
        self.diagnostics = self._initialize_diagnostics(
            raw_namelist['diagnostics'])
        self.tasks = self.initialize_tasks() if initialize_tasks else None

    @staticmethod
    def _need_ncl(raw_diagnostics):
        if not raw_diagnostics:
            return False
        for diagnostic in raw_diagnostics.values():
            if not diagnostic.get('scripts'):
                continue
            for script in diagnostic['scripts'].values():
                if script.get('script', '').lower().endswith('.ncl'):
                    logger.info("NCL script detected, checking NCL version")
                    check_ncl_version()
                    return True
        return False

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
            diagnostic['preprocessor_output'] = \
                self._initialize_preprocessor_output(
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

        for model in models:
            for key in model:
                MODEL_KEYS.add(key)

        check_duplicate_models(models)
        return models

    def _initialize_variables(self, raw_variable, models):
        """Define variables for all models."""
        # TODO: rename `variables` to `attributes` and store in dict
        # using filenames as keys?
        variables = []

        for model in models:
            variable = copy.deepcopy(raw_variable)
            variable.update(model)
            if ('cmor_table' not in variable
                    and variable.get('project') in CMOR_TABLES):
                variable['cmor_table'] = variable['project']
            variables.append(variable)

        required_keys = {
            'short_name', 'field', 'model', 'project', 'start_year',
            'end_year', 'preprocessor', 'diagnostic'
        }

        for variable in variables:
            _update_from_others(variable, ['cmor_table', 'mip'], models)
            check_variable(variable, required_keys)
            variable['filename'] = get_output_file(variable,
                                                   self._cfg['preproc_dir'])

        return variables

    def _initialize_preprocessor_output(self, diagnostic_name, raw_variables,
                                        models):
        """Define variables in diagnostic"""
        if not raw_variables:
            return {}

        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        preprocessor_output = {}

        for variable_name, raw_variable in raw_variables.items():
            if 'short_name' not in raw_variable:
                raw_variable['short_name'] = variable_name
            raw_variable['diagnostic'] = diagnostic_name
            raw_variable['preprocessor'] = str(raw_variable['preprocessor'])
            preprocessor_output[variable_name] = \
                self._initialize_variables(raw_variable, models)

        return preprocessor_output

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
            settings['namelist'] = self._namelist_file
            settings['version'] = __version__
            settings['script'] = script_name
            # Add output dirs to settings
            for dir_name in ('run_dir', 'plot_dir', 'work_dir'):
                settings[dir_name] = os.path.join(self._cfg[dir_name],
                                                  diagnostic_name, script_name)
            # Copy other settings
            if self._support_ncl:
                settings['exit_on_ncl_warning'] = self._cfg['exit_on_warning']
            for key in ('max_data_filesize', 'output_file_type', 'log_level',
                        'write_plots', 'write_netcdf'):
                settings[key] = self._cfg[key]

            scripts[script_name] = {
                'script': raw_script,
                'output_dir': settings['work_dir'],
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
            for variable_name in diagnostic['preprocessor_output']:
                logger.info("Creating preprocessor tasks for variable %s",
                            variable_name)
                task = _get_preprocessor_task(
                    variables=diagnostic['preprocessor_output'][variable_name],
                    profiles=self._preprocessors,
                    config_user=self._cfg,
                    write_ncl_interface=self._support_ncl)
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
