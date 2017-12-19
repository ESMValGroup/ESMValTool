"""Namelist parser"""
from __future__ import print_function

import copy
import logging
import os

import yaml

from .interface_scripts.data_finder import (
    get_input_filelist, get_input_filename, get_output_file,
    get_start_end_year)
from .interface_scripts.data_interface import (get_legacy_ncl_env,
                                               write_legacy_ncl_interface)
from .interface_scripts.preprocessing_tools import merge_callback
from .preprocessor import (DEFAULT_ORDER, PreprocessingTask,
                           select_multi_model_settings,
                           select_single_model_settings)
from .preprocessor._download import synda_search
from .preprocessor._reformat import CMOR_TABLES
from .task import DiagnosticTask

logger = logging.getLogger(__name__)


class NamelistError(Exception):
    """Namelist contains an error."""


def read_namelist_file(filename, config_user, initialize_tasks=True):
    """Read a namelist from file."""
    raw_namelist = check_namelist(filename)
    return Namelist(
        raw_namelist, config_user, initialize_tasks, namelist_file=filename)


def check_namelist(filename):
    """Check a namelist file and return it in raw form."""
    # TODO: use yaml schema for checking basic properties
    with open(filename, 'r') as file:
        raw_namelist = yaml.safe_load(file)

    # TODO: add more checks
    models = []
    for model in raw_namelist['models']:
        if model in models:
            raise NamelistError(
                "Duplicate model {} in models section".format(model))
        models.append(model)

    for name, settings in raw_namelist['preprocessors'].items():
        for step in settings:
            if step not in DEFAULT_ORDER:
                raise NamelistError(
                    "Unknown function {} in preprocessor {}, choose from: {}"
                    .format(step, name, ', '.join(DEFAULT_ORDER)))
            # TODO: check for correct keyword arguments using inspect module

    return raw_namelist


def _models_match(full_model, short_model):
    """Check if short_model description matches full_model."""
    if isinstance(short_model, str):
        return short_model in full_model.values()

    if isinstance(short_model, dict):
        for key, value in short_model.items():
            if key not in full_model or value != full_model[key]:
                return False
        return True

    raise ValueError("Cannot decide if short_model {} describes "
                     "full_model {}".format(short_model, full_model))


def _find_model(full_models, short_model, warn=False):

    matches = [m for m in full_models if _models_match(m, short_model)]
    if len(matches) == 1:
        return matches[0]

    if warn:
        if not matches:
            logger.warning("No model matching description %s found.",
                           short_model)
        else:
            logger.warning("Found multiple models %s matching description %s.",
                           matches, short_model)


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


def get_single_model_settings(variable):
    """Get default preprocessor settings."""
    cfg = {}

    # Configure loading
    cfg['load'] = {
        'mode': 'concatenate',
        'constraints': variable['standard_name'],
        'callback': merge_callback,
    }

    # Configure fixes
    fixcfg = {
        'project': variable['project'],
        'model': variable['model'],
        'short_name': variable['short_name'],
    }
    cfg['fix_file'] = dict(fixcfg)
    # Only supply mip if the CMOR check fixes are implemented.
    if variable['project'] in CMOR_TABLES:
        fixcfg['mip'] = variable['mip']
    cfg['fix_metadata'] = dict(fixcfg)
    cfg['fix_data'] = dict(fixcfg)

    # Configure time extraction
    cfg['extract_time'] = {
        'yr1': variable['start_year'],
        'yr2': variable['end_year'] + 1,
        'mo1': 1,
        'mo2': 1,
        'd1': 1,
        'd2': 1,
    }

    return cfg


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


def get_single_model_task(variable, settings, user_config):
    """Get a preprocessor task for a single model"""
    logger.info("Configuring single-model task for variable %s model %s",
                variable['short_name'], variable['model'])
    settings = select_single_model_settings(settings)
    # TODO: implement variable derivation task

    # Get default preprocessor configuration.
    cfg = get_single_model_settings(variable)

    rootpath = user_config['rootpath']
    drs = user_config['drs']

    # Set up input files.
    input_files = get_input_filelist(
        variable=variable, rootpath=rootpath, drs=drs)

    # Set up downloading using synda if requested.
    # Do not download if files are already available locally.
    if not input_files and user_config['synda_download']:
        input_files = synda_search(variable)
        local_dir = os.path.dirname(
            get_input_filename(variable=variable, rootpath=rootpath, drs=drs))
        cfg['download'] = {
            'dest_folder': local_dir,
        }

    check_data_availability(input_files, variable['start_year'],
                            variable['end_year'])

    logger.info("Using input files:\n%s", '\n'.join(input_files))

    # Configure saving to output files
    output_file = get_output_file(
        variable=variable, preproc_dir=user_config['preproc_dir'])
    logger.info("Output will be written to:\n%s", output_file)

    cfg['save'] = {
        'target': output_file,
    }

    # Use settings from preprocessor profile
    for step, args in settings.items():
        # Remove disabled preprocessor functions
        if args is False:
            if step in cfg:
                del cfg[step]
            continue
        # Enable/update functions without keywords
        if step not in cfg:
            cfg[step] = {}
        if isinstance(args, dict):
            cfg[step].update(args)

    if user_config['save_intermediary_cubes']:
        debug = {'paths': [os.path.splitext(output_file)[0]]}
    else:
        debug = None
    return PreprocessingTask(settings=cfg, input_data=input_files, debug=debug)


def get_multi_model_task(standard_name, settings, ancestors, user_config):
    """Get multi-model preprocessor tasks."""
    logger.info("Configuring multi-model task")
    settings = select_multi_model_settings(settings)

    # Remove disabled steps
    settings = {
        step: settings[step]
        for step in settings if settings[step] is not False
    }
    # Provide default keyword arguments for enabled steps
    for step in settings:
        if settings[step] is True:
            settings[step] = {}

    # Configure loading
    settings['load'] = {'constraint': standard_name, 'mode': 'ordered'}

    # Configure saving
    files = [task.settings['save']['target'] for task in ancestors]
    settings['save'] = {'targets': files}

    if user_config['save_intermediary_cubes']:
        debug = {'paths': [os.path.splitext(f)[0] for f in files]}
    else:
        debug = None
    task = PreprocessingTask(
        settings=settings, ancestors=ancestors, debug=debug)

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
                name, raw_diagnostic['additional_models'])
            diagnostic['models'] = models
            diagnostic['variables'] = self._initialize_variables(
                name, raw_diagnostic['variables'], models)
            diagnostic['scripts'] = self._initialize_scripts(
                name, raw_diagnostic['scripts'])
            diagnostics[name] = diagnostic

        return diagnostics

    def _initialize_models(self, diagnostic_name, raw_additional_models):
        """Define models in diagnostic"""
        logger.debug("Setting models for diagnostic %s", diagnostic_name)

        models = self.models + raw_additional_models

        return models

    @staticmethod
    def _initialize_variables(diagnostic_name, raw_variables, models):
        """Define variables in diagnostic"""
        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        variables = {}
        cmor_keys = ['standard_name', 'long_name', 'units']
        for variable_name, raw_variable in raw_variables.items():
            variables[variable_name] = []
            for model in models:
                variable = copy.deepcopy(raw_variable)
                variable.update(model)
                if 'short_name' not in variable:
                    variable['short_name'] = variable_name
                _add_cmor_info(variable, cmor_keys)
                variables[variable_name].append(variable)

            for variable in variables[variable_name]:
                for key in ['mip'] + cmor_keys:
                    if key not in variable:
                        value = _get_value(key, variables[variable_name])
                        if value is not None:
                            variable[key] = value

        return variables

    def _initialize_scripts(self, diagnostic_name, raw_scripts):
        """Define script in diagnostic"""
        logger.debug("Setting script for diagnostic %s", diagnostic_name)

        scripts = {}

        # Dummy diagnostic for running only preprocessors
        if raw_scripts == 'None':
            return scripts

        for script_name, raw_settings in raw_scripts.items():
            raw_script = raw_settings.pop('script')
            ancestors = raw_settings.pop('ancestors', None)
            settings = copy.deepcopy(raw_settings)
            # Add output dir to settings
            settings['output_dir'] = os.path.join(
                self._cfg['output_dir'],
                diagnostic_name,
                script_name,
            )
            settings['exit_on_ncl_warning'] = self._cfg['exit_on_warning']

            scripts[script_name] = {
                'script': raw_script,
                'settings': settings,
                'ancestors': ancestors,
            }

        return scripts

    def _update_target_grid(self, variable, settings, all_models):
        """Configure regrid target_grid in preprocessor settings."""

        def issubset(sub, ref):
            """Check if dict sub is a subset of dict ref."""
            return all(ref.get(key, object()) == sub[key] for key in sub)

        if ('regrid' in settings and settings['regrid'] is not False
                and 'target_grid' in settings['regrid']):
            # if the target grid is a model, replace with filename
            target_grid_model = _find_model(
                full_models=all_models,
                short_model=settings['regrid']['target_grid'])
            if target_grid_model:
                if issubset(target_grid_model, variable):
                    # No need to regrid model onto itself
                    del settings['regrid']
                else:
                    target_variable = dict(variable)
                    target_variable.update(target_grid_model)
                    files = get_input_filelist(
                        variable=target_variable,
                        rootpath=self._cfg['rootpath'],
                        drs=self._cfg['drs'])
                    settings['regrid']['target_grid'] = files[0]

    def initialize_tasks(self):
        """Define tasks in namelist"""
        logger.info("Creating tasks from namelist")

        tasks = []

        all_preproc_tasks = {}
        for diagnostic_name, diagnostic in self.diagnostics.items():
            logger.info("Creating tasks for diagnostic %s", diagnostic_name)

            # Create preprocessor tasks
            for variable_name in diagnostic['variables']:
                # Create single model tasks
                single_model_tasks = []
                for variable in diagnostic['variables'][variable_name]:
                    preproc_id = variable['preprocessor']
                    settings = copy.deepcopy(self._preprocessors[preproc_id])
                    self._update_target_grid(variable, settings,
                                             diagnostic['models'])

                    task_id = '_'.join(
                        str(variable[key]) for key in sorted(variable))
                    if task_id not in all_preproc_tasks:
                        task = get_single_model_task(
                            variable=variable,
                            settings=settings,
                            user_config=self._cfg,
                        )
                        all_preproc_tasks[task_id] = task
                    single_model_tasks.append(all_preproc_tasks[task_id])

                # Create multi model task
                task_id = diagnostic_name + '_' + variable_name
                if task_id not in all_preproc_tasks:
                    task = get_multi_model_task(
                        standard_name=variable['standard_name'],
                        settings=settings,
                        ancestors=list(single_model_tasks),
                        user_config=self._cfg,
                    )
                    all_preproc_tasks[task_id] = task
                multi_model_task = all_preproc_tasks[task_id]

            # Create diagnostic tasks
            diagnostic_tasks = {}
            for task_id, script_cfg in diagnostic['scripts'].items():
                logger.info("Creating diagnostic task %s", task_id)
                diagnostic_task = DiagnosticTask(
                    script=script_cfg['script'],
                    settings=script_cfg['settings'],
                )
                diagnostic_tasks[task_id] = diagnostic_task
                # TODO: remove code below once new interface implemented
                os.makedirs(diagnostic_task.settings['output_dir'])
                write_legacy_ncl_interface(
                    variables=diagnostic['variables'],
                    settings=diagnostic_task.settings,
                    config_user=self._cfg,
                    output_dir=diagnostic_task.settings['output_dir'],
                    namelist_file=self._namelist_file,
                    script=diagnostic_task.script)
                diagnostic_task.settings['env'] = get_legacy_ncl_env(
                    config_user=self._cfg,
                    output_dir=diagnostic_task.settings['output_dir'],
                    namelist_basename=os.path.basename(self._namelist_file))

            # Preprocessor run only
            if not diagnostic_tasks:
                tasks.append(multi_model_task)

            # Add diagnostic tasks
            tasks.extend(diagnostic_tasks.values())

            # Establish relations between diagnostic tasks and minimize `tasks`
            for task_id, script_cfg in diagnostic['scripts'].items():
                if script_cfg['ancestors']:
                    ancestors = []
                    for ancestor_id in script_cfg['ancestors']:
                        diagnostic_task = diagnostic_tasks[ancestor_id]
                        ancestors.append(diagnostic_task)
                        tasks.remove(diagnostic_task)
                else:
                    ancestors = [multi_model_task]
                diagnostic_tasks[task_id].ancestors = ancestors

        return tasks

    def __str__(self):
        """Get human readable summary."""
        return '\n\n'.join(str(task) for task in self.tasks)

    def run(self):
        """Run all tasks in the namelist."""
        for task in self.tasks:
            task.run()
