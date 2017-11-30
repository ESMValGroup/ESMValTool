"""Namelist parser"""
from __future__ import print_function

import copy
import logging
import os

import yaml

from .interface_scripts.data_finder import get_input_filelist, get_output_file
from .interface_scripts.preprocessing_tools import merge_callback
from .preprocessor.reformat import CMOR_TABLES
from .preprocessor.run import DEFAULT_ORDER, PreprocessingTask
from .task import DiagnosticTask

logger = logging.getLogger(__name__)


class NamelistError(Exception):
    """Namelist contains an error."""


def read_namelist_file(filename, config_user):
    """Read a namelist from file."""
    raw_namelist = check_namelist(filename)
    return Namelist(raw_namelist, config_user)


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
            # TODO: check for correct keyword arguments?

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


def _get_mip(variable, model, all_models):
    """Try to get mip.

    First look in variable, then in model, then in all models.
    Raises a NamelistError if it has to resort to all models,
    but the mip is not the same for all models.
    """
    mip = None
    if 'mip' in model:
        mip = model['mip']
    elif 'mip' in variable:
        mip = variable['mip']
    else:
        for _model in all_models:
            if 'mip' in _model:
                if mip is None:
                    mip = _model['mip']
                elif mip != _model['mip']:
                    raise NamelistError("Ambigous mip for variable {}"
                                        .format(variable))
    return mip


def _get_standard_name(variable, model, mip=None):
    """Get standard_name.

    First look in variable, then try to look it up in CMOR tables.
    """
    standard_name = None
    if 'standard_name' in variable:
        standard_name = variable['standard_name']
    elif mip and model['project'] in CMOR_TABLES:
        variable_info = CMOR_TABLES[model['project']].get_variable(
            mip, variable['short_name'])
        if variable_info.standard_name:
            standard_name = variable_info.standard_name

    return standard_name


def get_preprocessor_task(settings, all_models, model, variable, user_config):
    """Get a preprocessor task for a single model"""
    # Preprocessor configuration
    pre = {}

    # Configure loading and saving
    input_files = get_input_filelist(
        project_info={'GLOBAL': user_config}, model=model, var=variable)
    if not input_files:
        raise NamelistError("No input files found for model {} "
                            "variable {}".format(model, variable))
    output_file = get_output_file(
        project_info={'GLOBAL': user_config}, model=model, var=variable)

    pre['load'] = {
        'uris': input_files,
        'callback': merge_callback,
    }
    pre['save'] = {
        'target': output_file,
    }

    # Find mip for CMOR fixes and determining standard_name.
    mip = _get_mip(variable, model, all_models)

    # If possible, constrain loading to a cube with required standard_name.
    standard_name = _get_standard_name(variable, model, mip)
    if standard_name:
        pre['load']['constraint'] = standard_name

    # Configure fixes
    cfg = {
        'project': model['project'],
        'model': model['name'],
        'short_name': variable['short_name'],
    }
    pre['fix_file'] = dict(cfg)
    # Only supply mip if the CMOR check fixes are implemented.
    if model['project'] in CMOR_TABLES:
        cfg['mip'] = mip
    pre['fix_metadata'] = dict(cfg)
    pre['fix_data'] = dict(cfg)

    # Configure time extraction
    pre['extract_time'] = {
        'yr1': model['start_year'],
        'yr2': model['end_year'] + 1,
        'mo1': 1,
        'mo2': 1,
        'd1': 1,
        'd2': 1,
    }

    # Override with settings from preprocessor profile
    pre.update(copy.deepcopy(settings))

    # Remove disabled preprocessor functions
    pre = {step: args for step, args in pre.items() if args is not False}

    # Configure regrid if enabled
    if 'regrid' in pre and 'target_grid' in pre['regrid']:
        # if the target grid is a model, replace with filename
        target_grid_model = _find_model(
            full_models=all_models, short_model=pre['regrid']['target_grid'])
        if target_grid_model:
            if target_grid_model == model:
                del pre['regrid']
            else:
                files = get_input_filelist(
                    project_info={'GLOBAL': user_config},
                    model=target_grid_model,
                    var=variable)
                pre['regrid']['target_grid'] = files[0]

    return PreprocessingTask(settings=pre)


class Namelist(object):
    """Namelist object"""

    def __init__(self, raw_namelist, config_user):
        """Parse a namelist file into an object."""
        self._cfg = config_user
        self._preprocessors = copy.deepcopy(raw_namelist['preprocessors'])
        self._models = copy.deepcopy(raw_namelist['models'])
        self._diagnostics = self._initialize_diagnostics(
            raw_namelist['diagnostics'])
        self.tasks = self._initialize_tasks()

    def _initialize_diagnostics(self, raw_diagnostics):
        """Define diagnostics in namelist"""
        logger.debug("Retrieving diagnostics from namelist")

        diagnostics = {}

        for name, raw_diagnostic in raw_diagnostics.items():
            diagnostic = {}
            diagnostic['name'] = name
            diagnostic['models'] = self._initialize_models(
                name, raw_diagnostic['models'])
            diagnostic['variables'] = self._initialize_variables(
                name, raw_diagnostic['variables'])
            diagnostic['script'] = self._initialize_script(
                name, raw_diagnostic['script'])
            diagnostic['settings'] = copy.deepcopy(raw_diagnostic['settings'])
            diagnostics[name] = diagnostic

        return diagnostics

    def _initialize_models(self, diagnostic_name, raw_models):
        """Define models in diagnostic"""
        logger.debug("Resolving models for diagnostic %s", diagnostic_name)

        models = []
        for short_model in raw_models:
            model = _find_model(
                full_models=self._models, short_model=short_model, warn=True)
            if not model:
                raise NamelistError("Unable to find model matching {} in "
                                    "diagnostic {}".format(
                                        short_model, diagnostic_name))
            models.append(model)

        return models

    @staticmethod
    def _initialize_variables(diagnostic_name, raw_variables):
        """Define variables in diagnostic"""
        logger.debug("Populating list of variables for diagnostic %s",
                     diagnostic_name)

        variables = []
        for variable_name, raw_variable in raw_variables.items():
            variable = copy.deepcopy(raw_variable)
            if 'short_name' not in variable:
                variable['short_name'] = variable_name
            variables.append(variable)

        return variables

    @staticmethod
    def _initialize_script(diagnostic_name, raw_script):
        """Define script in diagnostic"""
        logger.debug("Setting script for diagnostic %s", diagnostic_name)

        # Dummy diagnostic for running only preprocessors
        if raw_script == 'None':
            return

        diagnostics_root = os.path.join(
            os.path.dirname(__file__), 'diag_scripts')
        script_file = os.path.abspath(
            os.path.join(diagnostics_root, raw_script))

        bad_script = NamelistError(
            "Cannot execute script {} ({}) of diagnostic {}".format(
                raw_script, script_file, diagnostic_name))

        if not os.path.isfile(script_file):
            raise bad_script

        script = []
        if not os.access(script_file, os.X_OK):  # if not executable
            extension = os.path.splitext(raw_script)[1].lower()[1:]
            executables = {
                'py': ['python'],
                'ncl': ['ncl'],
                'r': ['Rscript', '--slave', '--quiet'],
            }
            if extension not in executables:
                raise bad_script
            script = executables[extension]
        script.append(script_file)

        return script

    def _initialize_tasks(self):
        """Define tasks in namelist"""
        logger.debug("Creating tasks from namelist")

        tasks = []

        all_preproc_tasks = {}
        for diagnostic_name, diagnostic in self._diagnostics.items():
            logger.debug("Creating tasks for diagnostic %s", diagnostic_name)

            # Create preprocessor tasks
            preproc_tasks = []
            for variable in diagnostic['variables']:
                for model in diagnostic['models']:
                    task_id = '_'.join(
                        str(item[key])
                        for item in (variable, model) for key in sorted(item))
                    if task_id not in all_preproc_tasks:
                        preproc_id = variable['preprocessor']
                        task = get_preprocessor_task(
                            settings=self._preprocessors[preproc_id],
                            all_models=diagnostic['models'],
                            model=model,
                            variable=variable,
                            user_config=self._cfg,
                        )
                        all_preproc_tasks[task_id] = task
                    preproc_tasks.append(all_preproc_tasks[task_id])

            # Create diagnostic task
            task = DiagnosticTask(
                # TODO: enable once DiagnosticTask.run() implemented
                script=None,  # diagnostic['script'],
                settings=diagnostic['settings'],
                ancestors=preproc_tasks,
            )
            tasks.append(task)

        return tasks

    def __str__(self):
        """Get human readable summary."""
        return '\n\n'.join(str(task) for task in self.tasks)

    def run(self):
        """Run all tasks in the namelist."""
        for task in self.tasks:
            task.run()
