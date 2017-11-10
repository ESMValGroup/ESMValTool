"""Namelist parser"""
from __future__ import print_function

import copy
import logging
import pprint

import yaml

from .interface_scripts.data_finder import get_input_filelist

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


def configure_preprocessor(user_config, preprocessor_profile, all_models,
                           model, variable):
    """Write preprocessor configuration"""

    # Preprocessor configuration
    pre = {}

    # Configure cube loading
    files = get_input_filelist(
        project_info={'GLOBAL': user_config}, model=model, var=variable)
    if not files:
        raise NamelistError("No input files found for model {} "
                            "variable {}".format(model, variable))
    pre['load'] = {
        'uris': files,
        'contstraint': variable['short_name'],
    }

    # Configure fixes
    cfg = {
        'project': model['project'],
        'model': model['name'],
        'short_name': variable['short_name'],
    }
    pre['fix_file'] = dict(cfg)
    if 'mip' in model:
        cfg['mip'] = model['mip']
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
    pre.update(copy.deepcopy(preprocessor_profile))

    # Remove disabled preprocessor functions
    for function, args in pre.items():
        if args is False:
            del pre[function]

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

    return pre


class Namelist(object):
    """Namelist object"""

    def __init__(self, raw_namelist, config_user):
        """Parse a namelist file into an object."""

        self._cfg = config_user
        self._preprocessors = copy.deepcopy(raw_namelist['preprocessors'])
        self._models = copy.deepcopy(raw_namelist['models'])

        self.diagnostics = {}
        for name, raw_diagnostic in raw_namelist['diagnostics'].items():
            diag = self.diagnostics[name] = {}
            diag['name'] = name
            diag['script'] = raw_diagnostic['script']
            diag['settings'] = copy.deepcopy(raw_diagnostic['settings'])
            diag['models'] = self._initialize_models(name,
                                                     raw_diagnostic['models'])
            diag['variables'] = self._initialize_variables(
                name, raw_diagnostic['variables'])

    def _initialize_models(self, name, raw_models):
        """Add models to diagnostic"""

        logger.debug("Resolving models for diagnostic %s", name)

        result = []
        for short_model in raw_models:
            model = _find_model(
                full_models=self._models, short_model=short_model, warn=True)
            if not model:
                raise NamelistError("Unable to find model matching {} in "
                                    "diagnostic {}".format(short_model, name))
            result.append(model)

        return result

    def _initialize_variables(self, name, raw_variables):

        logger.debug("Populating list of variables for diagnostic %s", name)

        result = []
        for raw_variable in raw_variables:
            variable = copy.deepcopy(raw_variable)

            # Preprocessor settings
            variable['preprocessor'] = {}
            for model in self.diagnostics[name]['models']:
                # Add preprocessor settings for each model
                modelkey = '_'.join(str(model[key]) for key in sorted(model))
                variable['preprocessor'][modelkey] = configure_preprocessor(
                    user_config=self._cfg,
                    preprocessor_profile=self._preprocessors[raw_variable[
                        'preprocessor']],
                    all_models=self.diagnostics[name]['models'],
                    model=model,
                    variable=raw_variable,
                )

            result.append(variable)

        return result

    def __str__(self):
        return pprint.pformat(self.diagnostics, indent=2)
