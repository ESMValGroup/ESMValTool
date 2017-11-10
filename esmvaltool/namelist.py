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

    with open(filename, 'r') as file:
        raw_namelist = yaml.safe_load(file)

    check_namelist(raw_namelist)

    return Namelist(raw_namelist, config_user)


def check_namelist(raw_namelist):

    # TODO: use yaml schema for checking basic properties

    models = []
    for model in raw_namelist['models']:
        if model in models:
            raise NamelistError(
                "Duplicate model {} in models section".format(model))
        models.append(model)


class Namelist:
    """Namelist object"""

    def __init__(self, raw_namelist, config_user):
        """Parse a namelist file into an object."""

        self._cfg = config_user
        self.preprocessors = {}
        self.models = []
        self.diagnostics = {}

        for name, profile in raw_namelist['preprocessors'].items():
            # TODO: check profile
            self.preprocessors[name] = copy.deepcopy(profile)
            self.preprocessors[name]['name'] = name

        for model in raw_namelist['models']:
            # TODO: check model
            self.models.append(model)

        for name, diagnostic in raw_namelist['diagnostics'].items():
            diag = {}
            diag['name'] = name
            diag['script'] = diagnostic['script']
            diag['settings'] = copy.deepcopy(diagnostic['settings'])
            diag['models'] = \
                self._initialize_diagnostic_models(name, diagnostic['models'])
            diag['variables'] = \
                self._initialize_diagnostic_variables(name, diagnostic['variables'], diag['models'])
            self.diagnostics[name] = diag

    def _initialize_diagnostic_models(self, name, models):
        """Add models to diagnostic"""

        logger.debug("Populating list of models for diagnostic %s", name)

        def match(full_model, short_model):
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

        result = []
        for short_model in models:
            matches = [m for m in self.models if match(m, short_model)]
            if len(matches) == 1:
                result.append(matches[0])
            else:
                if not matches:
                    msg = ("No model matching description {} found in "
                           "diagnostic {}".format(short_model, name))
                else:
                    msg = ("Found multiple models {} matching description {} "
                           "in diagnostic {}".format(matches, short_model,
                                                     name))
                raise NamelistError(msg)

        return result

    def _initialize_diagnostic_variables(self, name, variables, models):

        logger.debug("Populating list of variables for diagnostic %s", name)

        result = []
        for variable in variables:
            var = copy.deepcopy(variable)

            # preprocessor settings
            var['preprocessor'] = {}
            for model in models:
                modelkey = '_'.join(str(model[key]) for key in sorted(model))
                var['preprocessor'][modelkey] = {}
                pre = var['preprocessor'][modelkey]

                # configure cube loading
                pre['load'] = {
                    'uris':
                    get_input_filelist(
                        project_info={'GLOBAL': self._cfg},
                        model=model,
                        var=variable),
                    'contstraint':
                    variable['short_name'],
                }

                # configure fixes
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

                # configure time extraction
                pre['extract_time'] = {
                    'start_year': model['start_year'],
                    'end_year': model['end_year'],
                }

            result.append(var)

        return result

    def __str__(self):
        return pprint.pformat(self.diagnostics, indent=2)
