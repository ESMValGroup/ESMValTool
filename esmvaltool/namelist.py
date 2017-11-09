"""Namelist parser"""
import copy
import logging

import yaml

logger = logging.getLogger(__name__)


class NamelistError(Exception):
    """Raise if a namelist contains an error."""


class Namelist:
    """Namelist object"""

    def __init__(self, filename):
        """Parse a namelist file into an object."""
        with open(filename, 'r') as file:
            self._raw_dict = yaml.safe_load(file)

        self.preprocessor_profiles = {}
        self.models = []
        self.diagnostics = {}

        for name, profile in self._raw_dict['PREPROCESS'].items():
            # TODO: check profile
            self.preprocessor_profiles[name] = copy.deepcopy(profile)
            self.preprocessor_profiles[name]['name'] = name

        for model in self._raw_dict['MODELS']:
            # TODO: check model
            self.models.append(model)

        for name, diagnostic in self._raw_dict['DIAGNOSTICS'].items():
            self.diagnostics[name] = {'name': name}
            self.diagnostics[name]['scripts'] = copy.deepcopy(
                diagnostic['scripts'])
            self.diagnostics[name][
                'models'] = self._initialize_diagnostic_models(
                    name, diagnostic['models'])
            self.diagnostics[name]['variables'] = self._initialize_variables(
                name, diagnostic['variables'])

    def _initialize_diagnostic_models(self, name, models):
        """Add models to diagnostic"""

        def match(full_model, short_model):
            """Check if short_model description matches full_model."""
            if not short_model:
                return False
            for key in short_model:
                if (key not in full_model
                        or short_model[key] != full_model[key]):
                    return False
            return True

        selection = []
        for short_model in models:
            matches = [m for m in self.models if match(m, short_model)]
            if len(matches) == 1:
                selection.append(copy.deepcopy(matches[0]))
            else:
                if not matches:
                    msg = ("No model matching description {} found in "
                           "diagnostic {}".format(short_model, name))
                else:
                    msg = ("Multiple models {} matching description {} found "
                           "in diagnostic {}".format(matches, short_model,
                                                     name))
                raise NamelistError(msg)

        return selection

    def _initialize_variables(self, name, variables):
        pass
