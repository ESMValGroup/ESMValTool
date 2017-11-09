import yaml
import logging

logger = logging.getLogger(__name__)


class Namelist:

    def __init__(self, filename):
        """Parses an namelist file into an object.

        Uses the pyyaml parser for building a raw dict, then
        builds the objects from this dict by hand. This can probably
        be done smarter, but this is simpler to experiment with for now.
        """
        with open(filename, 'r') as f:
            self.raw_dict = yaml.load(f)

        self.preprocess_profiles = {}

        self.models = []

        if not self.raw_dict['PREPROCESS']:
            raise Exception('Mandatory preprocess_presets section not found')

        for profile in self.raw_dict['PREPROCESS']:
            # TODO: check profile
            self.preprocess_profiles[profile['id']] = profile

        for model in self.raw_dict['MODELS']:
            # TODO: check model
            self.models.append(model)

        for name, diagnostic in self.raw_dict['DIAGNOSTICS'].items():
            print(name)
            print(diagnostic['description'])
