"""Convenience functions for running a diagnostic script."""
import collections
import contextlib
import logging
import os
import sys

import yaml

logger = logging.getLogger(__name__)


VAR_STANDARD_NAME_STR = 'standard_name'
VAR_SHORT_NAME_STR = 'short_name'
VAR_LONG_NAME_STR = 'long_name'
VAR_UNITS_STR = 'units'
INPUT_DATA_STR = 'input_data'


Variable = collections.namedtuple('Variable', [VAR_STANDARD_NAME_STR,
                                               VAR_LONG_NAME_STR,
                                               VAR_UNITS_STR])


class Variables(object):
    """
    Class containing all variables.

    Methods:
        short_names
        standard_names
    """

    def __init__(self, cfg=None, **names):
        self._dict = {}

        # Add variables from cfg file
        if (cfg is not None):
            if isinstance(cfg, dict):
                data = cfg.get(INPUT_DATA_STR)
                if isinstance(data, dict):
                    for info in data.values():
                        default = 'not_specified'
                        name = info.get(VAR_SHORT_NAME_STR, default)
                        attr = Variable(
                            info.get(VAR_STANDARD_NAME_STR, default),
                            info.get(VAR_LONG_NAME_STR, default),
                            info.get(VAR_UNITS_STR, default))
                        self._add_to_dict(name, attr)
                else:
                    logger.warning("{} is not a valid ".format(repr(cfg)) +
                                   "configuration file!")
            else:
                logger.warning("{} is not a valid ".format(repr(cfg)) +
                               "configuration file!")

        # Add costum variables
        for name in names:
            attr = Variable(*names[name])
            attr = Variable(getattr(attr, VAR_STANDARD_NAME_STR),
                            getattr(attr, VAR_LONG_NAME_STR),
                            getattr(attr, VAR_UNITS_STR))
            self._add_to_dict(name, attr)
        if (not self._dict):
            logger.warning("No variables found!")

    def __repr__(self):
        return 'Variables: ' + repr(self.short_names())

    def _add_to_dict(self, name, attr):
        if (name not in self._dict):
            logger.debug("Added variable '{}' to collection".format(name))
        setattr(self, name, attr)
        self._dict[name] = attr

    def short_names(self):
        return list(self._dict)

    def standard_names(self):
        return [getattr(getattr(self, name), VAR_STANDARD_NAME_STR) for
                name in self._dict]


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def _get_input_data_files(cfg):
    """Get a dictionary containing all data input files."""
    input_files = {}
    for filename in cfg['input_files']:
        with open(filename) as file:
            metadata = yaml.safe_load(file)
            input_files.update(metadata)

    return input_files


@contextlib.contextmanager
def run_diagnostic():
    """Run a diagnostic."""
    # Implemented as context manager so we can support clean up actions later
    cfg = get_cfg()

    # Set up logging
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    logging.getLogger().setLevel(cfg['log_level'].upper())

    # Read input metadata
    cfg['input_data'] = _get_input_data_files(cfg)

    # Create output directories
    if cfg['write_netcdf']:
        os.makedirs(cfg['work_dir'])
    if cfg['write_plots']:
        os.makedirs(cfg['plot_dir'])

    yield cfg
