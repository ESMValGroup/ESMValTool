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
EXPERIMENT_STR = 'exp'


Variable = collections.namedtuple('Variable', [VAR_STANDARD_NAME_STR,
                                               VAR_LONG_NAME_STR,
                                               VAR_UNITS_STR,
                                               EXPERIMENT_STR])


class Variables(object):
    """
    Class containing all variables.
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
                        exp = info.get(EXPERIMENT_STR, default)
                        attr = Variable(
                            info.get(VAR_STANDARD_NAME_STR, default),
                            info.get(VAR_LONG_NAME_STR, default),
                            info.get(VAR_UNITS_STR, default),
                            self._get_exps(name, exp))
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
            exp = getattr(attr, EXPERIMENT_STR)
            attr = Variable(getattr(attr, VAR_STANDARD_NAME_STR),
                            getattr(attr, VAR_LONG_NAME_STR),
                            getattr(attr, VAR_UNITS_STR),
                            self._get_exps(name, exp))
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

    def _get_exps(self, name, new_exp=None):
        try:
            exps = getattr(getattr(self, name), EXPERIMENT_STR)
        except AttributeError:
            exps = []
        if isinstance(new_exp, list):
            for e in new_exp:
                if (e not in exps):
                    exps.append(e)
                    logger.debug("Added experiment '{}' ".format(e) +
                                 "to variable '{}'".format(name))
        else:
            if (new_exp not in exps and new_exp is not None):
                exps.append(new_exp)
                logger.debug("Added experiment '{}' ".format(new_exp) +
                             "to variable '{}'".format(name))
        return exps

    def short_names(self):
        return list(self._dict)

    def standard_names(self):
        return [getattr(getattr(self, name), VAR_STANDARD_NAME_STR) for
                name in self._dict]


class Experiments(object):
    """
    Class containing all experiments.
    """

    def __init__(self, cfg=None, *names):
        self._dict = {}

        # Add experiments from cfg file
        if (cfg is not None):
            if isinstance(cfg, dict):
                data = cfg.get(INPUT_DATA_STR)
                if isinstance(data, dict):
                    for info in data.values():
                        default = 'not_specified'
                        name = info.get(EXPERIMENT_STR, default)
                        self._add_to_dict(name, name)
                else:
                    logger.warning("{} is not a valid ".format(repr(cfg)) +
                                   "configuration file!")
            else:
                logger.warning("{} is not a valid ".format(repr(cfg)) +
                               "configuration file!")

        # Add costum experiments
        for name in names:
            self._add_to_dict(name, name)
        if (not self._dict):
            logger.warning("No experiments found!")

    def __repr__(self):
        return 'Experiments: ' + repr(self.names())

    def _add_to_dict(self, name, attr):
        if (name not in self._dict):
            logger.debug("Added experiment '{}' to collection".format(name))
        self._dict[name] = attr

    def names(self):
        return list(self._dict)


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
