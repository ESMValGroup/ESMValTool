"""Convenience functions for running a diagnostic script."""
import collections
import contextlib
import logging
import os
import sys

import yaml

logger = logging.getLogger(__name__)


DEFAULT_INFO_STR = 'not_specified'
EXP_STR = 'exp'
INPUT_DATA_STR = 'input_data'
MODEL_NAME_STR = 'model'
VAR_LONG_NAME_STR = 'long_name'
VAR_SHORT_NAME_STR = 'short_name'
VAR_STD_NAME_STR = 'standard_name'
VAR_UNITS_STR = 'units'


Variable = collections.namedtuple('Variable', [VAR_SHORT_NAME_STR,
                                               VAR_STD_NAME_STR,
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
        """
        Create (private) dictionary containing all variable information and add
        attributes to the class with the variables' short names.
        """
        self._dict = {}

        # Add variables from cfg file
        if (cfg is not None):
            if isinstance(cfg, dict):
                data = cfg.get(INPUT_DATA_STR)
                if isinstance(data, dict):
                    for info in data.values():
                        name = info.get(VAR_SHORT_NAME_STR, DEFAULT_INFO_STR)
                        attr = Variable(
                            name,
                            info.get(VAR_STD_NAME_STR, DEFAULT_INFO_STR),
                            info.get(VAR_LONG_NAME_STR, DEFAULT_INFO_STR),
                            info.get(VAR_UNITS_STR, DEFAULT_INFO_STR))
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
            self._add_to_dict(name, attr)
        if (not self._dict):
            logger.warning("No variables found!")

    def __repr__(self):
        """
        Return all the short names of the included variables.
        """
        return 'Variables: ' + repr(self.short_names())

    def _add_to_dict(self, name, attr):
        """
        Add variables to the class member.
        """
        if (name not in self._dict):
            logger.debug("Added variable '{}' to collection".format(name))
        setattr(self, name, attr)
        self._dict[name] = attr

    def short_names(self):
        """
        Return a list of the variables' short names.
        """
        return list(self._dict)

    def standard_names(self):
        """
        Return a list of the variables' standard names.
        """
        return [getattr(getattr(self, name), VAR_STD_NAME_STR) for
                name in self._dict]


class Models(object):
    """
    Class to easy set and get model data.

    Methods:
    """

    def __init__(self, cfg):
        """
        Create (private) dictionary of the form
            {(path1, var1, exp1, model1): data1,
             (path2, var2, exp2, model2): data2,
             ...}
        to get easy access to the models' data.
        """
        self._dict = {}
        success = True
        if isinstance(cfg, dict):
            input_data = cfg.get(INPUT_DATA_STR)
            if isinstance(input_data, dict):
                for path in input_data:
                    model_info = input_data[path]
                    if (not isinstance(model_info, dict)):
                        success = False
                        break
                    var = Variable(
                        model_info.get(VAR_SHORT_NAME_STR, DEFAULT_INFO_STR),
                        model_info.get(VAR_STD_NAME_STR, DEFAULT_INFO_STR),
                        model_info.get(VAR_LONG_NAME_STR, DEFAULT_INFO_STR),
                        model_info.get(VAR_UNITS_STR, DEFAULT_INFO_STR))
                    dict_key = (path,
                                var,
                                model_info[EXP_STR],
                                model_info[MODEL_NAME_STR])
                    dict_value = None
                    self._dict[dict_key] = dict_value
            else:
                success = False
        else:
            success = False
        if (not success):
            raise TypeError("{} is not a valid ".format(repr(cfg)) +
                             "configuration file")

    def __repr__(self):
        """
        Return the representation of the class member's dictionary keys.
        """
        output=''
        for model_info in self._dict:
            output += repr(model_info[1:])
            output += '\n'
        return output


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
