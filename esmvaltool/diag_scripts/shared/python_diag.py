"""
Convenience classes and functions to implement python diagnostics.
"""

import collections
import logging

logger = logging.getLogger(__name__)


DEFAULT_INFO_STR = 'not_specified'
INPUT_DATA_STR = 'input_data'
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
            success = True
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
                    success = False
            else:
                success = False
            if (not success):
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
    Class to easy save and access model data.

    Methods:
        get_data
        get_paths
        set_data
    """

    def __init__(self, cfg):
        """
        Create (private) dictionaries of the form
            {path1: {model1=..., project1=..., ...},
             path2: {model2=..., project2=..., ...},
             ...}
             (identical to 'input_data')
        and
            {path1: data1,
             path2: data2,
             ...}
        to get easy access to the models' data.
        """
        self._data = {}
        success = True
        if isinstance(cfg, dict):
            input_data = cfg.get(INPUT_DATA_STR)
            if isinstance(input_data, dict):
                for path in input_data:
                    model_info = input_data[path]
                    if (not isinstance(model_info, dict)):
                        success = False
                        break
                    self._data[path] = None
                self._models = input_data
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
        output = ''
        for path in self._models:
            output += repr(self._models[path]) + '\n'
        return output

    def get_data(self, **model_info):
        """
        Return all data in a list that match the given model informations.
        """
        pass

    def get_paths(self, **model_info):
        """
        Return all paths in a list that match the given model informations.
        """
        pass

    def set_data(self, **model_info):
        """
        Set the data of a specified model (the description must be unique).
        """
        pass
