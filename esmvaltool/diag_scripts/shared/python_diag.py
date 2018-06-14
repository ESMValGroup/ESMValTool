"""
Convenience classes and functions to implement python diagnostics.
"""

import collections
import logging

logger = logging.getLogger(__name__)


# Global variables relevant for all diagnostics
TIME = 'time'
YEAR = 'year'
MONTH = 'month_number'
DAY_Y = 'day_of_year'
DAY_M = 'day_of_month'
LAT = 'latitude'
LON = 'longitude'
HEIGHT = 'height'


# Variables for the following classes
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
    Class containing all variables. Short names of variables can be called by
    member.var, all other information by member.VAR.

    Methods:
        add_var
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
        return repr(self.short_names())

    def _add_to_dict(self, name, attr):
        """
        Add variables to the class member.
        """
        if (name not in self._dict):
            logger.debug("Added variable '{}' to collection".format(name))
        setattr(self, name, name)
        setattr(self, name.upper(), attr)
        self._dict[name] = attr

    def add_var(self, name, *attr):
        """
        Add a costum variable to the class member.
        """
        attr = Variable(*attr)
        self._add_to_dict(name, attr)

    def short_names(self):
        """
        Return a list of the variables' short names.
        """
        return list(self._dict)

    def standard_names(self):
        """
        Return a list of the variables' standard names.
        """
        return [getattr(self._dict[name], VAR_STD_NAME_STR) for
                name in self._dict]


class Models(object):
    """
    Class to easily save and access model data.

    Methods:
        add_model
        get_data
        get_paths
        set_data
    """

    def __init__(self, cfg):
        """
        Create (private) lists and dictionaries of the form
            [path1, path2, path3, ...]
        and
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
        self._paths = []
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
                    self._paths.append(path)
                    self._data[path] = None
                self._models = input_data
            else:
                success = False
        else:
            success = False
        if (not success):
            raise TypeError("{} is not a valid ".format(repr(cfg)) +
                            "configuration file")
        self._n_models = len(self._paths)

    def __repr__(self):
        """
        Return the representation of the class member's dictionary keys.
        """
        output = ''
        for path in self._models:
            output += repr(self._models[path]) + '\n'
        return output

    def __iter__(self):
        """
        Allow iteration through class member.
        """
        self._iter_counter = 0
        return self

    def __next__(self):
        """
        Allow iteration through class member.
        """
        if (self._iter_counter >= self._n_models):
            raise StopIteration()
        else:
            next_element = self._paths[self._iter_counter]
            self._iter_counter += 1
            return next_element

    def add_model(self, path, data=None, **model_info):
        """
        Add costum model to the class member. Hint: 'path' can also be a unique
        identifier for the model if not used for file access.
        """
        if (path in self._paths):
            logger.warning("{} already exists! ".format(path) +
                           "Overwriting old data.")
            self._paths.remove(path)
        self._paths.append(path)
        self._data[path] = data
        self._models[path] = model_info

    def get_data(self, path_to_model=None, **model_info):
        """
        Return all data in a list that match the given model informations.
        """
        if (path_to_model is not None):
            if (path_to_model in self._paths):
                return self._data.get(path_to_model)
            else:
                logger.warning("{} is not a valid ".format(path_to_model) +
                               "model path")
                return None
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
        return [self._data[path] for path in paths]

    def get_paths(self, **model_info):
        """
        Return all paths in a list that match the given model informations.
        """
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
        return paths

    def set_data(self, data, path_to_model=None, **model_info):
        """
        Set the data of a specified model (the description must be unique).
        """
        if (path_to_model is not None):
            if (path_to_model in self._paths):
                self._data[path_to_model] = data
                return None
            else:
                logger.warning("{} is not a valid ".format(path_to_model) +
                               "model path")
                return None
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("Data could no be saved: " +
                           "{} does not match any model".format(model_info))
            return None
        if (len(paths) != 1):
            logger.warning("Data could no be saved: " +
                           "{} is ambiguous".format(model_info))
            return None
        self._data[path] = data
