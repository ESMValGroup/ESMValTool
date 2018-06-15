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

EXP_STR = 'exp'
MODEL_STR = 'model'
OBS_STR = 'OBS'
PROJECT_STR = 'project'
SHORT_NAME_STR = 'short_name'


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

    def add_var(self, **names):
        """
        Add a costum variable to the class member.
        """
        for name in names:
            attr = Variable(*names[name])
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
        add_to_data
        get_data
        get_data_list
        get_exp
        get_model
        get_model_info
        get_model_info_list
        get_path
        get_path_list
        get_short_name
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

    def add_to_data(self, data, model_path=None, **model_info):
        """
        Add data to already given data (the description must be unique).
        """
        if (model_path is not None):
            if (model_path in self._paths):
                self._data[model_path] += data
                return None
            else:
                logger.warning("{} is not a valid ".format(model_path) +
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
        self._data[path] += data

    def get_data(self, model_path=None, **model_info):
        """
        Return data that matches the given model information or path. Fails if
        ambiguous.
        """
        if (model_path is not None):
            if (model_path in self._paths):
                return self._data.get(model_path)
            else:
                logger.warning("{} is not a valid ".format(model_path) +
                               "model path")
                return None
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
            return None
        if (len(paths) > 1):
            msg = 'Given model information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return self._data[paths[0]]

    def get_data_list(self, **model_info):
        """
        Return all data in a list that match the given model information
        (sorted alphabetically by path).
        """
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
        paths = sorted(paths)
        return [self._data[path] for path in paths]

    def get_exp(self, model_path):
        """
        Return 'exp' entry of given model.
        """
        if (model_path in self._paths):
            output = self._models[model_path].get(EXP_STR)
            if (output is None):
                logger.warning("Model {} does not ".format(model_path) +
                               "contain 'exp' information")
            return output
        else:
            logger.warning("{} is not a valid ".format(model_path) +
                           "model path")
            return None

    def get_model(self, model_path):
        """
        Return 'model' entry of given model.
        """
        if (model_path in self._paths):
            output = self._models[model_path].get(MODEL_STR)
            if (output is None):
                logger.warning("Model {} does not ".format(model_path) +
                               "contain 'model' information")
            return output
        else:
            logger.warning("{} is not a valid ".format(model_path) +
                           "model path")
            return None

    def get_model_info(self, model_path=None, **model_info):
        """
        Return model info that matches the given model information or path.
        Fails if ambiguous.
        """
        if (model_path is not None):
            if (model_path in self._paths):
                return self._models.get(model_path)
            else:
                logger.warning("{} is not a valid ".format(model_path) +
                               "model path")
                return None
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
            return None
        if (len(paths) > 1):
            msg = 'Given model information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return self._models[paths[0]]

    def get_model_info_list(self, model_path=None, **model_info):
        """
        Return all model info in a list that match the given model
        information (sorted alphabetically by path).
        """
        if (model_path is not None):
            if (model_path in self._paths):
                return [self._models.get(model_path)]
            else:
                logger.warning("{} is not a valid ".format(model_path) +
                               "model path")
                return None
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No data could be returned: " +
                           "{} does not match any model".format(model_info))
        paths = sorted(paths)
        return [self._models[path] for path in paths]

    def get_path(self, **model_info):
        """
        Return a path that match the given model information. Fails if
        ambiguous.
        """
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No paths could be returned: " +
                           "{} does not match any model".format(model_info))
            return None
        if (len(paths) > 1):
            msg = 'Given model information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return paths[0]

    def get_path_list(self, **model_info):
        """
        Return all paths in a list that match the given model information
        (sorted alphabetically).
        """
        paths = list(self._models)
        for info in model_info:
            paths = [path for path in paths if
                     self._models[path].get(info) == model_info[info]]
        if (not paths):
            logger.warning("No paths could be returned: " +
                           "{} does not match any model".format(model_info))
        paths = sorted(paths)
        return paths

    def get_project(self, model_path):
        """
        Return 'project' entry of given model.
        """
        if (model_path in self._paths):
            output = self._models[model_path].get(PROJECT_STR)
            if (output is None):
                logger.warning("Model {} does not ".format(model_path) +
                               "contain 'project' information")
            return output
        else:
            logger.warning("{} is not a valid ".format(model_path) +
                           "model path")
            return None

    def get_short_name(self, model_path):
        """
        Return 'short_name' entry of given model.
        """
        if (model_path in self._paths):
            output = self._models[model_path].get(SHORT_NAME_STR)
            if (output is None):
                logger.warning("Model {} does not ".format(model_path) +
                               "contain 'short_name' information")
            return output
        else:
            logger.warning("{} is not a valid ".format(model_path) +
                           "model path")
            return None

    def set_data(self, data, model_path=None, **model_info):
        """
        Set the data of a specified model (the description must be unique).
        """
        if (model_path is not None):
            if (model_path in self._paths):
                self._data[model_path] = data
                return None
            else:
                logger.warning("{} is not a valid ".format(model_path) +
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
