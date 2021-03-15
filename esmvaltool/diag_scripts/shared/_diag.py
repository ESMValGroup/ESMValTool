"""Convenience classes and functions to implement python diagnostics.

Example
-------
Import and use these basic classes by e.g.::

    import esmvaltool.diag_scripts.shared as e
    datasets = e.Datasets(cfg)
    variables = e.Variables(cfg)

"""


import collections
import logging
import warnings

from esmvaltool import ESMValToolDeprecationWarning

from . import names as n

logger = logging.getLogger(__name__)


# Global variables
DEFAULT_INFO = 'not_specified'


# Variable class containing all relevant information
BaseVariable = collections.namedtuple('Variable', [n.SHORT_NAME,
                                                   n.STANDARD_NAME,
                                                   n.LONG_NAME,
                                                   n.UNITS])


class Variable(BaseVariable):
    """Variable class containing all relevant information.

    Note
    ----
    This class has been deprecated in version 2.2 and will be removed two minor
    releases later in version 2.4.

    """

    def __new__(cls, short_name, standard_name, long_name, units):
        """Deprecate this class."""
        warnings.warn(
            "'Variable' has been deprecated in version 2.2 and will be "
            "removed two minor releases later in version 2.4",
            ESMValToolDeprecationWarning)
        self = super().__new__(cls, short_name, standard_name, long_name,
                               units)
        return self


class Variables(object):
    """Class to easily access a recipe's variables in a diagnostic.

    Note
    ----
    This class has been deprecated in version 2.2 and will be removed two minor
    releases later in version 2.4.

    Examples
    --------
    Get all variables of a recipe configuration `cfg`::

        variables = Variables(cfg)

    Access information of a variable `tas`::

        variables.short_name('tas')
        variables.standard_name('tas')
        variables.long_name('tas')
        variables.units('tas')

    Access :mod:`iris`-suitable dictionary of a variable `tas`::

        variables.iris_dict('tas')

    Check if variables `tas` and `pr` are available::

        variables.vars_available('tas', 'pr')

    """

    def __init__(self, cfg=None, **names):
        """Load variables.

        Parameters
        ----------
        cfg : dict, optional
            Configuation dictionary of the recipe.
        **names : dict or Variable, optional
            Keyword arguments of the form `short_name=Variable_object` where
            `Variable_object` can be given as :obj:`dict` or :class:`Variable`.

        """
        warnings.warn(
            "'Variables' has been deprecated in version 2.2 and will be "
            "removed two minor releases later in version 2.4",
            ESMValToolDeprecationWarning)
        self._dict = {}

        # Add variables from cfg file
        if cfg is not None:
            success = True
            if isinstance(cfg, dict):
                data = cfg.get(n.INPUT_DATA)
                if isinstance(data, dict):
                    for info in data.values():
                        name = info.get(n.SHORT_NAME, DEFAULT_INFO)
                        attr = Variable(
                            name,
                            info.get(n.STANDARD_NAME, DEFAULT_INFO),
                            info.get(n.LONG_NAME, DEFAULT_INFO),
                            info.get(n.UNITS, DEFAULT_INFO))
                        self._add_to_dict(name, attr)
                else:
                    success = False
            else:
                success = False
            if not success:
                logger.warning("%s is not a valid configuration file!", cfg)
        if not self._dict:
            logger.warning("Empty recipe configuration: the automatic "
                           "import of variables does not work for chained "
                           "scripts (using 'ancestors' key)")

        # Add costum variables
        self.add_vars(**names)
        if not self._dict:
            logger.warning("No variables found!")

    def __repr__(self):
        """Representation of the class."""
        output = ''
        for (name, attr) in self._dict.items():
            output += '{}: {}\n'.format(name, attr)
        return output

    def _add_to_dict(self, name, attr):
        """Add variable to class dictionary.

        Parameters
        ----------
        name : str
            `short_name` of the variable.
        attr : Variable
            All other information of the variable.

        """
        if name not in self._dict:
            logger.debug("Added variable '%s' to collection", name)
        self._dict[name] = attr

    def add_vars(self, **names):
        """Add costum variables to the class.

        Parameters
        ----------
        **names : dict or Variable, optional
            Keyword arguments of the form `short_name=Variable_object` where
            `Variable_object` can be given as :obj:`dict` or :class:`Variable`.

        """
        for (name, attr) in names.items():
            if isinstance(attr, Variable):
                attr_var = attr
            else:
                attr_var = Variable(
                    name,
                    attr.get(n.STANDARD_NAME, DEFAULT_INFO),
                    attr.get(n.LONG_NAME, DEFAULT_INFO),
                    attr.get(n.UNITS, DEFAULT_INFO))
            self._add_to_dict(name, attr_var)

    def iris_dict(self, var):
        """Access :mod:`iris` dictionary of the variable.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        dict
            Dictionary containing all attributes of the variable which can be
            used directly in :mod:`iris` (`short_name` replaced by `var_name`).

        """
        iris_dict = self._dict[var]._asdict()
        iris_dict[n.VAR_NAME] = iris_dict.pop(n.SHORT_NAME)
        return iris_dict

    def long_name(self, var):
        """Access long name.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        str
            Long name of the variable.

        """
        return getattr(self._dict[var], n.LONG_NAME)

    def modify_var(self, var, **names):
        """Modify an already existing  variable of the class.

        Parameters
        ----------
        var : str
            (Short) name of the existing variable.
        **names
            Keyword arguments of the form `short_name=tas`.

        Raises
        ------
        ValueError
            If `var` is not an existing variable.
        TypeError
            If a non-valid keyword argument is given.

        """
        if var not in self._dict:
            raise ValueError("Variable '{}' does not exist yet and cannot be "
                             "modified".format(var))
        old_var = self._dict.pop(var)
        new_var = {}
        for name in Variable._fields:
            new_var[name] = names.pop(name, getattr(old_var, name))

        # Check if names is not empty (=non-valid keyword argument given)
        if names:
            raise TypeError("Non-valid keyword arguments "
                            "given: {}".format(names))
        new_var = Variable(**new_var)
        self._add_to_dict(var, new_var)

    def short_name(self, var):
        """Access short name.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        str
            Short name of the variable.

        """
        return getattr(self._dict[var], n.SHORT_NAME)

    def short_names(self):
        """Get list of all `short_names`.

        Returns
        -------
        list
            List of all `short_names`.

        """
        return list(self._dict)

    def standard_name(self, var):
        """Access standard name.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        str
            Standard name of the variable.

        """
        return getattr(self._dict[var], n.STANDARD_NAME)

    def standard_names(self):
        """Get list of all `standard_names`.

        Returns
        -------
        list
            List of all `standard_names`.

        """
        return [self.standard_name(name) for name in self._dict]

    def units(self, var):
        """Access units.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        str
            Units of the variable.

        """
        return getattr(self._dict[var], n.UNITS)

    def var_name(self, var):
        """Access var name.

        Parameters
        ----------
        var : str
            (Short) name of the variable.

        Returns
        -------
        str
            Var name (=short name) of the variable.

        """
        return getattr(self._dict[var], n.SHORT_NAME)

    def vars_available(self, *args):
        """Check if given variables are available.

        Parameters
        ----------
        *args
            Short names of the variables to be tested.

        Returns
        -------
        bool
            `True` if variables are available, `False` if not.

        """
        for var in args:
            if var not in self._dict:
                return False
        return True


class Datasets(object):
    """Class to easily access a recipe's datasets in a diagnostic script.

    Note
    ----
    This class has been deprecated in version 2.2 and will be removed two minor
    releases later in version 2.4.

    Examples
    --------
    Get all variables of a recipe configuration `cfg`::

        datasets = Datasets(cfg)

    Access data of a dataset with path `dataset_path`::

        datasets.get_data(path=dataset_path)

    Access dataset information of the dataset::

        datasets.get_dataset_info(path=dataset_path)

    Access the data of all datasets with `exp=piControl`::

        datasets.get_data_list(exp=piControl)

    """

    def __init__(self, cfg):
        """Load datasets.

        Load all datasets of the recipe and store them in three internal
        :obj:`dict`/:obj:`list` containers: `self._paths`, `self._data` and
        `self._datasets`.

        Parameters
        ----------
        cfg : dict, optional
            Configuation dictionary of the recipe.

        Raises
        ------
        TypeError
            If recipe configuration dictionary is not valid.

        """
        warnings.warn(
            "'Datasets' has been deprecated in version 2.2 and will be "
            "removed two minor releases later in version 2.4",
            ESMValToolDeprecationWarning)
        self._iter_counter = 0
        self._paths = []
        self._data = {}
        success = True
        if isinstance(cfg, dict):
            input_data = cfg.get(n.INPUT_DATA)
            if isinstance(input_data, dict):
                for path in input_data:
                    dataset_info = input_data[path]
                    if not isinstance(dataset_info, dict):
                        success = False
                        break
                    self._paths.append(path)
                    self._data[path] = None
                self._datasets = input_data
            else:
                success = False
        else:
            success = False
        if not success:
            raise TypeError("{} is not a valid configuration "
                            "file".format(repr(cfg)))
        self._n_datasets = len(self._paths)
        if not self._paths:
            logger.warning("No datasets found!")
            logger.warning("Note: the automatic import of datasets does not "
                           "work for chained scripts (using 'ancestors' key)")

    def __repr__(self):
        """Representation of the class."""
        output = ''
        for path in self._datasets:
            output += repr(self._datasets[path]) + '\n'
        return output

    def __iter__(self):
        """Allow iteration through class."""
        self._iter_counter = 0
        return self

    def __next__(self):
        """Allow iteration through class."""
        if self._iter_counter >= self._n_datasets:
            raise StopIteration()
        next_element = self._paths[self._iter_counter]
        self._iter_counter += 1
        return next_element

    def _is_valid_path(self, path):
        """Check if path is in class.

        Parameters
        ----------
        path : str
            Path to be tested.

        Returns
        -------
        bool
            `True` if valid path, `False` if not.

        """
        if path in self._paths:
            return True
        logger.warning("%s is not a valid dataset path", path)
        return False

    def _extract_paths(self, dataset_info, fail_when_ambiguous=False):
        """Get all paths matching a given `dataset_info`.

        Parameters
        ----------
        dataset_info : dict
            Description of the desired datasets.
        fail_when_ambiguous : bool, optional
            Raise an exception when retrieved paths are ambiguous.

        Returns
        -------
        list
            All matching paths.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous and
            `fail_when_ambiguous` is set to `True`.

        """
        paths = list(self._datasets)
        for info in dataset_info:
            paths = [path for path in paths if
                     self._datasets[path].get(info) == dataset_info[info]]
        if not paths:
            logger.warning("%s does not match any dataset", dataset_info)
            return paths
        if not fail_when_ambiguous:
            return sorted(paths)
        if len(paths) > 1:
            msg = 'Given dataset information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return sorted(paths)

    def add_dataset(self, path, data=None, **dataset_info):
        """Add dataset to class.

        Parameters
        ----------
        path : str
            (Unique) path to the dataset.
        data: optional
            Arbitrary object to be saved as data for the dataset.
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        """
        if path in self._paths:
            logger.warning("%s already exists! Overwriting old data", path)
            self._paths.remove(path)
        self._paths.append(path)
        self._data[path] = data
        self._datasets[path] = dataset_info

    def add_to_data(self, data, path=None, **dataset_info):
        """Add element to a dataset's data.

        Notes
        -----
        Either `path` or a unique `dataset_info` description have to be given.
        Fails when given information is ambiguous.

        Parameters
        ----------
        data
            Element to be added to the dataset's data.
        path : str, optional
            Path to the dataset
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if path is not None:
            if self._is_valid_path(path):
                self._data[path] += data
                return None
            return None
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if paths:
            self._data[paths[0]] += data
        return None

    def get_data(self, path=None, **dataset_info):
        """Access a dataset's data.

        Notes
        -----
        Either `path` or a unique `dataset_info` description have to be
        given. Fails when given information is ambiguous.

        Parameters
        ----------
        path : str, optional
            Path to the dataset
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        `data_object`
            Data of the selected dataset.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if path is not None:
            if self._is_valid_path(path):
                return self._data.get(path)
            return None
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if not paths:
            return None
        return self._data[paths[0]]

    def get_data_list(self, **dataset_info):
        """Access the datasets' data in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Data of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        return [self._data[path] for path in paths]

    def get_dataset_info(self, path=None, **dataset_info):
        """Access a dataset's information.

        Notes
        -----
        Either `path` or a unique `dataset_info` description have to be
        given. Fails when given information is ambiguous.

        Parameters
        ----------
        path : str, optional
            Path to the dataset.
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        dict
            All dataset information.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if path is not None:
            if self._is_valid_path(path):
                return self._datasets.get(path)
            return None
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if not paths:
            return None
        return self._datasets[paths[0]]

    def get_dataset_info_list(self, **dataset_info):
        """Access dataset's information in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Information dictionaries of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        return [self._datasets[path] for path in paths]

    def get_info(self, key, path=None, **dataset_info):
        """Access a 'dataset_info`'s `key`.

        Notes
        -----
        Either `path` or a unique `dataset_info` description have to be
        given. Fails when given information is ambiguous. If the `dataset_info`
        does not contain the `key`, returns None.

        Parameters
        ----------
        key : str
            Desired dictionary key.
        path : str
            Path to the dataset.
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        str
            `key` information of the given dataset.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if path is not None:
            if self._is_valid_path(path):
                output = self._datasets[path].get(key)
                if output is None:
                    logger.warning("Dataset %s does not contain '%s' "
                                   "information", path, key)
                return output
            return None
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if not paths:
            return None
        output = self._datasets[paths[0]].get(key)
        if output is None:
            logger.warning("Dataset %s does not contain '%s' information",
                           path, key)
        return output

    def get_info_list(self, key, **dataset_info):
        """Access `dataset_info`'s `key` values.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            `key` information of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        output = [self._datasets[path].get(key) for path in paths]
        if None in output:
            logger.warning("One or more datasets do not contain '%s' "
                           "information", key)
        return output

    def get_path(self, **dataset_info):
        """Access a dataset's path.

        Notes
        -----
        A unique `dataset_info` description has to be given. Fails when given
        information is ambiguous.

        Parameters
        ----------
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        str
            Path of the selected dataset.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if not paths:
            return None
        return paths[0]

    def get_path_list(self, **dataset_info):
        """Access dataset's paths in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Paths of the selected datasets.

        """
        return self._extract_paths(dataset_info)

    def set_data(self, data, path=None, **dataset_info):
        """Set element as a dataset's data.

        Notes
        -----
        Either `path` or a unique `dataset_info` description have to be
        given. Fails when if given information is ambiguous.

        Parameters
        ----------
        data
            Element to be set as the dataset's data.
        path : str, optional
            Path to the dataset.
        **dataset_info: optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if path is not None:
            if self._is_valid_path(path):
                self._data[path] = data
                return None
            return None
        paths = self._extract_paths(dataset_info, fail_when_ambiguous=True)
        if paths:
            self._data[paths[0]] = data
        return None
