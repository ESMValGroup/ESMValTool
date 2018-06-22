"""Convenience classes and functions to implement python diagnostics.

Example
-------
Import and use these basic classes by e.g.::

    import esmvaltool.diag_scripts.shared as e
    datasets = e.Datasets(cfg)

Notes
-----
An example diagnostic using these classes is given in
`diag_scripts/examples/diagnostic.py`

            if self._is_valid_path(dataset_path):

"""


import collections
import logging

from . import names as n

logger = logging.getLogger(__name__)


# Global variables
DEFAULT_INFO = 'not_specified'
INPUT_DATA = 'input_data'


# Variable class containing all relevant information
Variable = collections.namedtuple('Variable', [n.SHORT_NAME,
                                               n.STANDARD_NAME,
                                               n.LONG_NAME,
                                               n.UNITS])


class Variables(object):
    """Class to easily access a recipe's variables.

    This class is designed to easily access variables in the diagnostic script.

    Examples
    --------
    Get all variables of a recipe configuration `cfg`::

        vars = Variables(cfg)

    Access `short_name` (as str) of a variable `tas`::

        vars.tas

    Access all other information of a variable `tas`::

        vars.TAS.short_name
        vars.TAS.standard_name
        vars.TAS.long_name
        vars.TAS.units

    """

    def __init__(self, cfg=None, **names):
        """Load variables.

        Parameters
        ----------
        cfg : dict, optional
            Configuation dictionary of the recipe.
        **names : tuple or Variable, optional
            Keyword arguments of the form `short_name=Variable_object` where
            `Variable_object` can be given as tuple or Variable.

        """
        self._dict = {}

        # Add variables from cfg file
        if cfg is not None:
            success = True
            if isinstance(cfg, dict):
                data = cfg.get(INPUT_DATA)
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
        for name in names:
            attr = Variable(*names[name])
            self._add_to_dict(name, attr)
        if not self._dict:
            logger.warning("No variables found!")

    def __repr__(self):
        """Representation of the class."""
        return repr(self.short_names())

    def _add_to_dict(self, name, attr):
        """Add variable to class (internal method).

        Parameters
        ----------
        name : str
            `short_name` of the variable.
        attr : Variable
            All other information of the variable.

        """
        if name not in self._dict:
            logger.debug("Added variable '%s' to collection", name)
        setattr(self, name, name)
        setattr(self, name.upper(), attr)
        self._dict[name] = attr

    def add_var(self, **names):
        """Add a costum variable to the class member.

        Parameters
        ----------
        **names : tuple or Variable, optional
            Keyword arguments of the form `short_name=Variable_object` where
            `Variable_object` can be given as tuple or Variable.

        """
        for name in names:
            attr = Variable(*names[name])
            self._add_to_dict(name, attr)

    def short_names(self):
        """Get list of all `short_names`.

        Returns
        -------
        list
            List of all `short_names`.

        """
        return list(self._dict)

    def standard_names(self):
        """Get list of all `standard_names`.

        Returns
        -------
        list
            List of all `standard_names`.

        """
        return [getattr(self._dict[name], n.STANDARD_NAME) for
                name in self._dict]


class Datasets(object):
    """Class to easily access a recipe's datasets

    This class is designed to easily access datasets in the diagnostic script.

    Examples
    --------
    Get all variables of a recipe configuration `cfg`::

        datasets = Datasets(cfg)

    Access data of a dataset with path `path`::

        datasets.get_data(dataset_path=path)

    Access dataset information of the dataset::

        datasets.get_dataset_info(dataset_path=path)

    Access the data of all datasets with `exp=piControl'::

        datasets.get_data_list(exp=piControl)

    """

    def __init__(self, cfg):
        """Load datasets.

        Load all datasets of the recipe and store them in three internal
        dictionaries/lists (`self._paths`, `self._data` and `self._datasets`).

        Parameters
        ----------
        cfg : dict, optional
            Configuation dictionary of the recipe.

        """
        self._iter_counter = 0
        self._paths = []
        self._data = {}
        success = True
        if isinstance(cfg, dict):
            input_data = cfg.get(INPUT_DATA)
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
            raise TypeError("{} is not a valid ".format(repr(cfg)) +
                            "configuration file")
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
            True if valid path, False if not.

        """
        if path in self._paths:
            return True
        logger.warning("%s is not a valid dataset path", path)
        return False

    def _extract_paths(self, dataset_info):
        """Get all paths matching a given `dataset_info`.

        Parameters
        ----------
        dataset_info : dict
            Description of the desired datasets.

        Returns
        -------
        list
            All matching paths.

        """
        paths = list(self._datasets)
        for info in dataset_info:
            paths = [path for path in paths if
                     self._datasets[path].get(info) == dataset_info[info]]
        if not paths:
            logger.warning("%s does not match any dataset", dataset_info)
        return sorted(paths)

    def add_dataset(self, path, data=None, **dataset_info):
        """Add dataset to class.

        Parameters
        ----------
        path : str
            (Unique) path to the dataset.
        data, optional
            Arbitrary object to be saved as data for the dataset.
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        """
        if path in self._paths:
            logger.warning("%s already exists! Overwriting old data", path)
            self._paths.remove(path)
        self._paths.append(path)
        self._data[path] = data
        self._datasets[path] = dataset_info

    def add_to_data(self, data, dataset_path=None, **dataset_info):
        """Add element to a dataset's data.

        Notes
        -----
        Either `dataset_path` or a unique `dataset_info` description have to be
        given. Prints warning and does nothing if given information is
        ambiguous.

        Parameters
        ----------
        data
            Element to be added to the dataset's data.
        dataset_path : str, optional
            Path to the dataset
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        """
        if dataset_path is not None:
            if self._is_valid_path(dataset_path):
                self._data[dataset_path] += data
                return None
            return None
        paths = self._extract_paths(dataset_info)
        if not paths:
            return None
        if len(paths) > 1:
            logger.warning("Data could no be saved: %s is ambiguous",
                           dataset_info)
            return None
        self._data[paths[0]] += data
        return None

    def get_data(self, dataset_path=None, **dataset_info):
        """Access a dataset's data.

        Notes
        -----
        Either `dataset_path` or a unique `dataset_info` description have to be
        given. Fails when given information is ambiguous.

        Parameters
        ----------
        dataset_path : str, optional
            Path to the dataset
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        data_object
            Data of the selected dataset

        Raises
        ------
        RuntimeError
            If data given by `dataset_info` is ambiguous.

        """
        if dataset_path is not None:
            if self._is_valid_path(dataset_path):
                return self._data.get(dataset_path)
            return None
        paths = self._extract_paths(dataset_info)
        if not paths:
            return None
        if len(paths) > 1:
            msg = 'Given dataset information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return self._data[paths[0]]

    def get_data_list(self, **dataset_info):
        """Access the datasets' data in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Data of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        return [self._data[path] for path in paths]

    def get_exp(self, dataset_path):
        """Access a dataset's `exp`.

        Notes
        -----
        If the `dataset_info` does not contain an `exp` value, returns None.

        Parameters
        ----------
        dataset_path : str
            Path to the dataset

        Returns
        -------
        str
            `exp` information of the given dataset.

        """
        if self._is_valid_path(dataset_path):
            output = self._datasets[dataset_path].get(n.EXP)
            if output is None:
                logger.warning("Dataset %s does not contain '%s' information",
                               dataset_path, n.EXP)
            return output
        return None

    def get_dataset(self, dataset_path):
        """Access a dataset's `dataset`.

        Notes
        -----
        If the `dataset_info` does not contain a `dataset` value, returns None.

        Parameters
        ----------
        dataset_path : str
            Path to the dataset

        Returns
        -------
        str
            `dataset` information of the given dataset.

        """
        if self._is_valid_path(dataset_path):
            output = self._datasets[dataset_path].get(n.DATASET)
            if output is None:
                logger.warning("Dataset %s does not contain '%s' information",
                               dataset_path, n.DATASET)
            return output
        return None

    def get_dataset_info(self, dataset_path=None, **dataset_info):
        """Access a dataset's information.

        Notes
        -----
        Either `dataset_path` or a unique `dataset_info` description have to be
        given. Fails when given information is ambiguous.

        Parameters
        ----------
        dataset_path : str, optional
            Path to the dataset
        **dataset_info, optional
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
        if dataset_path is not None:
            if self._is_valid_path(dataset_path):
                return self._datasets.get(dataset_path)
            return None
        paths = self._extract_paths(dataset_info)
        if not paths:
            return None
        if len(paths) > 1:
            msg = 'Given dataset information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return self._datasets[paths[0]]

    def get_dataset_info_list(self, **dataset_info):
        """Access datasets information in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Information dictionaries of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        return [self._datasets[path] for path in paths]

    def get_path(self, **dataset_info):
        """Access a dataset's path

        Notes
        -----
        A unique `dataset_info` description has to be given. Fails when given
        information is ambiguous.

        Parameters
        ----------
        **dataset_info, optional
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
        paths = self._extract_paths(dataset_info)
        if not paths:
            return None
        if len(paths) > 1:
            msg = 'Given dataset information is ambiguous'
            logger.error(msg)
            raise RuntimeError(msg)
        return paths[0]

    def get_path_list(self, **dataset_info):
        """Access datasets paths in a list.

        Notes
        -----
        The returned data is sorted alphabetically respective to the `paths`.

        Parameters
        ----------
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        Returns
        -------
        list
            Paths of the selected datasets.

        """
        paths = self._extract_paths(dataset_info)
        return paths

    def get_project(self, dataset_path):
        """Access a dataset's `project`.

        Notes
        -----
        If the `dataset_info` does not contain a `project` value, returns None.

        Parameters
        ----------
        dataset_path : str
            Path to the dataset

        Returns
        -------
        str
            `project` information of the given dataset.

        """
        if self._is_valid_path(dataset_path):
            output = self._datasets[dataset_path].get(n.PROJECT)
            if output is None:
                logger.warning("Dataset %s does not contain '%s' information",
                               dataset_path, n.PROJECT)
            return output
        return None

    def get_short_name(self, dataset_path):
        """Access a dataset's `short_name`.

        Notes
        -----
        If the `dataset_info` does not contain a `short_name` value, returns
        None.

        Parameters
        ----------
        dataset_path : str
            Path to the dataset

        Returns
        -------
        str
            `short_name` information of the given dataset.

        """
        if self._is_valid_path(dataset_path):
            output = self._datasets[dataset_path].get(n.SHORT_NAME)
            if output is None:
                logger.warning("Dataset %s does not contain '%s' information",
                               dataset_path, n.SHORT_NAME)
            return output
        return None

    def get_standard_name(self, dataset_path):
        """Access a dataset's `standard_name`.

        Notes
        -----
        If the `dataset_info` does not contain a `standard_name` value, returns
        None.

        Parameters
        ----------
        dataset_path : str
            Path to the dataset

        Returns
        -------
        str
            `standard_name` information of the given dataset.

        """
        if self._is_valid_path(dataset_path):
            output = self._datasets[dataset_path].get(n.STANDARD_NAME)
            if output is None:
                logger.warning("Dataset %s does not contain '%s' information",
                               dataset_path, n.STANDARD_NAME)
            return output
        return None

    def set_data(self, data, dataset_path=None, **dataset_info):
        """Set element as a dataset's data.

        Notes
        -----
        Either `dataset_path` or a unique `dataset_info` description have to be
        given. Prints warning and does nothing if given information is
        ambiguous.

        Parameters
        ----------
        data
            Element to be set as the dataset's data.
        dataset_path : str, optional
            Path to the dataset
        **dataset_info, optional
            Keyword arguments describing the dataset, e.g. `dataset=CanESM2`,
            `exp=piControl` or `short_name=tas`.

        """
        if dataset_path is not None:
            if self._is_valid_path(dataset_path):
                self._data[dataset_path] = data
                return None
            return None
        paths = self._extract_paths(dataset_info)
        if not paths:
            return None
        if len(paths) != 1:
            logger.warning("Data could no be saved: %s is ambiguous",
                           dataset_info)
            return None
        self._data[paths[0]] = data
        return None
