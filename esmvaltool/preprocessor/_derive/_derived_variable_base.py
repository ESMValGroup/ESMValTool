"""Contains the base class for derived variables."""


import importlib
import logging
import os

logger = logging.getLogger(__name__)


ALL_DERIVED_VARIABLES = {}


class DerivedVariableBase(object):
    """Base class for derived variables."""

    def __init__(self, short_name='unknown_variable'):
        """Save `short_name` of derived variable."""
        self.short_name = short_name

    def get_required(self, frequency):
        """Get variable `short_name` and `field` pairs required for derivation.

        This method needs to be overridden in the child class belonging to the
        desired variable to derive.

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        list of tuples
            List of tuples (`short_name`, `field`) of all variables required
            for derivation.

        Raises
        ------
        NotImplementedError
            If the desired variable derivation is not implemented, i.e. if this
            method is called from this base class and not a child class.

        """
        raise NotImplementedError("Don't know how to derive variable "
                                  "'{}'".format(self.short_name))

    def calculate(self, cubes):
        """Compute desired derived variable.

        This method needs to be overridden in the child class belonging to the
        desired variable to derive.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            Includes all the needed variables for derivation defined in
            :func:`get_required`.

        Returns
        -------
        iris.cube.Cube
            New derived variable.

        Raises
        ------
        NotImplementedError
            If the desired variable derivation is not implemented, i.e. if this
            method is called from this base class and not a child class.

        """
        raise NotImplementedError("Don't know how to derive variable "
                                  "'{}'".format(self.short_name))

    @staticmethod
    def get_derived_variable(short_name):
        """Select correct python module for derived variable.

        Get derived variable by searching for a file `short_name.py` in the
        module esmvaltool.preprocessor._derived_variables.

        Parameters
        ----------
        short_name : str
            `short_name` of the requested variable.

        Returns
        -------
        DerivedVariableBase

        """
        derived_var = DerivedVariableBase(short_name)
        try:
            derived_var_module = importlib.import_module(
                'esmvaltool.preprocessor._derive.{0}'.format(short_name))
            try:
                derived_var = getattr(derived_var_module,
                                      'DerivedVariable')(short_name)
            except AttributeError:
                logger.warning("File ESMValTool/esmvaltool/preprocessor/"
                               "_derive/%s.py for variable derivation does "
                               "not contain required class 'DerivedVariable'",
                               short_name)
        except ImportError:
            logger.warning("No module named '%s' in ESMValTool/esmvaltool/"
                           "preprocessor/_derive/ for variable derivation",
                           short_name)
        return derived_var

    @staticmethod
    def get_all_derived_variables():
        """Get all possible derived variables.

        Returns
        -------
        dict
            All derived variables with `short_name` (keys) and the associated
            python classes (values).

        """
        if ALL_DERIVED_VARIABLES:
            return ALL_DERIVED_VARIABLES
        current_path = os.path.dirname(os.path.realpath(__file__))
        for var_file in os.listdir(current_path):
            var_name = os.path.splitext(var_file)[0]
            try:
                var_module = importlib.import_module(
                    'esmvaltool.preprocessor._derive.{}'.format(var_name))
                try:
                    derived_var = getattr(var_module,
                                          'DerivedVariable')(var_name)
                    ALL_DERIVED_VARIABLES[var_name] = derived_var
                except AttributeError:
                    pass
            except ImportError:
                pass

        return ALL_DERIVED_VARIABLES
