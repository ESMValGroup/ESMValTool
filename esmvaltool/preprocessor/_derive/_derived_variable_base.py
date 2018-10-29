"""Contains the base class for derived variables."""


import importlib
import logging

logger = logging.getLogger(__name__)


class DerivedVariableBase(object):
    """Base class for derived variables."""

    # Static class member containing the required variables for derivation,
    # may include '{frequency}' in the field string
    _required_variables = {}

    def __init__(self, short_name='unknown_variable'):
        """Save `short_name` of derived variable."""
        self.short_name = short_name

    def get_required(self, frequency):
        """Return all required variables for derivation.

        Get variable `short_name` and `field` pairs required for derivation
        and optionally a list of needed fx files from the static class member
        _required_variables (may include '{frequency}' keyword in the `field`
        string for dynamic frequency adaption).

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        dict
            Dictionary containing a list of tuples `(short_name, field)` with
            the `vars` key and optionally a list of fx files with the key
            `fx_files`.

        Raises
        ------
        NotImplementedError
            If the desired variable derivation is not implemented, i.e. if no
            required variables are declared.

        """
        required_variables = dict(self.__class__._required_variables)
        if not required_variables:
            raise NotImplementedError("Don't know how to derive variable "
                                      "'{}'".format(self.short_name))
        if 'vars' not in required_variables:
            raise NotImplementedError(
                "Don't know how to derive variable '{}', all required "
                "variables have to be specified in the 'vars' key of the "
                "_required_variables dictionary (derivation from fx files "
                "only is not supported yet)".format(self.short_name))
        new_vars = []
        for (short_name, field) in required_variables['vars']:
            new_vars.append((short_name, field.format(frequency=frequency)))
        required_variables['vars'] = new_vars

        return required_variables

    def calculate(self, cubes):
        """Compute desired derived variable.

        This method needs to be overridden in the child class belonging to the
        desired variable to derive.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            Includes all the needed variables (incl. fx variables) for
            derivation defined in the static class variable
            `_required_variables`.

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
        module esmvaltool.preprocessor._derive.

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
