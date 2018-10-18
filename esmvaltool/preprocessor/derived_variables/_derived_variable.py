"""Contains the base class for derived variables."""


import importlib


class DerivedVariable(object):
    """Base class for derived variables."""

    def __init__(self, variable=None):
        """Save `variable` dict and add `short_name` if necessary."""
        if variable is None:
            self.variable = {}
        else:
            self.variable = dict(variable)
        if 'short_name' not in self.variable:
            self.variable['short_name'] = self.__class__.__name__

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
                                  "'{}'".format(self.variable['short_name']))

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
                                  "'{}'".format(self.variable['short_name']))

    @staticmethod
    def get_derived_variable(variable):
        """Select correct python module for derived variable.

        Get derived variable by searching for a file `short_name.py` in the
        module esmvaltool.preprocessor._derived_variables.

        Parameters
        ----------
        variable : dict
            Contains all information of the requested variable.

        Returns
        -------
        DerivedVariable

        """
        derived_var = DerivedVariable(variable)
        try:
            derived_var_module = importlib.import_module(
                'esmvaltool.preprocessor.derived_variables.'
                '{0}'.format(variable['short_name']))
            try:
                derived_var = getattr(derived_var_module,
                                      variable['short_name'])(variable)
            except AttributeError:
                pass
        except ImportError:
            pass
        return derived_var
