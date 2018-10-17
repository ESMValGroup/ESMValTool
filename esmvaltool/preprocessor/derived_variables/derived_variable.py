"""Contains the base class for derived variables."""


import importlib


class DerivedVariable(object):
    """Base class for derived variables."""

    def __init__(self, short_name=None):
        """Save desired short_name."""
        if short_name:
            self.short_name = short_name
        else:
            self.short_name = self.__class__.__name__

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
            `short_name` of the variable to derive.

        Returns
        -------
        DerivedVariable

        """
        derived_var = DerivedVariable(short_name)
        try:
            derived_var_module = importlib.import_module(
                'esmvaltool.preprocessor._derived_variables.'
                '{0}'.format(short_name))
            try:
                derived_var = getattr(derived_var_module, short_name)()
            except AttributeError:
                pass
        except ImportError:
            pass
        return derived_var
