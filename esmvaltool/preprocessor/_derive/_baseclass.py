"""Contains the base class for derived variables."""
from abc import abstractmethod


class DerivedVariableBase:
    """Base class for derived variables."""

    @property
    @staticmethod
    @abstractmethod
    def required():
        """List of required variables for derivation."""

    @staticmethod
    @abstractmethod
    def calculate(cubes):
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
