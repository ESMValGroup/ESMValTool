# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for bcc-csm1-1-m."""
from . import add_irregular_bounds
from ..fix import Fix


class tos(Fix):
    """Fixes for tos."""

    def fix_data(self, cube):
        """Fix data."""
        return add_irregular_bounds(cube)


class sic(Fix):
    """Fixes for sic."""

    def fix_data(self, cube):
        """Fix data."""
        return add_irregular_bounds(cube)
