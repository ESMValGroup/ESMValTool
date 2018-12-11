"""Derivation of variable `nbp_global`.

Note
----
Should not be necessary.

"""

from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `nbp_global`."""

    # Required variables
    # DUMMY
    _required_variables = {'vars': [{'short_name': 'nbp_global'}]}
