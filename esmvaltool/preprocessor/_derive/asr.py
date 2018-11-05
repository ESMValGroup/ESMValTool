"""Derivation of variable `asr`."""

from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `asr`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'rsdt',
            'field': 'T2{frequency}s'
        }, {
            'short_name': 'rsut',
            'field': 'T2{frequency}s'
        }]
    }

    def calculate(self, cubes):
        """Compute absorbed shortwave radiation."""
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))

        asr_cube = rsdt_cube - rsut_cube
        asr_cube.units = rsdt_cube.units

        return asr_cube
