"""Derivation of variable `sm`."""

import cf_units
import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase

# Constants
WATER_DENSITY = 998.2
WATER_DENSITY_UNIT = cf_units.Unit('kg m^-3')


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `sm`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'mrsos',
            'field': 'T2{frequency}s'
        }]
    }

    def calculate(self, cubes):
        """Compute soil moisture.

        Note
        ----
        Convert moisture content of soil layer (kg/m2) into volumetric soil
        moisture (m3/m3), assuming density of water 998.2 kg/m2 (at temperature
        20 deg C).

        """
        mrsos_cube = cubes.extract_strict(
            Constraint(name='moisture_content_of_soil_layer'))

        depth = mrsos_cube.coord('depth').bounds
        height = depth[..., 1] - depth[..., 0]

        sm_cube = mrsos_cube / height / WATER_DENSITY
        sm_cube.convert_units(WATER_DENSITY_UNITS)

        return sm_cube
