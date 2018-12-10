"""Derivation of variable `sm`."""

import cf_units
import iris
from iris import Constraint
import numpy as np

from ._derived_variable_base import DerivedVariableBase

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

        sm_cube = mrsos_cube / height / 998.2
        sm_cube.units = cf_units.Unit('m3 m^-3')
        sm_cube.data = np.ma.array(sm_cube.data, dtype=np.dtype('float32'))

        return sm_cube
