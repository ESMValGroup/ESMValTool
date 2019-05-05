"""Derivation of variable `sm`."""

import cf_units
import numpy as np
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `sm`."""

    # Required variables
    required = [{'short_name': 'mrsos'}]

    @staticmethod
    def calculate(cubes):
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
        layer_thickness = depth[..., 1] - depth[..., 0]

        sm_cube = mrsos_cube / layer_thickness / 998.2
        sm_cube.units = cf_units.Unit('m3 m^-3')
        sm_cube.data = np.ma.array(sm_cube.data, dtype=np.dtype('float32'))

        return sm_cube
