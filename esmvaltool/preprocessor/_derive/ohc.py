"""Derivation of variable `ohc`."""
import iris
from iris import Constraint

import numpy as np
from cf_units import Unit

from ._derived_variable_base import DerivedVariableBase


RHO_CP = iris.cube.Cube(4.09169e+6, units=Unit('kg m-3 J kg-1 K-1'))


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `ohc`."""

    # Required variables
    _required_variables = {
        'vars': [
          {'short_name': 'thetao', 'field': 'TO3M'},
        ],
        'fx_files':['volcello']
    }

    def calculate(self, cubes,):
        """
        Compute ocean heat content.
        Use c_p*rho_0= 4.09169e+6 J m-3 K-1 (Kuhlbrodt et al., 2015, Clim. Dyn.)

        Arguments
        ---------
        cube: iris.cube.Cube
           input cube.

        Returns
        ---------
        iris.cube.Cube
              Output OHC cube.
        """
        # 1. Load the thetao and volcello cubes
        cube = cubes.extract_strict(Constraint(cube_func=lambda c: c.var_name=='thetao'))
        volume = cubes.extract_strict(Constraint(cube_func=lambda c: c.var_name=='volcello'))
        # 2. multiply with each other and with cprho0
        cube.data *= volume.data * RHO_CP.data
        cube.units *= volume.units * RHO_CP.units
        return cube
