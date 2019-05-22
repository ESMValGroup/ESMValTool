"""Derivation of variable `ohc`."""
import iris
from iris import Constraint

from cf_units import Unit

from ._baseclass import DerivedVariableBase

RHO_CP = iris.coords.AuxCoord(4.09169e+6, units=Unit('kg m-3 J kg-1 K-1'))


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `ohc`."""

    # Required variables
    required = [
        {
            'short_name': 'thetao',
            'field': 'TO3M',
            'fx_files': [
                'volcello',
            ],
        },
    ]

    def calculate(
            self,
            cubes,
    ):
        """
        Compute ocean heat content.

        Use c_p*rho_0= 4.09169e+6 J m-3 K-1
        (Kuhlbrodt et al., 2015, Clim. Dyn.)

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
        cube = cubes.extract_strict(
            Constraint(cube_func=lambda c: c.var_name == 'thetao'))
        volume = cubes.extract_strict(
            Constraint(cube_func=lambda c: c.var_name == 'volcello'))
        # 2. multiply with each other and with cprho0
        # some juggling with coordinates needed since Iris is very
        # restrictive in this regard
        cube.remove_coord('day_of_month')
        cube.remove_coord('day_of_year')
        cube.remove_coord('month_number')
        cube.remove_coord('year')
        t_coord = cube.coord('time')
        cube.remove_coord('time')
        new_cube = cube * volume
        new_cube *= RHO_CP
        new_cube.add_dim_coord(t_coord, 0)
        return new_cube
