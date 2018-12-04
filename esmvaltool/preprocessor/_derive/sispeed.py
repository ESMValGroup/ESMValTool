"""Derivation of variable `sispeed`."""

import numpy as np
from iris import Constraint
import iris.cube
from iris.coords import DimCoord

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `sispeed`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'usi',
            'field': 'T2{frequency}'
        }, {
            'short_name': 'vsi',
            'field': 'T2{frequency}'
        }]
    }

    def calculate(self, cubes):
        """
        Compute sispeed module from velocity components siu and siv.

        Arguments
        ----
            cubes: cubelist containing velocity components.

        Returns
        -------
            Cube containing sea ice speed.

        """
        siu = cubes.extract_strict(Constraint(name='sea_ice_x_velocity'))
        siu.remove_coord('year')
        siu.remove_coord('day_of_year')

        siv = cubes.extract_strict(Constraint(name='sea_ice_y_velocity'))
        siv.remove_coord('latitude')
        siv.remove_coord('longitude')
        siv.remove_coord('year')
        siv.remove_coord('day_of_year')

        if isinstance(siu.coord('latitude'), DimCoord):
            siv.add_dim_coord(
                siu.coord('latitude'), siu.coord_dims('latitude')
            )
            siv.add_dim_coord(
                siu.coord('longitude'), siu.coord_dims('longitude')
            )
        else:
            siv.add_aux_coord(
                siu.coord('latitude'), siu.coord_dims('latitude')
            )
            siv.add_aux_coord(
                siu.coord('longitude'), siu.coord_dims('longitude')
            )

        speed = iris.cube.CubeList()
        for time in siu.coord('time').points:
            # Casting to float64 to avoid overflow errors
            siu_slice = siu.extract(Constraint(time=time))
            siu_slice.data = siu_slice.data.astype(np.float64)
            siv_slice = siu.extract(Constraint(time=time))
            siv_slice.data = siv_slice.data.astype(np.float64)
            speed_slice = (siu_slice ** 2 + siv_slice ** 2) ** 0.5
            del siu_slice
            del siv_slice
            # 64 bit precission no longer needed, cast to 32 to save memory
            speed_slice.data = speed_slice.data.astype(np.float32)
            speed.append(speed_slice)
        del siu
        del siv

        speed = speed.merge_cube()
        speed.short_name = 'sispeed'
        speed.standard_name = 'sea_ice_speed'
        speed.long_name = 'Sea-ice speed'
        return speed
