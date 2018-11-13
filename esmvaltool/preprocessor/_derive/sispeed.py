"""Derivation of variable `sispeed`."""

from iris import Constraint
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
        siv = cubes.extract_strict(Constraint(name='sea_ice_y_velocity'))
        siv.remove_coord('latitude')
        siv.remove_coord('longitude')
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
        speed = ((siu ** 2 + siv ** 2) ** 0.5)
        del siu
        del siv

        speed.short_name = 'sispeed'
        speed.standard_name = 'sea_ice_speed'
        speed.long_name = 'Sea-ice speed'
        speed.convert_units('km day-1')

        return speed
