"""Derivation of variable `amoc`."""
import iris
from iris import Constraint

import numpy as np
import cf_units

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `amoc`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'msftmyz',
            'field': 'TO2M'
         }],
    }

    def calculate(self, cubes,):
        """Compute Atlantic meriodinal overturning circulation."""

        cube = cubes.extract_strict(
            Constraint(name='ocean_meridional_overturning_mass_'
                       'streamfunction'))

        # 1: find the relevant region
        atlantic_region = 'atlantic_arctic_ocean'
        atl_constraint = iris.Constraint(region=atlantic_region)
        cube = cube.extract(constraint=atl_constraint)

        # 2: Remove the shallowest 500m to avoid wind driven mixed layer.
        depth_constraint = iris.Constraint(depth=lambda d: d >= 500.)
        cube = cube.extract(constraint=depth_constraint)

        # 3: Find the latitude closest to 26N
        rapid_location = 26.5
        lats = cube.coord('latitude').points
        rapid_index = np.argmin(np.abs(lats - rapid_location))
        rapid_constraint = iris.Constraint(latitude=lats[rapid_index])
        cube = cube.extract(constraint=rapid_constraint)

        # 4: find the maximum in the water column along the time axis.
        cube = cube.collapsed(['depth', 'region'],
                              iris.analysis.MAX,)
        return cube
