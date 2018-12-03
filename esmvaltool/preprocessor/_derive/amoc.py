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
            Constraint(name='ocean_meridional_overturning_mass_streamfunction'))

        # 1: find the relevant region
        regions = cube.coord('region')
        atlantic_region = 'atlantic_arctic_ocean'
        atl_constraint = iris.Constraint(region=atlantic_region)
        amoc_cube = cube.extract(constraint=atl_constraint)

        # 2: Find the latitude closest to 26N
        rapid_location = 26.5
        lats = amoc_cube.coord('latitude').points
        rapid_index = np.argmin(np.abs(lats - rapid_location))
        rapid_constraint = iris.Constraint(latitude=lats[rapid_index])

        rapid_cube = amoc_cube.extract(constraint=rapid_constraint)

        # 3: find the maximum in the water column along the time axis.
        result_cube = rapid_cube.collapsed(['depth', 'region',], iris.analysis.MAX,)
        print(result_cube)
        print(result_cube.data)
        return result_cube
