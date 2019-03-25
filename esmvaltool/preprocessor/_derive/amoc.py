"""Derivation of variable `amoc`."""
import iris
import numpy as np

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `amoc`."""

    # Required variables
    required = [{'short_name': 'msftmyz', 'mip': 'Omon'}]

    @staticmethod
    def calculate(cubes):
        """Compute Atlantic meriodinal overturning circulation.

        Arguments
        ---------
        cube: iris.cube.Cube
           input cube.

        Returns
        ---------
        iris.cube.Cube
              Output AMOC cube.
        """
        # 0. Load the msftmyz cube.
        cube = cubes.extract_strict(
            iris.Constraint(
                name='ocean_meridional_overturning_mass_streamfunction'))

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
        cube = cube.collapsed(
            ['depth', 'region'],
            iris.analysis.MAX,
        )
        return cube
