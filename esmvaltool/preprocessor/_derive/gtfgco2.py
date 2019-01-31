"""Derivation of variable `gtfgco2`."""
import iris
from iris import Constraint

import numpy as np

from ._derived_variable_base import DerivedVariableBase


def calculate_total_flux(fgco2_cube, cube_area):
    """
    Calculate the area of unmasked cube cells.

    Requires a cube with two spacial dimensions. (no depth coordinate).

    Parameters
    ----------
    cube: iris.cube.Cube
        Data Cube
    cube_area: iris.cube.Cube
        Cell area Cube

    Returns
    -------
    numpy.array:
        An numpy array containing the total flux of CO2.

    """
    data = []
    times = fgco2_cube.coord('time')

    fgco2_cube.data = np.ma.array(fgco2_cube.data)
    for time_itr in np.arange(len(times.points)):

        total_flux = fgco2_cube[time_itr].data * cube_area.data

        total_flux = np.ma.masked_where(fgco2_cube[time_itr].data.mask,
                                        total_flux)
        data.append(total_flux.sum())

    ######
    # Create a small dummy output array
    data = np.array(data)
    return data


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `gtfgco2`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'fgco2',
            'field': 'TO2M', }],
        'fx_files': ['areacello', ]}

    def calculate(self, cubes):
        """Compute longwave cloud radiative effect."""
        fgco2_cube = cubes.extract_strict(
            Constraint(name='surface_downward_mass_flux_of_carbon_dioxide'
                            '_expressed_as_carbon'))

        try:
            cube_area = cubes.extract_strict(
                Constraint(name='cell_area'))
        except iris.exceptions.ConstraintMismatchError:
            pass

        total_flux = calculate_total_flux(fgco2_cube, cube_area)

        # Dummy result cube
        result = fgco2_cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN,)
        result.units = fgco2_cube.units * cube_area.units

        result.data = total_flux
        return result
