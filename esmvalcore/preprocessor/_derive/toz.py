"""Derivation of variable `toz`."""

import cf_units
import iris
import numba
import numpy as np
from scipy import constants

from ._baseclass import DerivedVariableBase

# Constants
AVOGADRO_CONST = constants.value('Avogadro constant')
AVOGADRO_CONST_UNIT = constants.unit('Avogadro constant')
STANDARD_GRAVITY = 9.81
STANDARD_GRAVITY_UNIT = cf_units.Unit('m s^-2')
MW_AIR = 29
MW_AIR_UNIT = cf_units.Unit('g mol^-1')
MW_O3 = 48
MW_O3_UNIT = cf_units.Unit('g mol^-1')
DOBSON_UNIT = cf_units.Unit('2.69e20 m^-2')


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `toz`."""

    # Required variables
    required = [
        {
            'short_name': 'tro3'
        },
        {
            'short_name': 'ps'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute total column ozone.

        Note
        ----
        The surface pressure is used as a lower integration bound. A fixed
        upper integration bound of 0 Pa is used.

        """
        tro3_cube = cubes.extract_strict(
            iris.Constraint(name='mole_fraction_of_ozone_in_air'))
        ps_cube = cubes.extract_strict(
            iris.Constraint(name='surface_air_pressure'))

        p_layer_widths = _pressure_level_widths(
            tro3_cube, ps_cube, top_limit=0.0)
        toz_cube = (
            tro3_cube * p_layer_widths / STANDARD_GRAVITY * MW_O3 / MW_AIR)
        toz_cube = toz_cube.collapsed('air_pressure', iris.analysis.SUM)
        toz_cube.units = (tro3_cube.units * p_layer_widths.units /
                          STANDARD_GRAVITY_UNIT * MW_O3_UNIT / MW_AIR_UNIT)

        # Convert from kg m^-2 to Dobson unit (2.69e20 m^-2 )
        toz_cube = toz_cube / MW_O3 * AVOGADRO_CONST
        toz_cube.units = toz_cube.units / MW_O3_UNIT * AVOGADRO_CONST_UNIT
        toz_cube.convert_units(DOBSON_UNIT)
        toz_cube.data = np.ma.array(toz_cube.data, dtype=np.dtype('float32'))

        return toz_cube


# Helper functions
def _pressure_level_widths(tro3_cube, ps_cube, top_limit=0.0):
    """Create a cube with pressure level widths.

    This is done by taking a 2D surface pressure field as lower bound.

    Parameters
    ----------
        tro3_cube : iris.cube.Cube
            `Cube` containing `mole_fraction_of_ozone_in_air`.
        ps_cube : iris.cube.Cube
            `Cube` containing `surface_air_pressure`.
        top_limit : double
            Pressure in Pa.

    Returns
    -------
    iris.cube.Cube
    `Cube` of same shape as `tro3_cube` containing pressure level widths.

    """
    pressure_array = _create_pressure_array(tro3_cube, ps_cube, top_limit)

    data = _apply_pressure_level_widths(pressure_array)
    p_level_widths_cube = tro3_cube.copy(data=data)
    p_level_widths_cube.rename('pressure level widths')
    p_level_widths_cube.units = ps_cube.units

    return p_level_widths_cube


def _create_pressure_array(tro3_cube, ps_cube, top_limit):
    """Create an array filled with the `air_pressure` coord values.

    The array is created from the `tro3_cube` with the same dimensions
    as `tro3_cube`. This array is then sandwiched with a 2D array
    containing the surface pressure and a 2D array containing the top
    pressure limit.
    """
    # Create 4D array filled with pressure level values
    p_levels = tro3_cube.coord('air_pressure').points
    p_4d_array = iris.util.broadcast_to_shape(p_levels, tro3_cube.shape, [1])

    # Create 4d array filled with surface pressure values
    shape = tro3_cube.shape
    ps_4d_array = iris.util.broadcast_to_shape(ps_cube.data, shape, [0, 2, 3])

    # Set pressure levels below the surface pressure to NaN
    pressure_4d = np.where((ps_4d_array - p_4d_array) < 0, np.NaN, p_4d_array)

    # Make top_limit last pressure level
    top_limit_array = np.ones(ps_cube.shape) * top_limit
    data = top_limit_array[:, np.newaxis, :, :]
    pressure_4d = np.concatenate((pressure_4d, data), axis=1)

    # Make surface pressure the first pressure level
    data = ps_cube.data[:, np.newaxis, :, :]
    pressure_4d = np.concatenate((data, pressure_4d), axis=1)

    return pressure_4d


def _apply_pressure_level_widths(array, air_pressure_axis=1):
    """Compute pressure level widths.

    For a 1D array with pressure level columns, return a 1D array with
    pressure level widths.
    """
    return np.apply_along_axis(_p_level_widths, air_pressure_axis, array)


@numba.jit()  # ~10x faster
def _p_level_widths(array):
    """Create pressure level widths.

    The array with pressure levels is assumed to be monotonic and the
    values are decreasing.

    The first element is the lower boundary (surface pressure), the last
    value is the upper boundary. Thicknesses are only calculated for the
    values between these boundaries, the returned array, therefore,
    contains two elements less.

    >>> _p_level_widths(np.array([1020, 1000, 700, 500, 5]))
    array([170., 250., 595.])

    >>> _p_level_widths(np.array([990, np.NaN, 700, 500, 5]))
    array([  0., 390., 595.])
    """
    surface_pressure = array[0]
    top_limit = array[-1]
    array = array[1:-1]

    p_level_widths = np.ones(array.shape) * np.NAN

    last_pressure_level = len(array) - 1
    for i, val in enumerate(array):
        # numba would otherwise initialize it to 0 and
        # hide bugs that would occur in raw Python
        bounds_width = np.NAN
        if np.isnan(val):
            bounds_width = 0
        else:
            # Distance to lower bound
            if i == 0 or np.isnan(array[i - 1]):
                # First pressure level with value
                dist_to_lower_bound = surface_pressure - val
            else:
                dist_to_lower_bound = 0.5 * (array[i - 1] - val)

            # Distance to upper bound
            if i == last_pressure_level:  # last pressure level
                dist_to_upper_bound = val - top_limit
            else:
                dist_to_upper_bound = 0.5 * (val - array[i + 1])

            # Check monotonicity - all distances must be >= 0
            if dist_to_lower_bound < 0.0 or dist_to_upper_bound < 0.0:
                raise ValueError("Pressure level value increased with "
                                 "height.")

            bounds_width = dist_to_lower_bound + dist_to_upper_bound

        p_level_widths[i] = bounds_width
    return p_level_widths
