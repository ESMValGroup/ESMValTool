"""Derivation of variable `xco2`."""

import iris
from iris import Constraint
import numba
import numpy as np

from ._derived_variable_base import DerivedVariableBase

# Constants

# free air correction constant [s-2]
FAIR_COR = 3.0825958e-6
# Pi [1]
PI = 3.1415926536
# Molecular weight of the atmosphere [kg mol-1]
MW_AIR = 28.9644e-3
# Avogadro number [mol-1]
N_AVO = 6.022140857e23


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `xco2`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'co2',
            'field': 'T3{frequency}'
        }, {
            'short_name': 'hus',
            'field': 'T3{frequency}'
        }, {
            'short_name': 'zg',
            'field': 'T3{frequency}'
        }, {
            'short_name': 'ps',
            'field': 'T2{frequency}'
        }]
    }

    def calculate(self, cubes):
        """Calculate the average-column atmospheric CO2 [1e-6].

        The calculation follows the method described in the obs4MIPs
        technical note "Merged SCIMACHY/ENVISAT and TANSO-FTS/GOSAT
        atmospheric column-average dry-air mole fraction of CO2 (XCO2)"
        by M. Buchwitz and M. Reuter, available at:

        http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_obs4mips/
              esacci/ghg/data/obs4mips/crdp_3/CO2/v100/
              TechNote_SCIAGOSAT_L3_CRDP3_001_XCO2_FINAL.pdf

        """
        co2_cube = cubes.extract_strict(
            Constraint(name='mole_fraction_of_carbon_dioxide_in_air'))
        hus_cube = cubes.extract_strict(Constraint(name='specific_humidity'))
        zg_cube = cubes.extract_strict(Constraint(name='geopotential_height'))
        ps_cube = cubes.extract_strict(Constraint(name='surface_air_pressure'))

        # level thickness (note: Buchwitz & Reuter use hPa but we use Pa;
        # in fact, this does not matter as units cancel out when calculating
        # xco2)
        p_layer_widths = _pressure_level_widths(
            co2_cube, ps_cube, top_limit=0.0)

        #iris.save(p_layer_widths, '/pf/b/b380103/workesm/esmvaltool_output/p_layer_widths.nc')

        # latitudes (1-dim array)
        lat = co2_cube.coord('latitude').points

        # gravitational acceleration g0 on the geoid approximated by the
        # international gravity formula depending only on the latitude
        g0 = np.array(lat)
        g0 = 9.780327 * (1. + 0.0053024 * (np.sin(lat / 180. * PI))**2
             - 0.0000058 * (np.sin(2. * lat / 180. * PI))**2)

        # approximation of the gravitational acceleration including the
        # free air correction
        g_4d_array = iris.util.broadcast_to_shape(g0**2, zg_cube.shape, [2])
        g_4d_array = np.sqrt(g_4d_array.data + 2. * FAIR_COR * zg_cube.data)

        #outcube = co2_cube.copy(g_4d_array)
        #iris.save(outcube, '/pf/b/b380103/workesm/esmvaltool_output/g_4d_array.nc')

        # number of dry air particles (air molecules excluding water vapor)
        # within each layer
        n_dry = (hus_cube * -1. + 1.) * N_AVO * p_layer_widths.data / (MW_AIR * g_4d_array)

        #iris.save(n_dry, '/pf/b/b380103/workesm/esmvaltool_output/n_dry.nc')

        # number of CO2 molecules per layer
        co2_cube = co2_cube * n_dry

        # column-average CO2
        xco2_cube = (
            co2_cube.collapsed('air_pressure', iris.analysis.SUM) /
            n_dry.collapsed('air_pressure', iris.analysis.SUM))
        xco2_cube.units = co2_cube.units

        #iris.save(xco2_cube, '/pf/b/b380103/workesm/esmvaltool_output/xco2_new.nc')

        return xco2_cube

# Helper functions
def _pressure_level_widths(co2_cube, ps_cube, top_limit=0.0):
    """Create a cube with pressure level widths.

    This is done by taking a 2D surface pressure field as lower bound.

    Parameters
    ----------
        co2_cube : iris.cube.Cube
            `Cube` containing `mole_fraction_of_carbon_dioxide_in_air`.
        ps_cube : iris.cube.Cube
            `Cube` containing `surface_air_pressure`.
        top_limit : double
            Pressure in Pa.

    Returns
    -------
    iris.cube.Cube
    `Cube` of same shape as `co2_cube` containing pressure level widths.

    """
    pressure_array = _create_pressure_array(co2_cube, ps_cube, top_limit)

    data = _apply_pressure_level_widths(pressure_array)
    p_level_widths_cube = co2_cube.copy(data=data)
    p_level_widths_cube.rename('pressure level widths')
    p_level_widths_cube.units = ps_cube.units

    return p_level_widths_cube


def _create_pressure_array(co2_cube, ps_cube, top_limit):
    """Create an array filled with the `air_pressure` coord values.

    The array is created from the `co2_cube` with the same dimensions
    as `co2_cube`. This array is then sandwiched with a 2D array
    containing the surface pressure and a 2D array containing the top
    pressure limit.
    """
    # Create 4D array filled with pressure level values
    p_levels = co2_cube.coord('air_pressure').points
    p_4d_array = iris.util.broadcast_to_shape(p_levels, co2_cube.shape, [1])

    # Create 4d array filled with surface pressure values
    shape = co2_cube.shape
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
