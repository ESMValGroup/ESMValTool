"""
Miscellaneous functions for deriving variables.

"""


import cf_units
import numba
import numpy as np
import iris
from iris import Constraint
from scipy import constants


Avogadro_const = constants.value('Avogadro constant')
Avogadro_const_unit = constants.unit('Avogadro constant')
g = 9.81
g_unit = cf_units.Unit('m s^-2')
mw_air = 29
mw_air_unit = cf_units.Unit('g mol^-1')
mw_O3 = 48
mw_O3_unit = cf_units.Unit('g mol^-1')
Dobson_unit = cf_units.Unit('2.69e20 m^-2')


def calc_lwp(cubes):
    """
    Liquid water path is calculated by subtracting clivi (ice water) from clwvi
    (condensed water path).
    Note: Some models output the variable "clwvi" which only contains lwp. In
    these cases, the input clwvi cube is just returned.

    Args:
        * cubes: cubelist containing clwvi_cube and clivi_cube

    Returns:
        Cube containing liquid water path.
    """
    clwvi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_mass_content_of_cloud_condensed_water'))
    clivi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_mass_content_of_cloud_ice'))

    model = clwvi_cube.attributes['model_id']
    project = clwvi_cube.attributes['project_id']
    # Should we check that the model/project_id are the same on both cubes?

    BAD_MODELS = ['CESM1-CAM5-1-FV2', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM',
                  'CMCC-CMS', 'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR',
                  'CCSM4', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM',
                  'MIROC-ESM', 'CSIRO-Mk3-6-0', 'MPI-ESM-MR', 'MPI-ESM-LR',
                  'MPI-ESM-P']
    if (project in ["CMIP5", "CMIP5_ETHZ"] and model in BAD_MODELS) or \
        (project == 'OBS' and model == 'UWisc'):
            print("INFO: assuming that variable clwvi from {} model {} "
                  "contains only liquid water".format(project, model))
            lwp_cube = clwvi_cube
    else:
        lwp_cube = clwvi_cube - clivi_cube

    # TODO: Rename cube lwp_cube.name('liquid_water_path') here?
    # TODO: Fix units? lwp_cube.units = cf_units.Unit('kg') here?
    return lwp_cube


def calc_toz(cubes):
    """
    Create total column ozone from ozone mol fraction on pressure levels.

    The surface pressure is used as a lower integration bound. A fixed upper
    integration bound of 100 Pa is used.

    Args:
        * cubes: cubelist containing tro3_cube (mole_fraction_of_ozone_in_air)
                 and ps_cube (surface_air_pressure).

    Returns:
        Cube containing total column ozone.
    """
    tro3_cube = cubes.extract_strict(
        Constraint(name='mole_fraction_of_ozone_in_air'))
    ps_cube = cubes.extract_strict(Constraint(name='surface_air_pressure'))

    assert tro3_cube.coord_dims('time') and ps_cube.coord_dims('time'), \
            'No time dimension found.'
    p_layer_widths = _pressure_level_widths(tro3_cube, ps_cube, top_limit=100)
    toz = (tro3_cube * p_layer_widths / g * mw_O3 / mw_air).collapsed('air_pressure', iris.analysis.SUM)
    toz.units = tro3_cube.units * p_layer_widths.units / g_unit * mw_O3_unit / mw_air_unit
    toz.rename('atmosphere mass content of ozone')

    # Convert from kg m^-2 to Dobson unit (2.69e20 m^-2 )
    toz = toz / mw_O3 * Avogadro_const
    toz.units = toz.units / mw_O3_unit * Avogadro_const_unit
    toz.convert_units(Dobson_unit)
    return toz


def _pressure_level_widths(tro3_cube, ps_cube, top_limit=100):
    """
    Create a cube with pressure level widths taking a 2D surface pressure field
    as lower bound.

    Args:
        tro3_cube: Cube containing mole_fraction_of_ozone_in_air
        ps_cube: Surface air pressure cube.
        top_limit: Pressure in Pa.

    returns: Cube of same shape as tro3_cube containing pressure level widths.
    """
    assert ps_cube.units == 'Pa'
    assert tro3_cube.coord('air_pressure').units == 'Pa'

    pressure_array = _create_pressure_array(tro3_cube, ps_cube, top_limit)

    p_level_widths_cube = tro3_cube.copy(data=_apply_pressure_level_widths(pressure_array))
    p_level_widths_cube.rename('pressure level widths')
    p_level_widths_cube.units = ps_cube.units

    return p_level_widths_cube


def _create_pressure_array(tro3_cube, ps_cube, top_limit):
    """
    Create an array filled with the 'air_pressure' coord values from the
    tro3_cube with the same dimensions as tro3_cube.
    This array is then sandwiched with a 2D array containing the surface
    pressure, and a 2D array containing the top pressure limit.
    """
    # create 4D array filled with pressure level values
    p_levels = tro3_cube.coord('air_pressure').points
    p_4D_array = iris.util.broadcast_to_shape(p_levels, tro3_cube.shape, [1])
    assert p_4D_array.shape == tro3_cube.shape

    # create 4d array filled with surface pressure values
    ps_4D_array = iris.util.broadcast_to_shape(ps_cube.data, tro3_cube.shape, [0, 2, 3])
    assert ps_4D_array.shape == tro3_cube.shape

    # set pressure levels below the surface pressure to NaN
    pressure_4d = np.where((ps_4D_array - p_4D_array) < 0, np.NaN, p_4D_array)

    # make top_limit last pressure level
    top_limit_array = np.ones(ps_cube.shape) * top_limit
    pressure_4d = np.concatenate((pressure_4d, top_limit_array[:, np.newaxis, :, :]), axis=1)
    assert (pressure_4d[:, -1, :, :] == top_limit).all()

    # make surface pressure the first pressure level
    pressure_4d = np.concatenate((ps_cube.data[:, np.newaxis, :, :], pressure_4d, ), axis=1)
    assert (pressure_4d[:, 0, :, :] == ps_cube.data).all()

    return pressure_4d


def _apply_pressure_level_widths(array, air_pressure_axis=1):
    """
    For a  1D array with pressure level columns, return a 1D  array with pressure
    level widths.
    """
    return np.apply_along_axis(_p_level_widths, air_pressure_axis, array)


@numba.jit()  # ~10x faster
def _p_level_widths(array):
    """
    Creates pressure level widths from an array with pressure level values.
    The array is assumed to be monotonic and the values are decreasing.

    The first element is the lower boundary (surface pressure), the last value
    is the upper boundary. Thicknesses are only calculated for the values
    between these boundaries, the returned array, therefore, contains two
    elements less.

    >>> _p_level_widths(np.array([1020, 1000, 700, 500, 5]))
    [170 250 595]


    >>> _p_level_widths(np.array([990, np.NaN, 700, 500, 5]))
    [0 390 595]
    """
    surface_pressure = array[0]
    top_limit = array[-1]
    array = array[1:-1]

    p_level_widths = np.ones(array.shape) * np.NAN

    last_pressure_level = len(array) - 1
    for i, val in enumerate(array):
        bounds_width = np.NAN  # numba would otherwise initialise it to 0 and
                               # hide bugs that would occur in raw Python
        if np.isnan(val):
            bounds_width = 0
        else:
            # distance to lower bound
            if i == 0 or np.isnan(array[i-1]):  # first pressure level with value
                dist_to_lower_bound = surface_pressure - val
            else:
                dist_to_lower_bound = 0.5*(array[i-1] - val)

            # distance to upper bound
            if i == last_pressure_level:        # last pressure level
                dist_to_upper_bound = val - top_limit
            else:
                dist_to_upper_bound = 0.5*(val - array[i+1])

            # Check monotonicity - all distances must be >= 0
            if dist_to_lower_bound < 0.0 or dist_to_upper_bound < 0.0:
                raise ValueError('Pressure level value increased with height.')

            bounds_width = dist_to_lower_bound + dist_to_upper_bound

        p_level_widths[i] = bounds_width
    return p_level_widths

