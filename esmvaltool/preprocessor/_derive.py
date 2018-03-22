"""Miscellaneous functions for deriving variables."""

import logging

import cf_units
import iris
import numba
import numpy as np
from iris import Constraint
from scipy import constants

logger = logging.getLogger(__name__)

Avogadro_const = constants.value('Avogadro constant')
Avogadro_const_unit = constants.unit('Avogadro constant')
g = 9.81
g_unit = cf_units.Unit('m s^-2')
mw_air = 29
mw_air_unit = cf_units.Unit('g mol^-1')
mw_O3 = 48
mw_O3_unit = cf_units.Unit('g mol^-1')
Dobson_unit = cf_units.Unit('2.69e20 m^-2')


def get_required(short_name, field=None):
    """Get variable short_name and field pairs required to derive variable"""
    frequency = field[2] if field else 'M'
    required = {
        'lwp': [
            ('clwvi', 'T2' + frequency + 's'),
            ('clivi', 'T2' + frequency + 's'),
        ],
        'toz': [
            ('tro3', 'T3' + frequency),
            ('ps', 'T2' + frequency + 's'),
        ],
    }

    if short_name in required:
        return required[short_name]

    raise NotImplementedError("Don't know how to derive {}".format(short_name))


def derive(cubes, short_name):
    """Derive variable `short_name`"""
    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Derive
    functions = {
        'lwp': calc_lwp,
        'toz': calc_toz,
    }
    if short_name in functions:
        cubes = iris.cube.CubeList(cubes)
        cube = functions[short_name](cubes)
        cube.attributes['_filename'] = cubes[0].attributes['_filename']
        return cube

    raise NotImplementedError("Don't know how to derive {}".format(short_name))


def calc_lwp(cubes):
    """Compute liquid water path.

    Liquid water path is calculated by subtracting clivi (ice water) from clwvi
    (condensed water path).
    Note: Some models output the variable "clwvi" which only contains lwp. In
    these cases, the input clwvi cube is just returned.

    Arguments
    ---------
        cubes: cubelist containing clwvi_cube and clivi_cube

    Returns
    -------
        Cube containing liquid water path.

    """
    clwvi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_mass_content_of_cloud_condensed_water'))
    clivi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_mass_content_of_cloud_ice'))

    model = clwvi_cube.attributes.get('model_id')
    project = clwvi_cube.attributes.get('project_id')
    # Should we check that the model/project_id are the same on both cubes?

    bad_models = [
        'CESM1-CAM5-1-FV2', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS',
        'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'CCSM4',
        'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC-ESM',
        'CSIRO-Mk3-6-0', 'MPI-ESM-MR', 'MPI-ESM-LR', 'MPI-ESM-P'
    ]
    if ((project in ["CMIP5", "CMIP5_ETHZ"] and model in bad_models)
            or (project == 'OBS' and model == 'UWisc')):
        logger.info("Assuming that variable clwvi from %s model %s "
                    "contains only liquid water", project, model)
        lwp_cube = clwvi_cube
    else:
        lwp_cube = clwvi_cube - clivi_cube

    # TODO: Rename cube lwp_cube.name('liquid_water_path') here?
    # TODO: Fix units? lwp_cube.units = cf_units.Unit('kg') here?
    return lwp_cube


def calc_toz(cubes):
    """Compute total column ozone from ozone mol fraction on pressure levels.

    The surface pressure is used as a lower integration bound. A fixed upper
    integration bound of 100 Pa is used.

    Arguments
    ----
        cubes: cubelist containing tro3_cube (mole_fraction_of_ozone_in_air)
               and ps_cube (surface_air_pressure).

    Returns
    -------
        Cube containing total column ozone.

    """
    tro3_cube = cubes.extract_strict(
        Constraint(name='mole_fraction_of_ozone_in_air'))
    ps_cube = cubes.extract_strict(Constraint(name='surface_air_pressure'))

    assert tro3_cube.coord_dims('time') and ps_cube.coord_dims('time'), \
        'No time dimension found.'
    p_layer_widths = _pressure_level_widths(tro3_cube, ps_cube, top_limit=100)
    toz = tro3_cube * p_layer_widths / g * mw_O3 / mw_air
    toz = toz.collapsed('air_pressure', iris.analysis.SUM)
    toz.units = (tro3_cube.units * p_layer_widths.units / g_unit * mw_O3_unit /
                 mw_air_unit)

    # Convert from kg m^-2 to Dobson unit (2.69e20 m^-2 )
    toz = toz / mw_O3 * Avogadro_const
    toz.units = toz.units / mw_O3_unit * Avogadro_const_unit
    toz.convert_units(Dobson_unit)
    toz.units = cf_units.Unit('DU')

    # Set names
    toz.var_name = 'toz'
    toz.standard_name = (
        'equivalent_thickness_at_stp_of_atmosphere_ozone_content')
    toz.long_name = ' Total Ozone Column'
    return toz


def _pressure_level_widths(tro3_cube, ps_cube, top_limit=100):
    """Create a cube with pressure level widths.

    This is done by taking a 2D surface pressure field as lower bound.

    Arguments
    ---------
        tro3_cube: Cube containing mole_fraction_of_ozone_in_air
        ps_cube: Surface air pressure cube.
        top_limit: Pressure in Pa.

    Returns
    -------
        Cube of same shape as tro3_cube containing pressure level widths.

    """
    assert ps_cube.units == 'Pa'
    assert tro3_cube.coord('air_pressure').units == 'Pa'

    pressure_array = _create_pressure_array(tro3_cube, ps_cube, top_limit)

    data = _apply_pressure_level_widths(pressure_array)
    p_level_widths_cube = tro3_cube.copy(data=data)
    p_level_widths_cube.rename('pressure level widths')
    p_level_widths_cube.units = ps_cube.units

    return p_level_widths_cube


def _create_pressure_array(tro3_cube, ps_cube, top_limit):
    """Create an array filled with the 'air_pressure' coord values.

    The array is created from the tro3_cube with the same dimensions
    as tro3_cube. This array is then sandwiched with a 2D array containing
    the surface pressure, and a 2D array containing the top pressure limit.
    """
    # create 4D array filled with pressure level values
    p_levels = tro3_cube.coord('air_pressure').points
    p_4d_array = iris.util.broadcast_to_shape(p_levels, tro3_cube.shape, [1])
    assert p_4d_array.shape == tro3_cube.shape

    # create 4d array filled with surface pressure values
    shape = tro3_cube.shape
    ps_4d_array = iris.util.broadcast_to_shape(ps_cube.data, shape, [0, 2, 3])
    assert ps_4d_array.shape == tro3_cube.shape

    # set pressure levels below the surface pressure to NaN
    pressure_4d = np.where((ps_4d_array - p_4d_array) < 0, np.NaN, p_4d_array)

    # make top_limit last pressure level
    top_limit_array = np.ones(ps_cube.shape) * top_limit
    data = top_limit_array[:, np.newaxis, :, :]
    pressure_4d = np.concatenate((pressure_4d, data), axis=1)
    assert (pressure_4d[:, -1, :, :] == top_limit).all()

    # make surface pressure the first pressure level
    data = ps_cube.data[:, np.newaxis, :, :]
    pressure_4d = np.concatenate((data, pressure_4d), axis=1)
    assert (pressure_4d[:, 0, :, :] == ps_cube.data).all()

    return pressure_4d


def _apply_pressure_level_widths(array, air_pressure_axis=1):
    """Compute pressure level widths.

    For a  1D array with pressure level columns, return a 1D  array with
    pressure level widths.
    """
    return np.apply_along_axis(_p_level_widths, air_pressure_axis, array)


@numba.jit()  # ~10x faster
def _p_level_widths(array):
    """Create pressure level widths from an array with pressure level values.

    The array is assumed to be monotonic and the values are decreasing.

    The first element is the lower boundary (surface pressure), the last value
    is the upper boundary. Thicknesses are only calculated for the values
    between these boundaries, the returned array, therefore, contains two
    elements less.

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
        # numba would otherwise initialise it to 0 and
        # hide bugs that would occur in raw Python
        bounds_width = np.NAN
        if np.isnan(val):
            bounds_width = 0
        else:
            # distance to lower bound
            if i == 0 or np.isnan(array[i - 1]):
                # first pressure level with value
                dist_to_lower_bound = surface_pressure - val
            else:
                dist_to_lower_bound = 0.5 * (array[i - 1] - val)

            # distance to upper bound
            if i == last_pressure_level:  # last pressure level
                dist_to_upper_bound = val - top_limit
            else:
                dist_to_upper_bound = 0.5 * (val - array[i + 1])

            # Check monotonicity - all distances must be >= 0
            if dist_to_lower_bound < 0.0 or dist_to_upper_bound < 0.0:
                raise ValueError('Pressure level value increased with height.')

            bounds_width = dist_to_lower_bound + dist_to_upper_bound

        p_level_widths[i] = bounds_width
    return p_level_widths
