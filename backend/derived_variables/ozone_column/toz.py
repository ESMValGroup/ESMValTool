import cf_units
import numpy as np
import iris


g = 9.81
g_unit = cf_units.Unit('m s^-2')
mw_air = 29
mw_air_unit = cf_units.Unit('g mol^-1')
mw_O3 = 48
mw_O3_unit = cf_units.Unit('g mol^-1')


def total_column_ozone(tro3_cube, ps_cube):
    """
    Create total column ozone from ozone mol fraction on pressure levels.

    The surface pressure is used as a lower integration bound. A fixed upper
    integration bound of 100 Pa is used.

    tro3_cube: Cube containing mole_fraction_of_ozone_in_air
    ps_cube: Surface pressure cube.

    returns: Cube containing total column ozone.
    """
    p_layer_widths = pressure_layer_widths(tro3_cube, ps_cube, top_limit=100)
    toz = (tro3_cube * p_layer_widths / g * mw_O3 / mw_air).collapsed('air_pressure', iris.analysis.SUM)
    #toz.units = tro3_cube.units * p_layer_widths.units / g_unit * mw_O3.units / mw_air_unit
    # TODO check toz variable name
    # TODO check toz target unit
    return toz


def pressure_layer_widths(tro3_cube, ps_cube, top_limit=100):
    """
    Create a cube with pressure level thicknesses.
    Take a 2D surface pressure field as lower bound.

    tro3_cube: Cube containing mole_fraction_of_ozone_in_air
    ps_cube: Surface pressure cube.
    top_limit: Pressure in Pa.

    returns:
    """
    # TODO: assert all pressures in Pa
    pressure_cube = create_pressure_cube(tro3_cube, ps_cube, top_limit)
    return apply_pressure_level_widths(pressure_cube)


def create_pressure_cube(tro3_cube, ps_cube, top_limit):
    """
    Create a cube filled with the the 'air_pressure' coord values of the
    tro3_cube of the same dimensions as tro3_cube.
    This array is then sandwiched with a 2D array containing the surface
    pressure, and a 2D array containing the top pressure limit.
    """
    p = tro3_cube.coord('air_pressure').points
    pressure_4D_array = iris.util.broadcast_to_shape(p,            tro3_cube.shape, [1])
    ps_4D_array       = iris.util.broadcast_to_shape(ps_cube.data, tro3_cube.shape, [0, 2, 3])

    pressure_4d = np.where((ps_4D_array - pressure_4D_array) < 0, np.NaN, pressure_4D_array)

    # make top_limit last pressure level
    top_limit_array = np.ones(ps_cube.shape) * top_limit
    pressure_4d = np.concatenate((pressure_4d, top_limit_array[:, np.newaxis, :, :]), axis=1)

    # make surface pressure the first pressure level
    pressure_4d = np.concatenate((ps_cube.data[:, np.newaxis,: , :], pressure_4d, ), axis=1)

    # TODO return iris cube
    return pressure_4d


def apply_pressure_level_widths(array):
    """
    For an array with pressure level columns, return an array with pressure
    level widths.
    """
    # TODO assert monotonicity
    # TODO axis number -> work with cube here, then # cube.coord_dims('air_pressure')
    return np.apply_along_axis(pressure_level_widths, 1, array)


def pressure_level_widths(array):
    """
    Creates pressure level widths from an array with pressure level values.
    The array is assumed to be monotonic and the values are decreasing.

    The first element is the lower boundary (surface pressure), the last value
    is the upper boundary. Thicknesses are only calculated for the values
    between these boundaries, the returned array, therefore, contains two
    elements less.

    >>> pressure_level_widths(np.array([1020, 1000, 700, 500, 5]))
    [170 250 595]


    >>> pressure_level_widths(np.array([990, np.NaN, 700, 500, 5]))
    [0 390 595]
    """
    surface_pressure = array[0]
    top_limit = array[-1]
    array = array[1:-1]
    assert np.min(np.ma.fix_invalid(array, fill_value=np.Inf)) >= top_limit, \
            "Lowest value is below top_limit."

    pressure_level_widths = []
    num_array = len(array)
    for i, val in enumerate(array):
        if np.isnan(val):
            bounds_width = 0
        else:
            if i == 0 or np.isnan(array[i-1]):  # first pressure level with value
                dist_to_lower_bound = surface_pressure - val
            else:
                dist_to_lower_bound = 0.5*(array[i-1] - val)

            if i == num_array - 1:              # last pressure level
                dist_to_upper_bound = val - top_limit
            else:
                dist_to_upper_bound = 0.5*(val - array[i+1])

            bounds_width = dist_to_lower_bound + dist_to_upper_bound
        pressure_level_widths.append(bounds_width)
    return pressure_level_widths

