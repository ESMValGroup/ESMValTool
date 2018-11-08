"""Miscellaneous functions for deriving variables."""

import logging

import cf_units
import iris
import numba
import numpy as np
import yaml
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
        'lwcre': [
            ('rlut', 'T2' + frequency + 's'),
            ('rlutcs', 'T2' + frequency + 's'),
        ],
        'lwp': [
            ('clwvi', 'T2' + frequency + 's'),
            ('clivi', 'T2' + frequency + 's'),
        ],
        'netcre': [
            ('rlut', 'T2' + frequency + 's'),
            ('rlutcs', 'T2' + frequency + 's'),
            ('rsut', 'T2' + frequency + 's'),
            ('rsutcs', 'T2' + frequency + 's'),
        ],
        'swcre': [
            ('rsut', 'T2' + frequency + 's'),
            ('rsutcs', 'T2' + frequency + 's'),
        ],
        'toz': [
            ('tro3', 'T3' + frequency),
            ('ps', 'T2' + frequency + 's'),
        ],
        'rtnt': [('rsdt', 'T2' + frequency + 's'),
                 ('rsut', 'T2' + frequency + 's'), ('rlut',
                                                    'T2' + frequency + 's')],
        'rsnt': [
            ('rsdt', 'T2' + frequency + 's'),
            ('rsut', 'T2' + frequency + 's'),
        ],
        'rsns': [
            ('rsds', 'T2' + frequency + 's'),
            ('rsus', 'T2' + frequency + 's'),
        ],
        'rlns': [
            ('rlds', 'T2' + frequency + 's'),
            ('rlus', 'T2' + frequency + 's'),
        ],
        'cllmtisccp': [('clisccp', 'T4' + frequency)],
        'clltkisccp': [('clisccp', 'T4' + frequency)],
        'clmmtisccp': [('clisccp', 'T4' + frequency)],
        'clmtkisccp': [('clisccp', 'T4' + frequency)],
        'clhmtisccp': [('clisccp', 'T4' + frequency)],
        'clhtkisccp': [('clisccp', 'T4' + frequency)]
    }

    if short_name in required:
        return required[short_name]

    raise NotImplementedError("Don't know how to derive {}".format(short_name))


def derive(cubes, variable):
    """Derive variable"""
    short_name = variable['short_name']
    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Available derivation functions
    functions = {
        'lwcre': calc_lwcre,
        'lwp': calc_lwp,
        'netcre': calc_netcre,
        'swcre': calc_swcre,
        'toz': calc_toz,
        'rtnt': calc_rtnt,
        'rsnt': calc_rsnt,
        'rsns': calc_rsns,
        'rlns': calc_rlns,
        'cllmtisccp': calc_cllmtisccp,
        'clltkisccp': calc_clltkisccp,
        'clmmtisccp': calc_clmmtisccp,
        'clmtkisccp': calc_clmtkisccp,
        'clhmtisccp': calc_clhmtisccp,
        'clhtkisccp': calc_clhtkisccp
    }

    if short_name not in functions:
        raise NotImplementedError(
            "Don't know how to derive {}".format(short_name))

    # Preprare input cubes and derive
    cubes = iris.cube.CubeList(cubes)
    cube = functions[short_name](cubes)

    # Set standard attributes
    cube.var_name = short_name
    if variable['standard_name'] not in iris.std_names.STD_NAMES:
        iris.std_names.STD_NAMES[variable['standard_name']] = {
            'canonical_units': variable['units']
        }
    for attribute in ('standard_name', 'long_name', 'units'):
        setattr(cube, attribute, variable[attribute])

    # Set attributes required by preprocessor
    cube.attributes['_filename'] = variable['filename']
    cube.attributes['metadata'] = yaml.safe_dump(variable)

    return cube


def calc_lwcre(cubes):
    """Compute longwave cloud radiative effect from all-sky and clear-sky flux.

    Arguments
    ----
        cubes: cubelist containing rlut (toa_outgoing_longwave_flux) and rlutcs
               (toa_outgoing_longwave_flux_assuming_clear_sky).

    Returns
    -------
        Cube containing longwave cloud radiative effect.

    """
    rlut_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_longwave_flux'))
    rlutcs_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

    lwcre = rlutcs_cube - rlut_cube
    lwcre.units = rlut_cube.units

    return lwcre


def calc_lwp(cubes):
    """Compute liquid water path.

    Liquid water path is calculated by subtracting clivi (ice water) from clwvi
    (condensed water path).
    Note: Some datasets output the variable "clwvi" which only contains lwp. In
    these cases, the input clwvi cube is just returned.

    Arguments
    ---------
        cubes: cubelist containing clwvi_cube and clivi_cube

    Returns
    -------
        Cube containing liquid water path.

    """
    clwvi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_cloud_condensed_water_content'))
    clivi_cube = cubes.extract_strict(
        Constraint(name='atmosphere_cloud_ice_content'))

    dataset = clwvi_cube.attributes.get('model_id')
    project = clwvi_cube.attributes.get('project_id')
    # Should we check that the model_id/project_id are the same on both cubes?

    bad_datasets = [
        'CESM1-CAM5-1-FV2', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS',
        'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'CCSM4',
        'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC-ESM',
        'CSIRO-Mk3-6-0', 'MPI-ESM-MR', 'MPI-ESM-LR', 'MPI-ESM-P'
    ]
    if ((project in ["CMIP5", "CMIP5_ETHZ"] and dataset in bad_datasets)
            or (project == 'OBS' and dataset == 'UWisc')):
        logger.info(
            "Assuming that variable clwvi from %s dataset %s "
            "contains only liquid water", project, dataset)
        lwp_cube = clwvi_cube
    else:
        lwp_cube = clwvi_cube - clivi_cube

    return lwp_cube


def calc_netcre(cubes):
    """Compute net cloud radiative effect.

       Calculate net cre as sum of longwave and shortwave cloud
       radiative effects.

    Arguments
    ----
        cubes: cubelist containing rlut (toa_outgoing_longwave_flux), rlutcs
               (toa_outgoing_longwave_flux_assuming_clear_sky), rsut
               (toa_outgoing_shortwave_flux) and rsutcs
               (toa_outgoing_shortwave_flux_assuming_clear_sky).

    Returns
    -------
        Cube containing net cloud radiative effect.

    """
    lwcre = calc_lwcre(cubes)
    swcre = calc_swcre(cubes)

    netcre = lwcre + swcre
    netcre.units = lwcre.units

    return netcre


def calc_swcre(cubes):
    """Compute shortwave cloud radiative effect from all-sky and clear-sky

       flux.

    Arguments
    ----
        cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
               rsutcs (toa_outgoing_shortwave_flux_assuming_clear_sky).

    Returns
    -------
        Cube containing shortwave cloud radiative effect.

    """
    rsut_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_shortwave_flux'))
    rsutcs_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

    swcre = rsutcs_cube - rsut_cube

    return swcre


def calc_toz(cubes):
    """Compute total column ozone from ozone mol fraction on pressure levels.

    The surface pressure is used as a lower integration bound. A fixed upper
    integration bound of 0 Pa is used.

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

    p_layer_widths = _pressure_level_widths(tro3_cube, ps_cube, top_limit=0)
    toz = tro3_cube * p_layer_widths / g * mw_O3 / mw_air
    toz = toz.collapsed('air_pressure', iris.analysis.SUM)
    toz.units = (tro3_cube.units * p_layer_widths.units / g_unit * mw_O3_unit /
                 mw_air_unit)

    # Convert from kg m^-2 to Dobson unit (2.69e20 m^-2 )
    toz = toz / mw_O3 * Avogadro_const
    toz.units = toz.units / mw_O3_unit * Avogadro_const_unit
    toz.convert_units(Dobson_unit)
    toz.data = np.ma.array(toz.data, dtype=np.dtype('float32'))

    return toz


def calc_rtnt(cubes):
    """Compute rtnt: TOA Net downward Total Radiation.

    Arguments
    ----
        cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
               rsdt (toa_incoming_shortwave_flux) and
               rlut (toa_outgoing_longwave_flux).

    Returns
    -------
        Cube containing TOA Net downward Total Radiation.
        Units: W m-2

    """
    rsdt_cube = cubes.extract_strict(
        Constraint(name='toa_incoming_shortwave_flux'))
    rsut_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_shortwave_flux'))
    rlut_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_longwave_flux'))

    # rtnt = (rsdt - rsut) - rlut
    rtnt = rsdt_cube - rsut_cube - rlut_cube

    return rtnt


def calc_rsnt(cubes):
    """Compute rsnt: TOA Net downward Shortwave Radiation.

    Arguments
    ----
        cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
               rsdt (toa_incoming_shortwave_flux).

    Returns
    -------
        Cube containing TOA Net downward Shortwave Radiation.
        Units: W m-2

    """
    rsdt_cube = cubes.extract_strict(
        Constraint(name='toa_incoming_shortwave_flux'))
    rsut_cube = cubes.extract_strict(
        Constraint(name='toa_outgoing_shortwave_flux'))

    # rsnt = rsdt - rsut
    rsnt = rsdt_cube - rsut_cube

    return rsnt


def calc_rsns(cubes):
    """Compute rsns: Surface Net downward Shortwave Radiation.

    Arguments
    ----
        cubes: cubelist containing
               rsus (surface_upwelling_shortwave_flux_in_air) and
               rsds (surface_downwelling_shortwave_flux_in_air).

    Returns
    -------
        Cube containing Surface Net downward Shortwave Radiation.
        Units: W m-2

    """
    rsds_cube = cubes.extract_strict(
        Constraint(name='surface_downwelling_shortwave_flux_in_air'))
    rsus_cube = cubes.extract_strict(
        Constraint(name='surface_upwelling_shortwave_flux_in_air'))

    # rsns = rsds - rsus
    rsns = rsds_cube - rsus_cube

    return rsns


def calc_rlns(cubes):
    """Compute rlns: Surface Net downward Longwave Radiation.

    Arguments
    ----
        cubes: cubelist containing
               rlds (surface_downwelling_longwave_flux_in_air) and
               rlus (surface_upwelling_longwave_flux_in_air).

    Returns
    -------
        Cube containing Surface Net downward Longwave Radiation.
        Units: W m-2

    """
    rlds_cube = cubes.extract_strict(
        Constraint(name='surface_downwelling_longwave_flux_in_air'))
    rlus_cube = cubes.extract_strict(
        Constraint(name='surface_upwelling_longwave_flux_in_air'))

    # rlns = rlds - rlus
    rlns = rlds_cube - rlus_cube

    return rlns


def calc_cllmtisccp(cubes):
    """Compute cllmtisccp:

    long name: ISCCP Low Level Medium-Thickness Cloud Area Fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP Low Level Medium-Thickness Cloud Area Fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
    plev = iris.Constraint(air_pressure=lambda p: p > 68000.)
    cllmtisccp_cube = clisccp_cube
    cllmtisccp_cube = cllmtisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in cllmtisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        cllmtisccp_cube = cllmtisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        cllmtisccp_cube = cllmtisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return cllmtisccp_cube


def calc_clltkisccp(cubes):
    """Compute clltkisccp:

    long name: ISCCP low level thick cloud area fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP low level thick cloud area fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
    plev = iris.Constraint(air_pressure=lambda p: p > 68000.)
    clltkisccp_cube = clisccp_cube
    clltkisccp_cube = clltkisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in clltkisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        clltkisccp_cube = clltkisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        clltkisccp_cube = clltkisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return clltkisccp_cube


def calc_clmmtisccp(cubes):
    """Compute clmmtisccp:

    long name: ISCCP Middle Level Medium-Thickness Cloud Area Fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP Middle Level Medium-Thickness Cloud Area Fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
    plev = iris.Constraint(air_pressure=lambda p: 44000. < p <= 68000.)
    clmmtisccp_cube = clisccp_cube
    clmmtisccp_cube = clmmtisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in clmmtisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        clmmtisccp_cube = clmmtisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        clmmtisccp_cube = clmmtisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return clmmtisccp_cube


def calc_clmtkisccp(cubes):
    """Compute clmtkisccp:

    long name: ISCCP Middle Level Thick Cloud Area Fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP Middle Level Thick Cloud Area Fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
    plev = iris.Constraint(air_pressure=lambda p: 44000. < p <= 68000.)
    clmtkisccp_cube = clisccp_cube
    clmtkisccp_cube = clmtkisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in clmtkisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        clmtkisccp_cube = clmtkisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        clmtkisccp_cube = clmtkisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return clmtkisccp_cube


def calc_clhmtisccp(cubes):
    """Compute clhmtisccp:

    long name: ISCCP High Level Medium-Thickness Cloud Area Fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP High Level Medium-Thickness Cloud Area Fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
    plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
    clhmtisccp_cube = clisccp_cube
    clhmtisccp_cube = clhmtisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in clhmtisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        clhmtisccp_cube = clhmtisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        clhmtisccp_cube = clhmtisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return clhmtisccp_cube


def calc_clhtkisccp(cubes):
    """Compute clhtkisccp:

    long name: ISCCP high level thick cloud area fraction
    short name: same

    Arguments
    ----
        cubes: cubelist containing
               clisccp(isccp_cloud_area_fraction)

    Returns
    -------
        Cube: ISCCP high level thick cloud area fraction.
        Units: %

    """
    clisccp_cube = cubes.extract_strict(
        Constraint(name='isccp_cloud_area_fraction'))

    tau = iris.Constraint(
        atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
    plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
    clhtkisccp_cube = clisccp_cube
    clhtkisccp_cube = clhtkisccp_cube.extract(tau & plev)
    coord_names = [
        coord.standard_name for coord in clhtkisccp_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        clhtkisccp_cube = clhtkisccp_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        clhtkisccp_cube = clhtkisccp_cube.collapsed('air_pressure',
                                                    iris.analysis.SUM)

    return clhtkisccp_cube


def _pressure_level_widths(tro3_cube, ps_cube, top_limit=0):
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
