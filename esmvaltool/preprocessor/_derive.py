"""Automatically derive variables."""


import logging

import iris
import yaml
# from iris import Constraint

from .derived_variables.derived_variable import DerivedVariable

logger = logging.getLogger(__name__)


def get_required(short_name, field=None):
    """Get variable short_name and field pairs required to derive variable."""
    frequency = field[2] if field else 'M'
    derived_var = DerivedVariable.get_derived_variable(short_name)
    return derived_var.get_required(frequency)

    # required = {
    #     'lwcre': [
    #         ('rlut', 'T2' + frequency + 's'),
    #         ('rlutcs', 'T2' + frequency + 's'),
    #     ],
    #     'lwp': [
    #         ('clwvi', 'T2' + frequency + 's'),
    #         ('clivi', 'T2' + frequency + 's'),
    #     ],
    #     'netcre': [
    #         ('rlut', 'T2' + frequency + 's'),
    #         ('rlutcs', 'T2' + frequency + 's'),
    #         ('rsut', 'T2' + frequency + 's'),
    #         ('rsutcs', 'T2' + frequency + 's'),
    #     ],
    #     'swcre': [
    #         ('rsut', 'T2' + frequency + 's'),
    #         ('rsutcs', 'T2' + frequency + 's'),
    #     ],
    #     'toz': [
    #         ('tro3', 'T3' + frequency),
    #         ('ps', 'T2' + frequency + 's'),
    #     ],
    #     'rtnt': [('rsdt', 'T2' + frequency + 's'),
    #              ('rsut', 'T2' + frequency + 's'), ('rlut',
    #                                                 'T2' + frequency + 's')],
    #     'rsnt': [
    #         ('rsdt', 'T2' + frequency + 's'),
    #         ('rsut', 'T2' + frequency + 's'),
    #     ],
    #     'rsns': [
    #         ('rsds', 'T2' + frequency + 's'),
    #         ('rsus', 'T2' + frequency + 's'),
    #     ],
    #     'rlns': [
    #         ('rlds', 'T2' + frequency + 's'),
    #         ('rlus', 'T2' + frequency + 's'),
    #     ],
    #     'cllmtisccp': [('clisccp', 'T4' + frequency)],
    #     'clltkisccp': [('clisccp', 'T4' + frequency)],
    #     'clmmtisccp': [('clisccp', 'T4' + frequency)],
    #     'clmtkisccp': [('clisccp', 'T4' + frequency)],
    #     'clhmtisccp': [('clisccp', 'T4' + frequency)],
    #     'clhtkisccp': [('clisccp', 'T4' + frequency)]
    # }

    # if short_name in required:
    #     return required[short_name]


def derive(cubes, variable):
    """Derive variable."""
    short_name = variable['short_name']

    # Do nothing if variable is already available
    if short_name == cubes[0].var_name:
        return cubes[0]

    # Available derivation functions
    # functions = {
    #     'lwcre': calc_lwcre,
    #     'lwp': calc_lwp,
    #     'netcre': calc_netcre,
    #     'swcre': calc_swcre,
    #     'toz': calc_toz,
    #     'rtnt': calc_rtnt,
    #     'rsnt': calc_rsnt,
    #     'rsns': calc_rsns,
    #     'rlns': calc_rlns,
    #     'cllmtisccp': calc_cllmtisccp,
    #     'clltkisccp': calc_clltkisccp,
    #     'clmmtisccp': calc_clmmtisccp,
    #     'clmtkisccp': calc_clmtkisccp,
    #     'clhmtisccp': calc_clhmtisccp,
    #     'clhtkisccp': calc_clhtkisccp
    # }

    # Preprare input cubes and derive
    # cubes = iris.cube.CubeList(cubes)
    # cube = functions[short_name](cubes)

    # Preprare input cubes and derive correct variable
    cubes = iris.cube.CubeList(cubes)
    derived_var = DerivedVariable.get_derived_variable(short_name)
    cube = derived_var.calculate(cubes)

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


# def calc_lwcre(cubes):
#     """Compute longwave cloud radiative effect from all-sky and clear-sky flux.

#     Arguments
#     ----
#         cubes: cubelist containing rlut (toa_outgoing_longwave_flux) and rlutcs
#                (toa_outgoing_longwave_flux_assuming_clear_sky).

#     Returns
#     -------
#         Cube containing longwave cloud radiative effect.

#     """
#     rlut_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_longwave_flux'))
#     rlutcs_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

#     lwcre = rlutcs_cube - rlut_cube
#     lwcre.units = rlut_cube.units

#     return lwcre


# def calc_lwp(cubes):
#     """Compute liquid water path.

#     Liquid water path is calculated by subtracting clivi (ice water) from clwvi
#     (condensed water path).
#     Note: Some datasets output the variable "clwvi" which only contains lwp. In
#     these cases, the input clwvi cube is just returned.

#     Arguments
#     ---------
#         cubes: cubelist containing clwvi_cube and clivi_cube

#     Returns
#     -------
#         Cube containing liquid water path.

#     """
#     clwvi_cube = cubes.extract_strict(
#         Constraint(name='atmosphere_cloud_condensed_water_content'))
#     clivi_cube = cubes.extract_strict(
#         Constraint(name='atmosphere_cloud_ice_content'))

#     dataset = clwvi_cube.attributes.get('model_id')
#     project = clwvi_cube.attributes.get('project_id')
#     # Should we check that the model_id/project_id are the same on both cubes?

#     bad_datasets = [
#         'CESM1-CAM5-1-FV2', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS',
#         'IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'CCSM4',
#         'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MIROC-ESM',
#         'CSIRO-Mk3-6-0', 'MPI-ESM-MR', 'MPI-ESM-LR', 'MPI-ESM-P'
#     ]
#     if ((project in ["CMIP5", "CMIP5_ETHZ"] and dataset in bad_datasets)
#             or (project == 'OBS' and dataset == 'UWisc')):
#         logger.info(
#             "Assuming that variable clwvi from %s dataset %s "
#             "contains only liquid water", project, dataset)
#         lwp_cube = clwvi_cube
#     else:
#         lwp_cube = clwvi_cube - clivi_cube

#     return lwp_cube


# def calc_netcre(cubes):
#     """Compute net cloud radiative effect.

#        Calculate net cre as sum of longwave and shortwave cloud
#        radiative effects.

#     Arguments
#     ----
#         cubes: cubelist containing rlut (toa_outgoing_longwave_flux), rlutcs
#                (toa_outgoing_longwave_flux_assuming_clear_sky), rsut
#                (toa_outgoing_shortwave_flux) and rsutcs
#                (toa_outgoing_shortwave_flux_assuming_clear_sky).

#     Returns
#     -------
#         Cube containing net cloud radiative effect.

#     """
#     lwcre = calc_lwcre(cubes)
#     swcre = calc_swcre(cubes)

#     netcre = lwcre + swcre
#     netcre.units = lwcre.units

#     return netcre


# def calc_swcre(cubes):
#     """Compute shortwave cloud radiative effect from all-sky and clear-sky

#        flux.

#     Arguments
#     ----
#         cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
#                rsutcs (toa_outgoing_shortwave_flux_assuming_clear_sky).

#     Returns
#     -------
#         Cube containing shortwave cloud radiative effect.

#     """
#     rsut_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_shortwave_flux'))
#     rsutcs_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

#     swcre = rsutcs_cube - rsut_cube

#     return swcre


# def calc_rtnt(cubes):
#     """Compute rtnt: TOA Net downward Total Radiation.

#     Arguments
#     ----
#         cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
#                rsdt (toa_incoming_shortwave_flux) and
#                rlut (toa_outgoing_longwave_flux).

#     Returns
#     -------
#         Cube containing TOA Net downward Total Radiation.
#         Units: W m-2

#     """
#     rsdt_cube = cubes.extract_strict(
#         Constraint(name='toa_incoming_shortwave_flux'))
#     rsut_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_shortwave_flux'))
#     rlut_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_longwave_flux'))

#     # rtnt = (rsdt - rsut) - rlut
#     rtnt = rsdt_cube - rsut_cube - rlut_cube

#     return rtnt


# def calc_rsnt(cubes):
#     """Compute rsnt: TOA Net downward Shortwave Radiation.

#     Arguments
#     ----
#         cubes: cubelist containing rsut (toa_outgoing_shortwave_flux) and
#                rsdt (toa_incoming_shortwave_flux).

#     Returns
#     -------
#         Cube containing TOA Net downward Shortwave Radiation.
#         Units: W m-2

#     """
#     rsdt_cube = cubes.extract_strict(
#         Constraint(name='toa_incoming_shortwave_flux'))
#     rsut_cube = cubes.extract_strict(
#         Constraint(name='toa_outgoing_shortwave_flux'))

#     # rsnt = rsdt - rsut
#     rsnt = rsdt_cube - rsut_cube

#     return rsnt


# def calc_rsns(cubes):
#     """Compute rsns: Surface Net downward Shortwave Radiation.

#     Arguments
#     ----
#         cubes: cubelist containing
#                rsus (surface_upwelling_shortwave_flux_in_air) and
#                rsds (surface_downwelling_shortwave_flux_in_air).

#     Returns
#     -------
#         Cube containing Surface Net downward Shortwave Radiation.
#         Units: W m-2

#     """
#     rsds_cube = cubes.extract_strict(
#         Constraint(name='surface_downwelling_shortwave_flux_in_air'))
#     rsus_cube = cubes.extract_strict(
#         Constraint(name='surface_upwelling_shortwave_flux_in_air'))

#     # rsns = rsds - rsus
#     rsns = rsds_cube - rsus_cube

#     return rsns


# def calc_rlns(cubes):
#     """Compute rlns: Surface Net downward Longwave Radiation.

#     Arguments
#     ----
#         cubes: cubelist containing
#                rlds (surface_downwelling_longwave_flux_in_air) and
#                rlus (surface_upwelling_longwave_flux_in_air).

#     Returns
#     -------
#         Cube containing Surface Net downward Longwave Radiation.
#         Units: W m-2

#     """
#     rlds_cube = cubes.extract_strict(
#         Constraint(name='surface_downwelling_longwave_flux_in_air'))
#     rlus_cube = cubes.extract_strict(
#         Constraint(name='surface_upwelling_longwave_flux_in_air'))

#     # rlns = rlds - rlus
#     rlns = rlds_cube - rlus_cube

#     return rlns


# def calc_cllmtisccp(cubes):
#     """Compute cllmtisccp:

#     long name: ISCCP Low Level Medium-Thickness Cloud Area Fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP Low Level Medium-Thickness Cloud Area Fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
#     plev = iris.Constraint(air_pressure=lambda p: p > 68000.)
#     cllmtisccp_cube = clisccp_cube
#     cllmtisccp_cube = cllmtisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in cllmtisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         cllmtisccp_cube = cllmtisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         cllmtisccp_cube = cllmtisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return cllmtisccp_cube


# def calc_clltkisccp(cubes):
#     """Compute clltkisccp:

#     long name: ISCCP low level thick cloud area fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP low level thick cloud area fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
#     plev = iris.Constraint(air_pressure=lambda p: p > 68000.)
#     clltkisccp_cube = clisccp_cube
#     clltkisccp_cube = clltkisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in clltkisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         clltkisccp_cube = clltkisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         clltkisccp_cube = clltkisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return clltkisccp_cube


# def calc_clmmtisccp(cubes):
#     """Compute clmmtisccp:

#     long name: ISCCP Middle Level Medium-Thickness Cloud Area Fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP Middle Level Medium-Thickness Cloud Area Fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
#     plev = iris.Constraint(air_pressure=lambda p: 44000. < p <= 68000.)
#     clmmtisccp_cube = clisccp_cube
#     clmmtisccp_cube = clmmtisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in clmmtisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         clmmtisccp_cube = clmmtisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         clmmtisccp_cube = clmmtisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return clmmtisccp_cube


# def calc_clmtkisccp(cubes):
#     """Compute clmtkisccp:

#     long name: ISCCP Middle Level Thick Cloud Area Fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP Middle Level Thick Cloud Area Fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
#     plev = iris.Constraint(air_pressure=lambda p: 44000. < p <= 68000.)
#     clmtkisccp_cube = clisccp_cube
#     clmtkisccp_cube = clmtkisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in clmtkisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         clmtkisccp_cube = clmtkisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         clmtkisccp_cube = clmtkisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return clmtkisccp_cube


# def calc_clhmtisccp(cubes):
#     """Compute clhmtisccp:

#     long name: ISCCP High Level Medium-Thickness Cloud Area Fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP High Level Medium-Thickness Cloud Area Fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
#     plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
#     clhmtisccp_cube = clisccp_cube
#     clhmtisccp_cube = clhmtisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in clhmtisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         clhmtisccp_cube = clhmtisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         clhmtisccp_cube = clhmtisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return clhmtisccp_cube


# def calc_clhtkisccp(cubes):
#     """Compute clhtkisccp:

#     long name: ISCCP high level thick cloud area fraction
#     short name: same

#     Arguments
#     ----
#         cubes: cubelist containing
#                clisccp(isccp_cloud_area_fraction)

#     Returns
#     -------
#         Cube: ISCCP high level thick cloud area fraction.
#         Units: %

#     """
#     clisccp_cube = cubes.extract_strict(
#         Constraint(name='isccp_cloud_area_fraction'))

#     tau = iris.Constraint(
#         atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
#     plev = iris.Constraint(air_pressure=lambda p: p <= 44000.)
#     clhtkisccp_cube = clisccp_cube
#     clhtkisccp_cube = clhtkisccp_cube.extract(tau & plev)
#     coord_names = [
#         coord.standard_name for coord in clhtkisccp_cube.coords()
#         if len(coord.points) > 1
#     ]
#     if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
#         clhtkisccp_cube = clhtkisccp_cube.collapsed(
#             'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
#     if 'air_pressure' in coord_names:
#         clhtkisccp_cube = clhtkisccp_cube.collapsed('air_pressure',
#                                                     iris.analysis.SUM)

#     return clhtkisccp_cube
