"""Auxiliary derivation functions used for multiple variables."""

import logging

import iris
from iris import Constraint

logger = logging.getLogger(__name__)


def _get_land_fraction(cubes, standard_name, derive_from_ocean_fraction=False):
    """Extract land fraction as :mod:`dask.array`."""
    cube = cubes.extract_strict(Constraint(name=standard_name))
    if derive_from_ocean_fraction:
        fx_vars = ['sftof', 'sftlf']
    else:
        fx_vars = ['sftlf']
    land_fraction = None
    for fx_var in fx_vars:
        if land_fraction is not None:
            break
        try:
            fx_cube = cubes.extract_strict(_var_name_constraint(fx_var))
        except iris.exceptions.ConstraintMismatchError:
            logger.debug(
                "Cannot correct cube '%s' with '%s', fx file not found",
                standard_name, fx_var)
        else:
            if not _shape_is_broadcastable(fx_cube.shape, cube.shape):
                logger.debug("Cannot broadcast fx cube '%s' to cube '%s'",
                             fx_var, standard_name)
            else:
                if fx_var == 'sftof':
                    land_fraction = 1.0 - fx_cube.core_data() / 100.0
                else:
                    land_fraction = fx_cube.core_data() / 100.0
                logger.debug("Using fx cube '%s' to fix '%s'", fx_var,
                             standard_name)
    return land_fraction


def _shape_is_broadcastable(shape_1, shape_2):
    """Check if two :mod:`numpy.array' shapes are broadcastable."""
    return all((m == n) or (m == 1) or (n == 1)
               for (m, n) in zip(shape_1[::-1], shape_2[::-1]))


def _var_name_constraint(var_name):
    """:mod:`iris.Constraint` using `var_name` of a :mod:`iris.cube.Cube`."""
    return Constraint(cube_func=lambda c: c.var_name == var_name)


def cloud_area_fraction(cubes, tau_constraint, plev_constraint):
    """Calculate cloud area fraction for different parameters."""
    clisccp_cube = cubes.extract_strict(
        iris.Constraint(name='isccp_cloud_area_fraction'))
    new_cube = clisccp_cube
    new_cube = new_cube.extract(tau_constraint & plev_constraint)
    coord_names = [
        coord.standard_name for coord in new_cube.coords()
        if len(coord.points) > 1
    ]
    if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
        new_cube = new_cube.collapsed(
            'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
    if 'air_pressure' in coord_names:
        new_cube = new_cube.collapsed('air_pressure', iris.analysis.SUM)

    return new_cube


def grid_area_correction(cubes, standard_name, ocean_var=False):
    """Correct (flux) variable defined relative to land/sea area."""
    cube = cubes.extract_strict(Constraint(name=standard_name))
    core_data = cube.core_data()
    land_fraction = _get_land_fraction(
        cubes, standard_name, derive_from_ocean_fraction=ocean_var)
    if land_fraction is not None:
        if ocean_var:
            land_fraction = 1.0 - land_fraction
        cube.data = core_data * land_fraction
    return cube
