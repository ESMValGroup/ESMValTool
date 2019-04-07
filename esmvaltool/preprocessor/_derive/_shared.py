"""Auxiliary derivation functions used for multiple variables."""

import logging

import dask.array as da
import iris
from iris import Constraint

logger = logging.getLogger(__name__)


def _get_fx_cube(cubes, standard_name, fx_land_name, fx_sea_name=None):
    """Extract desired fx cubes from list of cubes."""
    cube = cubes.extract_strict(Constraint(name=standard_name))
    fx_cube = None
    invert = False
    if fx_sea_name is not None:
        try:
            fx_cube = cubes.extract_strict(_var_name_constraint(fx_sea_name))
        except iris.exceptions.ConstraintMismatchError:
            logger.debug(
                "Cannot correct cube '%s' with '%s', fx file not found",
                standard_name, fx_sea_name)
        else:
            if not _shape_is_broadcastable(fx_cube.shape, cube.shape):
                fx_cube = None
                invert = True
                logger.debug(
                    "Cannot broadcast fx cube '%s' to cube '%s', trying '%s'",
                    fx_sea_name, standard_name, fx_land_name)
            else:
                logger.debug("Using fx cube '%s' to fix '%s'", fx_sea_name,
                             standard_name)
    if fx_cube is None:
        try:
            fx_cube = cubes.extract_strict(_var_name_constraint(fx_land_name))
        except iris.exceptions.ConstraintMismatchError:
            logger.debug(
                "Cannot correct cube '%s' with '%s', fx file not found",
                standard_name, fx_land_name)
        else:
            if not _shape_is_broadcastable(fx_cube.shape, cube.shape):
                fx_cube = None
                invert = False
                logger.debug("Cannot broadcast fx cube '%s' to cube '%s'",
                             fx_land_name, standard_name)
            else:
                logger.debug("Using fx cube '%s' to fix '%s'", fx_land_name,
                             standard_name)
    return (fx_cube, invert)


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


def grid_area_correction(cubes, standard_name, sea_var=False):
    """Correct (flux) variable defined relative to land/sea area."""
    cube = cubes.extract_strict(Constraint(name=standard_name))
    core_data = cube.core_data()
    (fraction_cube, invert) = _get_fx_cube(
        cubes,
        standard_name,
        fx_land_name='sftlf',
        fx_sea_name='sftof' if sea_var else None)
    if fraction_cube is not None:
        if not invert:
            core_data *= fraction_cube.core_data() / 100.0
        else:
            fraction = (da.ones_like(fraction_cube.core_data()) -
                        fraction_cube.core_data() / 100.0)
            core_data *= fraction
    return cube.copy(core_data)
