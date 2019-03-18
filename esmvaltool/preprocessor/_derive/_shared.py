"""Auxiliary derivation functions used for multiple variables."""

import logging

import iris
from iris import Constraint

logger = logging.getLogger(__name__)


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


def grid_area_correction(cubes,
                         standard_name,
                         integrate=True,
                         fraction_var='land_area_fraction',
                         area_var='cell_area'):
    """Correct (flux) variable defined relative to land area and integrate."""
    cube = cubes.extract_strict(Constraint(name=standard_name))
    try:
        fraction_cube = cubes.extract_strict(Constraint(name=fraction_var))
        cube.data *= fraction_cube.data / 100.0
    except iris.exceptions.ConstraintMismatchError:
        logger.debug("Cannot correct cube '%s' with '%s', fx file not found",
                     standard_name, fraction_var)
        pass
    if integrate:
        try:
            areacella_cube = cubes.extract_strict(Constraint(name=area_var))
            cube.data *= areacella_cube.data
            cube.units *= areacella_cube.units
        except iris.exceptions.ConstraintMismatchError:
            logger.debug(
                "Area files for cube '%s' not found, trying to calculate area "
                "using bounds", standard_name)
            try:
                for coord_name in ('latitude', 'longitude'):
                    if not cube.coord(coord_name).has_bounds():
                        cube.coord(coord_name).guess_bounds()
                area = iris.analysis.cartography.area_weights(cube)
                cube.data *= area
            except iris.exceptions.CoordinateMultiDimError as exc:
                logger.error(
                    "Cannot integrate irregular cube '%s', necessary area "
                    "files not found", standard_name)
                logger.error(cube)
                raise exc
    return cube
