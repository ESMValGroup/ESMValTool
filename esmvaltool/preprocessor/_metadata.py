"""
Metadata operations on data cubes.

Allows for unit conversions.
"""
import logging

logger = logging.getLogger(__name__)


def convert_units(cube, units):
    """
    Convert the units of a cube to new ones.

    This converts units of a cube.

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube

        units: str
            new units in udunits form

    Returns
    -------
    iris.cube.Cube
        converted cube.
    """
    cube.convert_units(units)
    return cube
