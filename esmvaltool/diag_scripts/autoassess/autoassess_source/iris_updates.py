'''
Use this module to store code that should be added to the main Iris code
base.

This should be copies of iris routines with appropriate changes.
'''

import numpy as np

import iris.exceptions


def _guess_bounds(coord, bound_position=0.5, bound_min=None, bound_max=None):
    """
    A copy of the iris.coords._guess_bounds() method, but applying physical
    limits to coordinates.

    Return bounds for this coordinate based on its points.

    Kwargs:

    * bound_position - The desired position of the bounds relative to the
                        position of the points.

    * bound_min - A bound minimum beyond which a bound cannot be extrapolated

    * bound_max - A bound maximum beyond which a bound cannot be extrapolated

    Returns:
        A numpy array of shape (len(coord.points), 2).

    .. note::

        This method only works for coordinates with ``coord.ndim == 1``.

    """
    # XXX Consider moving into DimCoord
    # ensure we have monotonic points
    if not coord.is_monotonic():
        raise ValueError("Need monotonic points to generate bounds for %s"
                         % coord.name())

    if coord.ndim != 1:
        raise iris.exceptions.CoordinateMultiDimError(coord)

    if coord.shape[0] < 2:
        raise ValueError('Cannot guess bounds for a coordinate of length '
                         '1.')

    if coord.bounds is not None:
        raise ValueError('Coord already has bounds. Remove the bounds '
                         'before guessing new ones.')

    if getattr(coord, 'circular', False):
        points = np.empty(coord.points.shape[0] + 2)
        points[1:-1] = coord.points
        direction = 1 if coord.points[-1] > coord.points[0] else -1
        points[0] = coord.points[-1] - (coord.units.modulus * direction)
        points[-1] = coord.points[0] + (coord.units.modulus * direction)
        diffs = np.diff(points)
    else:
        diffs = np.diff(coord.points)
        diffs = np.insert(diffs, 0, diffs[0])
        diffs = np.append(diffs, diffs[-1])

    min_bounds = coord.points - diffs[:-1] * bound_position
    max_bounds = coord.points + diffs[1:] * (1 - bound_position)

    # Apply given minimum bound
    # Using explicit test for bound_min as bound_min=0.0 fails test
    if bound_min is not None:
        min_bounds = np.maximum(min_bounds, bound_min)
        max_bounds = np.maximum(max_bounds, bound_min)

    # Apply given maximum bound
    # Using explicit test for bound_max as bound_max=0.0 fails test
    if bound_max is not None:
        min_bounds = np.minimum(min_bounds, bound_max)
        max_bounds = np.minimum(max_bounds, bound_max)

    bounds = np.array([min_bounds, max_bounds]).transpose()

    return bounds


def guess_bounds(coord, bound_position=0.5, bound_min=None, bound_max=None):
    '''
    A copy of the iris.coords.guess_bounds() method, but applying physical
    limits to coordinates.

    Add contiguous bounds to a coordinate, calculated from its points.

    Puts a cell boundary at the specified fraction between each point and
    the next, plus extrapolated lowermost and uppermost bound points, so
    that each point lies within a cell.

    With regularly spaced points, the resulting bounds will also be
    regular, and all points lie at the same position within their cell.
    With irregular points, the first and last cells are given the same
    widths as the ones next to them.

    Kwargs:

    * bound_position - The desired position of the bounds relative to the
                        position of the points.

    * bound_min - A bound minimum beyond which a bound cannot be extrapolated

    * bound_max - A bound maximum beyond which a bound cannot be extrapolated

    .. note::

        An error is raised if the coordinate already has bounds, is not
        one-dimensional, or is not monotonic.

    .. note::

        Unevenly spaced values, such from a wrapped longitude range, can
        produce unexpected results :  In such cases you should assign
        suitable values directly to the bounds property, instead.

    '''
    coord.bounds = _guess_bounds(coord, bound_position, bound_min, bound_max)
