"""Constraints as classes

Iris constraints in the form of classes, which is useful for
e.g. multiprocessing, where lambda functions can't be pickled and used
in other processes.

An alternative is the use of functions, and then use functools.partial
to set up the criteria for each constraint. With classes, the critera
are given through the constructor.

"""

import iris


class RangeConstraint:
    """Class for a range constraint

    This class can be used with an iris.Constraint, in particular in a
    multiprocessing environment, where lambda functionality (often
    used in Iris examples) can't be used.

    The constructor takes two argmuments, `start` and `end`, that
    specify the boundary of the range, and are *inclusive*.

    Usage examples with iris.Constraint:

        constraint = iris.Constraint(year=RangeConstraint(1950, 2050))
        century_cube = constraint.extract(cube)

        constraint = iris.Constraint(time=RangeConstraint(cftime.datetime(2010, 1, 1),
                                                          cftime.datetime(2010, 6, 30))
        halfyear_cube = constraint.extract(cube)

    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __call__(self, cell):
        return self.start <= cell.point <= self.end


class CoordConstraint(RangeConstraint):
    """Class for a coordinate range constraint

    This class is mainly for clarity in code.

    The class is identical to the RangeConstraint, except for the
    constructor parameter names: `lower` and `upper` instead of
    `start` and `end`.
    See the documentation of `RangeConstraint` for details.

    Example:

        # Assumes west longitude has smaller values (lower) than east longitude
        constraint = iris.Constraint(longitude=CoordConstraint(5, 15))
        europe_africa = constraint.extract(cube)
    """


class EqualConstraint:

    """Class for a equality constraint

    This class can be used with an iris.Constraint, in particular in a
    multiprocessing environment, where lambda functionality (often
    used in Iris examples) can't be used.

    The constructor takes one argmument, `value`, that specify the
    value which the relevant coordinate (given as parameter name in
    `iris.Constraint`) should equal.

    Usage examples with iris.Constraint:

        constraint = iris.Constraint(year=2010)
        cube-2010 = constraint.extract(cube)

        # Use an auxiliary season coordinate that equals one of 'djf', 'mam', 'jja' or 'son',
        # and extract all winters
        constraint = iris.Constraint(season='djf')
        winters = constraint.extract(cube)

    """

    def __init__(self, value):
        self.value = value

    def __call__(self, value):
        return self.value == value


def make_date_constraint(start, end):
    """Factory fucntion to create a range constraint specifically for the `time` coordinate

    Arguments
    ---------

        start, end: lower and upper limits (inclusive) of the time
            boundaries. Use a `cftime.datetime` or variant as type.

    """

    constraint = RangeConstraint(start, end)
    return iris.Constraint(time=constraint)
