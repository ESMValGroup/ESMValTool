"""DUMMY DOCSTRING"""

from datetime import datetime
from calendar import monthrange
import numpy as np
import cftime
from cf_units import date2num
import iris
from iris.coords import DimCoord


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


def make_year_constraint_all_calendars(start, end):
    """Utility function to create a dict of time constraints on year-basis

    This create a dict of the same time constraint, but for different calendars.
    Since comparisons between different calendar types are not (always) possible,
    a calendar for a time coordinate should be compared to the specific constraint
    with the same calendar.

    The calendar type (as a string) can be obtained through the coordinate's `units`
    attribute: `cube.coord('time').units.calendar`; the resulting string is the key
    for the dict, which then as a value yields the correct constraint

    Arguments
    ---------

        start, end: integer
            Start and end year. Month and day are 1, 1 for the starting year, and
            31, 12 or 30, 12 (for a 360-day calendar) for the end year

    """

    dates = {
        'default': (cftime.datetime(start, 1, 1), cftime.datetime(end, 12, 31)),
        '360_day': (cftime.Datetime360Day(start, 1, 1), cftime.Datetime360Day(end, 12, 30)),
        '365_day': (cftime.DatetimeNoLeap(start, 1, 1), cftime.DatetimeNoLeap(end, 12, 31)),
        'proleptic_gregorian': (cftime.DatetimeProlepticGregorian(start, 1, 1),
                                cftime.DatetimeProlepticGregorian(end, 12, 31)),
        'gregorian': (cftime.DatetimeGregorian(start, 1, 1), cftime.DatetimeGregorian(end, 12, 31)),
        'julian': (cftime.DatetimeJulian(start, 1, 1), cftime.DatetimeJulian(end, 12, 31)),
    }
    constraints = {key: make_date_constraint(*value) for key, value in dates.items()}
    return constraints


def months_coord_to_days_coord(coord):
    """Convert a dimension coordinate from 'months since' to 'days since'

    This function uses the `calendar.monthrange` function to calculate
    the days per month, and sets the lower and upper bound for each month.
    Once the bounds have been set, `cf_units.date2num` is used to convert
    the bounds to numeric values, in days since the original offset.
    These bounds are averaged, to produce midpoints, which are the actual
    points for the new dimension coordinate. The new dimension coordinate
    also includes bounds, which the original may not have.

    """

    units = coord.units
    origin = units.origin
    # Assume an origin format of "<unit> since <start-date>", so we can split on 'since'
    step, startdate = map(str.strip, origin.split('since'))
    if step != 'months':
        raise ValueError('units step is not months')

    # Parse the starting date; assume it has a YYYY-MM-DD HH:MM:SS format,
    # or YYYY-MM-DD without the timestamp
    # Note: leading zeros for months, days, hours, minutes or seconds
    # may be safely ignored: 2010-1-1 or 2010-01-01 will both parse fine
    try:
        t0 = datetime.strptime(startdate, "%Y-%m-%d %H:%M:%S")  # pylint: disable=invalid-name
    except ValueError:
        t0 = datetime.strptime(startdate, "%Y-%m-%d")  # pylint: disable=invalid-name

    points = coord.points.astype(np.int)
    bounds = []
    # Remember that 'point's are in whole months
    for point in points:
        year = t0.year + point // 12
        month = t0.month + point % 12
        current = datetime(year, month, 1)
        # Get number of days for this year (ignore starting weekday number)
        _, ndays = monthrange(year, month)
        # And set the boundary dates for this month
        bounds.append([current, datetime(year, month, ndays)])

    # date2num accepts a two-dimensional numpy array
    bounds = np.array(bounds)
    boundpoints = date2num(bounds, unit=f'days since {startdate}', calendar='gregorian')

    midpoints = boundpoints.mean(axis=1)
    day_coord = DimCoord(midpoints, bounds=boundpoints, standard_name=coord.standard_name,
                         long_name=coord.long_name, units=f'days since {startdate}',
                         var_name=coord.var_name)

    return day_coord
