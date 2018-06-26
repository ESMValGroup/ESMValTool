"""
Time and area operations on data cubes

Allows for selecting data subsets using certain time bounds;
selecting geographical regions; constructing seasonal and area
averages; checks on data time frequencies (daily, monthly etc)
"""
from datetime import timedelta
import iris


# slice cube over a restricted time period
def time_slice(mycube, yr1, mo1, d1, yr2, mo2, d2):
    """
    Slice cube on time

    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or
    time_slice(cube,'2006','2','2','2010','1','1');

    Returns a cube
    """
    import datetime
    time_units = mycube.coord('time').units
    if time_units.calendar == '360_day':
        if d1 > 30:
            d1 = 30
        if d2 > 30:
            d2 = 30
    my_date1 = datetime.datetime(int(yr1), int(mo1), int(d1))
    my_date2 = datetime.datetime(int(yr2), int(mo2), int(d2))

    time1 = time_units.date2num(my_date1)
    time2 = time_units.date2num(my_date2)
    if iris.__version__ < '2.1':
        my_constraint = iris.Constraint(time=lambda t: (
            time1 < time_units.date2num(t.point) < time2))    
    else:
        my_constraint = iris.Constraint(time=lambda t: (time1 < t.point < time2))

    cube_slice = mycube.extract(my_constraint)
    return cube_slice


# get the time average
def time_average(mycube):
    """Get the time average over MEAN; returns a cube"""
    var_mean = mycube.collapsed('time', iris.analysis.MEAN)
    return var_mean


# get the probability a value is greater than a threshold
def proportion_greater(mycube, coord1, threshold):
    """
    Proportion greater

    Return the probability that a cetain variable coord1 (string)
    is greater than a threshold threshold (float or string),
    across a cube mycube; returns a cube
    """
    thr = float(threshold)
    result = mycube.collapsed(
        coord1, iris.analysis.PROPORTION, function=lambda values: values > thr)
    return result


# get the seasonal mean
def seasonal_mean(mycube):
    """
    Function to compute seasonal means with MEAN

    Chunks time in 3-month periods and computes means over them;
    Returns a cube
    """
    iris.coord_categorisation.add_season(mycube, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(
        mycube, 'time', name='season_year')
    annual_seasonal_mean = mycube.aggregated_by(['clim_season', 'season_year'],
                                                iris.analysis.MEAN)

    def spans_three_months(time):
        """Check for three months"""
        return (time.bound[1] - time.bound[0]) == 2160

    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)


# set of time axis checks
# funcs that perform checks on the time axis
# of data cubes and validates the type of data:
# daily, monthly, seasonal or yearly
class NoBoundsError(ValueError):
    """OBS files dont have bounds"""

    pass


def is_daily(cube):
    """Test whether the time coordinate contains only daily bound periods."""

    def is_day(bound):
        """Count days"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=1) == time_span

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_day(bound) for bound in cube.coord('time').bounds])


def is_monthly(cube):
    """A month is a period of at least 28 days, up to 31 days."""

    def is_month(bound):
        """Count months"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=31) >= time_span >= timedelta(days=28)

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_month(bound) for bound in cube.coord('time').bounds])


def is_seasonal(cube):
    """
    Check if data is seasonal

    A season is a period of 3 months, i.e.
    at least 89 days, and up to 92 days.
    """

    def is_season(bound):
        """Count seasons"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        is_seas = timedelta(days=31 + 30 + 31) >= time_span >= \
            timedelta(days=28 + 31 + 30)
        return is_seas

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_season(bound) for bound in cube.coord('time').bounds])


def is_yearly(cube):
    """A year is a period of at least 360 days, up to 366 days."""

    def is_year(bound):
        """Count years"""
        t_s = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=365) == t_s or timedelta(days=360) == t_s

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_year(bound) for bound in cube.coord('time').bounds])
