"""
Module to create 'Supermeans'.

Supermeans are averages over several years. For monthly averages each calendar
month is periodically averaged over several years, an example is the average
surface temperature in April between 2000 and 2010. For seasonal averages
seasons are averaged periodically.
Annual 'Supermeans' are averages over several full years.
"""

import os.path

import cf_units
import iris
import iris.coord_categorisation
from iris.coord_categorisation import _pt_date
import numpy as np


class NoBoundsError(ValueError):
    """Return error and pass."""


class InvalidPeriod(ValueError):
    """Return error and pass."""


def get_supermean(name, season, data_dir, obs_flag=None):
    """Calculated supermeans from retrieved data, which are pickled Iris cubes.

    :param name: Cube name. Should be CF-standard name. If no CF-standard name
                 exists the STASH code in msi format (for example m01s30i403)
                 is used as name.
    :param season: Supermean for a season (including annual).
                   ['ann', 'djf', 'mam', 'jja', 'son']
    :param data_dir: Directory containing cubes of model output data for
                     supermeans.
    :returns: Supermeaned cube.
    :rtype Cube:

    The monthly and seasonal supermeans are periodic averages, for example
    the seasonal supermean consists of the averaged season, where each
    season is averaged over several years.
    The annual supermean is a continuous mean over multiple years.

    Supermeans are only applied to full clima years (Starting Dec 1st).
    """
    name_constraint = iris.Constraint(name=name)

    if not obs_flag:
        cubes_path = os.path.join(data_dir, 'cubeList.nc')
    else:
        cubes_path = os.path.join(data_dir, obs_flag + '_cubeList.nc')
    cubes = iris.load(cubes_path)

    # use STASH if no standard name
    for cube in cubes:
        if cube.name() == 'unknown':
            cube.rename(str(cube.attributes['STASH']))

    cube = cubes.extract_cube(name_constraint)

    if season in ['djf', 'mam', 'jja', 'son']:
        supermeans_cube = periodic_mean(cube, period='season')
        return supermeans_cube.extract(iris.Constraint(season=season))
    elif season == 'ann':
        return periodic_mean(cube)
    else:
        raise ValueError(
            "Argument 'season' must be one of "
            "['ann', 'djf', 'mam', 'jja', 'son']. "
            "It is: " + str(season))


def contains_full_climate_years(cube):
    """Test whether cube covers full climate year(s).

    A climate year begins at YYYY-12-01 00:00:00,
    and ends at YYYY-12-01 00:00:00.

    In case of diurnal data, which is sampled at certain hours of the day, the
    climate year is shifted by up to 23 hours. The climate year boundaries of
    data sampled at 18:00 would be YYYY-12-01 18:00:00.

    :param Cube: Cube.
    :returns: True if first and last time bound
              in cube are at YYYY-12-01 00:00:00.
    :rtype: boolean
    """
    origin = cube.coord('time').units.origin
    calendar = cube.coord('time').units.calendar
    format_ = 'YYYY-%m-%d %H:%M:%S'

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()

    def _num2date(num):
        return cf_units.num2date(num, origin, calendar)

    if is_24h_sampled(cube):
        # find out number of sampling intervals (difference < 24 h)
        intervals = []
        for i in range(len(cube.coord('time').points) - 1):
            diff = cube.coord('time').points[i] - cube.coord('time').points[0]
            if diff < 24:
                intervals.append(round(diff))
        intervals = len(intervals)

        year_boundaries = [
            'YYYY-12-01 {:02d}:00:00'.format(hour) for hour in range(24)
        ]

        bounding_datetimes = []
        time_bounds = cube.coord('time').bounds
        for i in range(intervals):
            start = _num2date(time_bounds[i][0]).strftime(format_)
            end = _num2date(time_bounds[i - intervals][1]).strftime(format_)
            bounding_datetimes.append((start, end))
        return all(start == end and start in year_boundaries and
                   end in year_boundaries
                   for start, end in bounding_datetimes)
    else:
        start = _num2date(cube.coord('time').bounds[0][0]).strftime(format_)
        end = _num2date(cube.coord('time').bounds[-1][1]).strftime(format_)
        year_boundary = 'YYYY-12-01 00:00:00'
        return start == year_boundary and end == year_boundary


def is_24h_sampled(cube):
    """Check if cube data was sample once per day."""
    meaning_periods = []
    for c_m in cube.cell_methods:
        if c_m.method == 'mean' and 'time' in c_m.coord_names:
            meaning_periods.extend(c_m.intervals)
    return '24 hour' in meaning_periods


def periodic_mean(cube, period=None):
    """Return cube in which all identical periods are averaged into one.

    In case of months this would be averages over all Januaries, Februaries,
    etc. In case of season this would averages over all Winters, Springs,
    Summers and Autumns.
    If no period is specified the average of all data in `cube` is calculated.

    Averaging works with data sampled multiple times per day (diurnal data).

    The averaging takes the different lengths of periods in the Gregorian
    calendar into account.

    Requires cube with data for full Climate Years. Climate years start at the
    1st of December.

    :param cube: Cube with data for each calendar month.
    :param period: 'month', 'season'
    :returns: Cube with periodic monthly averages.
    :rtype: Cube

    Note: In the returned cube, the bounds for each
    period are the start boundary
    of the first period that is averaged over,
    and the end boundary of the last period that is averaged over.
    """
    if period not in [None, 'month', 'season']:
        raise InvalidPeriod('Invalid period: ' + str(period))

    _cube = cube.copy()

    if _cube.coord('time').has_bounds():
        add_start_hour(_cube, 'time', name='start_hour')
    else:
        iris.coord_categorisation.add_hour(_cube, 'time', name='start_hour')

    if period == 'month':
        iris.coord_categorisation.add_month(_cube, 'time', name='month')
    elif period == 'season':
        iris.coord_categorisation.add_season(_cube, 'time')
    elif period is None:
        pass
    else:
        raise InvalidPeriod('Invalid period: ' + str(period))

    time_points_per_day = len(set(_cube.coord('start_hour').points))
    if period is None:  # multi-annual mean
        if time_points_per_day > 1:
            _cube = time_average_by(_cube, 'start_hour')
        else:
            _cube.remove_coord('start_hour')
            _cube = time_average_by(_cube)
    else:
        if time_points_per_day > 1:
            _cube = time_average_by(_cube, [period, 'start_hour'])
        else:
            _cube.remove_coord('start_hour')
            _cube = time_average_by(_cube, period)

    return _cube


def add_start_hour(cube, coord, name='diurnal_sampling_hour'):
    """Add AuxCoord for diurnal data. Diurnal data is sampled every 24 hours.

    The hour value is taken from the first time bound, or the time point if no
    bounds exist.
    """
    _add_categorised_coord(cube, name, coord, start_hour_from_bounds)


def start_hour_from_bounds(coord, _, bounds):
    """Add hour from bounds."""
    return np.array([_pt_date(coord, _bounds[0]).hour for _bounds in bounds])


def _add_categorised_coord(cube,
                           name,
                           from_coord,
                           category_function,
                           units='1'):
    """
    Add categorized coordinate.

    This function creates a category from coordinate bounds. To derive the
    category from the points use:
        `iris.coord_categorisation.add_categorised_coord`

    This function has the same interface as
        `iris.coord_categorisation.add_categorised_coord`

    ######################################################################

    Add a new coordinate to a cube, by categorising an existing one.
    Make a new :class:`iris.coords.AuxCoord` from mapped values, and add it to
    the cube.


    Args:

    * cube (:class:`iris.cube.Cube`):
        the cube containing 'from_coord'.  The new coord will be added into it.
    * name (string):
        name of the created coordinate
    * from_coord (:class:`iris.coords.Coord` or string):
        coordinate in 'cube', or the name of one
    * category_function (callable):
        function(coordinate, value), returning a category value for a
        coordinate point-value

    Kwargs:

    * units:
        units of the category value, typically 'no_unit' or '1'.
    """
    # Interpret coord, if given as a name
    if isinstance(from_coord, str):
        from_coord = cube.coord(from_coord)

    if cube.coords(name):
        msg = 'A coordinate "%s" already exists in the cube.' % name
        raise ValueError(msg)

    new_coord = iris.coords.AuxCoord(
        category_function(from_coord, from_coord.points, from_coord.bounds),
        units=units,
        attributes=from_coord.attributes.copy())
    new_coord.rename(name)

    # Add into the cube
    cube.add_aux_coord(new_coord, cube.coord_dims(from_coord))


def time_average_by(cube, periods='time'):
    """Average cube over time or over periods.

    i. e. time-based categorical
    coordinates, with calendar dependent weighting.
    """
    if isinstance(periods, str):
        periods = [periods]

    # create new cube with time coord and orig duration as data
    durations_cube = iris.cube.Cube(
        # durations normalised to 1
        durations(cube.coord('time')) / np.max(durations(cube.coord('time'))),
        long_name='duration',
        units='1',
        attributes=None,
        dim_coords_and_dims=[(cube.coord('time').copy(), 0)])
    # there must be an AuxCoord for each period
    for period in periods:
        if period != 'time':
            durations_cube.add_aux_coord(cube.coord(period), 0)

    # calculate weighted sum
    orig_cell_methods = cube.cell_methods

    # multiply each time slice by its duration
    idx_obj = [None] * cube.data.ndim
    idx_obj[cube.coord_dims('time')[0]] = slice(
        None)  # [None, slice(None), None] == [np.newaxis, :, np.newaxis]
    cube.data *= durations_cube.data[idx_obj]

    if periods == ['time']:  # duration weighted averaging
        cube = cube.collapsed(periods, iris.analysis.SUM)
        durations_cube = durations_cube.collapsed(periods, iris.analysis.SUM)
    else:
        cube = cube.aggregated_by(periods, iris.analysis.SUM)
        durations_cube = durations_cube.aggregated_by(periods,
                                                      iris.analysis.SUM)

    # divide by aggregated weights
    if durations_cube.data.shape == ():
        cube.data /= durations_cube.data
    else:
        cube.data /= durations_cube.data[idx_obj]

    # correct cell methods
    cube.cell_methods = orig_cell_methods
    time_averaging_method = iris.coords.CellMethod(
        method='mean', coords=periods)
    cube.add_cell_method(time_averaging_method)

    return cube


def durations(time_coord):
    """Return durations of time periods."""
    assert time_coord.has_bounds(), 'No bounds. Do not guess.'
    durs = np.array(
        [bounds[1] - bounds[0] for bounds in time_coord.bounds])
    return durs
