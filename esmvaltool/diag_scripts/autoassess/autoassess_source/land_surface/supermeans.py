"""
(C) Crown Copyright 2017, the Met Office

Module to create 'Supermeans'.

Supermeans are averages over several years. For monthly averages each calendar
month is periodically averaged over several years, an example is the average
surface temperature in April between 2000 and 2010. For seasonal averages
seasons are averaged periodically. Annual 'Supermeans' are averages over
several full years.
"""

import os.path
import iris
import iris.coord_categorisation
import numpy as np


class InvalidPeriod(ValueError):
    """Type error of wrong period"""

    pass


def get_supermean(name, season, data_dir):
    """Calculated supermeans from retrieved data, which are Iris cubes.

    :param name: Cube name. Should be CF-standard name. If no CF-standard name
                 exists the STASH code in msi format (for example m01s30i403)
                 is used as name.
    :param season: Supermean for a season (including annual).
                   ['ann', 'djf', 'mam', 'jja', 'son']
    :param data_dir: Directory containing pickled cubes of model output data
                     for supermeans.
    :returns: Supermeaned cube.
    :rtype Cube:

    The monthly and seasonal supermeans are periodic averages, for example
    the seasonal supermean consists of the averaged season, where each
    season is averaged over several years.
    The annual supermean is a continuous mean over multiple years.
    Supermeans are only applied to full climat years (Starting Dec 1st).
    """
    name_constraint = iris.Constraint(name=name)

    cubes_path = os.path.join(data_dir, 'cubeList.nc')
    cubes = iris.load(cubes_path)

    # use STASH if no standard name
    for cube in cubes:
        if cube.name() == 'unknown':
            cube.rename(str(cube.attributes['STASH']))

    cube = cubes.extract_strict(name_constraint)

    if season in ['djf', 'mam', 'jja', 'son']:
        supermeans_cube = periodic_mean(cube, period='season')
        return supermeans_cube.extract(iris.Constraint(season=season))
    elif season == 'ann':
        return periodic_mean(cube)
    else:
        raise ValueError(
            "Argument 'season' must be one of"
            "['ann', 'djf', 'mam', 'jja', 'son']. "
            "It is: " + str(season))


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

    Note: Returned cube, the bounds for each period are the start boundary
    of the first period that is averaged over, and the end boundary of the last
    period that is averaged over.
    """
    if period not in [None, 'month', 'season']:
        raise InvalidPeriod('Invalid period: ' + str(period))

    _cube = cube.copy()

    iris.coord_categorisation.add_hour(_cube, 'time', name='start_hour')

    if period == 'month':
        iris.coord_categorisation.add_month(_cube, 'time', name='month')
        all_months = set([
            'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
            'Oct', 'Nov', 'Dec'
        ])
        assert set(_cube.coord('month').points) == all_months, \
            'Not all months in Cube.'
    elif period == 'season':
        iris.coord_categorisation.add_season(_cube, 'time')
        all_seasons = set(['djf', 'mam', 'jja', 'son'])
        assert set(_cube.coord('season').points) == all_seasons, \
            'Not all seasons in Cube.'
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


def time_average_by(cube, periods='time'):
    """Average cube over time or over periods

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
