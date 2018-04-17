#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Module to replicate data access API of the previous version of AutoAssess:

    - from loaddata import load_run_ss

For information on the PP-Header attributes, see:

Unified Model Documentation Paper F03: "Input and Output File Formats"
available here: https://code.metoffice.gov.uk/doc/um/vn10.5/umdp.html

"""

import cf_units
import os.path
import pickle
import re
import datetime
from datetime import timedelta as td

import iris
import iris.coord_categorisation as coord_cat
from iris.time import PartialDateTime as pdt


def is_daily(cube):
    """Test whether the time coordinate contains only daily bound periods."""
    def is_day(bound):
        time_span = td(hours=(bound[1] - bound[0]))
        return td(days=1) == time_span

    return all([is_day(bound) for bound in cube.coord('time').bounds])


def is_monthly(cube):
    """A month is a period of at least 28 days, up to 31 days."""
    def is_month(bound):
        # VPREDOI
        #time_span = td(hours=(bound[1] - bound[0]))
        time_span = td(days=(bound[1] - bound[0]))
        return td(days=31) >= time_span >= td(days=28)

    return all([is_month(bound) for bound in cube.coord('time').bounds])


def is_seasonal(cube):
    """
    A season is a period of 3 months, i.e. at least 89 days, and up to 92 days.
    """
    def is_season(bound):
        time_span = td(hours=(bound[1] - bound[0]))
        return td(days=31+30+31) >= time_span >= td(days=28+31+30)

    return all([is_season(bound) for bound in cube.coord('time').bounds])


def is_yearly(cube):
    """A year is a period of at least 360 days, up to 366 days."""
    def is_year(bound):
        time_span = td(hours=(bound[1] - bound[0]))
        return td(days=365) == time_span or td(days=360) == time_span

    return all([is_year(bound) for bound in cube.coord('time').bounds])


def is_time_mean(cube):
    """Check whether a cube was averaged over time."""
    for cell_method in cube.cell_methods:
        if cell_method.method == 'mean' and 'time' in cell_method.coord_names:
            return True


def is_zonal_mean(cube):
    """Check whether a cube was averaged over longitude."""
    for cell_method in cube.cell_methods:
        if cell_method.method == 'mean' and 'longitude' in cell_method.coord_names:
            return True


def is_minimum_in_period(cube):
    """Check if cube contains minimum values during time period."""
    for cell_method in cube.cell_methods:
        if cell_method.method == 'minimum' and 'time' in cell_method.coord_names:
            return True


def is_maximum_in_period(cube):
    """Check if cube contains maximum values during time period."""
    for cell_method in cube.cell_methods:
        if cell_method.method == 'maximum' and 'time' in cell_method.coord_names:
            return True


def select_by_variable_name(cubes, variable_name):
    """
    Select subset from CubeList matching a CF-name or STASH code.

    :param CubeList cubes: Iris CubeList.
    :param sting variable_name: CF-name or model-section-item STASH code, e.g.
        'm01s30i404'
    :returns: CubeList with Cubes that have a matching STASH attribute.
    :rtype: CubeList
    """
    regex = '^m01s[0-9]{2}i[0-9]{3}$'
    if re.match(regex, variable_name):
        constraint = iris.AttributeConstraint(STASH=variable_name)
    else:
        constraint = iris.Constraint(variable_name)
    return cubes.extract(constraint)


def select_by_averaging_period(cubes, averaging_period):
    """
    Select subset from CubeList depending on averaging period.

    :param CubeList cubes: Iris CubeList.
    :param string averaging period: Must be one of 'daily', 'monthly',
        'seasonal', and 'annual'.
    :returns: CubeList with Cubes that are averaged over a certain period.
    :rtype: CubeList
    :raises: `AssertionError` if `cubes` is not a `list`.
    """
    assert isinstance(cubes, list)
    select_period = {'daily': is_daily,
                     'monthly': is_monthly,
                     'seasonal': is_seasonal,
                     'annual': is_yearly}
    selected_cubes = [cube for cube in cubes if select_period[averaging_period](cube)]
    return iris.cube.CubeList(selected_cubes)


def select_by_pressure_level(cubes, lblev):
    """
    Select data from CubeList on the specified pressure levels.

    :param CubeList cubes: Iris CubeList.
    :param list lblev: List of pressure levels in hPa.
    :returns: CubeList with Cubes only containing specified pressure levels.
    :rtype: CubeList
    """
    pressure_level = iris.Constraint(pressure=lblev)
    return cubes.extract(pressure_level)  # CubeList.extract returns always CubeList


def select_by_processing(cubes, lbproc):
    """
    Select subset from CubeList by the processing that has been applied.

    :param CubeList cubes: Iris CubeList.
    :param list lbproc: List with PP-header attributes describing processing
        steps (lbproc). Currently, only 128 ('meaned') is implemented.
    :returns: CubeList with Cubes that were processed according to lbproc
        attribute.
    :rtype: CubeList
    """
    assert isinstance(cubes, list)
    assert lbproc != 0

    # bits are used to indicate processing
    select_processing = {64:  is_zonal_mean,
                         128: is_time_mean,
                         4096: is_minimum_in_period,
                         8192: is_maximum_in_period}
    selected_cubes = []
    for cube in cubes:
        missing_method = False
        _lbproc = lbproc
        for key, processed in select_processing.items():
            # check processing only for set bits, using bitwise AND
            # 192 & 32 = 0
            # 192 & 64 = 64
            # 192 & 128 = 128
            if lbproc & key == key:
                if not processed(cube):
                    missing_method = True
                # if processing set bit to zero using bitwise XOR
                _lbproc = _lbproc ^ key

        # _lbproc == 0 -> processing for each bit has been tested
        if not missing_method and _lbproc == 0:
            selected_cubes.append(cube)

        if _lbproc != 0:
            raise NotImplementedError('Lbproc ' + str(lbproc) + ' is not '
                                      'implemented.')

    return iris.cube.CubeList(selected_cubes)


def select_by_initial_meaning_period(cubes, lbtim):
    """
    Select subset from CubeList by matching the some of the information
    encoded in the 'Time indicator' `lbtim`. Namely, the initial meaning
    period and the used calendar.

    :param CubeList cubes: Iris CubeList.
    :param list lbtim: List with PP-Header attributes `lbtim`, a three digit
        number. Currently implemented: [121, 122, 621, 622]
        - First digit: Initial averaging time in hours. Must be either 1 or 6.
        - Second digit: Ignored.
        - Third digit: Calendar: 1 - Proleptic Gregorian Calendar
                                 2 - 360d Calendar
    :returns: CubeList of Cubes with the matching initial meaning period, and
        calendar.
    :rtype: Iris CubeList
    :raises: `NotImplementedError` for not implemented values in `lbtim`.
    """
    implemented_values = [121, 122, 621, 622]

    assert isinstance(cubes, list)
    assert isinstance(lbtim, list)

    lbtims = lbtim
    if any(lbtim not in implemented_values for lbtim in lbtims):
        msg = 'Implemented values:' + str(implemented_values) +\
              'Given:' + str(lbtims)
        raise NotImplementedError(msg)

    selected_cubes = iris.cube.CubeList()
    for lbtim in lbtims:
        IA, IB, IC = str(lbtim)[:]  # pylint: disable=unused-variable
        # IA - time interval in hours between the individual fields from which
        #      the mean was calculated
        # IB - = 2 if the field is a time mean between T1 and T2, or represents
        #          a sequence of times between T1 and T2.
        # IC - = 1 if the Proleptic Gregorian calendar is used for T1 and T2.
        #      = 2 if the '360-day' calendar (i.e. 12 30-day months) is used
        #          for T1 and T2.

        for cube in cubes:
            # select by original meaning interval (IA)
            select_meaning_interval = {1: ('1 hour',),
                                       6: ('6 hour',)}
            if not select_meaning_interval[int(IA)] == cube.cell_methods[0].intervals:
                continue

            # select by IB
            # Iris cubes have no T1 and T2 attributes, or equivalent
            # Unclear how to select Iris cubes on IB
            pass  # pylint: disable=unnecessary-pass

            # select calendar (IC)
            # see cf_units.CALENDARS for possible cube calendars
            select_calendar = {1: 'gregorian',  # TODO does iris distinguish between
                               2: '360_day'}    # proleptic_greorian and gregorian?
            if select_calendar[int(IC)] == cube.coord('time').units.calendar:
                selected_cubes.append(cube)
    return selected_cubes


def select_certain_months(cubes, lbmon):
    """
    Select data from CubeList that matches the specified months.

    :param CubeList cubes: Iris CubeList.
    :param list lbmon: List with month numbers, e.g. lbmon=[5,6,7] for Mai,
        June, and July.
    :returns: CubeList with Cubes containing only data for the specified months.
    :rtype: CubeList
    :raises: `AssertionError` if `cubes` is not an `iris.cube.CubeList`.
    """
    # add 'month number' coordinate
    add_time_coord = {
        'monthly':  lambda cube: coord_cat.add_month_number(cube, 'time', name='month_number'),
        'seasonal': lambda cube: coord_cat.add_season(cube, 'time', name='clim_season'),
        'annual':   lambda cube: coord_cat.add_season_year(cube, 'time', name='season_year')
    }
    assert isinstance(cubes, iris.cube.CubeList)

    for cube in cubes:
        add_time_coord['monthly'](cube)

    # filter by month number
    month_constraint = iris.Constraint(month_number=lbmon)
    return cubes.extract(month_constraint)  # CubeList.extract returns always CubeList


def extract_time_range(cubes, start, end):
    """
    For each cube in `cubes` keep only the data between `start` and `end`.

    It uses the time point of the data, or if available the beginning of a time
    period. This time point has to be at or after `start`, and
    before or at `end`.

    :param CubeList cubes: Iris CubeList.
    :param datetime.date start:
    :param datetime.date end:
    :returns: CubeList with Cubes that contain only data between `start` and
        `end`.
    :rtype: Iris CubeList
    """
    def make_time_range(unit, calendar):
        """Produces a function to be used in an Iris constraint."""
        start_time = pdt(year=start.year, month=start.month, day=start.day,
                         hour=0, minute=0, second=0)
        end_time = pdt(year=end.year, month=end.month, day=end.day,
                       hour=0, minute=0, second=0)
        #TODO it would be conceptually better to use datetime.datetime, but that
        # does not work due to a bug in IRIS. PartialDateTime is meant to be
        # used for periodical selections, such as 'each March'
        #start_time = datetime.datetime(start.year, start.month, start.day, 0, 0, 0)
        #end_time = datetime.datetime(end.year, end.month, end.day, 0, 0, 0)
        def time_range(cell):
            time_unit = cf_units.Unit(unit, calendar)
            if cell.bound:
                return start_time <= time_unit.num2date(cell.bound[0]) and \
                                     time_unit.num2date(cell.bound[0]) <= end_time
            else:
                cell_point_datetime = time_unit.num2date(cell.point)
                # set everything smaller than days to zero
                cell_point_datetime = datetime.datetime(cell_point_datetime.year,
                                                        cell_point_datetime.month,
                                                        cell_point_datetime.day)
                return start_time <= cell_point_datetime <= end_time
        return time_range


    time_ranged_cubes = iris.cube.CubeList()
    for cube in cubes:
        unit = cube.coord('time').units
        calendar = unit.calendar
        in_range = make_time_range(unit, calendar)
        time_range = iris.Constraint(coord_values={'time':in_range})

        time_ranged_cubes.append(cube.extract(time_range))
    return time_ranged_cubes



def load_run_ss(run_object, averaging_period, variable_name,
                lbmon=None, lbproc=None, lblev=None, lbtim=None,
                from_dt=None, to_dt=None):
    """
    DEPRECATED: Do not use for new Assessment Areas. Instead, read the
    pickled CubeList `cubes.pickle` in the directory with the retrieved data.

    Select a single Cube from the data that was retrieved for a single
    Assessment Area.

    This function replaces the function with the same name in the previous
    version of AutoAssess. It was the access interface to a common
    `STASH Split` directory, which all Assessment Areas used.
    This function required PP-header attributes to select data.

    For information on the PP-Header attributes, see:

    Unified Model Documentation Paper F03: "Input and Output File Formats"
    available here: https://code.metoffice.gov.uk/doc/um/vn10.5/umdp.html

    :param dict run_object: Dictionary specifying the assessment run. For its
        contents see `function: create_run_object` in the module
        `autoassess.run_area`.
    :param str averaging period: Must be one of 'instantaneous', 'daily',
        'monthly', 'seasonal', and 'annual'.
    :param str variable_name: CF-name or model-section-item STASH code, e.g.
        'm01s30i404'.
    :param list lbmon: List with month numbers, e.g. lbmon=[5,6,7] for Mai,
        June, and July.
    :param int lbproc: List with PP-header attributes describing processing
        steps (lbproc).
    :param list lblev: List of pressure levels in hPa.
    :param list lbtim: List with PP-Header attributes `lbtim`, a three digit
        number.
    :param datetime from_dt: Return data from this date onwards.
    :param datetime to_dt: Return data up to this date.
    :returns: A single Iris cube.
    :rtype: Iris cube
    :raises: `AssertionError` if not exactly one cube is selected.
    """
    pickle_file = 'cubes.pickle'
    pickle_path = os.path.join(run_object['data_root'], run_object['runid'],
                               run_object['_area'], pickle_file)

    with open(pickle_path, 'r') as fh:
        cubes = pickle.load(fh)
    cubes.sort(key=lambda c: c.standard_name)

    return _load_run_ss(cubes, run_object, averaging_period, variable_name,
                        lbmon=lbmon, lbproc=lbproc, lblev=lblev, lbtim=lbtim,
                        from_dt=from_dt, to_dt=to_dt)


def _load_run_ss(cubes, run_object, averaging_period, variable_name,
                 lbmon=None, lbproc=None, lblev=None, lbtim=None,
                 from_dt=None, to_dt=None):
    """
    Select a single Cube from the given cubes.

    See `load_run_ss` for explanation of the arguments.
    """
    arguments = 'cubes: ' + str(cubes) + '\n' + \
                'run_object: ' + str(run_object) + '\n' + \
                'averaging_period: ' + averaging_period + ', ' + \
                'variable_name' + variable_name + ', ' + \
                'lbmon=' + str(lbmon) + ', ' + 'lbproc=' + str(lbproc) + ', ' + \
                'lblev=' + str(lblev) + ', ' + 'lbtim=' + str(lbtim) + ', ' + \
                'from_dt=' + str(from_dt) + ', ' + 'to_dt=' + str(to_dt)

    selected_cubes = select_by_variable_name(cubes, variable_name)

    if averaging_period in ['daily', 'monthly', 'seasonal', 'annual']:
        selected_cubes = select_by_averaging_period(selected_cubes, averaging_period)
    if lblev:
        selected_cubes = select_by_pressure_level(selected_cubes, lblev)

    if lbtim:
        selected_cubes = select_by_initial_meaning_period(selected_cubes, lbtim)

    if lbmon:
        selected_cubes = select_certain_months(selected_cubes, lbmon)

    # time constraint
    if from_dt:
        start = from_dt.date()
    else:
        #start = datetime.date.min  # 1/01/01, i.e. no constraint
        start = run_object['from_' + averaging_period]  # Old AutoAssess behaviour
    assert isinstance(start, datetime.date)

    if to_dt:
        end = to_dt.date()
    else:
        #end = datetime.date.max  # 9999/12/31
        end = run_object['to_' + averaging_period]  # Old AutoAssess behaviour
    assert isinstance(end, datetime.date)

    selected_cubes = extract_time_range(selected_cubes, start, end)

    assert len(selected_cubes) > 0, 'No cube found.' + arguments
    assert selected_cubes[0] != None, 'No cube found.' + arguments
    assert len(selected_cubes) < 2, 'More than one cube found.' + arguments

    return selected_cubes[0]

