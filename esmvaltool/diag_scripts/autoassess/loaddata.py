"""
(C) Crown Copyright 2017, the Met Office.

Module to replicate data access API of the previous version of AutoAssess:

    - from loaddata import load_run_ss

For information on the PP-Header attributes, see:

Unified Model Documentation Paper F03: "Input and Output File Formats"
available here: https://code.metoffice.gov.uk/doc/um/vn10.5/umdp.html.
"""

import os.path
import re
import datetime
from datetime import timedelta as td
from datetime import datetime as dd

import cf_units
import iris
import iris.coord_categorisation as coord_cat


def is_daily(cube):
    """Test whether the time coordinate contains only daily bound periods."""
    def is_day(bound):
        """Check if day."""
        time_span = td(hours=(bound[1] - bound[0]))
        return td(days=1) == time_span

    return all([is_day(bound) for bound in cube.coord('time').bounds])


def is_monthly(cube):
    """A month is a period of at least 28 days, up to 31 days."""
    def is_month(bound):
        """Check if month."""
        time_span = td(days=(bound[1] - bound[0]))
        return td(days=31) >= time_span >= td(days=28)

    return all([is_month(bound) for bound in cube.coord('time').bounds])


def is_seasonal(cube):
    """Season is 3 months, i.e. at least 89 days, and up to 92 days."""
    def is_season(bound):
        """Check if season."""
        time_span = td(days=(bound[1] - bound[0]))
        return td(days=31 + 30 + 31) >= time_span >= td(days=28 + 31 + 30)

    return all([is_season(bound) for bound in cube.coord('time').bounds])


def is_yearly(cube):
    """A year is a period of at least 360 days, up to 366 days."""
    def is_year(bound):
        """Check if year."""
        time_span = td(days=(bound[1] - bound[0]))
        return td(days=365) == time_span or td(days=360) == time_span

    return all([is_year(bound) for bound in cube.coord('time').bounds])


def is_time_mean(cube):
    """Check whether a cube was averaged over time."""
    for cell_method in cube.cell_methods:
        if cell_method.method == 'mean' and 'time' in cell_method.coord_names:
            return True


def is_zonal_mean(cube):
    """Check whether a cube was averaged over longitude."""
    for cell_met in cube.cell_methods:
        if cell_met.method == 'mean' and 'longitude' in cell_met.coord_names:
            return True


def is_minimum_in_period(cube):
    """Check if cube contains minimum values during time period."""
    for cell_met in cube.cell_methods:
        if cell_met.method == 'minimum' and 'time' in cell_met.coord_names:
            return True


def is_maximum_in_period(cube):
    """Check if cube contains maximum values during time period."""
    for cell_met in cube.cell_methods:
        if cell_met.method == 'maximum' and 'time' in cell_met.coord_names:
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


# get the seasonal mean
def seasonal_mean(mycube):
    """
    Function to compute seasonal means with MEAN.

    Chunks time in 3-month periods and computes means over them;
    Returns a cube.
    """
    if 'clim_season' not in mycube.coords():
        coord_cat.add_season(mycube, 'time', name='clim_season')
    if 'season_year' not in mycube.coords():
        coord_cat.add_season_year(mycube, 'time', name='season_year')
    annual_seasonal_mean = mycube.aggregated_by(['clim_season', 'season_year'],
                                                iris.analysis.MEAN)

    def spans_three_months(time):
        """Check for three months."""
        return (time.bound[1] - time.bound[0]) == 90  # days

    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)


# get annual mean
def annual_mean(mycube):
    """
    Function to compute annual mean with MEAN.

    Chunks time in 365-day periods and computes means over them;
    Returns a cube.
    """
    coord_cat.add_year(mycube, 'time')
    yr_mean = mycube.aggregated_by('year', iris.analysis.MEAN)

    def spans_year(time):
        """Check for 12 months."""
        return (time.bound[1] - time.bound[0]) == 365

    t_bound = iris.Constraint(time=spans_year)
    return yr_mean.extract(t_bound)


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
    select_period = {
        'daily': is_daily,
        'monthly': is_monthly,
        'seasonal': is_seasonal,
        'annual': is_yearly
    }
    if averaging_period == 'seasonal':
        selected_cubes = [
            cube for cube in cubes
            if select_period[averaging_period](seasonal_mean(cube))
        ]
    elif averaging_period == 'annual':
        selected_cubes = [
            cube for cube in cubes
            if select_period[averaging_period](annual_mean(cube))
        ]
    else:
        selected_cubes = [
            cube for cube in cubes if select_period[averaging_period](cube)
        ]
    return iris.cube.CubeList(selected_cubes)


def select_by_pressure_level(cubes, lblev):
    """
    Select data from CubeList on the specified pressure levels.

    :param CubeList cubes: Iris CubeList.
    :param list lblev: List of pressure levels in hPa.
    :returns: CubeList with Cubes only containing specified pressure levels.
    :rtype: CubeList.
    """
    pressure_level = iris.Constraint(pressure=lblev)
    return cubes.extract(
        pressure_level)  # CubeList.extract returns always CubeList


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
    select_processing = {
        64: is_zonal_mean,
        128: is_time_mean,
        4096: is_minimum_in_period,
        8192: is_maximum_in_period
    }
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
    Select cube.

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
        i_a, i_c = str(lbtim)[:][0], str(lbtim)[:][2]

        for cube in cubes:
            # select by original meaning interval (IA)
            select_meaning_interval = {1: ('1 hour', ), 6: ('6 hour', )}
            if select_meaning_interval[int(
                    i_a)] != cube.cell_methods[0].intervals:
                continue

            # select calendar (I_C)
            # see cf_units.CALENDARS for possible cube calendars
            select_calendar = {1: 'gregorian', 2: '360_day'}
            if select_calendar[int(i_c)] == cube.coord('time').units.calendar:
                selected_cubes.append(cube)
    return selected_cubes


def select_certain_months(cubes, lbmon):
    """
    Select data from CubeList that matches the specified months.

    :param CubeList cubes: Iris CubeList.
    :param list lbmon: List with month numbers, e.g. lbmon=[5,6,7] for Mai,
        June, and July.
    :returns: CubeList with Cubes containing only data for the specified mnth.
    :rtype: CubeList
    :raises: `AssertionError` if `cubes` is not an `iris.cube.CubeList`.
    """
    # add 'month number' coordinate
    add_time_coord = {
        'monthly': lambda cube: coord_cat.add_month_number(
            cube, 'time', name='month_number'),
        'seasonal': lambda cube: coord_cat.add_season(cube,
                                                      'time',
                                                      name='clim_season'),
        'annual': lambda cube: coord_cat.add_season_year(cube,
                                                         'time',
                                                         name='season_year')
    }
    assert isinstance(cubes, iris.cube.CubeList)

    for cube in cubes:
        add_time_coord['monthly'](cube)

    # filter by month number
    month_constraint = iris.Constraint(month_number=lbmon)
    return cubes.extract(
        month_constraint)  # CubeList.extract returns always CubeList


def get_time_offset(time_unit):
    """Return a datetime object equivalent to tunit."""
    # tunit e.g. 'day since 1950-01-01 00:00:00.0000000 UTC'
    cfunit = cf_units.Unit(time_unit, calendar=cf_units.CALENDAR_STANDARD)
    time_offset = cfunit.num2date(0)
    return time_offset


def datetime_to_int_days(date_obj, tunit):
    """Return time point converted from cube datetime cell."""
    if float(iris.__version__.split('.')[0]) >= 2.0:
        time_offset = get_time_offset(tunit)
        real_date = dd(date_obj.year, date_obj.month, date_obj.day, 0, 0, 0)
        days = (real_date - time_offset).days
    else:
        days = date_obj
    return days


def extract_time_range(cubes, start, end):
    """Extract time ranged data."""
    time_ranged_cubes = []
    iris.util.unify_time_units(cubes)
    time_unit = cubes[0].coord('time').units.name
    dd_start = dd(start.year, start.month, start.day, 0, 0, 0)
    t_1 = cf_units.date2num(dd_start, time_unit, cf_units.CALENDAR_STANDARD)
    dd_end = dd(end.year, end.month, end.day, 0, 0, 0)
    t_2 = cf_units.date2num(dd_end, time_unit, cf_units.CALENDAR_STANDARD)
    for cube in cubes:
        time_constraint = iris.Constraint(
            time=lambda t: (t_1 <= datetime_to_int_days(t.point,
                                                        time_unit) <= t_2))
        cube_slice = cube.extract(time_constraint)
        time_ranged_cubes.append(cube_slice)
    return time_ranged_cubes


def load_run_ss(run_object,
                averaging_period,
                variable_name,
                lbmon=None,
                lbproc=None,
                lblev=None,
                lbtim=None,
                from_dt=None,
                to_dt=None):
    """
    Use - this is still used.

    DEPRECATED: Do not use for new Assessment Areas. Instead, read the
    CubeList `cubeList.nc` in the directory with the retrieved data.

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
    cubelist_file = 'cubeList.nc'
    cubelist_path = os.path.join(run_object['data_root'], run_object['runid'],
                                 run_object['_area'], cubelist_file)

    cubes = iris.load(cubelist_path)
    cubes.sort(key=lambda c: c.standard_name)

    return _load_run_ss(
        cubes,
        run_object,
        averaging_period,
        variable_name,
        lbmon=lbmon,
        lbproc=lbproc,
        lblev=lblev,
        lbtim=lbtim,
        from_dt=from_dt,
        to_dt=to_dt)


def _load_run_ss(cubes,
                 run_object,
                 averaging_period,
                 variable_name,
                 lbmon=None,
                 lbproc=None,
                 lblev=None,
                 lbtim=None,
                 from_dt=None,
                 to_dt=None):
    """
    Select a single Cube from the given cubes.

    See `load_run_ss` for explanation of the arguments.
    """
    arguments = 'cubes: ' + str(cubes) + '\n' + \
                'run_object: ' + str(run_object) + '\n' + \
                'averaging_period: ' + averaging_period + ', ' + \
                'variable_name' + variable_name + ', ' + \
                'lbmon=' + str(lbmon) + ', ' + \
                'lbproc=' + str(lbproc) + ', ' + \
                'lblev=' + str(lblev) + ', ' + \
                'lbtim=' + str(lbtim) + ', ' + \
                'from_dt=' + str(from_dt) + ', ' + 'to_dt=' + str(to_dt)

    selected_cubes = select_by_variable_name(cubes, variable_name)

    if averaging_period in ['daily', 'monthly', 'seasonal', 'annual']:
        selected_cubes = select_by_averaging_period(selected_cubes,
                                                    averaging_period)
    if lblev:
        selected_cubes = select_by_pressure_level(selected_cubes, lblev)

    if lbtim:
        selected_cubes = select_by_initial_meaning_period(
            selected_cubes, lbtim)

    if lbmon:
        selected_cubes = select_certain_months(selected_cubes, lbmon)

    # time constraint
    if from_dt:
        start = from_dt.date()
    else:
        start = run_object['from_' +
                           averaging_period]  # Old AutoAssess behaviour
    assert isinstance(start, datetime.date)

    if to_dt:
        end = to_dt.date()
    else:
        end = run_object['to_' + averaging_period]  # Old AutoAssess behaviour
    assert isinstance(end, datetime.date)

    selected_cubes = extract_time_range(selected_cubes, start, end)

    assert len(selected_cubes) > 0, 'No cube found.' + arguments
    assert selected_cubes[0] is not None, 'No cube found.' + arguments
    assert len(selected_cubes) < 2, 'More than one cube found.' + arguments

    return selected_cubes[0]
