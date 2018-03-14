"""multimodel statistics

Functions for multi-model operations
supports a multitude of multimodel statistics
computations; the only requisite is the ingested
cubes have (TIME-LAT-LON) or (TIME-PLEV-LAT-LON)
dimensions; and obviously consistent units.

It operates on different (time) spans:
- full: computes stats on full model time;
- overlap: computes common time overlap between models;

"""

import logging
from functools import reduce

import os
from datetime import datetime
from datetime import timedelta
import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _time_pts_by_calendar(cube):
    """Extract time points depending on calendar"""
    # a_p is a prefactor that brings time points
    # to 365-day calendar
    time_units = cube.coord('time').units
    if time_units.calendar is not None:
        if time_units.calendar == '360_day':
            a_p = 365. / 360.
        else:
            a_p = 1.
    else:
        logger.debug("No calendar type: assuming 365-day")
        a_p = 1.
    # return tpoints
    t_points = [c * a_p for c in cube.coord('time').points]
    return t_points


def _plev_fix(dataset, pl_idx):
    """Extract valid plev data

    this function takes care of situations
    in which certain plevs are completely
    masked due to unavailable interpolation
    boundaries.
    """
    if np.ma.is_masked(dataset):
        # keep only the valid plevs
        if not np.all(dataset.mask[pl_idx]):
            statj = np.ma.array(dataset[pl_idx],
                                mask=dataset.mask[pl_idx])
        else:
            logger.debug('All vals in plev are masked, ignoring.')
            statj = None
    else:
        mask = np.zeros_like(dataset[pl_idx], bool)
        statj = np.ma.array(dataset[pl_idx], mask=mask)

    return statj


def _compute_statistic(datas, name):
    """Compute multimodel statistic"""
    datas = np.ma.array(datas)
    statistic = datas[0]

    if name == 'median':
        statistic_function = np.ma.median
    elif name == 'mean':
        statistic_function = np.ma.mean
    else:
        raise NotImplementedError

    # no plevs
    if len(datas[0].shape) < 3:
        statistic = statistic_function(datas, axis=0)
        return statistic

    # plevs
    for j in range(statistic.shape[0]):
        plev_check = [_plev_fix(cdata, j)
                      for cdata in datas
                      if _plev_fix(cdata, j) is not None]
        len_stat_j = sum(1 for _ in plev_check)
        stat_all = np.ma.zeros((len_stat_j,
                                statistic.shape[1],
                                statistic.shape[2]))
        for i, e_l in enumerate(plev_check):
            stat_all[i] = e_l

        # check for nr models
        if len_stat_j >= 2:
            statistic[j] = statistic_function(stat_all, axis=0)
        else:
            mask = np.ones(statistic[j].shape, bool)
            statistic[j] = np.ma.array(statistic[j],
                                       mask=mask)

    return statistic


def _put_in_cube(template_cube,
                 cube_data,
                 stat_name,
                 file_name,
                 t_axis):
    """Quick cube building and saving"""
    # grab coordinates from any cube
    times = template_cube.coord('time')
    # or get the FULL time axis
    if t_axis is not None:
        times = iris.coords.DimCoord(
            t_axis, standard_name='time',
            units=template_cube.coord('time').units)
    lats = template_cube.coord('latitude')
    lons = template_cube.coord('longitude')

    # no plevs
    if len(template_cube.shape) == 3:
        cspec = [(times, 0), (lats, 1), (lons, 2)]
    # plevs
    elif len(template_cube.shape) == 4:
        plev = template_cube.coord('air_pressure')
        cspec = [(times, 0), (plev, 1), (lats, 2), (lons, 3)]

    # correct dspec if necessary
    fixed_dspec = np.ma.fix_invalid(cube_data,
                                    copy=False,
                                    fill_value=1e+20)
    # put in cube
    stats_cube = iris.cube.Cube(fixed_dspec,
                                dim_coords_and_dims=cspec,
                                long_name=stat_name)
    coord_names = [coord.name() for coord in template_cube.coords()]
    if 'air_pressure' in coord_names:
        if len(template_cube.shape) == 3:
            stats_cube.add_aux_coord(template_cube.
                                     coord('air_pressure'))
    stats_cube.attributes['_filename'] = file_name
    # complete metadata
    stats_cube.var_name = template_cube.var_name
    stats_cube.long_name = template_cube.long_name
    stats_cube.standard_name = template_cube.standard_name
    stats_cube.units = template_cube.units
    return stats_cube


def _to_datetime(delta_t, unit_type):
    """Convert to a datetime point"""
    if unit_type == 'day since 1950-01-01 00:00:00.0000000':
        new_date = datetime(1950, 1, 1, 0) + timedelta(np.int(delta_t))
    elif unit_type == 'day since 1850-01-01 00:00:00.0000000':
        new_date = datetime(1850, 1, 1, 0) + timedelta(np.int(delta_t))
    # add more supported units here
    return new_date


def _get_overlap(cubes):
    """
    Get discrete time overlaps.

    This method gets the bounds of coord time
    from the cube and assembles a continuous time
    axis with smallest unit 1; then it finds the
    overlaps by doing a 1-dim intersect;
    takes the floor of first date and
    ceil of last date.
    """
    all_times = []
    for cube in cubes:
        # monthly data ONLY
        # converts all to 365-days calendars
        time_units = cube.coord('time').units
        bnd1 = float(cube.coord('time').points[0])
        bnd2 = float(cube.coord('time').points[-1])
        if time_units.calendar is not None:
            if time_units.calendar == '360_day':
                bnd1 = int(bnd1 / 360.) * 365.
                bnd2 = (int(bnd2 / 360.) + 1) * 365.
            else:
                bnd1 = int(bnd1 / 365.) * 365.
                bnd2 = (int(bnd2 / 365.) + 1) * 365.
        else:
            logger.debug("No calendar type: assuming 365-day")
            bnd1 = int(bnd1 / 365.) * 365.
            bnd2 = (int(bnd2 / 365.) + 1) * 365.
        all_times.append([bnd1, bnd2])
    bounds = [range(int(b[0]), int(b[-1]) + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, bounds)
    if len(time_pts) > 1:
        time_bounds_list = [time_pts[0], time_pts[-1]]
        return time_bounds_list


def _slice_cube(cube, min_t, max_t):
    """slice cube on time"""
    fmt = '%Y-%m-%d-%H'
    ctr = iris.Constraint(time=lambda x:
                          min_t <=
                          datetime.strptime(x.point.strftime(fmt),
                                            fmt) <= max_t)
    cube_slice = cube.extract(ctr)
    return cube_slice


def _slice_cube2(cube, t_1, t_2):
    """
    Efficient slicer

    Slice cube on time more memory-efficiently than
    _slice_cube(); using this function adds virtually
    no extra memory to the process
    """
    time_pts = [t for t in cube.coord('time').points]
    converted_t = _time_pts_by_calendar(cube)
    idxs = sorted([time_pts.index(ii)
                   for ii, jj in zip(time_pts, converted_t)
                   if t_1 <= jj <= t_2])
    cube_t_slice = cube.data[idxs[0]:idxs[-1] + 1]
    return cube_t_slice


def _monthly_t(cubes):
    """Rearrange time points for monthly data"""
    # get original cubes tpoints
    tpts = []
    # make all 365-day calendards
    for cube in cubes:
        tpts.append([int(t) for t in _time_pts_by_calendar(cube)])

    # convert to months for MONTHLY data
    t_x = list(set().union(*tpts))
    t_x = list({int((a / 365.) * 12.) for a in t_x})
    t_x.sort()
    # remake the time axis for roughly the 15th of the month
    t_x0 = list({t * 365. / 12. + 15. for t in t_x})
    t_x0.sort()

    return t_x, t_x0


def _full_time(cubes):
    """Construct a contiguous collection over time"""
    datas = []
    # get rearranged time points
    t_x, t_x0 = _monthly_t(cubes)
    # loop through cubes and apply masks
    for cube in cubes:
        # recast time points
        time_redone = [int(t) for t in _time_pts_by_calendar(cube)]
        # construct new shape
        fine_shape = tuple([len(t_x)] + list(cube.data.shape[1:]))
        # find indices of present time points
        oidx = [
            t_x.index(int((s / 365.) * 12.))
            for s in time_redone]
        # reshape data to include all possible times
        ndat = np.ma.resize(cube.data, fine_shape)
        # build the time mask
        c_ones = np.ones(fine_shape, bool)
        c_ones[oidx] = False
        ndat.mask |= c_ones

        # stitch new datas
        datas.append(ndat)

    # return datas
    return datas, t_x0


def _apply_overlap(cube, tx1, tx2):
    """
    Slice cubes

    This operation is memory-intensive;
    do it ONLY IF you need to do it;
    otherwise use span == constant.
    """
    logger.debug("Bounds: %s and %s", str(tx1), str(tx2))
    # slice cube on time
    with iris.FUTURE.context(cell_datetime_objects=True):
        cube = _slice_cube(cube, tx1, tx2)

    return cube


def multi_model_stats(cubes, span, filename, exclude):
    """Compute multi-model mean and median."""
    logger.debug('Multi model statistics: excluding files: %s', str(exclude))

    iris.util.unify_time_units(cubes)
    selection = [
        cube for cube in cubes
        if not all(cube.attributes.get(k) in exclude[k] for k in exclude)
    ]

    if len(selection) < 2:
        logger.info("Single model in list: will not compute statistics.")
        return cubes

    # check if we have any time overlap
    ovlp = _get_overlap(selection)
    if ovlp is None:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes

    start, stop = ovlp
    utype = str(selection[0].coord('time').units)
    start_dtime, stop_dtime = [_to_datetime(t, utype) for t in (start, stop)]

    # cases
    if span == 'overlap':
        logger.debug("Using common time overlap between "
                     "models to compute statistics.")

        # assemble data
        mean_dats = np.ma.zeros(_slice_cube2(selection[0],
                                             start,
                                             stop).shape)
        med_dats = np.ma.zeros(_slice_cube2(selection[0],
                                            start,
                                            stop).shape)

        for i in range(mean_dats.shape[0]):
            time_data = [_slice_cube2(cube, start, stop)[i]
                         for cube in selection]
            mean_dats[i] = _compute_statistic(time_data,
                                              'mean')
            med_dats[i] = _compute_statistic(time_data,
                                             'median')
        c_mean = _put_in_cube(_apply_overlap(selection[0],
                                             start_dtime,
                                             stop_dtime),
                              mean_dats,
                              'means',
                              filename['file_mean'],
                              t_axis=None)
        c_med = _put_in_cube(_apply_overlap(selection[0],
                                            start_dtime,
                                            stop_dtime),
                             med_dats,
                             'medians',
                             filename['file_median'],
                             t_axis=None)

    elif span == 'full':
        logger.debug("Using full time spans "
                     "to compute statistics.")
        # assemble data
        slices, time_axis = _full_time(selection)
        mean_dats = np.ma.zeros(slices[0].shape)
        med_dats = np.ma.zeros(slices[0].shape)

        for i in range(slices[0].shape[0]):
            time_data = [data[i] for data in slices]
            mean_dats[i] = _compute_statistic(time_data,
                                              'mean')
            med_dats[i] = _compute_statistic(time_data,
                                             'median')
        c_mean = _put_in_cube(selection[0],
                              mean_dats,
                              'means',
                              filename['file_mean'],
                              t_axis=time_axis)
        c_med = _put_in_cube(selection[0],
                             med_dats,
                             'medians',
                             filename['file_median'],
                             t_axis=time_axis)
    # save up
    save_cubes([c_mean])
    save_cubes([c_med])

    return cubes
