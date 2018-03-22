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

from datetime import datetime
from datetime import timedelta
from dateutil import parser
import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _parse_time_unit(tunit):
    """Return a datetime object equivalent to tunit"""
    # tunit e.g. 'day since 1950-01-01 00:00:00.0000000 UTC'
    # FIXME this needs more work to account for all
    # CF_UNIT UDUNIT cases
    unit_datetime = parser.parse(' '.join([tunit.split()[2],
                                           tunit.split()[3]]))
    return unit_datetime


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
            statj = np.ma.array(dataset[pl_idx], mask=dataset.mask[pl_idx])
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
        # get all NOT fully masked data - u_data
        # datas is per time point
        # so we can safely NOT compute stats for single points
        u_datas = [data for data in datas if not np.all(data.mask)]
        if len(u_datas) > 1:
            statistic = statistic_function(datas, axis=0)
        else:
            statistic.mask = True
        return statistic

    # plevs
    for j in range(statistic.shape[0]):
        plev_check = []
        for cdata in datas:
            fixed_data = _plev_fix(cdata, j)
            if fixed_data is not None:
                plev_check.append(fixed_data)

        # check for nr models
        if len(plev_check) > 1:
            plev_check = np.ma.array(plev_check)
            statistic[j] = statistic_function(plev_check, axis=0)
        else:
            statistic.mask[j] = True

    return statistic


def _put_in_cube(template_cube, cube_data, stat_name, file_name, t_axis):
    """Quick cube building and saving"""
    # grab coordinates from any cube
    times = template_cube.coord('time')
    # or get the FULL time axis
    if t_axis is not None:
        times = iris.coords.DimCoord(
            t_axis,
            standard_name='time',
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
    fixed_dspec = np.ma.fix_invalid(cube_data, copy=False, fill_value=1e+20)
    # put in cube
    stats_cube = iris.cube.Cube(
        fixed_dspec, dim_coords_and_dims=cspec, long_name=stat_name)
    coord_names = [coord.name() for coord in template_cube.coords()]
    if 'air_pressure' in coord_names:
        if len(template_cube.shape) == 3:
            stats_cube.add_aux_coord(template_cube.coord('air_pressure'))
    stats_cube.attributes['_filename'] = file_name
    # complete metadata
    stats_cube.var_name = template_cube.var_name
    stats_cube.long_name = template_cube.long_name
    stats_cube.standard_name = template_cube.standard_name
    stats_cube.units = template_cube.units
    return stats_cube


def _datetime_to_int_days(cube):
    """Return list of int(days) converted from cube datetime cells"""
    time_cells = [cell.point
                  for cell in cube.coord('time').cells()]
    time_units = cube.coord('time').units
    unit_type = time_units.name
    # extract date info
    real_dates = []
    for date_obj in time_cells:
        # real_date resets the actual data point day
        # to the 1st of the month so that there are no
        # wrong overlap indeces
        # NOTE: this workaround is good only
        # for monthly data
        real_date = datetime(date_obj.year,
                             date_obj.month,
                             1,
                             0,
                             0,
                             0)
        real_dates.append(real_date)
    days = [(date_obj - _parse_time_unit(unit_type)).days
            for date_obj in real_dates]
    return days


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
        bnd1 = _datetime_to_int_days(cube)[0]
        bnd2 = _datetime_to_int_days(cube)[-1]
        all_times.append([bnd1, bnd2])
    bounds = [range(b[0], b[-1] + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, bounds)
    if len(time_pts) > 1:
        time_bounds_list = [time_pts[0], time_pts[-1]]
        return time_bounds_list


def _slice_cube(cube, t_1, t_2):
    """
    Efficient slicer

    Cube indexing has a large memory leak!
    using this function adds virtually
    no extra memory to the process
    """
    time_pts = [t for t in cube.coord('time').points]
    converted_t = _datetime_to_int_days(cube)
    idxs = sorted([
        time_pts.index(ii) for ii, jj in zip(time_pts, converted_t)
        if t_1 <= jj <= t_2
    ])
    return [idxs[0], idxs[-1]]


def _monthly_t(cubes):
    """Rearrange time points for monthly data"""
    # get original cubes tpoints
    tpts = []
    for cube in cubes:
        tpts.append(_datetime_to_int_days(cube))
    t_x = list(set().union(*tpts))
    t_x.sort()
    t_x0 = list({float(t) for t in t_x})
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
        time_redone = _datetime_to_int_days(cube)
        # construct new shape
        fine_shape = tuple([len(t_x)] + list(cube.data.shape[1:]))
        # find indices of present time points
        oidx = [t_x.index(s) for s in time_redone]
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


def _assemble_overlap_data(selection, ovlp, stat_type, fname):
    """Get statistical data in iris cubes for OVERLAP"""
    start, stop = ovlp
    sl_1, sl_2 = _slice_cube(selection[0], start, stop)
    stats_dats = np.ma.zeros(selection[0].data[sl_1:sl_2 + 1].shape)

    for i in range(stats_dats.shape[0]):
        indices = [_slice_cube(cube, start, stop)
                   for cube in selection]
        time_data = [cube.data[indx[0]:indx[1] + 1][i]
                     for cube, indx in zip(selection, indices)]
        stats_dats[i] = _compute_statistic(time_data, stat_type)
    stats_cube = _put_in_cube(
        selection[0][sl_1:sl_2 + 1],
        stats_dats,
        stat_type,
        fname,
        t_axis=None)
    return stats_cube


def _assemble_full_data(selection, stat_type, fname):
    """Get statistical data in iris cubes for FULL"""
    slices, time_axis = _full_time(selection)
    stats_dats = np.ma.zeros(slices[0].shape)

    for i in range(slices[0].shape[0]):
        time_data = [data[i] for data in slices]
        stats_dats[i] = _compute_statistic(time_data, stat_type)
    stats_cube = _put_in_cube(selection[0], stats_dats, stat_type, fname,
                              time_axis)
    return stats_cube


def _update_fname(curr_filename, ovlp_interval, unit_type):
    """Update netCDF file names based on time properties"""
    froot = "_".join(curr_filename.split('_')[0:-1])
    start, stop = ovlp_interval
    unit_type = unit_type.name
    yr_1 = str((_parse_time_unit(unit_type)
                + timedelta(np.int(start))).year)
    yr_2 = str((_parse_time_unit(unit_type)
                + timedelta(np.int(stop))).year)
    new_file = "{}_{}-{}.nc".format(froot, yr_1, yr_2)
    return new_file


def multi_model_statistics(cubes, span, filenames, exclude, statistics):
    """Compute multi-model mean and median."""
    logger.debug('Multi model statistics: excluding files: %s', exclude)

    logger.debug('Multimodel statistics: computing: %s', statistics)
    selection = [
        cube for cube in cubes
        if not all(cube.attributes.get(k) in exclude[k] for k in exclude)
    ]

    if len(selection) < 2:
        logger.info("Single model in list: will not compute statistics.")
        return cubes

    # unify units
    iris.util.unify_time_units(selection)

    # check if we have any time overlap
    ovlp = _get_overlap(selection)
    if ovlp is None:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes
    u_type = selection[0].coord('time').units

    # cases
    if span == 'overlap':
        logger.debug("Using common time overlap between "
                     "models to compute statistics.")

        # assemble data
        for stat_name in statistics:
            updated_fname = _update_fname(filenames[stat_name], ovlp, u_type)
            cube_of_stats = _assemble_overlap_data(selection, ovlp, stat_name,
                                                   updated_fname)
            save_cubes([cube_of_stats])

    elif span == 'full':
        logger.debug("Using full time spans " "to compute statistics.")
        # assemble data
        for stat_name in statistics:
            fovlp = [
                min(_monthly_t(selection)[1]),
                max(_monthly_t(selection)[1])
            ]
            updated_fname = _update_fname(filenames[stat_name], fovlp, u_type)
            cube_of_stats = _assemble_full_data(selection, stat_name,
                                                updated_fname)
            save_cubes([cube_of_stats])

    return cubes
