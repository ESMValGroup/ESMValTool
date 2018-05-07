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
from datetime import datetime, timedelta
from functools import reduce

import cf_units
import iris
import numpy as np
import yaml

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _get_time_offset(time_unit):
    """Return a datetime object equivalent to tunit"""
    # tunit e.g. 'day since 1950-01-01 00:00:00.0000000 UTC'
    cfunit = cf_units.Unit(time_unit, calendar=cf_units.CALENDAR_STANDARD)
    time_offset = cfunit.num2date(0)
    return time_offset


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


def _put_in_cube(template_cube, cube_data, stat_name,
                 file_name, time_bounds, t_axis):
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

    metadata = {'model': 'MultiModel' + stat_name.title(),
                'filename': file_name}
    metadata_template = yaml.safe_load(template_cube.attributes['metadata'])
    for attr in ('short_name', 'standard_name', 'long_name', 'units', 'field',
                 'start_year', 'end_year', 'diagnostic', 'preprocessor'):
        if attr in metadata_template:
            metadata[attr] = metadata_template[attr]
            metadata['start_year'] = time_bounds[0]
            metadata['end_year'] = time_bounds[1]
    stats_cube.attributes['metadata'] = yaml.safe_dump(metadata)
    # complete metadata
    stats_cube.var_name = template_cube.var_name
    stats_cube.long_name = template_cube.long_name
    stats_cube.standard_name = template_cube.standard_name
    stats_cube.units = template_cube.units
    return stats_cube


def _datetime_to_int_days(cube):
    """Return list of int(days) converted from cube datetime cells"""
    # TODO replace the block when using iris 2.0
    # time_cells = [cell.point for cell in cube.coord('time').cells()]
    time_cells = [cube.coord('time').units.num2date(cell.point)
                  for cell in cube.coord('time').cells()]
    time_unit = cube.coord('time').units.name
    time_offset = _get_time_offset(time_unit)

    # extract date info
    real_dates = []
    for date_obj in time_cells:
        # real_date resets the actual data point day
        # to the 1st of the month so that there are no
        # wrong overlap indeces
        # NOTE: this workaround is good only
        # for monthly data
        real_date = datetime(date_obj.year, date_obj.month, 1, 0, 0, 0)
        real_dates.append(real_date)

    days = [(date_obj - time_offset).days for date_obj in real_dates]
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
        span = _datetime_to_int_days(cube)
        start, stop = span[0], span[-1]
        all_times.append([start, stop])
    bounds = [range(b[0], b[-1] + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, bounds)
    if len(time_pts) > 1:
        time_bounds_list = [time_pts[0], time_pts[-1]]
        return time_bounds_list


def _slice_cube(cube, t_1, t_2):
    """
    Efficient slicer

    Simple cube data slicer on indices
    of common time-data elements
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
    days = {day for cube in cubes for day in _datetime_to_int_days(cube)}
    return sorted(days)


def _full_time_slice(cubes, ndat, indices, ndatarr, t_idx):
    """Construct a contiguous collection over time"""
    for idx_cube, cube in enumerate(cubes):
        # reset mask
        ndat.mask = True
        ndat[indices[idx_cube]] = cube.data
        ndat.mask[indices[idx_cube]] = cube.data.mask
        ndatarr[idx_cube] = ndat[t_idx]

    # return time slice
    return ndatarr


def _assemble_overlap_data(cubes, ovlp, stat_type, filename, time_bounds):
    """Get statistical data in iris cubes for OVERLAP"""
    start, stop = ovlp
    sl_1, sl_2 = _slice_cube(cubes[0], start, stop)
    stats_dats = np.ma.zeros(cubes[0].data[sl_1:sl_2 + 1].shape)

    for i in range(stats_dats.shape[0]):
        indices = [_slice_cube(cube, start, stop) for cube in cubes]
        time_data = [
            cube.data[indx[0]:indx[1] + 1][i]
            for cube, indx in zip(cubes, indices)
        ]
        stats_dats[i] = _compute_statistic(time_data, stat_type)
    stats_cube = _put_in_cube(
        cubes[0][sl_1:sl_2 + 1], stats_dats, stat_type, filename,
        time_bounds, t_axis=None)
    return stats_cube


def _assemble_full_data(cubes, stat_type, filename, time_bounds):
    """Get statistical data in iris cubes for FULL"""
    # all times, new MONTHLY data time axis
    time_axis = [float(fl) for fl in _monthly_t(cubes)]

    # new big time-slice array shape
    new_shape = [len(time_axis)] + list(cubes[0].shape[1:])

    # assemble an array to hold all time data
    # for all cubes; shape is (ncubes,(plev), lat, lon)
    new_arr = np.ma.empty([len(cubes)] + list(new_shape[1:]))

    # data array for stats computation
    stats_dats = np.ma.zeros(new_shape)

    # assemble indices list to chop new_arr on
    indices_list = []

    # empty data array to hold time slices
    empty_arr = np.ma.empty(new_shape)

    # loop through cubes and populate empty_arr with points
    for cube in cubes:
        time_redone = _datetime_to_int_days(cube)
        oidx = [time_axis.index(s) for s in time_redone]
        indices_list.append(oidx)
    for i in range(new_shape[0]):
        # hold time slices only
        new_datas_array = _full_time_slice(cubes, empty_arr, indices_list,
                                           new_arr, i)
        # list to hold time slices
        time_data = []
        for j in range(len(cubes)):
            time_data.append(new_datas_array[j])
        stats_dats[i] = _compute_statistic(time_data, stat_type)
    stats_cube = _put_in_cube(cubes[0], stats_dats, stat_type, filename,
                              time_bounds, time_axis)
    return stats_cube


def _update_filename(filename, interval, time_unit):
    """Update netCDF file names based on time properties"""
    start, stop = [(_get_time_offset(time_unit) + timedelta(int(ts))).year
                   for ts in interval]
    filename = "{}_{}-{}.nc".format(filename.rpartition('_')[0], start, stop)
    return filename, start, stop


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
    interval = _get_overlap(selection)
    if interval is None:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes

    time_unit = selection[0].coord('time').units.name

    # cases
    files = []
    if span == 'overlap':
        logger.debug("Using common time overlap between "
                     "models to compute statistics.")

        # assemble data
        for stat_name in statistics:
            filename, startT, stopT = _update_filename(filenames[stat_name],
                                                       interval,
                                                       time_unit)
            time_bounds = [startT, stopT]
            cube_of_stats = _assemble_overlap_data(selection, interval,
                                                   stat_name, filename,
                                                   time_bounds)
            cube_of_stats.data = np.ma.array(cube_of_stats.data,
                                             dtype=np.dtype('float32'))
            save_cubes([cube_of_stats])
            files.append(filename)

    elif span == 'full':
        logger.debug("Using full time spans " "to compute statistics.")
        # assemble data
        time_points = _monthly_t(selection)
        interval = [min(time_points), max(time_points)]
        for stat_name in statistics:
            filename, startT, stopT = _update_filename(filenames[stat_name],
                                                       interval,
                                                       time_unit)
            time_bounds = [startT, stopT]
            cube_of_stats = _assemble_full_data(selection, stat_name, filename,
                                                time_bounds)
            cube_of_stats.data = np.ma.array(cube_of_stats.data,
                                             dtype=np.dtype('float32'))
            save_cubes([cube_of_stats])
            files.append(filename)

    cubes.extend(files)
    return cubes
