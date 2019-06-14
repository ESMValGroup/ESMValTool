"""multimodel statistics.

Functions for multi-model operations
supports a multitude of multimodel statistics
computations; the only requisite is the ingested
cubes have (TIME-LAT-LON) or (TIME-PLEV-LAT-LON)
dimensions; and obviously consistent units.

It operates on different (time) spans:
- full: computes stats on full dataset time;
- overlap: computes common time overlap between datasets;

"""

import logging
from datetime import datetime
from functools import reduce

import cf_units
import iris
import numpy as np

from .._config import use_legacy_iris

logger = logging.getLogger(__name__)


def _get_time_offset(time_unit):
    """Return a datetime object equivalent to tunit."""
    # tunit e.g. 'day since 1950-01-01 00:00:00.0000000 UTC'
    cfunit = cf_units.Unit(time_unit, calendar=cf_units.CALENDAR_STANDARD)
    time_offset = cfunit.num2date(0)
    return time_offset


def _plev_fix(dataset, pl_idx):
    """Extract valid plev data.

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


def _compute_statistic(datas, statistic_name):
    """Compute multimodel statistic."""
    datas = np.ma.array(datas)
    statistic = datas[0]

    if statistic_name == 'median':
        statistic_function = np.ma.median
    elif statistic_name == 'mean':
        statistic_function = np.ma.mean
    else:
        raise NotImplementedError

    # no plevs
    if len(datas[0].shape) < 3:
        # get all NOT fully masked data - u_data
        # datas is per time point
        # so we can safely NOT compute stats for single points
        if datas.ndim == 1:
            u_datas = [data for data in datas]
        else:
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

        # check for nr datasets
        if len(plev_check) > 1:
            plev_check = np.ma.array(plev_check)
            statistic[j] = statistic_function(plev_check, axis=0)
        else:
            statistic.mask[j] = True

    return statistic


def _put_in_cube(template_cube, cube_data, statistic, t_axis):
    """Quick cube building and saving."""
    if t_axis is None:
        times = template_cube.coord('time')
    else:
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
    elif len(template_cube.shape) == 1:
        cspec = [
            (times, 0),
        ]
    elif len(template_cube.shape) == 2:
        # If you're going to hardwire air_pressure into this,
        # might as well have depth here too.
        plev = template_cube.coord('depth')
        cspec = [
            (times, 0),
            (plev, 1),
        ]

    # correct dspec if necessary
    fixed_dspec = np.ma.fix_invalid(cube_data, copy=False, fill_value=1e+20)
    # put in cube
    stats_cube = iris.cube.Cube(
        fixed_dspec, dim_coords_and_dims=cspec, long_name=statistic)
    coord_names = [coord.name() for coord in template_cube.coords()]
    if 'air_pressure' in coord_names:
        if len(template_cube.shape) == 3:
            stats_cube.add_aux_coord(template_cube.coord('air_pressure'))

    stats_cube.var_name = template_cube.var_name
    stats_cube.long_name = template_cube.long_name
    stats_cube.standard_name = template_cube.standard_name
    stats_cube.units = template_cube.units
    return stats_cube


def _datetime_to_int_days(cube):
    """Return list of int(days) converted from cube datetime cells."""
    if use_legacy_iris():
        time_cells = [
            cube.coord('time').units.num2date(cell.point)
            for cell in cube.coord('time').cells()
        ]
    else:
        time_cells = [cell.point for cell in cube.coord('time').cells()]

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
    Efficient slicer.

    Simple cube data slicer on indices
    of common time-data elements.
    """
    time_pts = [t for t in cube.coord('time').points]
    converted_t = _datetime_to_int_days(cube)
    idxs = sorted([
        time_pts.index(ii) for ii, jj in zip(time_pts, converted_t)
        if t_1 <= jj <= t_2
    ])
    return [idxs[0], idxs[-1]]


def _monthly_t(cubes):
    """Rearrange time points for monthly data."""
    # get original cubes tpoints
    days = {day for cube in cubes for day in _datetime_to_int_days(cube)}
    return sorted(days)


def _full_time_slice(cubes, ndat, indices, ndatarr, t_idx):
    """Construct a contiguous collection over time."""
    for idx_cube, cube in enumerate(cubes):
        # reset mask
        ndat.mask = True
        ndat[indices[idx_cube]] = cube.data
        if np.ma.is_masked(cube.data):
            ndat.mask[indices[idx_cube]] = cube.data.mask
        else:
            ndat.mask[indices[idx_cube]] = False
        ndatarr[idx_cube] = ndat[t_idx]

    # return time slice
    return ndatarr


def _assemble_overlap_data(cubes, interval, statistic):
    """Get statistical data in iris cubes for OVERLAP."""
    start, stop = interval
    sl_1, sl_2 = _slice_cube(cubes[0], start, stop)
    stats_dats = np.ma.zeros(cubes[0].data[sl_1:sl_2 + 1].shape)

    # keep this outside the following loop
    # this speeds up the code by a factor of 15
    indices = [_slice_cube(cube, start, stop) for cube in cubes]

    for i in range(stats_dats.shape[0]):
        time_data = [
            cube.data[indx[0]:indx[1] + 1][i]
            for cube, indx in zip(cubes, indices)
        ]
        stats_dats[i] = _compute_statistic(time_data, statistic)
    stats_cube = _put_in_cube(
        cubes[0][sl_1:sl_2 + 1], stats_dats, statistic, t_axis=None)
    return stats_cube


def _assemble_full_data(cubes, statistic):
    """Get statistical data in iris cubes for FULL."""
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
        stats_dats[i] = _compute_statistic(time_data, statistic)
    stats_cube = _put_in_cube(cubes[0], stats_dats, statistic, time_axis)
    return stats_cube


def multi_model_statistics(products, span, output_products, statistics):
    """
    Compute multi-model statistics.

    Multimodel statistics computed along the time axis. Can be
    computed across a common overlap in time (set span: overlap)
    or across the full length in time of each model (set span: full).
    Restrictive compuation is also available by excluding any set of
    models that the user will not want to include in the statistics
    (set exclude: [excluded models list]).

    Restrictions needed by the input data:
    - model datasets must have consistent shapes,
    - higher dimesnional data is not supported (ie dims higher than four:
    time, vertical axis, two horizontal axes).

    Parameters
    ----------
    products: list
        list of data products to be used in multimodel stat computation;
        cube attribute of product is the data cube for computing the stats.
    span: str
        overlap or full; if overlap stas are computed on common time-span;
        if full stats are computed on full time spans.
    output_products: dict
        dictionary of output products.
    statistics: str
        statistical measure to be computed (mean or median).
    Returns
    -------
    list
        list of data products containing the multimodel stats computed.
    Raises
    ------
    ValueError
        If span is neither overlap nor full.

    """
    logger.debug('Multimodel statistics: computing: %s', statistics)
    if len(products) < 2:
        logger.info("Single dataset in list: will not compute statistics.")
        return products

    cubes = [cube for product in products for cube in product.cubes]
    # check if we have any time overlap
    interval = _get_overlap(cubes)
    if interval is None:
        logger.info("Time overlap between cubes is none or a single point."
                    "check datasets: will not compute statistics.")
        return products

    if span == 'overlap':
        logger.debug("Using common time overlap between "
                     "datasets to compute statistics.")
    elif span == 'full':
        logger.debug("Using full time spans to compute statistics.")
    else:
        raise ValueError(
            "Unexpected value for span {}, choose from 'overlap', 'full'"
            .format(span))

    statistic_products = set()
    for statistic in statistics:
        # Compute statistic
        if span == 'overlap':
            statistic_cube = _assemble_overlap_data(cubes, interval, statistic)
        elif span == 'full':
            statistic_cube = _assemble_full_data(cubes, statistic)
        statistic_cube.data = np.ma.array(
            statistic_cube.data, dtype=np.dtype('float32'))

        # Add to output product and log provenance
        statistic_product = output_products[statistic]
        statistic_product.cubes = [statistic_cube]
        for product in products:
            statistic_product.wasderivedfrom(product)
        logger.info("Generated %s", statistic_product)
        statistic_products.add(statistic_product)

    products |= statistic_products

    return products
