"""multimodel statistics

Functions for multi-model operations
supports a multitude of multimodel statistics
computations; the only requisite is the ingested
cubes have (TIME-LAT-LON) or (TIME-PLEV-LAT-LON)
dimensions; and obviously consistent units.

"""

import logging
from functools import reduce

import os
from datetime import datetime as dd
from datetime import timedelta as td
import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _plev_fix(dataset, pl_idx):
    """Extract valid plev data

    this function takes care of situations
    in which certain plevs are completely
    masked due to unavailable interpolation
    boundaries.
    """
    if np.ma.is_masked(dataset) is True:
        # keep only the valid plevs
        if not np.all(dataset.mask[:, pl_idx]):
            statj = np.ma.array(dataset[:, pl_idx],
                                mask=dataset.mask[:, pl_idx])
            return statj
    else:
        mask = np.zeros_like(dataset[:, pl_idx])
        statj = np.ma.array(dataset[:, pl_idx], mask=mask)
        return statj


def _compute_means(cubes):
    """Compute multimodel means"""
    datas = [c.data for c in cubes]

    # plevs
    if len(datas[0].shape) == 4:
        statistic = np.ma.zeros(datas[0].shape)
        for j in range(datas[0].shape[1]):
            stat_j = []
            for dataset in datas:
                if _plev_fix(dataset, j) is not None:
                    stat_j.append(_plev_fix(dataset, j))

            # check for nr models
            if len(stat_j) >= 2:
                statistic[:, j] = np.ma.mean(stat_j, axis=0)
            else:
                mask = np.ones_like(statistic[:, j])
                statistic[:, j] = np.ma.array(statistic[:, j],
                                              mask=mask)
    # no plevs
    else:
        statistic = np.ma.mean(datas, axis=0)

    return statistic


def _compute_medians(cubes):
    """Compute multimodel medians"""
    datas = [c.data for c in cubes]

    # plevs
    if len(datas[0].shape) == 4:
        statistic = np.ma.zeros(datas[0].shape)
        for j in range(datas[0].shape[1]):
            stat_j = []
            for dataset in datas:
                if _plev_fix(dataset, j) is not None:
                    stat_j.append(_plev_fix(dataset, j))
                # need to force nan's
                # otherwise np.ma.median is not good
                stat_l = []
                for c_dat in stat_j:
                    c_dat[c_dat.mask] = np.nan
                    stat_l.append(c_dat)

                # check for nr models
                if len(stat_l) >= 2:
                    statistic[:, j] = np.ma.median(stat_l, axis=0)
                else:
                    mask = np.ones_like(statistic[:, j])
                    statistic[:, j] = np.ma.array(statistic[:, j],
                                                  mask=mask)
                    statistic[:, j].mask[:] = np.nan

                # mask NaNs
                statmask = np.isnan(statistic[:, j])
                if np.any(statmask):
                    statistic[:, j] = np.ma.array(statistic[:, j],
                                                  mask=statmask,
                                                  fill_value=1e+20)

    # no plevs
    else:
        ncubes = []
        for cub in cubes:
            if np.ma.is_masked(cub.data) is True:
                cub.data[cub.data.mask] = np.nan
                ncubes.append(cub)
        statistic = np.ma.median([s.data for s in ncubes], axis=0)
        statmask = np.isnan(statistic)
        if np.any(statmask):
            statistic = np.ma.array(statistic,
                                    mask=statmask,
                                    fill_value=1e+20)

    return statistic


def _put_in_cube(cubes, stats_name, ncfiles, fname):
    """Quick cube building and saving"""
    # grab coordinates from any cube
    times = cubes[0].coord('time')
    lats = cubes[0].coord('latitude')
    lons = cubes[0].coord('longitude')
    if len(cubes[0].shape) == 3:
        cspec = [(times, 0), (lats, 1), (lons, 2)]
    elif len(cubes[0].shape) == 4:
        plev = cubes[0].coord('air_pressure')
        cspec = [(times, 0), (plev, 1), (lats, 2), (lons, 3)]
    # compute stats and put in cube
    if stats_name == 'means':
        dspec = _compute_means(cubes)
    elif stats_name == 'medians':
        dspec = _compute_medians(cubes)
    fixed_dspec = np.ma.fix_invalid(dspec, copy=False, fill_value=1e+20)
    stats_cube = iris.cube.Cube(fixed_dspec,
                                dim_coords_and_dims=cspec,
                                long_name=stats_name)
    coord_names = [coord.name() for coord in cubes[0].coords()]
    if 'air_pressure' in coord_names:
        if len(cubes[0].shape) == 3:
            stats_cube.add_aux_coord(cubes[0].coord('air_pressure'))
    stats_cube.attributes['_filename'] = fname
    stats_cube.attributes['NCfiles'] = str(ncfiles)
    return stats_cube


def _sdat(srl_no, unit_type):
    """Convert to a datetime point"""
    if unit_type == 'day since 1950-01-01 00:00:00.0000000':
        new_date = dd(1950, 1, 1, 0) + td(srl_no)
    elif unit_type == 'day since 1850-01-01 00:00:00.0000000':
        new_date = dd(1850, 1, 1, 0) + td(srl_no)
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
    utype = str(cubes[0].coord('time').units)
    all_times = []
    for cube in cubes:
        # monthly data
        # 1 month = 30 days ~= 0.082
        # assume time gating is already done to
        # order years
        bnd1 = float(cube.coord('time').points[0])
        bnd2 = float(cube.coord('time').points[-1])
        bnd1 = int(bnd1 / 365.) * 365.
        bnd2 = (int(bnd2 / 365.) + 1) * 365.
        all_times.append(np.array([[bnd1, bnd2]]))
    bounds = [range(int(b[0][0]), int(b[-1][-1]) + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, (i for i in bounds))
    if len(time_pts) > 1:
        return _sdat(time_pts[0], utype), _sdat(time_pts[-1], utype)


def _slice_cube(cube, min_t, max_t):
    """slice cube on time"""
    fmt = '%Y-%m-%d-%H'
    ctr = iris.Constraint(time=lambda x:
                          min_t <= dd.strptime(x.point.strftime(fmt),
                                               fmt) <= max_t)
    cube_slice = cube.extract(ctr)
    return cube_slice


def multi_model_mean(cubes, span, filename, exclude):
    """Compute multi-model mean and median."""
    logger.debug('Multi model statistics: excluding files: %s', str(exclude))
    selection = [
        cube for cube in cubes
        if not all(cube.attributes.get(k) in exclude[k] for k in exclude)
    ]

    if len(selection) < 2:
        logger.info("Single model in list: will not compute statistics.")
        return cubes

    # check if we have any time overlap
    if _get_overlap(cubes) is None:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes
    else:
        # add file name info
        file_names = [
            os.path.basename(cube.attributes.get('_filename'))
            for cube in cubes
        ]

        # look at overlap type
        if span == 'overlap':
            tmeans = []
            tx1, tx2 = _get_overlap(cubes)
            for cube in selection:
                logger.debug("Using common time overlap between "
                             "models to compute statistics.")
                logger.debug("Bounds: %s and %s", str(tx1), str(tx2))
                # slice cube on time
                with iris.FUTURE.context(cell_datetime_objects=True):
                    cube = _slice_cube(cube, tx1, tx2)

                # record data
                tmeans.append(cube)

#        elif span == 'full':
#            logger.debug("Using full time spans "
#                         "to compute statistics.")
#            tmeans = []
#            # lay down time points
#            tpts = [c.coord('time').points for c in selection]
#            t_x = list(set().union(*tpts))
#            t_x.sort()
#
#            # loop through cubes and apply masks
#            for cube in selection:
#                # construct new shape
#                fine_shape = tuple([len(t_x)] + list(cube.data.shape[1:]))
#                # find indices of present time points
#                oidx = [t_x.index(s) for s in cube.coord('time').points]
#                # reshape data to include all possible times
#                ndat = np.ma.resize(cube.data, fine_shape)
#
#                # build the time mask
#                c_ones = np.ones(fine_shape, bool)
#                for t_i in oidx:
#                    c_ones[t_i] = False
#                ndat.mask |= c_ones
#
#                # build the new coords
#                # preserve units
#                time_c = iris.coords.DimCoord(
#                    t_x, standard_name='time',
#                    units=cube.coord('time').units)
#                lat_c = cube.coord('latitude')
#                lon_c = cube.coord('longitude')
#
#                # time-lat-lon
#                if len(fine_shape) == 3:
#                    cspec = [(time_c, 0), (lat_c, 1), (lon_c, 2)]
#                # time-plev-lat-lon
#                elif len(fine_shape) == 4:
#                    pl_c = cube.coord('air_pressure')
#                    cspec = [(time_c, 0), (pl_c, 1), (lat_c, 2), (lon_c, 3)]
#
#                # build cube
#                ncube = iris.cube.Cube(ndat, dim_coords_and_dims=cspec)
#                coord_names = [coord.name() for coord in cube.coords()]
#                if 'air_pressure' in coord_names:
#                    if len(fine_shape) == 3:
#                        ncube.add_aux_coord(cube.coord('air_pressure'))
#                ncube.attributes = cube.attributes
#                tmeans.append(ncube)

        else:
            logger.debug("No type of time overlap specified "
                         "- will not compute cubes statistics")
            return cubes

    c_mean = _put_in_cube(tmeans, 'means', file_names, filename)

    c_med = _put_in_cube(tmeans, 'medians', file_names, filename)

    save_cubes([c_mean, c_med])

    return cubes
