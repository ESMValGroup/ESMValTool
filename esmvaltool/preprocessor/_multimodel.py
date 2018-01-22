"""Functions for multi-model operations"""
import logging
from functools import reduce

import os
import datetime
import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _put_in_cube(stats, stats_name, ncfiles, fname):
    """quick cube building"""
    stats_cube = iris.cube.Cube(stats, long_name=stats_name)
    stats_cube.attributes['_filename'] = fname
    stats_cube.attributes['NCfiles'] = str(ncfiles)
    return stats_cube


def _sdat(srl_no):
    """convert to a datatime point"""
    new_date = datetime.datetime(1950, 1, 1, 0) +\
        datetime.timedelta(srl_no)
    return new_date


def _get_overlap(cubes):
    """Get discrete time overlaps."""
    all_times = [cube.coord('time').points for cube in cubes]
    bounds = [range(int(b[0]), int(b[-1]) + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, (i for i in bounds))
    if len(time_pts) > 1:
        print(_sdat(time_pts[0]), _sdat(time_pts[-1]))
        return _sdat(time_pts[0]), _sdat(time_pts[-1])


def _slice_cube(cube, min_t, max_t):
    """slice cube on time"""
    irisConstr = 
        iris.Constraint(time=lambda x:
                        min_t <= datetime.datetime.strptime(x.point.strftime('%Y-%m-%d-%H'),
                                                            '%Y-%m-%d-%H') <= max_t)
    cube_slice = cube.extract(irisConstr)
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
        # empty lists to hold data
        tmeans = []
        data_size = []

        # add file name info
        file_names = [
            os.path.basename(cube.attributes.get('_filename'))
            for cube in cubes
        ]

        # look at overlap type
        if span == 'overlap':
            tx1, tx2 = _get_overlap(cubes)
            print(tx1, tx2)
            for cube in selection:
                logger.debug("Using common time overlap between "
                             "models to compute statistics.")
                # slice cube on time
                with iris.FUTURE.context(cell_datetime_objects=True):
                    cube = _slice_cube(cube, tx1, tx2)

                # store nr time points
                data_size.append(len(list(cube.coord('time').points)))

                # compute stats and append to lists
                tmeans.append(cube.collapsed('time', iris.analysis.MEAN))

        elif span == 'full':
            for cube in selection:
                logger.debug("Using full time spans "
                             "to compute statistics.")

                # store nr time points
                data_size.append(len(list(cube.coord('time').points)))

                # compute stats and append to lists
                tmeans.append(cube.collapsed('time', iris.analysis.MEAN))

        else:
            logger.debug("No type of time overlap specified "
                         "- will not compute cubes statistics")
            return cubes

    c_mean = _put_in_cube(np.mean([s.data for s in tmeans], axis=0),
                          'means', file_names, filename)

    c_med = _put_in_cube(np.median([s.data for s in tmeans], axis=0),
                         'medians', file_names, filename)

    c_dat = _put_in_cube(data_size, 'data_size', file_names, filename)

    save_cubes([c_mean, c_med, c_dat])

    return cubes
