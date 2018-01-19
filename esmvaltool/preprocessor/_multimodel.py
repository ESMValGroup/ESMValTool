"""Functions for multi-model operations"""
import logging
from functools import reduce

import os
import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _find_index(flist, idm):
    """Find nearest element to idm in flist"""
    max_t = min(flist, key=lambda x: abs(x - idm))
    id_min = flist.index(max_t)
    return id_min


def _put_in_cube(stats, stats_name, ncfiles, fname):
    """quick cube building"""
    stats_cube = iris.cube.Cube(stats, long_name=stats_name)
    stats_cube.attributes['_filename'] = fname
    stats_cube.attributes['NCfiles'] = str(ncfiles)
    return stats_cube


def _get_overlap(cubes):
    """Get discrete time overlaps."""
    all_times = [cube.coord('time').points for cube in cubes]
    bounds = [range(int(np.min(b)), int(np.max(b)) + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, (i for i in bounds))
    if len(time_pts) > 1:
        return time_pts[0], time_pts[-1]


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
        means = []
        medians = []
        data_size = []

        # add file name info
        file_names = [
            os.path.basename(cube.attributes.get('_filename'))
            for cube in cubes
        ]

        # look at overlap type
        if span == 'overlap':
            tx1, tx2 = _get_overlap(cubes)
            for cube in selection:
                logger.debug("Using common time overlap between "
                             "models to compute statistics.")
                # find the nearest points to the overlap region
                id_min = _find_index(list(cube.coord('time').points), tx1)
                id_max = _find_index(list(cube.coord('time').points), tx2)
                logger.debug("Indexing time axis for overlap on "
                             "indices %i and %i "
                             "of %i", id_min, id_max,
                             len(list(cube.coord('time').points)) - 1)

                # index data on these points and
                # get rid of masked values
                flat_data = cube.data[id_min:id_max + 1, :, :]
                data_size.append(flat_data[~flat_data.mask].shape[0])

                # compute stats and append to lists
                means.append(np.mean(flat_data[~flat_data.mask]))
                medians.append(np.median(flat_data[~flat_data.mask]))

        elif span == 'full':
            for cube in selection:
                logger.debug("Using full time spans "
                             "to compute statistics.")
                # get rid of masked
                data_size.append(cube.data[~cube.data.mask].shape[0])

                # compute stats and append to lists
                means.append(np.mean(cube.data[~cube.data.mask]))
                medians.append(np.median(cube.data[~cube.data.mask]))

        else:
            logger.debug("No type of time overlap specified "
                         "- will not compute cubes statistics")
            return cubes

    logger.debug("Global means: %s", means)
    c_mean = _put_in_cube(means, 'means', file_names, filename)

    logger.debug("Global medians: %s", medians)
    c_med = _put_in_cube(medians, 'medians', file_names, filename)

    logger.debug("Data sizes: %s", data_size)
    c_dat = _put_in_cube(data_size, 'data_size', file_names, filename)

    save_cubes([c_mean, c_med, c_dat])

    return cubes
