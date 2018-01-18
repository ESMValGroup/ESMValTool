"""Functions for multi-model operations"""
import logging
from functools import reduce

import iris
import numpy as np

from ._io import save_cubes

logger = logging.getLogger(__name__)


def _get_overlap(cubes):
    """Get discrete time overlaps."""
    all_times = [cube.coord('time').points for cube in cubes]
    bounds = [range(int(np.min(b)), int(np.max(b)) + 1) for b in all_times]
    time_pts = reduce(np.intersect1d, (i for i in bounds))
    if len(time_pts) > 1:
        time1 = time_pts[0]
        time2 = time_pts[-1]
        return time1, time2
    else:
        return 0, 0


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

    means = []
    medians = []
    data_size = []
    file_names = []
    tx1, tx2 = _get_overlap(cubes)

    # check if we have any time overlap
    if tx1 == tx2:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes

    for cube in selection:
        file_name = cube.attributes.get('_filename').split('/')[-1]
        file_names.append(file_name)
        for coord in cube.coords():
            if coord.standard_name == 'time':
                if span == 'overlap':

                    logger.debug("Using common time overlap between \
                                 models to compute statistics.")
                    # find the nearest points to the overlap region
                    all_times = list(cube.coord('time').points)
                    min_t = min(all_times, key=lambda x: abs(x - tx1))
                    max_t = min(all_times, key=lambda x: abs(x - tx2))
                    id_min = all_times.index(min_t)
                    id_max = all_times.index(max_t)
                    tot = len(all_times) - 1
                    logger.debug("Indexing time axis for overlap on \
                                 indices %i and %i \
                                 of %i", id_min, id_max, tot)

                    # index data on these points and
                    # get rid of masked values
                    flat_data = cube.data[id_min:id_max + 1, :, :]
                    cdata = flat_data[~flat_data.mask]
                    data_size.append(cdata.shape[0])

                    # compute stats and append to lists
                    cmean = np.mean(cdata)
                    cmedian = np.median(cdata)
                    means.append(cmean)
                    medians.append(cmedian)

                elif span == 'full':

                    logger.debug("Using full time spans \
                                 to compute statistics.")
                    # get rid of masked
                    flat_data = cube.data
                    cdata = flat_data[~flat_data.mask]
                    data_size.append(cdata.shape[0])

                    # compute stats and append to lists
                    cmean = np.mean(cdata)
                    cmedian = np.median(cdata)
                    means.append(cmean)
                    medians.append(cmedian)

    logger.debug("Global means: %s", means)
    means_cube = iris.cube.Cube(means, long_name='means')
    means_cube.attributes['_filename'] = filename
    means_cube.attributes['NCfiles'] = str(file_names)

    logger.debug("Global medians: %s", medians)
    medians_cube = iris.cube.Cube(medians, long_name='medians')
    medians_cube.attributes['_filename'] = filename
    medians_cube.attributes['NCfiles'] = str(file_names)

    logger.debug("Data sizes: %s", data_size)
    datasize_cube = iris.cube.Cube(data_size, long_name='data_size')
    datasize_cube.attributes['_filename'] = filename
    datasize_cube.attributes['NCfiles'] = str(file_names)

    save_cubes([means_cube, medians_cube, datasize_cube])

    return cubes
