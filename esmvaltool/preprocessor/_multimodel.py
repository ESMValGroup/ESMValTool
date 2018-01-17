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
        t1 = time_pts[0]
        t2 = time_pts[-1]
        return t1, t2
    else:
        return 0, 0


def multi_model_mean(cubes, span, filename, exclude):
    """Compute multi-model mean and median."""

    logger.debug('Multi model statistics: excluding: %s' % str(exclude))
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
    tx1, tx2 = _get_overlap(cubes)

    # check if we have any time overlap
    if tx1 == tx2:
        logger.info("Time overlap between cubes is none or a single point.")
        logger.info("check models: will not compute statistics.")
        return cubes

    for cube in selection:
        for coord in cube.coords():
            if coord.standard_name == 'time':
                if span == 'overlap':

                    logger.debug("Using common time overlap between models to compute statistics.")
                    # find the nearest points to the overlap region
                    all_times = list(cube.coord('time').points)
                    a1 = min(all_times, key=lambda x: abs(x - tx1))
                    a2 = min(all_times, key=lambda x: abs(x - tx2))
                    i1 = all_times.index(a1)
                    i2 = all_times.index(a2)
                    logger.debug("Indexing time axis for overlap on indices %i and %i of %i" % (i1,i2,len(all_times)-1))

                    # index data on these points and
                    # get rid of masked values
                    flat_data = cube.data[i1:i2+1,:,:]
                    cdata = flat_data[~flat_data.mask]
                    data_size.append(cdata.shape[0])

                    # compute stats and append to lists
                    cmean = np.mean(cdata)
                    cmedian = np.median(cdata)
                    means.append(cmean)
                    medians.append(cmedian)

                elif span == 'full':

                    logger.debug("Using full time spans to compute statistics.")
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

    logger.debug("Global medians: %s", medians)
    medians_cube = iris.cube.Cube(medians, long_name='medians')
    medians_cube.attributes['_filename'] = filename

    logger.debug("Data sizes: %s", data_size)
    datasize_cube = iris.cube.Cube(data_size, long_name='data_size')
    datasize_cube.attributes['_filename'] = filename

    save_cubes([means_cube, medians_cube, datasize_cube])

    return cubes
