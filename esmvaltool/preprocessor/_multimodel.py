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
    t1 = time_pts[0]
    t2 = time_pts[-1]
    return t1, t2


def multi_model_mean(cubes, span, filename, exclude=None):
    """Compute multi-model mean and median."""
    selection = [
        cube for cube in cubes
        if not all(cube.attributes.get(k) in exclude[k] for k in exclude)
    ]

    if len(selection) < 2:
        logger.info("Single model in list: will not compute statistics.")
        return cubes

    means = []
    medians = []
    tx1, tx2 = _get_overlap(cubes)

    for cube in selection:
        for coord in cube.coords():
            if coord.standard_name == 'time':
                if span == 'overlap':
                    all_times = list(cube.coord('time').points)
                    a1 = min(all_times, key=lambda x: abs(x - tx1))
                    a2 = min(all_times, key=lambda x: abs(x - tx2))
                    i1 = all_times.index(a1)
                    i2 = all_times.index(a2)
                    cmean = np.mean(cube.data[i1:i2, :, :])
                    cmedian = np.median(cube.data[i1:i2, :, :])
                    means.append(cmean)
                    medians.append(cmedian)
                elif span == 'full':
                    cmean = np.mean(cube.data)
                    cmedian = np.median(cube.data)
                    means.append(cmean)
                    medians.append(cmedian)

    logger.info("Global means: %s", means)
    means_cube = iris.cube.Cube(means, long_name='means')
    means_cube.attributes['_filename'] = filename

    logger.info("Global medians: %s", medians)
    medians_cube = iris.cube.Cube(medians, long_name='medians')
    medians_cube.attributes['_filename'] = filename

    save_cubes([means_cube, medians_cube])

    return cubes
