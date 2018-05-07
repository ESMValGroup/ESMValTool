"""Common plot functions."""
import logging

import iris.quickplot as qplt
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def example_map_plot(cube, filename):
    """Plot a cube on a world map."""
    logger.debug("Creating %s", filename)
    fig = plt.figure()
    qplt.pcolormesh(cube)
    plt.gca().coastlines()
    fig.savefig(filename)
