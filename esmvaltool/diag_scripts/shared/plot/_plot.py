"""Common plot functions."""
import logging

import iris.quickplot
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def quickplot(cube, filename, plot_type, **kwargs):
    """Plot a cube using one of the iris.quickplot functions."""
    logger.debug("Creating '%s' plot %s", plot_type, filename)
    plot_function = getattr(iris.quickplot, plot_type)
    fig = plt.figure()
    plot_function(cube, **kwargs)
    plt.gca().coastlines()
    fig.savefig(filename)
