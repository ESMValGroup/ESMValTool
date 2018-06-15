"""Common plot functions."""
import logging
import os

import iris.quickplot
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def get_path_to_mpl_style(style_file):
    """
    Get path to matplotlib style file.
    """
    if (not isinstance(style_file, str)):
        raise TypeError("Invalid input: {} is not ".format(style_file) +
                        "a string")
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'styles_python', 'matplotlib', style_file)


def quickplot(cube, filename, plot_type, **kwargs):
    """Plot a cube using one of the iris.quickplot functions."""
    logger.debug("Creating '%s' plot %s", plot_type, filename)
    plot_function = getattr(iris.quickplot, plot_type)
    fig = plt.figure()
    plot_function(cube, **kwargs)
    plt.gca().coastlines()
    fig.savefig(filename)
