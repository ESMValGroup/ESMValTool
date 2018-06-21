"""Common plot functions."""
import logging
import os
import yaml

import iris.quickplot
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def get_path_to_mpl_style(style_file):
    """Get path to matplotlib style file."""
    if not isinstance(style_file, str):
        raise TypeError("Invalid input: {} is not ".format(style_file) +
                        "a string")
    base_dir = os.path.dirname(__file__)
    filepath = os.path.join(base_dir, 'styles_python', 'matplotlib',
                            style_file)
    logger.debug("Using matplotlib style: %s", filepath)
    return filepath


def get_dataset_style(dataset, style_file='cmip5.yml'):
    """Retrieve the style information for the given dataset."""
    # Default path
    base_dir = os.path.dirname(__file__)
    default_dir = os.path.join(base_dir, 'styles_python')

    # Check if style_file is valid
    filepath = os.path.join(default_dir, style_file)
    if os.path.isfile(filepath):
        with open(filepath, 'r') as infile:
            style = yaml.safe_load(infile)
    else:
        raise IOError("Invalid input: could not open style file " +
                      "'{}'".format(filepath))
    logger.debug("Using style file %s for dataset %s", filepath, dataset)

    # Check if file has entry for unknown dataset
    default_dataset = 'default'
    options = ['color', 'dash', 'thick', 'mark', 'avgstd', 'facecolor']
    if default_dataset not in style:
        raise IOError("Style file '{}' does not ".format(filepath) +
                      "contain default information for unknown datasets")
    for option in options:
        if option not in style[default_dataset]:
            raise IOError("Style file '{}' ".format(filepath) +
                          "does not contain '{}' ".format(option) +
                          "default information for unknown datasets")

    # Check if dataset is available
    if not style.get(dataset):
        logger.warning("Dataset '%s' not found in style file, using default " +
                       "entry", dataset)
        return style[default_dataset]

    # Get compulsory information
    for option in options:
        if option not in style[dataset]:
            logger.warning("No style information '%s' found for dataset " +
                           "'%s', using default value for unknown datasets",
                           option, dataset)
            style[dataset].update({option: style[default_dataset][option]})

    return style[dataset]


def quickplot(cube, filename, plot_type, **kwargs):
    """Plot a cube using one of the iris.quickplot functions."""
    logger.debug("Creating '%s' plot %s", plot_type, filename)
    plot_function = getattr(iris.quickplot, plot_type)
    fig = plt.figure()
    plot_function(cube, **kwargs)
    # plt.gca().coastlines()
    fig.savefig(filename)
