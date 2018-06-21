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


def get_model_style(model, style_file='cmip5.yml'):
    """Retrieve the style information for the given model."""
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
    logger.debug("Using style file %s for model %s", filepath, model)

    # Check if file has entry for unknown model
    default_model = 'default'
    options = ['color', 'dash', 'thick', 'mark', 'avgstd', 'facecolor']
    if default_model not in style:
        raise IOError("Style file '{}' does not ".format(filepath) +
                      "contain default information for unknown models")
    for option in options:
        if option not in style[default_model]:
            raise IOError("Style file '{}' ".format(filepath) +
                          "does not contain '{}' ".format(option) +
                          "default information for unknown models")

    # Check if model is available
    if not style.get(model):
        logger.warning("Model '%s' not found in style file, using default " +
                       "entry", model)
        return style[default_model]

    # Get compulsory information
    for option in options:
        if option not in style[model]:
            logger.warning("No style information '%s' found for model '%s', " +
                           "using default value for unknown models",
                           option, model)
            style[model].update({option: style[default_model][option]})

    return style[model]


def quickplot(cube, filename, plot_type, **kwargs):
    """Plot a cube using one of the iris.quickplot functions."""
    logger.debug("Creating '%s' plot %s", plot_type, filename)
    plot_function = getattr(iris.quickplot, plot_type)
    fig = plt.figure()
    plot_function(cube, **kwargs)
    # plt.gca().coastlines()
    fig.savefig(filename)
