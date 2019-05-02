"""Common plot functions."""
import logging
import os

import iris.quickplot
import matplotlib.pyplot as plt
import yaml

logger = logging.getLogger(__name__)


def _process_axes_functions(axes, axes_functions):
    """Process axes functions of the form `axes.functions(*args, **kwargs)."""
    if axes_functions is None:
        return None
    output = None
    for (func, attr) in axes_functions.items():
        axes_function = getattr(axes, func)

        # Simple functions (argument directly given)
        if not isinstance(attr, dict):
            try:
                out = axes_function(*attr)
            except TypeError:
                out = axes_function(attr)

        # More complicated functions (args and kwargs given)
        else:
            args = attr.get('args', [])
            kwargs = attr.get('kwargs', {})

            # Process 'transform' kwargs
            if 'transform' in kwargs:
                kwargs['transform'] = getattr(axes, kwargs['transform'])
            out = axes_function(*args, **kwargs)

        # Return legend if possible
        if func == 'legend':
            output = out
    return output


def _check_size_of_parameters(*args):
    """Check if the size of (array-like) args is identical."""
    if len(args) < 2:
        logger.warning("Less than two arguments given, comparing not possible")
        return
    arg_0 = args[0]
    for arg in args:
        try:
            if len(arg_0) != len(arg):
                raise ValueError("Invalid input: array-like parameters need "
                                 "to have the same size")
        except TypeError:
            raise TypeError("Invalid input: some parameters are not "
                            "array-like")
    return


def get_path_to_mpl_style(style_file=None):
    """Get path to matplotlib style file."""
    if style_file is None:
        style_file = 'default.mplstyle'
    if not style_file.endswith('.mplstyle'):
        style_file += '.mplstyle'
    base_dir = os.path.dirname(os.path.realpath(__file__))
    filepath = os.path.join(base_dir, 'styles_python', 'matplotlib',
                            style_file)
    logger.debug("Using matplotlib style: %s", filepath)
    return filepath


def get_dataset_style(dataset, style_file=None):
    """Retrieve the style information for the given dataset."""
    if style_file is None:
        style_file = 'cmip5.yml'
    if not style_file.endswith('.yml'):
        style_file += '.yml'
    base_dir = os.path.dirname(os.path.realpath(__file__))
    default_dir = os.path.join(base_dir, 'styles_python')

    # Check if style_file is valid
    filepath = os.path.join(default_dir, style_file)
    if os.path.isfile(filepath):
        with open(filepath, 'r') as infile:
            style = yaml.safe_load(infile)
    else:
        raise IOError("Invalid input: could not open style file "
                      "'{}'".format(filepath))
    logger.debug("Using style file %s for dataset %s", filepath, dataset)

    # Check if file has entry for unknown dataset
    default_dataset = 'default'
    options = ['color', 'dash', 'thick', 'mark', 'avgstd', 'facecolor']
    if default_dataset not in style:
        raise IOError("Style file '{}' does not contain default information "
                      "for unknown datasets".format(filepath))
    for option in options:
        if option not in style[default_dataset]:
            raise IOError("Style file '{}' does not contain '{}' default "
                          "information for unknown "
                          "datasets".format(filepath, option))

    # Check if dataset is available
    if not style.get(dataset):
        logger.warning(
            "Dataset '%s' not found in style file, using default "
            "entry", dataset)
        return style[default_dataset]

    # Get compulsory information
    for option in options:
        if option not in style[dataset]:
            logger.warning(
                "No style information '%s' found for dataset "
                "'%s', using default value for unknown datasets", option,
                dataset)
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
    plt.close(fig)


def multi_dataset_scatterplot(x_data, y_data, datasets, filepath, **kwargs):
    """Plot a multi dataset scatterplot.

    Notes
    -----
    Allowed keyword arguments:

    * `mpl_style_file` (:obj:`str`):  Path to the matplotlib style file.

    * `dataset_style_file` (:obj:`str`): Path to the dataset style file.

    * `plot_kwargs` (`array-like`): Keyword arguments for the plot (e.g.
      `label`, `makersize`, etc.).

    * `save_kwargs` (:obj:`dict`): Keyword arguments for saving the plot.

    * `axes_functions` (:obj:`dict`): Arbitrary functions for axes, i.e.
      `axes.set_title('title')`.

    Parameters
    ----------
    x_data : array-like
        x data of each dataset.
    y_data : array-like
        y data of each dataset.
    datasets : array-like
        Names of the datasets.
    filepath : str
        Path to which plot is written.
    **kwargs
        Keyword arguments.

    Raises
    ------
    TypeError
        A non-valid keyword argument is given or `x_data`, `y_data`, `datasets`
        or (if given) `plot_kwargs` is not array-like.
    ValueError
        `x_data`, `y_data`, `datasets` or `plot_kwargs` do not have the same
         size.

    """
    # Allowed kwargs
    allowed_kwargs = [
        'mpl_style_file',
        'dataset_style_file',
        'plot_kwargs',
        'save_kwargs',
        'axes_functions',
    ]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            raise TypeError("{} is not a valid keyword argument".format(kwarg))

    # Check parameters
    _check_size_of_parameters(x_data, y_data, datasets,
                              kwargs.get('plot_kwargs', x_data))
    empty_dict = [{} for _ in x_data]

    # Create matplotlib instances
    plt.style.use(get_path_to_mpl_style(kwargs.get('mpl_style_file')))
    (fig, axes) = plt.subplots()

    # Plot data
    for (idx, dataset) in enumerate(datasets):
        style = get_dataset_style(dataset, kwargs.get('dataset_style_file'))

        # Fix problem when plotting ps file
        facecolor = style['color'] if filepath.endswith('ps') else \
            style['facecolor']

        # Plot
        axes.plot(
            x_data[idx],
            y_data[idx],
            markeredgecolor=style['color'],
            markerfacecolor=facecolor,
            marker=style['mark'],
            **(kwargs.get('plot_kwargs', empty_dict)[idx]))

    # Costumize plot
    legend = _process_axes_functions(axes, kwargs.get('axes_functions'))

    # Save plot
    fig.savefig(
        filepath, additional_artists=[legend], **kwargs.get('save_kwargs', {}))
    logger.info("Wrote %s", filepath)
    plt.close()


def scatterplot(x_data, y_data, filepath, **kwargs):
    """Plot a scatterplot.

    Notes
    -----
    Allowed keyword arguments:

    * `mpl_style_file` (:obj:`str`):  Path to the matplotlib style file.

    * `plot_kwargs` (`array-like`): Keyword arguments for the plot (e.g.
      `label`, `makersize`, etc.).

    * `save_kwargs` (:obj:`dict`): Keyword arguments for saving the plot.

    * `axes_functions` (:obj:`dict`): Arbitrary functions for axes, i.e.
      `axes.set_title('title')`.

    Parameters
    ----------
    x_data : array-like
        x data of each dataset.
    y_data : array-like
        y data of each dataset.
    filepath : str
        Path to which plot is written.
    **kwargs
        Keyword arguments.

    Raises
    ------
    TypeError
        A non-valid keyword argument is given or `x_data`, `y_data` or (if
        given) `plot_kwargs` is not array-like.
    ValueError
        `x_data`, `y_data` or `plot_kwargs` do not have the same size.

    """
    # Allowed kwargs
    allowed_kwargs = [
        'mpl_style_file',
        'plot_kwargs',
        'save_kwargs',
        'axes_functions',
    ]
    for kwarg in kwargs:
        if kwarg not in allowed_kwargs:
            raise TypeError("{} is not a valid keyword argument".format(kwarg))

    # Check parameters
    _check_size_of_parameters(x_data, y_data, kwargs.get(
        'plot_kwargs', x_data))
    empty_dict = [{} for _ in x_data]

    # Create matplotlib instances
    plt.style.use(get_path_to_mpl_style(kwargs.get('mpl_style_file')))
    (fig, axes) = plt.subplots()

    # Plot data
    for (idx, x_vals) in enumerate(x_data):
        plot_kwargs = kwargs.get('plot_kwargs', empty_dict)[idx]

        # Fix problem when plotting ps file
        if 'markerfacecolor' in plot_kwargs and filepath.endswith('ps'):
            plot_kwargs.pop('markerfacecolor')

        # Plot
        axes.plot(x_vals, y_data[idx],
                  **(kwargs.get('plot_kwargs', empty_dict)[idx]))

    # Costumize plot
    legend = _process_axes_functions(axes, kwargs.get('axes_functions'))

    # Save plot
    fig.savefig(
        filepath, additional_artists=[legend], **kwargs.get('save_kwargs', {}))
    logger.info("Wrote %s", filepath)
    plt.close()
