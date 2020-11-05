"""A collection of utility functions for dealing with weights."""
from collections import defaultdict

import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds['weight']


def read_metadata(cfg: dict, groupby: str = 'variable_group') -> dict:
    """Read the metadata from the config file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item[groupby]

        datasets[variable].append(item)

    return datasets


def weighted_quantile(values: list,
                      quantiles: list,
                      weights: list = None) -> 'np.array':
    """Calculate weighted quantiles.

    Analogous to np.quantile, but supports weights.

    Based on: https://stackoverflow.com/a/29677616/6012085

    Parameters
    ----------
    values: array_like
        List of input values.
    quantiles: array_like
        List of quantiles between 0.0 and 1.0.
    weights: array_like
        List with same length as `values` containing the weights.

    Returns
    -------
    np.array
        Numpy array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if weights is None:
        weights = np.ones(len(values))
    weights = np.array(weights)

    if not np.all((quantiles >= 0) & (quantiles <= 1)):
        raise ValueError('Quantiles should be between 0.0 and 1.0')

    idx = np.argsort(values)
    values = values[idx]
    weights = weights[idx]

    weighted_quantiles = np.cumsum(weights) - 0.5 * weights

    # Cast weighted quantiles to 0-1 To be consistent with np.quantile
    min_val = weighted_quantiles.min()
    max_val = weighted_quantiles.max()
    weighted_quantiles = (weighted_quantiles - min_val) / max_val

    return np.interp(quantiles, weighted_quantiles, values)


def calculate_percentiles(data: 'xr.DataArray',
                          percentiles: list,
                          weights: dict = None) -> 'xr.DataArray':
    """Calculate (weighted) percentiles.

    Calculate the (weighted) percentiles for the given data.

    Percentiles is a list of values between 0 and 100.

    The `model_ensemble` dimension in weights has to contain at
    least the same elements as in data.
    If `weights` is not specified, the non-weighted percentiles are calculated.

    Returns a DataArray with 'percentiles' as the dimension.
    """
    if weights is not None:
        weights = weights.sel(model_ensemble=data.model_ensemble)

    output = xr.apply_ufunc(weighted_quantile,
                            data,
                            input_core_dims=[['model_ensemble']],
                            output_core_dims=[['percentiles']],
                            kwargs={
                                'weights': weights,
                                'quantiles': percentiles / 100
                            },
                            vectorize=True)

    output['percentiles'] = percentiles

    return output

def boxplot(axes,
            pos=None,
            data=None,
            median=None,
            mean=None,
            box=None,
            whisk=None,
            dots=None,
            weights=None,
            width=.8,
            color=sns.xkcd_rgb['greyish'],
            alpha=1.,
            showcaps=True,
            median_kwargs=None,
            mean_kwargs=None,
            whisk_kwargs=None,
            dots_kwargs=None,
            dots_sizes=None):
    """Custom-made boxplot routine based on user set statistics.

    Parameters
    ----------
    axes : plt.axis object
    pos : float, optional
        Location (center) of the box.
    data: array-like, optional
        An array of data to calculate the statistics from. If this is
        not None, median, mean, box, and whisk will be overwritten.
    median : float or array-like, optional
        Location of the median or data to calculate the median.
    mean : float or array-like, optional
        Location of the mean or data to calculate the mean.
    box : tuple of float or array-like, shape (2), optional
        Extend of the box or data to calculate the box.
    whisk : tuple of float or array-like, shape (2), optional
        Extend of the whiskers or data to calculate the whiskers.
    dots : array-like, optional
        An array of dots to plot.
    weights : list, optional
        weights to calculate weighted quantiles if applies
    width : float, optional
        Width of box, median, mean, and caps (caps have .4*width).
    color : string, optional
        Box color and default color for median, mean, and whiskers.
    alpha : float, optional
        A number between 0. and 1. to determine opacity of drawn objects
    showcaps : bool, optional
        Whether to draw caps at the end of the whiskers.
    median_kwargs : dict, optional
        Keyword arguments passed on to axes.hlines for the median.
    mean_kwargs : dict, optional
        Keyword arguments passed on to axes.hlines for the mean.
    whisk_kwargs : dict, optional
        Keyword arguments passed on to axes.hlines and axes.vlines for whiskers
        and caps.
    dots_kwargs : dict, optional
        Keyword arguments passed on to axes.scatter for the dots.
    dots_sizes : tuple of (min, max), optional
    """
    if data is not None:
        if mean is None:
            mean = data
        if median is None:
            median = data
        if box is None:
            box = data
        if whisk is None:
            whisk = data
        if dots is None:
            dots = data

    if mean is not None and not isinstance(mean, (int, float)):
        mean = np.average(mean, weights=weights)
    if median is not None and not isinstance(median, (int, float)):
        median = weighted_quantile(median, .5, weights)
    if box is not None and len(box) != 2 and not isinstance(box, tuple):
        box = weighted_quantile(box, (.25, .75), weights)
    elif tuple(box) == (None, None):
        box = None
    if whisk is not None and len(whisk) != 2 and not isinstance(whisk, tuple):
        whisk = weighted_quantile(whisk, (.05, .95), weights)
    elif whisk is not None and tuple(whisk) == (None, None):
        whisk = None

    if pos is None:
        pos = 0
    if median_kwargs is None:
        median_kwargs = {}
    if mean_kwargs is None:
        mean_kwargs = {}
    if whisk_kwargs is None:
        whisk_kwargs = {}
    if dots_kwargs is None:
        dots_kwargs = {}
    if 'colors' not in median_kwargs.keys():
        median_kwargs['colors'] = 'k'
    if 'colors' not in mean_kwargs.keys():
        mean_kwargs['colors'] = 'k'
    if 'colors' not in whisk_kwargs.keys():
        whisk_kwargs['colors'] = color
    if 'alpha' not in median_kwargs.keys():
        median_kwargs['alpha'] = alpha
    if 'alpha' not in mean_kwargs.keys():
        mean_kwargs['alpha'] = alpha
    if 'alpha' not in whisk_kwargs.keys():
        whisk_kwargs['alpha'] = alpha
    if 'caps_width' in whisk_kwargs.keys():
        caps_width = whisk_kwargs.pop('caps_width')
    else:
        caps_width = .4
    if 'width' in mean_kwargs.keys():
        mean_width = mean_kwargs.pop('width')
    else:
        mean_width = 1.
    if 'width' in median_kwargs.keys():
        median_width = median_kwargs.pop('width')
    else:
        median_width = 1.
    if 'linestyle' not in mean_kwargs.keys():
        if median is not None:
            mean_kwargs['linestyle'] = '--'
    if 'color' not in dots_kwargs.keys() and 'c' not in dots_kwargs.keys():
        dots_kwargs['color'] = 'k'
    if 's' not in dots_kwargs.keys():
        dots_kwargs['s'] = 1

    zorder = 100

    handle = [mpatches.Patch(color=color, alpha=alpha)]
    if median is not None:
        handle.append(axes.hlines([], [], [], **median_kwargs))
    if mean is not None:
        handle.append(axes.hlines([], [], [], **mean_kwargs))
    if dots is not None:
        handle.append(axes.scatter([], [], **dots_kwargs))

    x_0, x_1 = pos - .5 * width, pos + .5 * width
    if box is not None:  # plot box
        patch = PatchCollection(
            [Rectangle((x_0, box[0]), width, box[1] - box[0])],
            facecolor=color,
            alpha=alpha,
            zorder=zorder)
        axes.add_collection(patch)

    if median is not None:
        x0_median = x_0 + (1 - median_width) * width * .5
        x1_median = x_1 - (1 - median_width) * width * .5
        axes.hlines(median,
                    x0_median,
                    x1_median,
                    zorder=zorder,
                    **median_kwargs)

    if mean is not None:  # plot mean
        x0_mean = x_0 + (1 - mean_width) * width * .5
        x1_mean = x_1 - (1 - mean_width) * width * .5
        axes.hlines(mean, x0_mean, x1_mean, zorder=zorder, **mean_kwargs)

    if dots is not None:
        if len(dots) == 1:
            x_pos = pos
        else:
            x_pos = np.random.RandomState(0).uniform(pos - .4 * width,
                                                     pos + .4 * width,
                                                     len(dots))
        if weights is not None and dots_sizes is not None:
            sizes = np.interp(
                weights, [np.min(weights), np.max(weights)], dots_sizes)
            dots_kwargs['s'] = sizes**2
        patch = axes.scatter(x_pos, dots, zorder=zorder, **dots_kwargs)

    if whisk is not None:  # plot whiskers
        if box is None:
            box = (whisk[0], whisk[0])
        axes.vlines((pos, pos), (whisk[0], box[1]), (box[0], whisk[1]),
                    zorder=zorder,
                    **whisk_kwargs)
        if showcaps:  # plot caps
            x_0 = pos - .5 * caps_width * width
            x_1 = pos + .5 * caps_width * width
            axes.hlines(whisk, (x_0, x_0), (x_1, x_1),
                        zorder=zorder,
                        **whisk_kwargs)

    return tuple(handle)
