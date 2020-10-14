#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-stamp: <2020-02-05 15:59:30 lukbrunn>

(c) 2019 under a MIT License (https://mit-license.org)

Authors:
- Lukas Brunner || lukas.brunner@env.ethz.ch

Abstract:

"""
import numpy as np
import seaborn as sns
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

from weighted_temperature_graph import weighted_quantile

def boxplot(ax,
            pos=None,
            data=None,
            median=None,
            mean=None,
            box=None,
            whis=None,
            dots=None,
            weights=None,
            width=.8,
            color=sns.xkcd_rgb['greyish'],
            alpha=1.,
            showcaps=True,
            fancy_legend=False,
            return_handle=False,
            zorder=100,
            box_quantiles=(.25, .75),
            whis_quantiles=(.05, .95),
            median_kwargs=None,
            mean_kwargs=None,
            whis_kwargs=None,
            dots_kwargs=None,
            dots_sizes=None):
    """
    A custom-mad boxplot routine based on user set statistics.

    Parameters
    ----------
    ax : plt.axis object
    pos : float, optional
        Location (center) of the box.
    data: array-like, optional
        An array of data to calculate the statistics from. If this is
        not None, median, mean, box, and whis will be overwritten.
    median : float or array-like, optional
        Location of the median or data to calculate the median.
    mean : float or array-like, optional
        Location of the mean or data to calculate the mean.
    box : tuple of float or array-like, shape (2), optional
        Extend of the box or data to calculate the box.
    whis : tuple of float or array-like, shape (2), optional
        Extend of the whiskers or data to calculate the whiskers.
    dots : array-like, optional
        An array of dots to plot.
    weights : list, optional
        weights to calculate weighted quantiles if applies
    width : float, optional
        Width of box, median, mean, and caps (caps have .4*width).
    color : string, optional
        Box color and default color for median, mean, and whiskers.
    showcaps : bool, optional
        Whether to draw caps at the end of the whiskers.
    zorder : int, optional
        zorder of the drawn objects.
    median_kwargs : dict, optional
        Keyword arguments passed on to ax.hlines for the median.
    mean_kwargs : dict, optional
        Keyword arguments passed on to ax.hlines for the mean.
    whis_kwargs : dict, optional
        Keyword arguments passed on to ax.hlines and ax.vlines for whiskers
        and caps.
    dots_kwargs : dict, optional
        Keyword arguments passed on to ax.scatter for the dots.
    dots_sizes : tuple of (min, max), optional
    """
    if data is not None:
        if mean is None:
            mean = data
        if median is None:
            median = data
        if box is None:
            box = data
        if whis is None:
            whis = data
        if dots is None:
            dots = data

    if mean is not None and not isinstance(mean, (int, float)):
        mean = np.average(mean, weights=weights)
    if median is not None and not isinstance(median, (int, float)):
        median = weighted_quantile(median, .5, weights)
    if box is not None and len(box) != 2 and not isinstance(box, tuple):
        box = weighted_quantile(box, box_quantiles, weights)
    elif tuple(box) == (None, None):
        box = None
    if whis is not None and len(whis) != 2 and not isinstance(whis, tuple):
        whis = weighted_quantile(whis, whis_quantiles, weights)
    elif whis is not None and tuple(whis) == (None, None):
        whis = None

    if pos is None:
        pos = 0
    if median_kwargs is None:
        median_kwargs = {}
    if mean_kwargs is None:
        mean_kwargs = {}
    if whis_kwargs is None:
        whis_kwargs = {}
    if dots_kwargs is None:
        dots_kwargs = {}
    if 'colors' not in median_kwargs.keys():
        median_kwargs['colors'] = 'k'
    if 'colors' not in mean_kwargs.keys():
        mean_kwargs['colors'] = 'k'
    if 'colors' not in whis_kwargs.keys():
        whis_kwargs['colors'] = color
    if 'alpha' not in median_kwargs.keys():
        median_kwargs['alpha'] = alpha
    if 'alpha' not in mean_kwargs.keys():
        mean_kwargs['alpha'] = alpha
    if 'alpha' not in whis_kwargs.keys():
        whis_kwargs['alpha'] = alpha
    if 'caps_width' in whis_kwargs.keys():
        caps_width = whis_kwargs.pop('caps_width')
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
        if weights is not None and dots_sizes is not None:
            sizes = np.interp(weights, [np.min(weights), np.max(weights)], dots_sizes)
            dots_kwargs['s'] = sizes**2

    handle = [mpatches.Patch(color=color, alpha=alpha)]
    if median is not None:
        handle.append(ax.hlines([], [], [], **median_kwargs))
    if mean is not None:
        handle.append(ax.hlines([], [], [], **mean_kwargs))
    if dots is not None:
        handle.append(ax.scatter([], [], **dots_kwargs))
    if return_handle:
        return tuple(handle)

    x0, x1 = pos - .5*width, pos + .5*width
    if box is not None:  # plot box
        patch = PatchCollection(
            [Rectangle((x0, box[0]), width, box[1] - box[0])],
            facecolor=color,
            alpha=alpha,
            zorder=zorder)
        ax.add_collection(patch)

    if median is not None:
        x0_median = x0 + (1 - median_width)*width*.5
        x1_median = x1 - (1 - median_width)*width*.5
        ax.hlines(median, x0_median, x1_median, zorder=zorder, **median_kwargs)

    if mean is not None:  # plot mean
        x0_mean = x0 + (1 - mean_width)*width*.5
        x1_mean = x1 - (1 - mean_width)*width*.5
        ax.hlines(mean, x0_mean, x1_mean, zorder=zorder, **mean_kwargs)

    if dots is not None:
        if len(dots) == 1:
            xx = pos
        else:
            xx = np.random.RandomState(0).uniform(pos - .4*width, pos + .4*width, len(dots))
        patch = ax.scatter(xx, dots, zorder=zorder, **dots_kwargs)

    if whis is not None:  # plot whiskers
        if box is None:
            box = (whis[0], whis[0])
        ax.vlines((pos, pos), (whis[0], box[1]), (box[0], whis[1]),
                  zorder=zorder, **whis_kwargs)
        if showcaps:  # plot caps
            x0, x1 = pos - .5*caps_width*width, pos + .5*caps_width*width
            ax.hlines(whis, (x0, x0), (x1, x1), zorder=zorder, **whis_kwargs)

    return tuple(handle)
