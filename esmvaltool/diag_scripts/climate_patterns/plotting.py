# (C) Crown Copyright 2022-2024, Met Office.
"""Script containing plotting functions for driving scripts.

Author
------
Gregory Munday (Met Office, UK)
"""

import os

import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import area_statistics


def subplot_positions(j):
    """Determine sub-plot positions in a 3x3 figure.

    Parameters
    ----------
    j : int
        index of cube position in cubelist

    Returns
    -------
    x_pos : int
        x subplot position
    y_pos : int
        y subplot position
    """
    if j <= 2:
        y_pos = j
        x_pos = 0
    elif 2 < j <= 5:
        y_pos = j - 3
        x_pos = 1
    else:
        y_pos = j - 6
        x_pos = 2

    return x_pos, y_pos


def plot_patterns(cube_list, plot_path):
    """Plot climate patterns for jules_mode: off.

    Parameters
    ----------
    cube_list : cubelist
        input cubelist for plotting patterns per variable
    plot_path : path
        path to plot_dir

    Returns
    -------
    None
    """
    fig, axis = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig.suptitle("Patterns from a random grid-cell", fontsize=18, y=0.98)

    plt.figure(figsize=(14, 12))
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("Global Patterns, January", fontsize=18, y=0.95)

    for j, cube in enumerate(cube_list):
        # determining plot positions
        x_pos, y_pos = subplot_positions(j)
        months = np.arange(1, 13)
        # plots patterns for an arbitrary grid cell
        axis[x_pos, y_pos].plot(months, cube[:, 50, 50].data)
        axis[x_pos, y_pos].set_ylabel(
            str(cube.var_name) + " / " + str(cube.units)
        )
        if j > 5:
            axis[x_pos, y_pos].set_xlabel("Time")

        # January patterns
        plt.subplot(3, 3, j + 1)
        qplt.pcolormesh(cube[0])

    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, "Patterns"), dpi=300)
    plt.close()

    fig.tight_layout()
    fig.savefig(os.path.join(plot_path, "Patterns Timeseries"), dpi=300)


def plot_timeseries(cubelist, plot_path, title, save_name):
    """Plot timeseries and maps of climatologies, anomalies and patterns.

    Parameters
    ----------
    cubelist : cubelist
        input cubelist for plotting per variable
    plot_path : path
        path to plot_dir
    title: str
        title for the figure
    save_name: str
        name for the saved figure

    Returns
    -------
    None
    """
    fig, axs = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig.suptitle(f"{title}", fontsize=18, y=0.98)

    for j, cube in enumerate(cubelist):
        # determining plot positions
        x_pos, y_pos = subplot_positions(j)
        yrs = (1850 + np.arange(cube.shape[0])).astype("float")
        months = np.arange(1, 13)

        # anomaly timeseries
        avg_cube = area_statistics(cube, "mean").data
        if save_name == "Climatologies":
            axs[x_pos, y_pos].plot(months, avg_cube)
        else:
            axs[x_pos, y_pos].plot(yrs, avg_cube)
        axs[x_pos, y_pos].set_ylabel(cube.long_name + " / " + str(cube.units))
        if j > 5:
            axs[x_pos, y_pos].set_xlabel("Time")

    fig.tight_layout()
    fig.savefig(os.path.join(plot_path, f"{save_name}"), dpi=300)
