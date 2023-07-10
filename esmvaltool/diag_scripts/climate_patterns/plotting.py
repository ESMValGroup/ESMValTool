"""Script containing plotting functions for driving scripts.

Author
------
Gregory Munday (Met Office, UK)
"""

import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import sub_functions as sf


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
    """Plot climate patterns for imogen_mode: off.

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
        axis[x_pos,
             y_pos].set_ylabel(str(cube.var_name) + " / " + str(cube.units))
        if j > 5:
            axis[x_pos, y_pos].set_xlabel("Time")

        # January patterns
        plt.subplot(3, 3, j + 1)
        qplt.pcolormesh(cube[0])

    plt.tight_layout()
    plt.savefig(plot_path + "Patterns")
    plt.close()

    fig.tight_layout()
    fig.savefig(plot_path + "Patterns Timeseries")


def plot_cp_timeseries(list_cubelists, plot_path):
    """Plot timeseries and maps of climatologies, anomalies and patterns.

    Parameters
    ----------
    list_cubelists : cubelist
        input cubelist for plotting per variable
    plot_path : path
        path to plot_dir

    Returns
    -------
    None
    """
    fig1, ax1 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig1.suptitle("40 Year Climatologies, 1850-1889", fontsize=18, y=0.98)

    fig2, ax2 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig2.suptitle("Anomaly Timeseries, 1850-2100", fontsize=18, y=0.98)

    fig3, ax3 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig3.suptitle("Patterns from a random grid-cell", fontsize=18, y=0.98)

    plt.figure(figsize=(14, 12))
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("Global Patterns, January", fontsize=18, y=0.95)

    for i, cube_list in enumerate(list_cubelists):
        for j, cube in enumerate(cube_list):
            # determining plot positions
            x_pos, y_pos = subplot_positions(j)
            yrs = (1850 + np.arange(cube.shape[0])).astype("float")
            if i == 0:
                # climatology
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax1[x_pos, y_pos].plot(yrs, avg_cube)
                ax1[x_pos,
                    y_pos].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax1[x_pos, y_pos].set_xlabel("Time")
            if i == 1:
                # anomaly timeseries
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax2[x_pos, y_pos].plot(yrs, avg_cube)
                ax2[x_pos,
                    y_pos].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax2[x_pos, y_pos].set_xlabel("Time")
            if i == 2:
                months = np.arange(1, 13)
                # avg_cube = sf.area_avg(cube, return_cube=False)
                ax3[x_pos, y_pos].plot(months, cube[:, 50, 50].data)
                ax3[x_pos, y_pos].set_ylabel(
                    str(cube.var_name) + " / " + str(cube.units))
                if j > 5:
                    ax3[x_pos, y_pos].set_xlabel("Time")
            if i == 2:
                # January patterns
                plt.subplot(3, 3, j + 1)
                qplt.pcolormesh(cube[0])

    plt.tight_layout()
    plt.savefig(plot_path + "Patterns")
    plt.close()

    fig1.tight_layout()
    fig1.savefig(plot_path + "Climatologies")

    fig2.tight_layout()
    fig2.savefig(plot_path + "Anomalies")

    fig3.tight_layout()
    fig3.savefig(plot_path + "Patterns Timeseries")


def plot_scores(cube_list, plot_path):
    """Plot color mesh of scores per variable per month.

    Parameters
    ----------
    cube_list : cube
        input cubelist for plotting
    plot_path : path
        path to plot_dir

    Returns
    -------
    None
    """
    for cube in cube_list:
        plt.figure(figsize=(14, 12))
        plt.subplots_adjust(hspace=0.5)
        plt.suptitle("Scores " + cube.var_name, fontsize=18, y=0.98)
        for j in range(0, 12):
            plt.subplot(4, 3, j + 1)
            qplt.pcolormesh(cube[j])

        plt.tight_layout()
        plt.savefig(plot_path + "R2_Scores_" + str(cube.var_name))
        plt.close()

    # plot global scores timeseries per variable
    plt.figure(figsize=(5, 8))
    for cube in cube_list:
        score = sf.area_avg(cube, return_cube=True)
        iplt.plot(score, label=cube.var_name)
        plt.xlabel("Time")
        plt.ylabel("R2 Score")
        plt.legend(loc="center left")

    plt.savefig(plot_path + "score_timeseries")
    plt.close()
