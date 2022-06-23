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
    """Function to manually determine sub-plot positions in a 3x3 figure.

    Parameters
    ----------
    j (int): index of cube position in cubelist

    Returns
    -------
    x, y (int): subplot positions
    """
    if j <= 2:
        y = j
        x = 0
    elif 2 < j <= 5:
        y = j - 3
        x = 1
    else:
        y = j - 6
        x = 2

    return x, y


def plot_patterns(cube_list, plot_path):
    """Plots climate patterns for imogen_mode: off.

    Parameters
    ----------
    cube_list (cubelist): input cubelist for plotting patterns per variable
    plot_path (path): path to plot_dir

    Returns
    -------
    None
    """
    fig, ax = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig.suptitle("Patterns from a random grid-cell", fontsize=18, y=0.98)

    plt.figure(figsize=(14, 12))
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("Global Patterns, January", fontsize=18, y=0.95)

    for j, cube in enumerate(cube_list):
        # determining plot positions
        x, y = subplot_positions(j)
        months = np.arange(1, 13)
        # plots patterns for a random grid cell
        ax[x, y].plot(months, cube[:, 50, 50].data)
        ax[x, y].set_ylabel(str(cube.var_name) + " / " + str(cube.units))
        if j > 5:
            ax[x, y].set_xlabel('Time')

        # January patterns
        plt.subplot(3, 3, j + 1)
        qplt.pcolormesh(cube[0])

    plt.tight_layout()
    plt.savefig(plot_path + 'Patterns')
    plt.close()

    fig.tight_layout()
    fig.savefig(plot_path + 'Patterns Timeseries')


def plot_cp_timeseries(list_cubelists, plot_path):
    """Plots timeseries and maps of climatologies, anomalies and patterns.

    Parameters
    ----------
    list_cubelists (cubelist): input cubelist for plotting per variable
    plot_path (path): path to plot_dir

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

    for i, _ in enumerate(list_cubelists):
        cube_list = list_cubelists[i]
        for j, cube in enumerate(cube_list):
            # determining plot positions
            x, y = subplot_positions(j)
            yrs = (1850 + np.arange(cube.shape[0])).astype('float')
            if i == 0:
                # climatology
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax1[x, y].plot(yrs, avg_cube)
                ax1[x, y].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax1[x, y].set_xlabel('Time')
            if i == 1:
                # anomaly timeseries
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax2[x, y].plot(yrs, avg_cube)
                ax2[x, y].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax2[x, y].set_xlabel('Time')
            if i == 2:
                months = np.arange(1, 13)
                # avg_cube = sf.area_avg(cube, return_cube=False)
                ax3[x, y].plot(months, cube[:, 50, 50].data)
                ax3[x, y].set_ylabel(
                    str(cube.var_name) + " / " + str(cube.units))
                if j > 5:
                    ax3[x, y].set_xlabel('Time')
            if i == 2:
                # January patterns
                plt.subplot(3, 3, j + 1)
                qplt.pcolormesh(cube[0])

    plt.tight_layout()
    plt.savefig(plot_path + 'Patterns')
    plt.close()

    fig1.tight_layout()
    fig1.savefig(plot_path + 'Climatologies')

    fig2.tight_layout()
    fig2.savefig(plot_path + 'Anomalies')

    fig3.tight_layout()
    fig3.savefig(plot_path + 'Patterns Timeseries')


def plot_scores(cube_list, plot_path):
    """Plots color mesh of scores per variable per month.

    Parameters
    ----------
    cube_list (cube): input cubelist for plotting
    plot_path (path): path to plot_dir

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
        plt.savefig(plot_path + 'R2_Scores_' + str(cube.var_name))
        plt.close()

    # plot global scores timeseries per variable
    plt.figure(figsize=(5, 8))
    for cube in cube_list:
        score = sf.area_avg(cube, return_cube=True)
        iplt.plot(score, label=cube.var_name)
        plt.xlabel('Time')
        plt.ylabel('R2 Score')
        plt.legend(loc='center left')

    plt.savefig(plot_path + 'score_timeseries')
    plt.close()


def ebm_plots(reg, avg_list, yrs, forcing, forcing_126, lambda_c, temp_global,
              tas_delta, temp_global_126, tas_126_delta, plot_path):
    """Plots the regression between tas and rtmt, as well as the resultant
    forcing timeseries.

    Parameters
    ----------
    reg (obj): linear regression outputs with attributes of tas_4x_cube and
        rtmt_4x_cube
    avg_list (arr): array of variable anomaly timeseries
    yrs (arr): array of timeseries years for plotting
    forcing (arr): derived effective forcing
    forcing_126 (arr): effective forcing timeseries, SSP126
    lambda_c (float): derived climate sensitivity, from reg.slope
    temp_global (arr): EBM predicted global near-surface air temperature
    tas_delta (cube): near-surface air temperature anomaly
    temp_global_126 (cube): EBM predicted global near-surface
        air temperature, SSP126
    tas_126_delta (cube): near-surface air temperature anomaly, SSP126
    plot_path (path): path to plot_dir

    Returns
    -------
    None
    """
    yrs = (1850 + np.arange(temp_global.shape[0])).astype('float')

    # regression line
    x_reg = np.linspace(-3.0, 4.0, 2)
    y_reg = reg.slope * x_reg + reg.intercept

    fig, ax = plt.subplots(3, 1, figsize=(6, 12))

    # plotting regression line
    ax[0].scatter(avg_list[2].data, avg_list[3].data, c='b', label='tas, rtmt')
    ax[0].plot(x_reg, y_reg, c='r', label='regression')
    ax[0].legend(loc='upper left')
    ax[0].set_ylabel(
        "Net Downward Radiative Flux at Top of Atmosphere / $Wm^{-2}$")
    ax[0].set_xlabel("Near Surface Air Temperature / $K$")

    # plotting forcing timeseries
    ax[1].plot(yrs, forcing, color='red', label='SSP585')
    ax[1].plot(yrs, forcing_126, color='blue', label='SSP126')
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("Radiative forcing (Wm-2)")
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    ld = r"$\lambda$"
    ax[1].text(0.2,
               0.9,
               f"{ld} = {lambda_c:.2f}",
               transform=ax[1].transAxes,
               ha="center",
               va="center",
               size=11,
               bbox=bbox_props)

    # plotting ebm prediction
    ax[2].plot(yrs,
               temp_global,
               color='black',
               zorder=10,
               linewidth=1.5,
               label='EBM Prediction SSP585')
    ax[2].plot(yrs,
               temp_global_126,
               color='purple',
               zorder=10,
               linewidth=1.5,
               label='EBM Prediction SSP126')
    ax[2].plot(yrs, tas_delta, color='red', label='Model SSP585')
    ax[2].plot(yrs, tas_126_delta, color='blue', label='Model SSP126')
    ax[2].legend(loc='upper left')
    ax[2].set_xlabel("Time")
    ax[2].set_ylabel("Air Surface Temperature (K)")

    fig.savefig(plot_path + 'ebm_plots', dpi=200)
    plt.close(fig)
