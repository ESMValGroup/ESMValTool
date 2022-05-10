import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import sub_functions as sf


def subplot_positions(j):
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


def plot_cp_timeseries(list_cubelists, plot_path):
    """Plots timeseries of aggregated cubes, across all scenarios."""
    fig1, ax1 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig1.suptitle("40 Year Climatologies, 1850-1889", fontsize=18, y=0.98)

    fig2, ax2 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig2.suptitle("Anomaly Timeseries, 1850-2100", fontsize=18, y=0.98)

    fig3, ax3 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig3.suptitle("Variable Timeseries, 1850-2100", fontsize=18, y=0.98)

    fig4, ax4 = plt.subplots(3, 3, figsize=(14, 12), sharex=True)
    fig4.suptitle("Aggregated Global Patterns", fontsize=18, y=0.98)

    plt.figure(figsize=(14, 12))
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("Aggregated Global Patterns, January", fontsize=18, y=0.95)

    for i in range(len(list_cubelists)):
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
                # full timeseries
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax2[x, y].plot(yrs, avg_cube)
                ax2[x, y].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax2[x, y].set_xlabel('Time')
            if i == 2:
                # anomalies
                avg_cube = sf.area_avg(cube, return_cube=False)
                ax3[x, y].plot(yrs, avg_cube)
                ax3[x, y].set_ylabel(cube.long_name + " / " + str(cube.units))
                if j > 5:
                    ax3[x, y].set_xlabel('Time')
            if i == 3:
                months = np.arange(1, 13)
                # avg_cube = sf.area_avg(cube, return_cube=False)
                ax4[x, y].plot(months, cube[:, 50, 50].data)
                ax4[x, y].set_ylabel(
                    str(cube.var_name) + " / " + str(cube.units))
                if j > 5:
                    ax3[x, y].set_xlabel('Time')
            if i == 3:
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
    fig3.savefig(plot_path + 'Standard Timeseries')

    fig4.tight_layout()
    fig4.savefig(plot_path + 'Patterns Timeseries')


def plot_scores(list_cubelists, plot_path):
    # plot color mesh of individual months
    cube_list = list_cubelists[4]
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


def plot_ebm_timeseries(list_cubes, plot_path, ocean_frac, land_frac):
    """Plots timeseries of aggregated cubes, across all scenarios."""
    fig, ax = plt.subplots(2, 3, figsize=(9, 8), sharey=True, sharex=True)
    for i, cube in enumerate(list_cubes):
        yrs = (1850 + np.arange(cube.shape[0])).astype('float')
        if i == 0:
            avg_cube = sf.area_avg(cube, return_cube=True)
            avg_cube.data -= np.mean(avg_cube.data[0:40])
            ax[0, 0].plot(yrs, avg_cube.data)
            ax[0, 0].set_title('Global')
            ax[0, 0].set_ylabel('Net TOA Radiation (Wm-2)')

            avg_cube2 = sf.area_avg_landsea(cube,
                                            ocean_frac,
                                            land_frac,
                                            land=True,
                                            return_cube=True)
            avg_cube2.data -= np.mean(avg_cube2.data[0:40])
            ax[0, 1].plot(yrs, avg_cube2.data)
            ax[0, 1].set_title('Land-Weighted')

            avg_cube3 = sf.area_avg_landsea(cube,
                                            ocean_frac,
                                            land_frac,
                                            land=False,
                                            return_cube=True)
            avg_cube3.data -= np.mean(avg_cube3.data[0:40])
            ax[0, 2].plot(yrs, avg_cube3.data)
            ax[0, 2].set_title('Ocean-Weighted')

        if i == 1:
            avg_cube = sf.area_avg(cube, return_cube=True)
            avg_cube.data -= np.mean(avg_cube.data[0:40])
            ax[1, 0].plot(yrs, avg_cube.data)
            ax[1, 0].set_ylabel('Air Temperature $(K)$')
            ax[1, 0].set_xlabel('Time')

            avg_cube2 = sf.area_avg_landsea(cube,
                                            ocean_frac,
                                            land_frac,
                                            land=True,
                                            return_cube=True)
            avg_cube2.data -= np.mean(avg_cube2.data[0:40])
            ax[1, 1].plot(yrs, avg_cube2.data)
            ax[1, 1].set_xlabel('Time')

            avg_cube3 = sf.area_avg_landsea(cube,
                                            ocean_frac,
                                            land_frac,
                                            land=False,
                                            return_cube=True)
            avg_cube3.data -= np.mean(avg_cube3.data[0:40])
            ax[1, 2].plot(yrs, avg_cube3.data)
            ax[1, 2].set_xlabel('Time')

            fig.savefig(plot_path + "anomaly_plots")
            plt.close(fig)

    plt.figure(figsize=(8, 4))
    for i, cube in enumerate(list_cubes):
        if i == 0:
            plt.subplot(121)
            qplt.pcolormesh(cube[0])
            plt.gca().coastlines()

        if i == 1:
            plt.subplot(122)
            qplt.pcolormesh(cube[0])
            plt.gca().coastlines()

    plt.tight_layout()
    plt.savefig(plot_path + cube.var_name + "_mesh_global_jan")
    plt.close()

    return fig


def plot_ebm_prediction(temp_global, tas_delta, plot_path):
    """Plots the EBM's tas prediction vs actual model output."""
    yrs = (1850 + np.arange(temp_global.shape[0])).astype('float')

    plt.plot(yrs,
             temp_global,
             color='black',
             zorder=10,
             linewidth=1.5,
             label='EBM Prediction')
    plt.plot(yrs, tas_delta, color='red', label='Model')
    plt.legend(loc='upper left')
    plt.xlabel("Time")
    plt.ylabel("Air Surface Temperature (K)")

    plt.savefig(plot_path + "_tas_check")
    plt.close()


def forcing_plot(reg, avg_list, yrs, forcing, plot_path):
    """Plots the regression between tas and rtmt, as well as the resultant
    forcing timeseries."""
    # regression line
    x_reg = np.linspace(-3.0, 4.0, 2)
    y_reg = reg.slope * x_reg + reg.intercept

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))

    # plotting regression line
    ax[0].scatter(avg_list[2].data, avg_list[3].data, c='b', label='tas, rtmt')
    ax[0].plot(x_reg, y_reg, c='r', label='regression')
    ax[0].legend(loc='upper left')

    # plotting forcing timeseries
    ax[1].plot(yrs, forcing)
    ax[1].set_xlabel("Time")
    ax[1].set_ylabel("Radiative forcing (Wm-2)")

    fig.savefig(plot_path + '_regression+forcing')
    plt.close(fig)
