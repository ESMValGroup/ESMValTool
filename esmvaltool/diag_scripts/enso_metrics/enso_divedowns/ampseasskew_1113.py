"""Script for enso metrics amplitude, seasonality, skewness dive downs."""

import logging
import os
from pprint import pformat

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    climate_statistics,
    extract_region,
    extract_season,
    meridional_statistics,
)
from matplotlib.lines import Line2D
from scipy.stats import skew

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def seasonality_level2(mod_obs_cubes, dt_ls):
    """Plot the for monthly standard deviation as line."""
    fig = plt.figure(figsize=(6, 6), dpi=150)
    pltcolors = ["black", "tab:blue"]
    for i, cube in enumerate(mod_obs_cubes):
        qplt.plot(cube, label=dt_ls[i], color=pltcolors[i], linewidth=4)

    months = ["Jan", "May", "Sep"]
    plt.xticks(range(1, 13, 4), labels=months)
    plt.xlim(1, 12)
    # Set the x and y axis labels
    plt.xlabel("Months")
    plt.ylabel("SSTA std (°C)")
    plt.title("SSTA standard deviation")
    plt.grid(linestyle="--")
    plt.legend()
    return fig


def seasonality_level3(mod_obs_cubes, dt_ls):
    """Plot seasonality along longitude as filled contour."""
    fig = plt.figure(figsize=(16, 8), dpi=300)
    i = 121
    # Define month labels for the y-axis
    months = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]

    # Plot data for each label and cube in the dictionary
    for label, cube in dict(zip(dt_ls, mod_obs_cubes, strict=False)).items():
        ax = plt.subplot(i)  # Add projection for geographic plots
        c1 = iplt.contourf(
            cube,
            cmap="Reds",
            levels=np.arange(0.0, 2, 0.1),
            extend="both",
        )  # Set levels to include <0 and >2
        ax.set_yticks(range(1, 13))
        ax.set_yticklabels(months)
        ax.set_title(label)
        ax.set_ylabel("Months")
        ax.set_xlabel("Longitude")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
        i += 1

    fig.subplots_adjust(bottom=0.15)  # Adjust bottom margin to fit colorbar
    cax = plt.axes([0.15, -0.08, 0.7, 0.05])  # Position for colorbar
    cbar = fig.colorbar(
        c1,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(0, 2.2, 0.5),
    )
    cbar.set_label("SSTA std (°C)")
    return fig


def seasonality_level4(mod_obs_cubes, dt_ls):
    """Plot seasonality along longitude as line for NDJ and MAM."""
    fig = plt.figure(figsize=(10, 8), dpi=300)
    lines = {f"{dt_ls[1]}": "solid", f"{dt_ls[0]}": "dashdot"}
    for label, cube in dict(zip(dt_ls, mod_obs_cubes, strict=False)).items():
        for seas, col in {"NDJ": "red", "MAM": "blue"}.items():
            cube_plot = extract_season(
                cube,
                seas,
            )  # get NDJ #{ NDJ, MAM, red, blue
            cube_plot = climate_statistics(
                cube_plot,
                operator="std_dev",
                period="full",
            )
            qplt.plot(
                cube_plot,
                color=col,
                label=f"{label} {seas}",
                linestyle=lines[label],
                linewidth=2,
            )

    plt.xlabel("longitude")
    plt.ylabel("SSTA std (°C)")
    plt.title("SSTA standard deviation")
    plt.xlim(150, 270)
    # plt.legend()
    create_legend_seas(dt_ls)
    plt.grid(linestyle="--")
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
    return fig


def create_legend_seas(dt_ls):
    """Create a legend for the season composite plot."""
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="blue",
            markersize=8,
            label="MAM",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=8,
            label="NDJ",
        ),
        Line2D([0], [0], linestyle="solid", color="k", lw=1.5, label=dt_ls[1]),
        Line2D(
            [0],
            [0],
            linestyle="dashdot",
            color="k",
            lw=1.5,
            label=f"Ref: {dt_ls[0]}",
        ),
    ]
    plt.legend(handles=legend_elements)


def amplitude_level2(mod_obs_cubes, dt_ls):
    """Plot SSTA std along longitude as line."""
    fig = plt.figure(figsize=(6, 6), dpi=300)
    pltcolors = ["black", "tab:blue"]

    for i, cube in enumerate(mod_obs_cubes):
        cube = climate_statistics(cube, operator="std_dev", period="full")
        qplt.plot(cube, label=dt_ls[i], color=pltcolors[i], linewidth=4)

    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))

    # Adding labels and title
    plt.xlabel("Longitude")
    plt.ylabel("SSTA std (°C)")
    plt.title("SSTA standard deviation")
    plt.grid(linestyle="--")
    plt.legend()
    return fig


def map_plots(mod_obs_cubes, dt_ls, i, valrange, cmap):
    """Plot maps for model and observation."""
    fig = plt.figure(figsize=(20, 7))
    if i < 200:
        extent = [120, 290, -20, 20]
        proj = ccrs.Orthographic(central_longitude=210.0)
        cbarext = [0.15, 0.08, 0.7, 0.05]
    else:
        extent = [120, 290, -15, 15]
        proj = ccrs.PlateCarree(central_longitude=180)
        cbarext = [0.15, 0.05, 0.7, 0.03]
    # i = 121
    for label, cube in dict(zip(dt_ls, mod_obs_cubes, strict=False)).items():
        ax1 = plt.subplot(i, projection=proj)
        ax1.add_feature(
            cfeature.LAND,
            facecolor="gray",
        )  # Add land feature with gray color
        ax1.coastlines()
        cf1 = iplt.contourf(
            cube,
            levels=np.arange(valrange[0], valrange[1], 0.1),
            extend="both",
            cmap=cmap,
        )
        # if i < 200:  # 121 or 122
        ax1.set_extent(extent, crs=ccrs.PlateCarree())
        ax1.set_title(label)

        # Add gridlines for latitude and longitude
        gl1 = ax1.gridlines(draw_labels=True, linestyle="--")
        gl1.top_labels = False
        gl1.right_labels = False
        i += 1

    # Add a single colorbar at the bottom
    cax = plt.axes(cbarext)
    cbar = fig.colorbar(
        cf1,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(valrange[0], valrange[1] + 0.2, 0.5),
    )
    cbar.set_label("SSTA std (°C)")
    return fig


def asym_level2(mod_obs_ls, dt_ls):
    """Plot SSTA skewness along longitude as line."""
    fig = plt.figure(figsize=(7, 6), dpi=150)
    pltcolors = ["black", "tab:blue"]

    for i, tuplepoints in enumerate(mod_obs_ls):
        plt.plot(*tuplepoints, label=dt_ls[i], color=pltcolors[i], linewidth=4)
    plt.axhline(y=0, color="black", linewidth=2)
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
    plt.xlim(150, 275)  # Set x-axis limits to match the region of interest
    # Adding labels and title
    plt.title("SSTA skewness")
    plt.xlabel("Longitude")
    plt.ylabel("SSTA skew (°C)")
    plt.grid(linestyle="--")
    plt.legend()
    return fig


def format_longitude(x, _pos):
    """Format longitude values for plotting."""
    if x > 180:
        return f"{int(360 - x)}°W"
    if x == 180:
        return f"{int(x)}°"
    return f"{int(x)}°E"


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute ENSO metrics and return figures."""
    model_obs = [input_pair[0][var_group[0]], input_pair[1][var_group[0]]]
    # input_pair: obs first
    if metric == "11amplitude":
        fig2 = amplitude_level2(model_obs, dt_ls)
        fig3 = map_plots(
            [input_pair[0][var_group[1]], input_pair[1][var_group[1]]],
            dt_ls,
            i=121,
            valrange=[0, 2],
            cmap="Reds",
        )
        return fig2, fig3

    if metric == "12seasonality":
        fig2 = seasonality_level2(model_obs, dt_ls)
        fig3 = seasonality_level3(
            [input_pair[0][var_group[1]], input_pair[1][var_group[1]]],
            dt_ls,
        )
        fig4 = seasonality_level4(
            [input_pair[0][var_group[2]], input_pair[1][var_group[2]]],
            dt_ls,
        )
        process = {}
        for seas in ["NDJ", "MAM"]:
            for label, cube in dict(
                zip(
                    dt_ls,
                    [input_pair[0][var_group[3]], input_pair[1][var_group[3]]],
                    strict=False,
                ),
            ).items():
                cube_plot = extract_season(cube, seas)
                cube_plot = climate_statistics(
                    cube_plot,
                    operator="std_dev",
                    period="full",
                )
                process[f"{label} {seas}"] = cube_plot
        fig5 = map_plots(
            list(process.values()),
            list(process.keys()),
            i=221,
            valrange=[0, 2],
            cmap="Reds",
        )
        return fig2, fig3, fig4, fig5

    if metric == "13asymmetry":
        asymls, asym_cubes = [], []
        region = {
            "start_longitude": 150.0,
            "end_longitude": 275.0,
            "start_latitude": -5.0,
            "end_latitude": 5.0,
        }
        for cube in model_obs:
            # lvl2
            cube2 = extract_region(cube, **region)
            cube2 = meridional_statistics(cube2, "mean")
            asymls.append(
                (cube2.coord("longitude").points, skew(cube2.data, axis=0)),
            )
            # lvl3
            skew_array = skew(cube.data, axis=0)
            new_cube = cube.collapsed("time", iris.analysis.MAX)
            new_cube.data = skew_array
            asym_cubes.append(new_cube)

        fig2 = asym_level2(asymls, dt_ls)
        fig3 = map_plots(
            asym_cubes,
            dt_ls,
            i=121,
            valrange=[-1.5, 1.5],
            cmap="RdBu_r",
        )
        return fig2, fig3


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "statistics": ["anomaly"],
        "domains": ["eq"],
        "plot_types": ["line"],
        "authors": [
            "chun_felicity",
            "sullivan_arnold",
            "beucher_romain",
        ],
        "references": [
            "access-nri",
        ],
        "ancestors": ancestor_files,
    }
    return record


def main(cfg):
    """Run ENSO metrics."""
    input_data = cfg["input_data"].values()
    mpl.rc("font", size=12)
    # iterate through each metric and get variable group, select_metadata, map to function call
    metrics = {
        "11amplitude": ["tos_seas_split", "tos_amp"],
        "12seasonality": [
            "tos_seas_area",
            "tos_seas_meri",
            "tos_seas_split",
            "tos_asym",
        ],
        "13asymmetry": ["tos_asym"],
    }

    # select twice with project to get obs, iterate through model selection
    for metric, var_preproc in metrics.items():
        logger.info("%s,%s", metric, var_preproc)
        obs, models = [], []
        for var_prep in var_preproc:
            obs += select_metadata(
                input_data,
                variable_group=var_prep,
                project="OBS",
            )
            obs += select_metadata(
                input_data,
                variable_group=var_prep,
                project="OBS6",
            )
            models += select_metadata(
                input_data,
                variable_group=var_prep,
                project="CMIP6",
            )

        # log
        msg = f"{metric} : observation datasets {len(obs)}, models {pformat(models)}"
        logger.info(msg)

        # list dt_files
        obs_files = [ds["filename"] for ds in obs]  # and models separate

        # obs datasets for each model
        obs_datasets = {
            dataset["variable_group"]: iris.load_cube(dataset["filename"])
            for dataset in obs
        }

        # group models by dataset
        model_ds = group_metadata(models, "dataset", sort="project")

        for dataset, mod_ds in model_ds.items():
            logger.info(
                "%s, preprocessed cubes:%d, dataset:%s",
                metric,
                len(mod_ds),
                dataset,
            )

            model_datasets = {
                attributes["variable_group"]: iris.load_cube(
                    attributes["filename"],
                )
                for attributes in mod_ds
            }
            input_pair = [obs_datasets, model_datasets]
            logger.info(pformat(model_datasets))

            # compute metric, get figure
            figs = compute_enso_metrics(
                input_pair,
                [obs[0]["dataset"], dataset],
                var_preproc,
                metric,
            )

            dt_files = obs_files + [ds["filename"] for ds in models]
            for i, fig in enumerate(figs):
                prov_record = get_provenance_record(
                    f"ENSO metrics {metric} level {i + 2}",
                    dt_files,
                )

                save_figure(
                    f"{dataset}_{metric}_level_{i + 2}",
                    prov_record,
                    cfg,
                    figure=fig,
                    dpi=300,
                )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
