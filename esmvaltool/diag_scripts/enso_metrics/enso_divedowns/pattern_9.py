"""diagnostic script to plot ENSO pattern dive downs."""

import logging
import os
from pprint import pformat

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    climate_statistics,
    mask_above_threshold,
    mask_below_threshold,
)
from matplotlib.lines import Line2D

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_ensofig(i, cube, label):
    """Subplot map with the given cube and label."""
    ax1 = plt.subplot(i, projection=ccrs.Orthographic(central_longitude=210.0))
    ax1.add_feature(
        cfeature.LAND,
        facecolor="gray",
    )  # Add land feature with gray color
    ax1.coastlines()
    cf2 = iplt.contourf(
        cube,
        levels=np.arange(-2, 2, 0.1),
        extend="both",
        cmap="RdBu_r",
    )

    ax1.set_extent([130, 290, -20, 20], crs=ccrs.PlateCarree())
    ax1.set_title(label)
    # Add gridlines for latitude and longitude
    gl1 = ax1.gridlines(draw_labels=True, linestyle="--")
    gl1.top_labels = False
    gl1.right_labels = False
    return i + 2, cf2


def plot_maps_pattern(processed_data):
    """Plot maps for the processed model and obs data."""
    fig = plt.figure(figsize=(20, 7))
    i = 121
    for label, cube in processed_data.items():
        cf1 = plot_ensofig(i, cube, label)[1]
        i += 1
    # Add a single colorbar at the bottom
    cax = plt.axes([0.15, 0.08, 0.7, 0.05])
    cbar = fig.colorbar(
        cf1,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(-2, 2.2, 0.5),
    )
    cbar.set_label("regression(ENSO SSTA, SSTA) (°C/°C)")
    return fig


def enso_regression(prep_datasets, i, label):
    """Plot composite maps."""
    events = enso_events(prep_datasets[1])
    for enso, years in events.items():
        year_enso = iris.Constraint(
            time=lambda cell, years=years: cell.point.year in years,
        )
        cube_2 = prep_datasets[0].extract(year_enso)
        cube = climate_statistics(cube_2, operator="mean")
        if enso == "nina":  # plot separate
            i, cf = plot_ensofig(i, cube, f"La Nina {label}")
        else:
            i, cf = plot_ensofig(i, cube, f"El Nino {label}")
    return cf


def plot_maps_pattern4(model_ds, obs_ds, dt_ls):
    """Plot maps for the model and obs data."""
    fig = plt.figure(figsize=(20, 10))
    i = 221  # nina top row
    enso_regression(obs_ds, i, dt_ls[1])
    i = 222
    cf = enso_regression(model_ds, i, dt_ls[0])

    # Add a single colorbar at the bottom
    cax = plt.axes([0.15, 0.08, 0.7, 0.05])
    cbar = fig.colorbar(
        cf,
        cax=cax,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(-2, 2.2, 0.5),
    )
    cbar.set_label("SSTA (°C)")
    return fig


def lin_regress_matrix(cubea, cubeb):
    """
    Calculate the linear regression of cubea on cubeb using matrix operations.

    Parameters
    ----------
    cubea: iris.cube.Cube
        The 2D input cube for which the regression is calculated.

    cubeb: iris.cube.Cube
        The cube used as the independent variable in the regression.

    Returns
    -------
    iris.cube.Cube
        A new cube containing the slope of the regression for each spatial point.
    """
    # Get data as flattened arrays
    a_data = cubea.data.reshape(
        cubea.shape[0],
        -1,
    )  # Shape (time, spatial_points)
    b_data = cubeb.data.flatten()  # Shape (time,)

    # Add intercept term by stacking a column of ones with cubeb
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T

    # Solve the linear equations using least squares method
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)

    # Extract slopes from coefficients #coefs 1
    slopes = coefs[0].reshape(cubea.shape[1], cubea.shape[2])

    # Create a new Iris Cube for the regression results
    result_cube = iris.cube.Cube(
        slopes,
        long_name="regression ENSO SSTA",
        dim_coords_and_dims=[
            (cubea.coord("latitude"), 0),
            (cubea.coord("longitude"), 1),
        ],
    )

    return result_cube


def mask_to_years(events):
    """Convert masked array of events to years."""
    maskedtime = np.ma.masked_array(
        events.coord("time").points,
        mask=events.data.mask,
    )
    # return years
    return [
        events.coord("time").units.num2date(time).year
        for time in maskedtime.compressed()
    ]


def enso_events(cube):
    """Identify ENSO events."""
    a_events = mask_to_years(mask_above_threshold(cube.copy(), -0.75))
    o_events = mask_to_years(mask_below_threshold(cube.copy(), 0.75))
    return {"nina": a_events, "nino": o_events}


def plot_enso_ssta(prep_datasets, line, label):
    """Plot composite SSTA patterns."""
    events = enso_events(prep_datasets[1])
    for enso, years in events.items():
        year_enso = iris.Constraint(
            time=lambda cell, years=years: cell.point.year in years,
        )
        cube_2 = prep_datasets[0].extract(year_enso)
        cube = climate_statistics(cube_2, operator="mean")

        if enso == "nina":  # plot separate
            qplt.plot(
                cube,
                color="blue",
                linestyle=line,
                lw=3,
                label=f"{label} La Nina",
            )
        else:
            qplt.plot(
                cube,
                color="red",
                linestyle=line,
                lw=3,
                label=f"{label} El Nino",
            )


def create_legend(dt_ls):
    """Create a legend for the ENSO composite plot."""
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="blue",
            markersize=8,
            label="La Niña",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=8,
            label="El Niño",
        ),
        Line2D([0], [0], linestyle="solid", color="k", lw=1.5, label=dt_ls[0]),
        Line2D(
            [0],
            [0],
            linestyle="dashdot",
            color="k",
            lw=1.5,
            label=f"Ref: {dt_ls[1]}",
        ),
    ]
    plt.legend(handles=legend_elements)


def plot_lvl3_pattern(model_ds, obs_ds, labels):
    """Create figure for the composite SSTA pattern."""
    fig = plt.figure(figsize=(10, 8), dpi=300)

    plot_enso_ssta(model_ds, "solid", labels[0])
    plot_enso_ssta(obs_ds, "dashdot", labels[1])

    plt.title("ENSO's SSTA pattern")
    plt.ylabel("ENSO SSTA (°C)")
    plt.xlabel("longitude")
    plt.grid(linestyle="--")
    plt.axhline(y=0, color="black", linewidth=1.5)
    plt.xlim(150, 270)
    # plt.legend()
    create_legend(labels)
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
    return fig


def format_longitude(x, _pos):
    """Format longitude values for plotting."""
    if x > 180:
        return f"{int(360 - x)}°W"
    if x == 180:
        return f"{int(x)}°"
    return f"{int(x)}°E"


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute metrics and return figures."""
    model_ssta = input_pair[1][var_group[2]]  # input_pair: obs first
    model_nino34 = input_pair[1][var_group[0]]
    reg_mod = lin_regress_matrix(model_ssta, model_nino34)
    reg_obs = lin_regress_matrix(
        input_pair[0][var_group[2]],
        input_pair[0][var_group[0]],
    )
    processed = dict(zip(dt_ls, [reg_mod, reg_obs], strict=False))
    fig2 = plot_maps_pattern(processed)

    fig3 = plot_lvl3_pattern(
        [input_pair[1][var_group[1]], model_nino34],
        [input_pair[0][var_group[1]], input_pair[0][var_group[0]]],
        dt_ls,
    )
    fig4 = plot_maps_pattern4(
        [model_ssta, model_nino34],
        [input_pair[0][var_group[2]], input_pair[0][var_group[0]]],
        dt_ls,
    )

    return fig2, fig3, fig4


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "statistics": ["anomaly"],
        "domains": ["eq"],
        "plot_types": ["line"],
        "authors": [
            "chun_felicity",
            "beucher_romain",
            "sullivan_arnold",
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

    # iterate through each metric and get variable groups
    metrics = {"09pattern": ["tos_patdiv1", "tos_pat2", "tos_pat_map"]}

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

        msg = f"{metric} : observation datasets {len(obs)}, models {pformat(models)}"
        logger.info(msg)

        # list dt_files
        obs_files = [ds["filename"] for ds in obs]

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
                len(model_ds),
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
                [dataset, obs[0]["dataset"]],
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
