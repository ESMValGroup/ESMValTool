"""diagnostic script to plot ENSO lifecycle dive down."""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    extract_month,
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


def dec_and_years(n34_cube):
    """Extract December data and years."""
    n34_dec = extract_month(n34_cube, 12)
    leadlagyr = 3
    n34_dec_years = [
        n34_dec.coord("time").units.num2date(time).year
        for time in n34_dec.coord("time").points
    ]
    event_years = n34_dec_years[leadlagyr:-leadlagyr]
    # Ensure that the selected years are not the last years
    return n34_dec, event_years


def sst_regressed_2d(event_years, n34_area, n34_dec):
    """Regression function for sst_time_series area on sst_enso."""
    n34_area_selected = []
    for yr in event_years:
        enso_epoch = [yr - 2, yr - 1, yr, yr + 1, yr + 2, yr + 3]
        year_enso = iris.Constraint(
            time=lambda cell, enso_epoch=enso_epoch: cell.point.year
            in enso_epoch,
        )

        n34_area_selected.append(n34_area.extract(year_enso).data)

    event_constr = iris.Constraint(
        time=lambda cell: cell.point.year in event_years,
    )
    n34_dec_ct = n34_dec.extract(event_constr)

    # 2 area linear regression
    b_data = n34_dec_ct.data
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T
    a_arr = np.array(n34_area_selected)
    a_data = a_arr.reshape(a_arr.shape[0], -1)

    coefs_area, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)
    slope_area = coefs_area[0].reshape(a_arr.shape[1], a_arr.shape[2])

    return slope_area


def plot_level2_contour(model_ds, obs_ds, dt_ls):
    """Plot the level 2 contour for ENSO SSTA regression."""
    # Create the subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 10), sharey=True)

    # First subplot: HadISST
    c1 = plot_ssta(
        ax1,
        obs_ds[0],
        obs_ds[1],
        f"{dt_ls[1]}:\nreg(ENSO SSTA, SSTA)",
        i=None,
    )

    # Second subplot: Model
    plot_ssta(
        ax2,
        model_ds[0],
        model_ds[1],
        f"{dt_ls[0]}:\nreg(ENSO SSTA, SSTA)",
        i=None,
    )

    # Adjust the layout to add space for the colorbar
    fig.subplots_adjust(bottom=0.15)

    # Add a horizontal colorbar beneath the subplots
    ticks = np.arange(-1.2, 1.3, 0.2)
    fig.colorbar(
        c1,
        ax=[ax1, ax2],
        label="Regression (°C / °C)",
        orientation="horizontal",
        ticks=ticks,
        fraction=0.05,
        pad=0.08,
    )

    return fig


def plot_ssta(ax, area_coordlon, ssta, label, i=None):
    """Plot the SSTA for lifecycle along longitude on a contour plot."""
    if i is not None:
        axpos = ax[i]
    else:
        axpos = ax
    c1 = axpos.contourf(
        area_coordlon,
        range(1, 73),
        ssta,
        levels=np.linspace(-1.2, 1.2, 14),
        cmap="RdBu_r",
        extend="both",
    )
    if i == 0:
        axpos.set_ylabel("Months")
    else:
        axpos.set_ylabel("")
    axpos.set_xlim([160, 280])
    axpos.set_title(label)  #: reg SSTA
    axpos.set_yticks(range(1, 73, 6))
    axpos.set_yticklabels(["Jan", "Jul"] * (len(range(1, 73, 6)) // 2))
    axpos.set_xticks([180, 220, 260])
    if label.split()[0] == "El":
        axpos.xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
        axpos.set_xlabel("Longitude")
    else:
        axpos.xaxis.set_major_formatter(plt.NullFormatter())
    return c1


def enso_composite_plot(model_n34, line):
    """Plot the ENSO composite lifecycle lines."""
    n34_dec = extract_month(model_n34, 12)
    events = enso_events_lc(n34_dec)  # check years not in first/last 3
    colors = {"la nina": "blue", "el nino": "red"}
    months = np.arange(1, 73) - 36
    logger.info(events)
    for enso, years in events.items():
        cube_data = []
        for yr in years:
            enso_epoch = [yr - 2, yr - 1, yr, yr + 1, yr + 2, yr + 3]
            year_enso = iris.Constraint(
                time=lambda cell, enso_epoch=enso_epoch: cell.point.year
                in enso_epoch,
            )
            cube_2 = model_n34.extract(year_enso)  # extract rolling 6yr
            cube_data.append(cube_2.data.data)

        # No regression, mean enso epoch
        mean = np.mean(np.array(cube_data), axis=0)
        plt.plot(months, mean, linestyle=line, color=colors[enso], lw=3)


def sst_2d(event_years, n34_area):
    """Extract SST area data for the event years."""
    n34_area_selected = []
    for yr in event_years:
        enso_epoch = [yr - 2, yr - 1, yr, yr + 1, yr + 2, yr + 3]
        year_enso = iris.Constraint(
            time=lambda cell, enso_epoch=enso_epoch: cell.point.year
            in enso_epoch,
        )
        n34_area_selected.append(n34_area.extract(year_enso).data)

    arr = np.array(n34_area_selected)
    a_data = arr.reshape(arr.shape[0], -1)
    means = np.mean(a_data, axis=0)
    return means.reshape(arr.shape[1], arr.shape[2])


def enso_composite_plot4(n34_cube, n34_area, label, i, axes):
    """Plot the ENSO composite lifecycle dive down 4."""
    n34_dec = extract_month(n34_cube, 12)
    events = enso_events_lc(n34_dec)

    for enso, years in events.items():
        # plot area
        means = sst_2d(years, n34_area)
        cplot = plot_ssta(
            axes[enso],
            n34_area.coord("longitude").points,
            means,
            f"{enso.title()} {label}",
            i,
        )

    return cplot


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


def plot_l3_enso_lifecycle(model_ds, obs_ds, dt_ls):
    """Plot the ENSO lifecycle composite lines."""
    fig = plt.figure(figsize=(10, 6), dpi=300)
    enso_composite_plot(model_ds, "solid")
    enso_composite_plot(obs_ds, "dashdot")

    plt.axhline(y=0, color="black", linewidth=2)

    xticks = np.arange(1, 73, 6) - 36  # Adjust for lead/lag months
    xtick_labels = ["Jan", "Jul"] * (len(xticks) // 2)
    plt.xticks(xticks, xtick_labels)
    plt.xlim(-35, 36)
    plt.grid(linestyle="--")
    plt.yticks(np.arange(-2, 2.5, step=1))
    # plt.legend()
    create_legend(dt_ls)
    return fig


def plot4_enso_contour(model_ds, obs_ds, dt_ls):
    """Plot the ENSO lifecycle dive down 4 with subplots."""
    fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(8, 20), sharey=True)
    # # layout rows la nina ax1, el nino ax2 for datasets 0,1
    axes_rows = {"la nina": ax1, "el nino": ax2}

    c1 = enso_composite_plot4(model_ds[0], model_ds[1], dt_ls[0], 1, axes_rows)
    # obs
    c1 = enso_composite_plot4(obs_ds[0], obs_ds[1], dt_ls[1], 0, axes_rows)
    plt.subplots_adjust(hspace=0.08, wspace=0.08)
    fig.colorbar(
        c1,
        ax=[ax2[0], ax2[1]],
        label="ENSO SSTA °C",
        orientation="horizontal",
        ticks=np.arange(-1.2, 1.3, 0.4),
        fraction=0.05,
        pad=0.08,
    )

    return fig


def compute_enso_metrics(input_pair, dt_ls, var_group):
    """Compute the ENSO lifecycle dive downs."""
    mod_dec, mod_years = dec_and_years(input_pair[1][var_group[0]])
    model_area = sst_regressed_2d(
        mod_years,
        input_pair[1][var_group[1]],
        mod_dec,
    )
    # level 2, input_pair: obs first
    obs_dec, obs_years = dec_and_years(input_pair[0][var_group[0]])
    obs_area = sst_regressed_2d(
        obs_years,
        input_pair[0][var_group[1]],
        obs_dec,
    )

    # plot function #need xticks, labels as dict/ls
    modds = (input_pair[1][var_group[1]].coord("longitude").points, model_area)
    obsds = (input_pair[0][var_group[1]].coord("longitude").points, obs_area)
    fig2 = plot_level2_contour(modds, obsds, dt_ls)

    fig3 = plot_l3_enso_lifecycle(
        input_pair[1][var_group[0]],
        input_pair[0][var_group[0]],
        dt_ls,
    )

    fig4 = plot4_enso_contour(
        list(input_pair[1].values()),
        list(input_pair[0].values()),
        dt_ls,
    )
    return fig2, fig3, fig4


def mask_to_years(events):
    """Convert masked array to years."""
    maskedtime = np.ma.masked_array(
        events.coord("time").points,
        mask=events.data.mask,
    )
    # return years
    return [
        events.coord("time").units.num2date(time).year
        for time in maskedtime.compressed()
    ]


def enso_events_lc(cube):
    """Get ENSO events from the cube."""
    datayears = [
        cube.coord("time").units.num2date(time).year
        for time in cube.coord("time").points
    ]
    leadlagyrs = datayears[:3] + datayears[-3:]
    # get cube years min/max, remove 3:-3
    cb_std = cube.data.std()
    a_events = mask_to_years(mask_above_threshold(cube.copy(), -0.5 * cb_std))
    o_events = mask_to_years(mask_below_threshold(cube.copy(), 0.5 * cb_std))
    events = {"la nina": a_events, "el nino": o_events}
    for key, yrls in events.items():
        events[key] = [yr for yr in yrls if yr not in leadlagyrs]

    return events


def format_longitude(x, _pos):
    """Format longitude values for plotting."""
    if x > 180:
        return f"{int(360 - x)}°W"
    if x == 180:
        return f"{int(x)}°"
    return f"{int(x)}°E"


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

    # iterate through each metric and get variable group
    metrics = {
        "10lifecycle": ["tos_lifdur1", "tos_lifdurdiv2"],
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
        msg = (
            f"{metric} : observation datasets {len(obs)}, models {len(models)}"
        )
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

            # compute metric, get figure
            figs = compute_enso_metrics(
                input_pair,
                [dataset, obs[0]["dataset"]],
                var_preproc,
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
