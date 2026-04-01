"""Diagnostic script to plot ENSO metrics duration, diversity dive downs."""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    extract_month,
    mask_above_threshold,
    mask_below_threshold,
    rolling_window_statistics,
)

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def mask_to_years(events):
    """Build list of years with mask."""
    maskedtime = np.ma.masked_array(
        events.coord("time").points,
        mask=events.data.mask,
    )
    return [
        events.coord("time").units.num2date(time).year
        for time in maskedtime.compressed()
    ]


def enso_events_lc(cube, metric):
    """Get ENSO events for lifecycle or duration."""
    datayears = [
        cube.coord("time").units.num2date(time).year
        for time in cube.coord("time").points
    ]
    leadlagyrs = datayears[:3] + datayears[-3:]

    cb_std = cube.data.std()
    a_events = mask_to_years(mask_above_threshold(cube.copy(), -0.5 * cb_std))
    o_events = mask_to_years(mask_below_threshold(cube.copy(), 0.5 * cb_std))
    events = {"la nina": a_events, "el nino": o_events}
    if metric == "14duration":  # check years not in first/last 3
        for key, yrls in events.items():
            events[key] = [yr for yr in yrls if yr not in leadlagyrs]

    return events


def enso_composite(n34):
    """Get ENSO composite for duration."""
    n34_dec = extract_month(n34, 12)
    events = enso_events_lc(
        n34_dec,
        "14duration",
    )  # check years not in first/last 3

    enso_res = {}
    for enso, years in events.items():
        years_of_interest = []
        for yr in years:
            years_of_interest.append(
                [yr - 2, yr - 1, yr, yr + 1, yr + 2, yr + 3],
            )

        cube_data = {}
        for enso_epoch in years_of_interest:
            year_enso = iris.Constraint(
                time=lambda cell, enso_epoch=enso_epoch: cell.point.year
                in enso_epoch,
            )
            cube_2 = n34.extract(year_enso)  # extract rolling 6
            yr = enso_epoch[2]
            cube_data[yr] = cube_2.data

        durations = [
            threshold_duration(line, 0.5, enso)
            for yr, line in cube_data.items()
        ]
        enso_res[enso] = durations  # ls of durations for each event

    return enso_res


def threshold_duration(line, value, enso):
    """Count duration in months for each dataset and enso composite.

    Creates a boolean list indicating the months above/below the threshold.
    Then iterates through the list to count the maximum number of consecutive
    months that meet the condition to give duration.

    Args:
    ----
    line (np.ndarray) : ENSO lifecycle data line
    value (float) : Threshold value
    enso (str) : ENSO phase ("el nino" or "la nina")
    """
    cnt_month = line > value if enso == "el nino" else line < -value

    cnt = 0
    durations = []
    # count number of consecutive true months
    for a in cnt_month:
        if a:
            cnt += 1
        else:
            if cnt != 0:
                durations.append(cnt)
            cnt = 0
    if not durations:
        durations.append(cnt)  # if no events, append 0
    return max(durations)


def duration_composite_plot(obs_model, dt_ls):
    """Plot ENSO duration for composites."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    symbol = {"el nino": "> ", "la nina": "< -"}
    for ninanino, ax in zip(["la nina", "el nino"], [ax1, ax2], strict=False):
        bplt = ax.boxplot(
            [obs_model[0][ninanino], obs_model[1][ninanino]],
            tick_labels=dt_ls,
            showmeans=True,
        )

        ax.set_title(f"{ninanino.title()} duration")
        ax.set_ylabel(f"SSTA{symbol[ninanino]}0.5 (months)")
        ax.grid(linestyle="--", axis="y")
        colour_boxplots(bplt)

    return fig


def colour_boxplots(bplt):
    """Colour boxplots for model and observations."""
    colour = ["black", "tab:blue"]
    for k in bplt.keys():  # colour separately
        logger.info("%s: %s", k, bplt[k])
        for j, line in enumerate(bplt[k]):
            if len(bplt[k]) > 2:  # caps,whiskers first 2
                n = 0 if j < len(bplt[k]) / 2 else 1
                line.set_color(colour[n])
            else:
                line.set_color(colour[j])
            line.set_linewidth(2)
            if k == "means":
                line.set_markerfacecolor(colour[j])
                line.set_markeredgecolor(colour[j])


def diversity_plots3(div_data, dt_ls):
    """Plot ENSO diversity for ENSO and composites."""
    fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(16, 10), sharey=True)
    symbol = {"el nino": "max", "la nina": "min", "enso": "max/min"}
    for ninanino, ax in zip(
        ["la nina", "el nino", "enso"],
        [ax1[0], ax1[1], ax2[0]],
        strict=False,
    ):
        bplt = ax.boxplot(
            [div_data[0][ninanino], div_data[1][ninanino]],
            tick_labels=dt_ls,
            showmeans=True,
        )

        ax.yaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
        ax.set_ylabel(f"longitude of {symbol[ninanino]} SSTA")
        ax.set_title(f"{ninanino.title()} Diversity")
        ax.grid(linestyle="--", axis="y")
        colour_boxplots(bplt)

    ax2[0].set_title("ENSO Diversity")
    ax2[1].remove()

    return fig


def diversity(ssta_cube, events_dict):
    """Get ENSO diversity for longitude of max/min SSTA."""
    ssta_cube = extract_month(ssta_cube, 12)  # diversity preprocessing
    ssta_cube = rolling_window_statistics(
        ssta_cube,
        coordinate="longitude",
        operator="mean",
        window_length=5,
    )

    res_lon = {}
    for (
        enso,
        events,
    ) in events_dict.items():  # each enso year, max/min SSTA, get longitude
        year_enso = iris.Constraint(
            time=lambda cell, events=events: cell.point.year in events,
        )
        cube = ssta_cube.extract(year_enso)
        if enso == "nina":
            cube = cube * -1
        # iterate through cube, each time get max/min value and return lon
        loc_ls = []
        for yr_slice in cube.slices(["longitude"]):
            indx = np.argmax(yr_slice.data)  # if nina multiply by -1 or min
            loc_ls.append(cube.coord("longitude").points[indx])

        res_lon[enso] = loc_ls
    return res_lon


def compute_enso_metrics(input_pair, dt_ls, var_group, metric):
    """Compute ENSO metrics for duration and diversity."""
    if metric == "14duration":
        # level 2, input_pair: obs first
        mod = enso_composite(input_pair[1][var_group[0]])
        obs = enso_composite(input_pair[0][var_group[0]])

        fig = duration_composite_plot([obs, mod], dt_ls)
    elif metric == "15diversity":
        data_box = []
        for ds in input_pair:  # obs first
            events = enso_events_lc(ds[var_group[0]], metric)
            results_lon = diversity(ds[var_group[1]], events)
            results_lon["enso"] = (
                results_lon["el nino"] + results_lon["la nina"]
            )
            data_box.append(results_lon)

        fig = diversity_plots3(data_box, dt_ls)
    else:
        fig = None
    return fig


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
        "14duration": ["tos_lifdur1"],
        "15diversity": ["tos_patdiv1", "tos_lifdurdiv2"],
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

        msg = (
            f"{metric} : observation datasets {len(obs)}, models {len(models)}"
        )
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

            # compute metric, get figure
            figs = compute_enso_metrics(
                input_pair,
                [obs[0]["dataset"], dataset],
                var_preproc,
                metric,
            )

            dt_files = obs_files + [ds["filename"] for ds in models]
            prov_record = get_provenance_record(
                f"ENSO metrics {metric} dive down",
                dt_files,
            )

            save_figure(
                f"{dataset}_{metric}_divedown",
                prov_record,
                cfg,
                figure=figs,
                dpi=300,
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
