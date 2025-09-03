"""diagnostic script to plot ENSO nhf feedback metrics."""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_figure,
    select_metadata,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(obs_ds, model_ds, title, metric_varls):
    """Plot level 1 diagnostics for ENSO feedback metrics."""
    var_set = {"sst": 5, "nhf": 100, "ssh": 40}
    xvar = metric_varls[0].split("_")[0]
    figure = plt.figure(figsize=(10, 8), dpi=300)
    xseq = np.linspace(var_set[xvar] * -1, var_set[xvar], num=100)

    mod_slope, intcpt = linreg_1d(
        model_ds[metric_varls[1]],
        model_ds[metric_varls[0]],
    )
    plt.plot(xseq, intcpt + mod_slope * xseq)
    # time overlap
    obs1, obs2 = obs_extract_overlap(
        obs_ds[metric_varls[2]],
        obs_ds[metric_varls[0]],
    )
    obs_slope, intcpt = linreg_1d(obs1, obs2)
    plt.plot(xseq, intcpt + obs_slope * xseq, color="black")

    metric_val = abs((mod_slope - obs_slope) / obs_slope) * 100

    plt.scatter(
        model_ds[metric_varls[0]].data,
        model_ds[metric_varls[1]].data,
        s=10,
    )
    plt.scatter(obs2.data, obs1.data, s=20, c="black", marker="D")

    plt_settings(
        [mod_slope, obs_slope, metric_val],
        lvl1=True,
        metric_varls=metric_varls,
    )
    plt.title(title)

    return metric_val, figure


def plot_level2(obs_ds, model_ds, metric_varls, ds_labels):
    """Plot level 2 diagnostics for ENSO feedback metrics."""
    fig = plt.figure(figsize=(13, 5))
    plt.subplot(121)
    plt_lvl2_subplot(
        model_ds[metric_varls[0]],
        model_ds[metric_varls[1]],
        ds_labels[1],
        metric_varls,
    )

    plt.subplot(122)
    obs1, obs2 = obs_extract_overlap(
        obs_ds[metric_varls[0]],
        obs_ds[metric_varls[2]],
    )

    plt_lvl2_subplot(obs1, obs2, ds_labels[0], metric_varls)
    return fig


def plt_lvl2_subplot(ts_cube, tauu_cube, dataset_label, metric_varls):
    """Level 2 subplots for ENSO feedback metrics."""
    df = pd.DataFrame({"tos": ts_cube.data, "tauu": tauu_cube.data})
    slopes = []
    logger.info("%s, shape: %s", dataset_label, df.shape)
    var_set = {"sst": 5, "taux": 100, "nhf": 100, "ssh": 40}
    xvar = metric_varls[0].split("_")[0]

    plt.scatter(ts_cube.data, tauu_cube.data, c="k", s=10)
    xseq = np.linspace(var_set[xvar] * -1, var_set[xvar], num=50)
    slope, intcpt = linreg_1d(df["tauu"], df["tos"])
    plt.plot(xseq, intcpt + slope * xseq, c="black")
    slopes.append(slope)

    xseq = np.linspace(var_set[xvar] * -1, 0, num=50)
    slope, intcpt = linreg_1d(
        df.loc[df["tos"] < 0, "tauu"],
        df.loc[df["tos"] < 0, "tos"],
    )

    plt.plot(xseq, intcpt + slope * xseq, linewidth=3)
    slopes.append(slope)

    xseq = np.linspace(0, var_set[xvar], num=50)
    slope, intcpt = linreg_1d(
        df.loc[df["tos"] > 0, "tauu"],
        df.loc[df["tos"] > 0, "tos"],
    )
    plt.plot(xseq, intcpt + slope * xseq, color="red", linewidth=3)
    slopes.append(slope)

    plt.title(dataset_label)
    plt_settings(slopes, lvl1=False, metric_varls=metric_varls)


def linreg_1d(tauu, ts):  # cube or df
    """Perform linear regression on 1D data."""
    logger.info("linreg_1d shapes %s, %s", tauu.shape, ts.shape)
    if isinstance(tauu, iris.cube.Cube):
        b_data = ts.data
        a_data = tauu.data
    else:
        b_data = np.array(ts)
        a_data = np.array(tauu)
    b_with_intercept = np.vstack([b_data, np.ones_like(b_data)]).T
    coefs, _, _, _ = np.linalg.lstsq(b_with_intercept, a_data, rcond=None)
    slope, intercept = coefs[0], coefs[1]
    return slope, intercept


def plt_settings(slopes, lvl1, metric_varls):
    """Set plot settings for ENSO feedback metrics."""
    var_set = {
        "sst": (5, "°C"),
        "taux": (100, "1e-3 N/m2"),
        "nhf": (100, "W m-2"),
        "ssh": (40, "cm"),
        "east": "nino3",
        "west": "nino4",
    }
    # variable, limit, units, region
    xvar = metric_varls[0].split("_")[0]
    yvar = metric_varls[1].split("_")[0]
    xregion = var_set[metric_varls[0].split("_")[1]]  # east, west, eqp
    yregion = var_set[metric_varls[1].split("_")[1]]
    plt.xlim(var_set[xvar][0] * -1, var_set[xvar][0])
    plt.xticks(
        np.arange(
            var_set[xvar][0] * -1,
            var_set[xvar][0] + 1,
            var_set[xvar][0] / 2,
        ),
    )
    plt.ylim(var_set[yvar][0] * -1, var_set[yvar][0])
    plt.yticks(
        np.arange(
            var_set[yvar][0] * -1,
            var_set[yvar][0] + 1,
            var_set[yvar][0] / 2,
        ),
    )
    plt.grid(linestyle="--")
    plt.ylabel(
        f"{yregion} {yvar.upper()}A ({var_set[yvar][1]})",
    )  # parse labels/units
    plt.xlabel(f"{xregion} {xvar.upper()}A ({var_set[xvar][1]})")

    if lvl1:
        # for lvl1 #slopes is mod 0, obs 1, metric_val 2
        plt.text(
            0.05,
            0.95,
            f"model slope: {slopes[0]:.2f}",
            color="C0",
            transform=plt.gca().transAxes,
        )
        plt.text(
            0.05,
            0.9,
            f"ref slope: {slopes[1]:.2f}",
            color="black",
            transform=plt.gca().transAxes,
        )

        plt.text(
            0.99,
            0.03,
            f"metric: {slopes[2]:.2f}%",
            fontsize=12,
            ha="right",
            transform=plt.gca().transAxes,
            backgroundcolor="white",
        )

    else:
        plt.text(
            0.02,
            0.85,
            f"slope(all): {slopes[0]:.2f}\nslope(x<0): {slopes[1]:.2f}\nslope(x>0): {slopes[2]:.2f}",
            fontsize=12,
            ha="left",
            transform=plt.gca().transAxes,
            backgroundcolor="white",
        )


def obs_extract_overlap(obs_1, obs_2):
    """Extract overlapping time range from two observation datasets."""
    start_1 = obs_1.coord("time").cell(0).point
    end_1 = obs_1.coord("time").cell(-1).point
    start_2 = obs_2.coord("time").cell(0).point
    end_2 = obs_2.coord("time").cell(-1).point

    start_overlap = max(start_1, start_2)
    end_overlap = min(end_1, end_2)
    logger.info(
        "%s, %s obs time overlap: %s to %s",
        obs_1.standard_name,
        obs_2.standard_name,
        start_overlap,
        end_overlap,
    )
    obs1 = obs_1.extract(
        iris.Constraint(
            time=lambda t: start_overlap <= t.point <= end_overlap,
        ),
    )
    obs2 = obs_2.extract(
        iris.Constraint(
            time=lambda t: start_overlap <= t.point <= end_overlap,
        ),
    )

    return obs1, obs2


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
    """Run ENSO NHF feedback metrics."""
    input_data = cfg["input_data"].values()

    sst_nhf = ["sst_east", "nhf_east_mod", "nhf_east_obs"]
    metric = "sst_nhf"
    obs, models = [], []

    obs += select_metadata(input_data, project="OBS")
    obs += select_metadata(input_data, project="OBS6")
    models += select_metadata(input_data, project="CMIP6")

    # log
    msg = f"{metric} : observation datasets {len(obs)}, models {len(models)}"
    logger.info(msg)

    # group models by dataset
    model_ds = group_metadata(models, "dataset", sort="project")

    # dataset name
    for dataset, mod_ds in model_ds.items():
        logger.info(
            "%s, preprocessed cubes:%d, dataset:%s",
            metric,
            len(model_ds),
            dataset,
        )
        dt_files = [ds["filename"] for ds in obs] + [
            ds["filename"] for ds in mod_ds
        ]

        obs_ds = {
            dataset["variable_group"]: iris.load_cube(dataset["filename"])
            for dataset in obs
        }
        model = {
            attributes["variable_group"]: iris.load_cube(
                attributes["filename"],
            )
            for attributes in mod_ds
        }
        model["nhf_east_mod"] = -model["nhf_east_mod"]  # make negative

        title = "net heat flux feedback"
        value, fig = plot_level1(obs_ds, model, title, sst_nhf)

        metricfile = get_diagnostic_filename("matrix", cfg, extension="csv")
        with open(metricfile, "a+", encoding="utf-8") as f:
            f.write(f"{dataset},{metric},{value}\n")

        prov_record = get_provenance_record(
            f"ENSO metrics {metric} feedback",
            dt_files,
        )
        save_figure(
            f"{dataset}_{metric}",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
        )

        ds_labels = [
            f"{obs[0]['dataset']}_{obs[1]['dataset']}",
            mod_ds[0]["dataset"],
        ]
        fig = plot_level2(obs_ds, model, sst_nhf, ds_labels)
        save_figure(
            f"{dataset}_{metric}_lvl2",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
        )
        # plot level 3 separate diagnostic script


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
