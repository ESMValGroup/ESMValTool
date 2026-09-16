import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from basic_functions import (
    get_provenance_record,
    load_and_update_dict,
)

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

MIN_VALID_POINTS = 2  # Minimum number of valid points for interpolation

def plot_ts(cfg, plot_dict, title, y_title, output_basename):
    """
    Plot monthly values for all datasets in a single figure.
    """

    logger.info("Plotting monthly timeseries")
    month_labels = [
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
    obs_datasets = {"ERA5", "NCEP", "HadISST", "EN4"}
    common_x = np.arange(12)

    plt.figure(figsize=(10, 5))
    highlight_colors = [
        "tab:orange",
        "tab:red",
        "tab:green",
        "tab:brown",
        "tab:pink",
        "tab:olive",
        "tab:gray",
    ]
    highlight_idx = 0
    input_filenames = set()
    max_months = 12
    multimodel_profiles = []

    for dataset, dict_info in plot_dict.items():
        cube = dict_info["cube"]
        file = dict_info["filename"]
        input_filenames.update(file if isinstance(file, list) else [file])

        values = np.asarray(cube.data, dtype=float).squeeze()
        if values.ndim != 1:
            logger.warning(
                "Skipping %s: expected 1D monthly data, got shape %s.",
                dataset,
                values.shape,
            )
            continue

        n_months = values.shape[0]
        x = np.arange(n_months)
        max_months = max(max_months, n_months)

        valid = np.isfinite(values)
        if dataset not in obs_datasets \
            and np.count_nonzero(valid) >= MIN_VALID_POINTS:
            interp_vals = np.interp(
                common_x,
                x[valid],
                values[valid],
                left=np.nan,
                right=np.nan,
            )
            multimodel_profiles.append(interp_vals)

        if dataset in obs_datasets:
            color = "black"
            linewidth = 2.5
            alpha = 1.0
        elif dataset in cfg.get("highlight_datasets", []):
            color = highlight_colors[highlight_idx % len(highlight_colors)]
            highlight_idx += 1
            linewidth = 1.2
            alpha = 0.6
        else:
            # Skip plotting non-highlight models while 
            # still using them for multimodel stats.
            continue

        plt.plot(
            x,
            values,
            label=dataset,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
        )

    if multimodel_profiles:
        profiles = np.asarray(multimodel_profiles, dtype=float)
        multimodel_median = np.nanmedian(profiles, axis=0)
        multimodel_std = np.nanstd(profiles, axis=0)

        lower = multimodel_median - multimodel_std
        upper = multimodel_median + multimodel_std

        plt.fill_between(
            common_x,
            lower,
            upper,
            color="tab:blue",
            alpha=0.2,
            linewidth=0,
            label="Multimodel median ±1 std",
        )

        plt.plot(
            common_x,
            multimodel_median,
            label="Multimodel median",
            color="tab:blue",
            linewidth=3,
        )

    plt.xlim(-0.5, max_months - 0.5)
    if max_months == len(month_labels):
        tick_labels = month_labels
    else:
        tick_labels = [str(i + 1) for i in range(max_months)]
    plt.xticks(np.arange(max_months), tick_labels)
    plt.xlabel("Month", fontsize=12)
    plt.ylabel(y_title, fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True)
    plt.tight_layout()

    provenance_record = get_provenance_record(
        output_basename, list(input_filenames)
    )
    save_figure(output_basename, provenance_record, cfg)
    logger.info(f"Monthly plot saved: {output_basename}")
    plt.close()


def plot_ts_gradient(
    cfg, plot_dict_west, plot_dict_east, title, y_title, output_basename
):
    """
    Plot monthly gradient (west - east) for all datasets in a single figure.
    """

    logger.info("Plotting monthly gradient timeseries")
    month_labels = [
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
    obs_datasets = {"ERA5", "NCEP", "HadISST", "EN4"}
    common_x = np.arange(12)

    plt.figure(figsize=(10, 5))
    highlight_colors = [
        "tab:orange",
        "tab:red",
        "tab:green",
        "tab:brown",
        "tab:pink",
        "tab:olive",
        "tab:gray",
    ]
    highlight_idx = 0
    input_filenames = set()
    max_months = 12
    multimodel_profiles = []

    # Compute gradient for each dataset
    for dataset in plot_dict_west:
        if dataset not in plot_dict_east:
            logger.warning(
                "Dataset %s not found in east dict, skipping gradient.",
                dataset,
            )
            continue

        dict_info_west = plot_dict_west[dataset]
        dict_info_east = plot_dict_east[dataset]

        cube_west = dict_info_west["cube"]
        cube_east = dict_info_east["cube"]
        file_west = dict_info_west["filename"]
        file_east = dict_info_east["filename"]

        input_filenames.update(
            file_west if isinstance(file_west, list) else [file_west]
        )
        input_filenames.update(
            file_east if isinstance(file_east, list) else [file_east]
        )

        values_west = np.asarray(cube_west.data, dtype=float).squeeze()
        values_east = np.asarray(cube_east.data, dtype=float).squeeze()

        if values_west.ndim != 1 or values_east.ndim != 1:
            logger.warning(
                "Skipping %s: expected 1D monthly data.",
                dataset,
            )
            continue

        # Ensure same length, use shorter length
        min_len = min(len(values_west), len(values_east))
        values = values_west[:min_len] - values_east[:min_len]

        n_months = values.shape[0]
        x = np.arange(n_months)
        max_months = max(max_months, n_months)

        valid = np.isfinite(values)
        if dataset not in obs_datasets \
            and np.count_nonzero(valid) >= MIN_VALID_POINTS:
            interp_vals = np.interp(
                common_x,
                x[valid],
                values[valid],
                left=np.nan,
                right=np.nan,
            )
            multimodel_profiles.append(interp_vals)

        if dataset in obs_datasets:
            color = "black"
            linewidth = 2.5
            alpha = 1.0
        elif dataset in cfg.get("highlight_datasets", []):
            color = highlight_colors[highlight_idx % len(highlight_colors)]
            highlight_idx += 1
            linewidth = 1.2
            alpha = 0.6
        else:
            # Skip plotting non-highlight models while 
            # still using them for multimodel stats.
            continue

        plt.plot(
            x,
            values,
            label=dataset,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
        )

    if multimodel_profiles:
        profiles = np.asarray(multimodel_profiles, dtype=float)
        multimodel_median = np.nanmedian(profiles, axis=0)
        multimodel_std = np.nanstd(profiles, axis=0)

        lower = multimodel_median - multimodel_std
        upper = multimodel_median + multimodel_std

        plt.fill_between(
            common_x,
            lower,
            upper,
            color="tab:blue",
            alpha=0.2,
            linewidth=0,
            label="Multimodel median ±1 std",
        )

        plt.plot(
            common_x,
            multimodel_median,
            label="Multimodel median",
            color="tab:blue",
            linewidth=3,
        )

    plt.xlim(-0.5, max_months - 0.5)
    if max_months == len(month_labels):
        tick_labels = month_labels
    else:
        tick_labels = [str(i + 1) for i in range(max_months)]
    plt.xticks(np.arange(max_months), tick_labels)
    plt.xlabel("Month", fontsize=12)
    plt.ylabel(y_title, fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True)
    plt.tight_layout()

    provenance_record = get_provenance_record(
        output_basename, list(input_filenames)
    )
    save_figure(output_basename, provenance_record, cfg)
    logger.info("Monthly gradient plot saved: %s", output_basename)
    plt.close()


def main(cfg):
    """Plot monthly climatologies for multiple datasets and observations."""
    input_data = cfg["input_data"].values()
    grouped_data = group_metadata(input_data, "dataset")
    (
        ceio_wind_monthly,
        east_sst_monthly,
        west_sst_monthly,
        east_pr_monthly,
        west_pr_monthly,
    ) = {}, {}, {}, {}, {}
    for _group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, "ceio_wind_monthly", ceio_wind_monthly)
        load_and_update_dict(group_md, "east_sst_monthly", east_sst_monthly)
        load_and_update_dict(group_md, "west_sst_monthly", west_sst_monthly)
        load_and_update_dict(group_md, "east_pr_monthly", east_pr_monthly)
        load_and_update_dict(group_md, "west_pr_monthly", west_pr_monthly)

    # Plot results for all datasets
    plot_ts(
        cfg,
        ceio_wind_monthly,
        "Monthly CEIO zonal wind speed",
        "Zonal wind speed / m $\\mathregular{s^{-1}}$",
        "ceio_winds_monthly",
    )
    plot_ts(
        cfg,
        east_sst_monthly,
        "Monthly SST in the EEIO",
        "SST / $^\\circ$C",
        "east_sst_monthly",
    )
    plot_ts(
        cfg,
        west_sst_monthly,
        "Monthly SST in the WEIO",
        "SST / $^\\circ$C",
        "west_sst_monthly",
    )
    plot_ts(
        cfg,
        east_pr_monthly,
        "Monthly precipitation in the EEIO",
        "Precipitation / mm day$^{-1}$",
        "east_pr_monthly",
    )
    plot_ts(
        cfg,
        west_pr_monthly,
        "Monthly precipitation in the WEIO",
        "Precipitation / mm day$^{-1}$",
        "west_pr_monthly",
    )

    # Plot gradients (west - east)
    plot_ts_gradient(
        cfg,
        west_sst_monthly,
        east_sst_monthly,
        "Monthly SST gradient (WEIO - EEIO)",
        "SST gradient / $^\\circ$C",
        "sst_gradient_monthly",
    )
    plot_ts_gradient(
        cfg,
        west_pr_monthly,
        east_pr_monthly,
        "Monthly precipitation gradient (WEIO - EEIO)",
        "Precipitation gradient / mm day$^{-1}$",
        "pr_gradient_monthly",
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
