from fileinput import filename
import logging
from pathlib import Path
from pprint import pformat
import matplotlib.pyplot as plt
import iris # type: ignore
import numpy as np
import numpy.ma as ma

from basic_functions import (get_provenance_record, iso_depth_4d, 
                             load_and_update_dict,
                             load_data)

from esmvaltool.diag_scripts.shared import ( # type: ignore
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot # type: ignore

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

def calculate_taylor_stats(model_cube, obs_cube):
    """
    Calculate Taylor diagram statistics for a given model and observation cube.
    Returns a dictionary with keys 'correlation', 'stddev', and 'rmse'.
    """

    model_data = model_cube.data.flatten()
    obs_data = obs_cube.data.flatten()
    model_mask = np.isfinite(model_data)
    obs_mask = np.isfinite(obs_data)
    model_data = model_data[model_mask]
    obs_data = obs_data[obs_mask]

    # Remove mean from model and obs data
    model_anom = model_data - np.mean(model_data)
    obs_anom = obs_data - np.mean(obs_data)

    # Calculate statistics
    correlation = np.corrcoef(model_anom, obs_anom)[0, 1]
    std_model = np.std(model_anom)
    std_obs = np.std(obs_anom)
    
    # Calculate centred RMS difference
    crmsd = np.sqrt(np.mean((model_anom - obs_anom) ** 2))

    return {
        'correlation': correlation,
        'stddev': std_model,
        'rmse': crmsd,
    }

def plot_taylor(cfg, plot_dict, title, output_basename):
    """Plot Taylor diagram statistics for all datasets in a single figure."""
    logger.info("Plotting Taylor diagram")
    obs_datasets = {'ERA5', 'NCEP', 'HadISST', 'EN4'}

    input_filenames = set()
    obs_entries = {}
    model_entries = {}

    for dataset, dict_info in plot_dict.items():
        file = dict_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])
        if dataset in obs_datasets:
            obs_entries[dataset] = dict_info
        else:
            model_entries[dataset] = dict_info

    if not obs_entries:
        logger.warning("No observations found for Taylor diagram, skipping %s", output_basename)
        return
    if not model_entries:
        logger.warning("No model datasets found for Taylor diagram, skipping %s", output_basename)
        return

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, polar=True)

    model_colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173',
        '#3182bd', '#31a354', '#756bb1', '#636363', '#e6550d',
    ]
    obs_markers = ['o', 's', '^', 'D', 'P', 'X', '*']

    max_std_ratio = 1.2
    obs_names = list(obs_entries.keys())
    for obs_name, obs_info in obs_entries.items():
        obs_cube = obs_info['cube']
        obs_data = obs_cube.data.flatten()
        obs_mask = np.isfinite(obs_data)
        obs_data = obs_data[obs_mask]
        obs_anom = obs_data - np.mean(obs_data)
        obs_std = np.std(obs_anom)

        # obs_std = np.std(np.asarray(obs_cube.data, dtype=float).flatten())
        logger.debug("Observation %s: stddev = %g", obs_name, obs_std)
        if not np.isfinite(obs_std) or obs_std == 0:
            logger.warning("Skipping observation %s due to invalid standard deviation", obs_name)
            continue

        for model_idx, (model_name, model_info) in enumerate(model_entries.items()):
            stats = calculate_taylor_stats(model_info['cube'], obs_cube)
            corr = stats['correlation']
            #logger.info("Model %s vs Observation %s: correlation = %g", model_name, obs_name, corr)

            if not np.isfinite(corr):
                continue
            corr = np.clip(corr, -1.0, 1.0)
            if corr < 0.0:
                logger.debug(
                    "Skipping %s vs %s due to negative correlation (%g) outside first quadrant",
                    model_name,
                    obs_name,
                    corr,
                )
                continue

            std_ratio = stats['stddev'] / obs_std
            logger.info("Model %s vs Observation %s: model stddev = %g, obs stddev = %g, stddev ratio = %g",
                model_name,
                obs_name,
                stats['stddev'],
                obs_std,
                std_ratio,
            )
            
            if not np.isfinite(std_ratio):
                continue

            theta = np.arccos(corr)
            marker = obs_markers[obs_names.index(obs_name) % len(obs_markers)]
            color = model_colors[model_idx % len(model_colors)]
            ax.plot(
                theta,
                std_ratio,
                linestyle='None',
                marker=marker,
                markersize=7,
                color=color,
                label=f"{model_name}",
                alpha=0.85,
            )
            max_std_ratio = max(max_std_ratio, std_ratio)

        obs_marker = obs_markers[obs_names.index(obs_name) % len(obs_markers)]
        ax.plot(
            0.0,
            1.0,
            linestyle='None',
            marker=obs_marker,
            markersize=10,
            color='black',
            label=f"{obs_name} reference",
        )

    max_std_ratio *= 1.1
    ax.set_xlim(0, np.pi / 2)
    ax.set_rlim(0, np.min([max_std_ratio, 2.0]))
    corr_ticks = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
    theta_ticks = np.arccos(corr_ticks[::-1])
    ax.set_xticks(theta_ticks)
    ax.set_xticklabels([f"{c:.2g}" for c in corr_ticks[::-1]])

    # Place correlation label above the top-right arc.
    ax.text(
        np.deg2rad(45),
        max_std_ratio * 1.10,
        "Correlation",
        ha='center',
        va='center',
    )

    # Centered-RMSD contours in normalized Taylor space.
    contour_levels = [0.5, 1.0, 1.5, 2.0]
    theta_grid = np.linspace(0.0, np.pi / 2, 400)
    for level in contour_levels:
        arc_term = level**2 - np.sin(theta_grid) ** 2
        valid = arc_term >= 0
        if not np.any(valid):
            continue

        radius = np.full_like(theta_grid, np.nan, dtype=float)
        radius[valid] = np.cos(theta_grid[valid]) + np.sqrt(arc_term[valid])
        ax.plot(
            theta_grid,
            radius,
            color='0.7',
            linestyle='--',
            linewidth=1.0,
            alpha=0.8,
        )

    ax.set_ylim(0.0, max_std_ratio)
    ax.set_ylabel("Normalised standard deviation ($\\sigma / \\sigma_{obs}$)", labelpad=30)
    ax.grid(True, alpha=0.4)

    ax.set_title(title)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.35, 1.0), loc='upper left')

    plt.tight_layout()

    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info("Taylor diagram saved: %s", output_basename)
    plt.close(fig)


def main(cfg):
    """Plot monthly climatologies for multiple datasets and observations."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    io_wind_son, io_wind_annual, io_sst_son, io_sst_annual = {}, {}, {}, {}
    for group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, 'IO_wind_son', io_wind_son)
        load_and_update_dict(group_md, 'IO_wind_annual', io_wind_annual)
        load_and_update_dict(group_md, 'IO_sst_son', io_sst_son)
        load_and_update_dict(group_md, 'IO_sst_annual', io_sst_annual)

    print(io_wind_son)
    logger.info("Data loaded, now plotting.")
    # Plot results for all datasets
    plot_taylor(
        cfg,
        io_wind_son,
        'Indian Ocean SON zonal winds',
        'taylor_io_wind_son',
    )
    plot_taylor(
        cfg,
        io_sst_son,
        'Indian Ocean SON SST',
        'taylor_io_sst_son',
    )
    plot_taylor(
        cfg,
        io_wind_annual,
        'Indian Ocean Annual zonal winds',
        'taylor_io_wind_annual',
    )
    plot_taylor(
        cfg,
        io_sst_annual,
        'Indian Ocean Annual SST',
        'taylor_io_sst_annual',
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

