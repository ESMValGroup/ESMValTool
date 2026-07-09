import logging
from pathlib import Path
import matplotlib.pyplot as plt
import iris # type: ignore
import numpy as np
import numpy.ma as ma

from basic_functions import get_provenance_record, load_and_update_dict, iso_depth_3d, iso_depth_4d

from esmvaltool.diag_scripts.shared import ( # type: ignore
    group_metadata,
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

OBS_DATASETS = {'ERA5', 'NCEP', 'HadISST', 'EN4'}
SEASON_LABELS = {0: 'DJF', 1: 'MAM', 2: 'JJA', 3: 'SON'}


def calculate_taylor_stats(model_cube, obs_cube):
    """
    Calculate Taylor diagram statistics for a given model and observation cube.
    Returns a dictionary with keys 'correlation', 'stddev', 'std_obs', and 'rmse'.
    """

    model_data = np.asarray(model_cube.data, dtype=float).flatten()
    obs_data = np.asarray(obs_cube.data, dtype=float).flatten()
    valid_mask = (
        np.isfinite(model_data)
        & np.isfinite(obs_data)
        & (model_data <= 1e10)
        & (obs_data <= 1e10)
    )


    if not np.any(valid_mask):
        return {
            'correlation': np.nan,
            'stddev': np.nan,
            'std_obs': np.nan,
            'rmse': np.nan,
        }

    model_data = model_data[valid_mask]
    obs_data = obs_data[valid_mask]
    if model_data.size < 2:
        return {
            'correlation': np.nan,
            'stddev': np.nan,
            'std_obs': np.nan,
            'rmse': np.nan,
        }

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
        'std_obs': std_obs,
        'rmse': crmsd,
    }


def _separate_obs_and_models(plot_dict):
    """Split plot_dict metadata into obs entries, model entries and file list."""
    input_filenames = set()
    obs_entries = {}
    model_entries = {}

    for dataset, dict_info in plot_dict.items():
        file = dict_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])
        if dataset in OBS_DATASETS:
            obs_entries[dataset] = dict_info
        else:
            model_entries[dataset] = dict_info

    return input_filenames, obs_entries, model_entries


def _slice_cube_for_season(cube, season_number):
    """Return cube sliced to one season; return None if season not present."""
    try:
        coord = cube.coord('season_number')
    except iris.exceptions.CoordinateNotFoundError:
        return cube

    season_points = np.atleast_1d(coord.points).astype(int)
    if season_points.size == 1:
        return cube if int(season_points[0]) == int(season_number) else None

    matches = np.where(season_points == int(season_number))[0]
    if matches.size == 0:
        return None

    coord_dims = cube.coord_dims('season_number')
    if not coord_dims:
        return cube if int(season_points[0]) == int(season_number) else None

    dim = coord_dims[0]
    slicer = [slice(None)] * cube.ndim
    slicer[dim] = int(matches[0])
    return cube[tuple(slicer)]


def _get_available_seasons(obs_entries, model_entries):
    """Return sorted season numbers found in obs/model cubes."""
    seasons = set()
    for entries in (obs_entries, model_entries):
        for info in entries.values():
            cube = info['cube']
            try:
                points = np.atleast_1d(cube.coord('season_number').points).astype(int)
                seasons.update(points.tolist())
            except iris.exceptions.CoordinateNotFoundError:
                continue
    return sorted(seasons)


def _plot_taylor_panel(ax, obs_entries, model_entries, panel_title, season_number=None):
    """Plot one Taylor diagram panel on the provided axis."""
    model_colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(model_entries)))  # type: ignore[attr-defined]
    obs_markers = ['o', 's', '^', 'D', 'P', 'X', '*']
    obs_names = list(obs_entries.keys())

    max_std_ratio = 1.2
    plotted_points = 0

    for obs_name, obs_info in obs_entries.items():
        obs_cube = obs_info['cube']
        if season_number is not None:
            obs_cube = _slice_cube_for_season(obs_cube, season_number)
        if obs_cube is None:
            continue

        for model_idx, (model_name, model_info) in enumerate(model_entries.items()):
            model_cube = model_info['cube']
            if season_number is not None:
                model_cube = _slice_cube_for_season(model_cube, season_number)
            if model_cube is None:
                continue

            stats = calculate_taylor_stats(model_cube, obs_cube)
            corr = stats['correlation']
            obs_std = stats['std_obs']

            if not np.isfinite(obs_std) or obs_std == 0 or not np.isfinite(corr):
                continue

            corr = np.clip(corr, -1.0, 1.0)
            if corr < 0.0:
                continue

            std_ratio = stats['stddev'] / obs_std
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
                label=model_name,
                alpha=0.85,
            )
            max_std_ratio = max(max_std_ratio, std_ratio)
            plotted_points += 1

        # Plot obs reference point at (corr=1, std_ratio=1).
        obs_marker = obs_markers[obs_names.index(obs_name) % len(obs_markers)]
        ax.plot(
            0.0,
            1.0,
            linestyle='None',
            marker=obs_marker,
            markersize=10,
            color='black',
            label=f'{obs_name} reference',
        )

    max_std_ratio *= 1.1
    ax.set_xlim(0, np.pi / 2)
    ax.set_rlim(0, min(max_std_ratio, 2.0))

    corr_ticks = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0])
    theta_ticks = np.arccos(corr_ticks[::-1])
    ax.set_xticks(theta_ticks)
    ax.set_xticklabels([f'{c:.2g}' for c in corr_ticks[::-1]])

    ax.text(np.deg2rad(45), max_std_ratio * 1.10, 'Correlation', ha='center', va='center')

    contour_levels = [0.5, 1.0, 1.5, 2.0]
    theta_grid = np.linspace(0.0, np.pi / 2, 400)
    for level in contour_levels:
        arc_term = level**2 - np.sin(theta_grid) ** 2
        valid = arc_term >= 0
        if not np.any(valid):
            continue
        radius = np.full_like(theta_grid, np.nan, dtype=float)
        radius[valid] = np.cos(theta_grid[valid]) + np.sqrt(arc_term[valid])
        ax.plot(theta_grid, radius, color='0.7', linestyle='--', linewidth=1.0, alpha=0.8)

    ax.set_ylim(0.0, max_std_ratio)
    ax.set_ylabel('Normalised standard deviation ($\\sigma / \\sigma_{obs}$)', labelpad=20)
    ax.grid(True, alpha=0.4)
    ax.set_title(panel_title)

    if plotted_points == 0:
        ax.text(np.deg2rad(30), 0.9, 'No valid points', ha='center', va='center')

    return max_std_ratio


def plot_taylor(cfg, plot_dict, title, output_basename):
    """Plot Taylor diagram(s); seasonal data becomes multi-panel figure."""
    logger.info('Plotting Taylor diagram: %s', output_basename)

    input_filenames, obs_entries, model_entries = _separate_obs_and_models(plot_dict)

    if not obs_entries:
        logger.warning('No observations found for Taylor diagram, skipping %s', output_basename)
        return
    if not model_entries:
        logger.warning('No model datasets found for Taylor diagram, skipping %s', output_basename)
        return

    seasons = _get_available_seasons(obs_entries, model_entries)

    if seasons:
        n_panels = len(seasons)
        n_cols = 2 if n_panels > 1 else 1
        n_rows = int(np.ceil(n_panels / n_cols))
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(7 * n_cols, 6 * n_rows),
            subplot_kw={'projection': 'polar'},
            squeeze=False,
        )
        flat_axes = axes.ravel()

        for idx, season_number in enumerate(seasons):
            season_label = SEASON_LABELS.get(int(season_number), f'Season {int(season_number)}')
            _plot_taylor_panel(
                flat_axes[idx],
                obs_entries,
                model_entries,
                panel_title=season_label,
                season_number=int(season_number),
            )

        for idx in range(n_panels, len(flat_axes)):
            flat_axes[idx].set_visible(False)

        handles, labels = flat_axes[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        fig.legend(
            by_label.values(),
            by_label.keys(),
            loc='center right',
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
        )
        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=(0.0, 0.0, 0.88, 0.96))
    else:
        fig = plt.figure(figsize=(12, 7))
        ax = fig.add_subplot(111, polar=True)
        _plot_taylor_panel(ax, obs_entries, model_entries, panel_title=title)

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.35, 1.0), loc='upper left')
        fig.tight_layout()

    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info('Taylor diagram saved: %s', output_basename)
    plt.close(fig)

def _replace_fill_values(cube, fill_value=1e20):
    """Replace values >= fill_value with NaN in a cube's data array."""
    data = np.asarray(cube.data, dtype=float)
    data[data >= fill_value] = np.nan
    cube.data = ma.masked_invalid(data)
    return cube

def _create_iso_depth_dict(
    cfg,
    plot_dict,
    iso_level=20.0,
    time_measure=None,
):
    """Create a plot_dict with 4D cubes converted to isotherm depth.

    Optionally applies a robust local-neighborhood coastal mask to remove
    isolated shallow outliers that can skew zonal means and multi-model stats.
    """
    new_plot_dict = {}
    for dataset, info in plot_dict.items():
        cube = info['cube']
        if time_measure is not None:
            iso_cube = iso_depth_4d(cube, iso_level, time_measure=time_measure)
        else:
            iso_cube = iso_depth_3d(cube, iso_level)
        
        iso_cube = _replace_fill_values(iso_cube)
        new_plot_dict[dataset] = {'cube': iso_cube, 'filename': info['filename']}

    return new_plot_dict



def main(cfg):
    """Create Taylor diagrams for seasonal and annual Indian Ocean diagnostics."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    io_wind_seas, io_pr_seas, io_sst_seas, io_theta_seas = {}, {}, {}, {}
    io_wind_annual, io_pr_annual, io_sst_annual, io_theta_annual = {}, {}, {}, {}

    for group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, 'IO_wind_seas', io_wind_seas)
        load_and_update_dict(group_md, 'IO_wind_annual', io_wind_annual)
        load_and_update_dict(group_md, 'IO_sst_seas', io_sst_seas)
        load_and_update_dict(group_md, 'IO_sst_annual', io_sst_annual)
        load_and_update_dict(group_md, 'IO_pr_seas', io_pr_seas)
        load_and_update_dict(group_md, 'IO_pr_annual', io_pr_annual)
        load_and_update_dict(group_md, 'IO_theta_seas', io_theta_seas)
        load_and_update_dict(group_md, 'IO_theta_annual', io_theta_annual)

    logger.info("Data loaded, now plotting.")

    io_t20d_seas = _create_iso_depth_dict(cfg, io_theta_seas, iso_level=20.0, time_measure='season_number')
    io_t20d_annual = _create_iso_depth_dict(cfg, io_theta_annual, iso_level=20.0, time_measure=None)

    seasonal_plots = {
        'wind': io_wind_seas,
        'pr': io_pr_seas,
        'sst': io_sst_seas,
        't20d': io_t20d_seas,
    }
    annual_plots = {
        'wind': io_wind_annual,
        'pr': io_pr_annual,
        'sst': io_sst_annual,
        't20d': io_t20d_annual,
    }

    for var_key, plot_dict in seasonal_plots.items():
        plot_taylor(
            cfg,
            plot_dict,
            f'Indian Ocean Seasonal {var_key.upper()}',
            f'taylor_io_{var_key}_seas',
        )

    for var_key, plot_dict in annual_plots.items():
        plot_taylor(
            cfg,
            plot_dict,
            f'Indian Ocean Annual {var_key.upper()}',
            f'taylor_io_{var_key}_annual',
        )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

