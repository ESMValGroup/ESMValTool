import logging
from pathlib import Path
from pprint import pformat
import matplotlib.pyplot as plt
import iris # type: ignore
import numpy as np
import numpy.ma as ma
import cmocean
import os

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


OBS_DATASETS = {'NCEP', 'HadISST', 'EN4'}


def is_mohc_dataset(dataset):
    """Return True if dataset likely belongs to MOHC."""
    dataset_upper = dataset.upper()
    mohc_markers = ("MOHC", "HADGEM", "HADCM", "UKESM")
    return any(marker in dataset_upper for marker in mohc_markers)


def _replace_fill_values(cube, fill_value=1e20):
    """Replace values >= fill_value with NaN in a cube's data array."""
    data = np.asarray(cube.data, dtype=float)
    data[data >= fill_value] = np.nan
    cube.data = ma.masked_invalid(data)
    return cube


def _collapse_latitude(cube):
    """Collapse latitude by mean while ignoring masked/invalid values."""
    try:
        cube.coord('latitude')
        masked = ma.masked_invalid(np.asarray(cube.data, dtype=float))
        return cube.copy(data=masked).collapsed('latitude', iris.analysis.MEAN)
    except iris.exceptions.CoordinateNotFoundError:
        return cube


def _coord_edges(centres):
    """Convert coordinate centre points to bin edges for use with pcolormesh."""
    centres = np.asarray(centres)
    if centres.size < 2:
        return np.array([centres[0] - 0.5, centres[0] + 0.5])
    d = np.diff(centres)
    return np.concatenate([
        [centres[0] - d[0] / 2],
        centres[:-1] + d / 2,
        [centres[-1] + d[-1] / 2],
    ])


def _compute_multimodel_bias_and_std(model_cubes, obs_data):
    """Compute multi-model median bias and inter-model std dev against obs_data.

    Parameters
    ----------
    model_cubes : list of iris.cube.Cube
        Model cubes with fill values already replaced by NaN.
    obs_data : np.ndarray
        Observation array with the same spatial/temporal shape as each model cube.

    Returns
    -------
    mm_median_bias : np.ndarray
        Multi-model median of (model − obs).
    mm_std : np.ndarray
        Inter-model standard deviation of model values.
    """
    bias_list = []
    model_list = []
    for mc in model_cubes:
        model_data = np.ma.filled(np.ma.asarray(mc.data, dtype=float), np.nan)
        bias_list.append(model_data - obs_data)
        model_list.append(model_data)
    bias_stack = np.array(bias_list)   # (n_models, ...)
    model_stack = np.array(model_list) # (n_models, ...)
    mm_median_bias = np.nanmedian(bias_stack, axis=0)
    mm_std = np.nanstd(model_stack, axis=0)
    return mm_median_bias, mm_std


def _separate_obs_and_models(plot_dict):
    """Split a plot_dict into (obs_name, obs_cube, obs_dataset, model_cubes, input_filenames)."""
    obs_entry = None
    model_cubes = []
    input_filenames = set()
    for dataset, info in plot_dict.items():
        file = info['filename']
        print(f"Processing dataset {dataset} with file(s): {file}")
        input_filenames.update(file if isinstance(file, list) else [file])
        cube = info['cube']
        if os.path.basename(file).startswith("OBS"):
            obs_entry = (dataset, cube)
        else:
            model_cubes.append(cube)
    return obs_entry, model_cubes, input_filenames


def _plot_lontime_panel(ax, fig, lon_edges, month_edges, data, cmap,
                        colorbar_label, panel_title, vmin=None, vmax=None):
    """Render a single lon-time pcolormesh panel with a horizontal colorbar.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw onto.
    fig : matplotlib.figure.Figure
        Parent figure (needed for colorbar).
    lon_edges, month_edges : array-like
        Bin edges for the x (longitude) and y (month) axes.
    data : np.ndarray
        2-D array shaped (month, longitude).
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : float
        Colour scale limits.
    colorbar_label : str
        Label for the horizontal colorbar.
    panel_title : str
        Title shown above the panel.
    """
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    n_months = data.shape[0]

    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)

    im = ax.pcolormesh(lon_edges, month_edges, data,
                       cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
    fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, label=colorbar_label)
    ax.set_title(panel_title, fontsize=12)
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Month')
    ax.set_yticks(list(range(1, n_months + 1))[:n_months])
    ax.set_yticklabels(month_labels[:n_months], fontsize=9)
    return im


def plot_lon_time_multimodel(cfg, plot_dict, cmap_list, title, output_basename,
                              variable=None, obs_vmin=None, obs_vmax=None,
                              bias_vlim=None,
                              std_vmax=None):
    """Plot lon-time climatology: obs, MM-median bias, and inter-model std dev.

    Creates a three-panel figure:
      Panel 1 — Observed climatology (lon x month).
      Panel 2 — Multi-model median bias (model - obs).
      Panel 3 — Inter-model standard deviation.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    cmap_list : list of str
        Three colourmap names for [obs, bias, std dev] panels.
    title : str
        Figure suptitle.
    output_basename : str
        Stem used when saving the figure via save_figure().
    variable : str, optional
        Variable name ('tos', 'ua', etc.) for determining default colour limits.
    obs_vmin, obs_vmax : float, optional
        Colour scale limits for the obs panel.  Auto-derived when omitted.
    bias_vlim : float, optional
        Symmetric colour limit ±bias_vlim for the bias panel.  Auto-derived when omitted.
    std_vmax : float, optional
        Upper colour limit for the std dev panel.  Auto-derived when omitted.
    """
    logger.info("Plotting lon-time plots: %s", output_basename)

    obs_entry, model_cubes, input_filenames = _separate_obs_and_models(plot_dict)
    if obs_entry is None:
        logger.warning("No obs found for lon-time plot %s, skipping.", output_basename)
        return
    if not model_cubes:
        logger.warning("No model data found for lon-time plot %s, skipping.", output_basename)
        return

    obs_dataset, raw_obs_cube = obs_entry

    obs_cube = _collapse_latitude(_replace_fill_values(raw_obs_cube.copy()))
    model_cubes = [_collapse_latitude(_replace_fill_values(mc.copy())) for mc in model_cubes]
    obs_data = np.ma.filled(np.ma.asarray(obs_cube.data, dtype=float), np.nan)  # (month, lon)
    mm_median_bias, mm_std = _compute_multimodel_bias_and_std(model_cubes, obs_data)

    # Coordinate edges for pcolormesh.
    try:
        lon_centres = obs_cube.coord('longitude').points
    except iris.exceptions.CoordinateNotFoundError:
        lon_centres = np.arange(obs_data.shape[1])
    try:
        month_centres = obs_cube.coord('month_number').points
    except iris.exceptions.CoordinateNotFoundError:
        month_centres = np.arange(1, obs_data.shape[0] + 1)

    lon_edges = _coord_edges(lon_centres)
    month_edges = _coord_edges(month_centres)

    # Derive colour limits, with optional overrides.
    if obs_vmin is None or obs_vmax is None:
        if variable == 'tos':
            _obs_vmin, _obs_vmax = 25.0, 30.0
        elif variable in ('ua', 'va', 'wind'):
            vabs = np.nanmax(np.abs(obs_data))
            _obs_vmin, _obs_vmax = -vabs, vabs
        else:
            _obs_vmin, _obs_vmax = np.nanmin(obs_data), np.nanmax(obs_data)
    obs_vmin = obs_vmin if obs_vmin is not None else _obs_vmin
    obs_vmax = obs_vmax if obs_vmax is not None else _obs_vmax

    if bias_vlim is None:
        # if variable == 'pr':
        #     bias_vlim = 5
        # else:
        bias_vlim = round(np.nanmax(np.abs(mm_median_bias)), 1)

    if std_vmax is None:
        std_vmax = np.nanmax(mm_std)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    _plot_lontime_panel(
        axes[0], fig, lon_edges, month_edges, obs_data,
        cmap=cmap_list[0], vmin=obs_vmin, vmax=obs_vmax,
        colorbar_label=obs_dataset,
        panel_title=f'Obs ({obs_dataset})',
    )
    _plot_lontime_panel(
        axes[1], fig, lon_edges, month_edges, mm_median_bias,
        cmap=cmap_list[1], vmin=-bias_vlim, vmax=bias_vlim,
        colorbar_label='Bias (model − obs)',
        panel_title='MM-median bias',
    )
    _plot_lontime_panel(
        axes[2], fig, lon_edges, month_edges, mm_std,
        cmap=cmap_list[2], vmin=0, vmax=std_vmax,
        colorbar_label='Std dev (models)',
        panel_title='Inter-model std dev',
    )

    fig.suptitle(title, fontsize=14)
    provenance_record = get_provenance_record(output_basename, sorted(list(input_filenames)))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info("Lon-time plot saved: %s", output_basename)
    plt.close(fig)


def plot_monthly_maps(cfg, data, lon_centres, lat_centres, cmap, title,
                      output_basename, input_filenames, vmin=None, vmax=None):
    """Plot 12 monthly map panels (3 rows × 4 columns) with a shared colorbar.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    data : np.ndarray
        3-D array shaped (month, latitude, longitude).
    lon_centres, lat_centres : array-like
        Coordinate centre points.
    cmap : str
        Matplotlib colormap name.
    title : str
        Figure suptitle and colorbar label.
    output_basename : str
        Stem used when saving the figure.
    input_filenames : set or list
        Source filenames for provenance tracking.
    vmin, vmax : float, optional
        Colour scale limits.  Derived from data when omitted.
    """
    data = np.asarray(data)
    if data.ndim != 3:
        logger.warning("Expected 3D data (month, lat, lon) for %s, got shape %s",
                       output_basename, data.shape)
        return

    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)

    n_months = min(12, data.shape[0])
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    lon_edges = _coord_edges(lon_centres)
    lat_edges = _coord_edges(lat_centres)

    fig, axes = plt.subplots(3, 4, figsize=(16, 10), constrained_layout=True)
    axes = axes.flatten()

    im = None
    for m in range(12):
        ax = axes[m]
        if m < n_months:
            im = ax.pcolormesh(lon_edges, lat_edges, data[m, :, :],
                               cmap=cmap, vmin=vmin, vmax=vmax, shading='flat')
            ax.set_title(month_labels[m], fontsize=10)
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
        else:
            ax.axis('off')

    if im is not None:
        cbar = fig.colorbar(im, ax=axes.tolist(), orientation='horizontal',
                            pad=0.06, shrink=0.9)
        cbar.set_label(title)

    fig.suptitle(title, fontsize=14)
    provenance_record = get_provenance_record(output_basename, sorted(list(input_filenames)))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info("Monthly map saved: %s", output_basename)
    plt.close(fig)


def plot_map_multimodel(cfg, plot_dict, cmap_list, title, output_basename):
    """Plot monthly maps for obs, MM-median bias, and inter-model std dev.

    Calls plot_monthly_maps three times — once per panel type — saving
    separate figures for obs, bias, and std dev.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary.
    plot_dict : dict
        Mapping of dataset name → {'cube': iris.cube.Cube, 'filename': ...}.
    cmap_list : list of str
        Three colourmap names for [obs, bias, std dev].
    title : str
        Base title string appended with the panel type for each figure.
    output_basename : str
        Base stem for output filenames; suffixes '_obs', '_bias', '_stddev' are appended.
    """
    logger.info("Plotting monthly map plots: %s", output_basename)

    obs_entry, model_cubes, input_filenames = _separate_obs_and_models(plot_dict)
    if obs_entry is None:
        logger.warning("No obs found for map plot %s, skipping.", output_basename)
        return
    if not model_cubes:
        logger.warning("No model data found for map plot %s, skipping.", output_basename)
        return

    obs_name, obs_cube = obs_entry
    obs_cube = _replace_fill_values(obs_cube.copy())
    model_cubes = [_replace_fill_values(mc.copy()) for mc in model_cubes]
    obs_data = np.ma.filled(np.ma.asarray(obs_cube.data, dtype=float), np.nan)

    mm_median_bias, mm_std = _compute_multimodel_bias_and_std(model_cubes, obs_data)

    try:
        lon_centres = obs_cube.coord('longitude').points
    except iris.exceptions.CoordinateNotFoundError:
        lon_centres = np.arange(obs_data.shape[-1])
    try:
        lat_centres = obs_cube.coord('latitude').points
    except iris.exceptions.CoordinateNotFoundError:
        lat_centres = np.arange(obs_data.shape[-2])

    plot_monthly_maps(
        cfg, obs_data, lon_centres, lat_centres,
        vmin=np.nanmin(obs_data), vmax=np.nanmax(obs_data),
        cmap=cmap_list[0],
        title=f'Obs ({obs_name}) - {title}',
        output_basename=output_basename + '_clim_obs',
        input_filenames=input_filenames,
    )
    bias_vlim = np.nanmax(np.abs(mm_median_bias))
    plot_monthly_maps(
        cfg, mm_median_bias, lon_centres, lat_centres,
        vmin=-bias_vlim, vmax=bias_vlim,
        cmap=cmap_list[1],
        title=f'MM-median bias - {title}',
        output_basename=output_basename + '_model_bias',
        input_filenames=input_filenames,
    )
    plot_monthly_maps(
        cfg, mm_std, lon_centres, lat_centres,
        vmin=0, vmax=np.nanmax(mm_std),
        cmap=cmap_list[2],
        title=f'Inter-model std dev - {title}',
        output_basename=output_basename + '_stddev',
        input_filenames=input_filenames,
    )

def _create_iso_depth_dict(
    cfg,
    plot_dict,
    iso_level=20.0,
    time_measure='month_number',
    apply_coastal_mask=True,
):
    """Create a plot_dict with 4D cubes converted to isotherm depth.

    Optionally applies a robust local-neighborhood coastal mask to remove
    isolated shallow outliers that can skew zonal means and multi-model stats.
    """
    new_plot_dict = {}
    for dataset, info in plot_dict.items():
        cube = info['cube']
        if cube.ndim == 4:
            iso_cube = iso_depth_4d(cube, iso_level, time_measure=time_measure)
            iso_cube = _replace_fill_values(iso_cube)
            if apply_coastal_mask:
                lat_all = iso_cube.coord('latitude').points
                lon_all = iso_cube.coord('longitude').points
                lat_mask = lat_all <= -4
                lon_mask = (lon_all >= 54) & (lon_all <= 56)

                box_mask = lat_mask[:, None] & lon_mask[None, :]  # (lat, lon)

                # Mask out the box
                iso_cube.data[:, box_mask] = np.nan

            new_plot_dict[dataset] = {'cube': iso_cube, 'filename': info['filename']}
        else:
            new_plot_dict[dataset] = info
    return new_plot_dict


def main(cfg):
    """Plot monthly climatologies for multiple datasets and observations."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    eio_wind_monthly, eio_sst_monthly, eio_theta_monthly, eio_pr_monthly = {}, {}, {}, {}
    for group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, 'eio_wind_monthly', eio_wind_monthly)
        load_and_update_dict(group_md, 'eio_sst_monthly', eio_sst_monthly)
        load_and_update_dict(group_md, 'eio_theta_monthly', eio_theta_monthly)
        load_and_update_dict(group_md, 'eio_pr_monthly', eio_pr_monthly)

    logger.info("Data loaded, now plotting.")

    eio_t20d_monthly = _create_iso_depth_dict(
        cfg,
        eio_theta_monthly,
        iso_level=20.0,
        time_measure='month_number',
        apply_coastal_mask=True,
    )

    plot_lon_time_multimodel(
        cfg,
        eio_wind_monthly,
        cmap_list=['cmo.delta', 'BrBG', 'RdPu'],
        title='Indian Ocean equatorial zonal wind — monthly climatology',
        output_basename='lon_time_eio_wind',
        variable='ua',
    )
    plot_lon_time_multimodel(
        cfg,
        eio_sst_monthly,
        cmap_list=['RdYlBu_r', 'RdBu_r', 'RdPu'],
        title='Indian Ocean equatorial SST — monthly climatology',
        output_basename='lon_time_eio_sst',
        variable='tos',
    )
    plot_lon_time_multimodel(
        cfg,
        eio_t20d_monthly,
        cmap_list=['cmo.deep', 'cmo.tarn', 'RdPu'],
        title='Indian Ocean equatorial 20°C isotherm depth — monthly climatology',
        output_basename='lon_time_eio_t20d',
        variable='t20d',
    )
    plot_lon_time_multimodel(
        cfg,
        eio_pr_monthly,
        cmap_list=['cmo.rain', 'BrBG', 'RdPu'],
        title='Indian Ocean equatorial precipitation — monthly climatology',
        output_basename='lon_time_eio_pr',
        variable='pr',
        bias_vlim=0.00008
    )

    plot_map_multimodel(
        cfg,
        eio_wind_monthly,
        cmap_list=['BrBG', 'BrBG', 'RdPu'],
        title='Indian Ocean equatorial zonal wind — monthly climatology',
        output_basename='map_eio_wind',
    )
    plot_map_multimodel(
        cfg,
        eio_sst_monthly,
        cmap_list=['RdYlBu_r', 'RdBu_r', 'RdPu'],
        title='Indian Ocean equatorial SST — monthly climatology',
        output_basename='map_eio_sst',
    )
    plot_map_multimodel(
        cfg,
        eio_t20d_monthly,
        cmap_list=['cmo.deep', 'cmo.tarn', 'RdPu'],
        title='Indian Ocean equatorial 20°C isotherm depth — monthly climatology',
        output_basename='map_eio_t20d',
    )
    plot_map_multimodel(
        cfg,
        eio_pr_monthly,
        cmap_list=['cmo.rain', 'BrBG', 'RdPu'],
        title='Indian Ocean equatorial precipitation — monthly climatology',
        output_basename='map_eio_pr',
    )




if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

