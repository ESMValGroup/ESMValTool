"""ML-based Downscaling Evaluation Diagnostic.

This diagnostic evaluates ML-based downscaling methods (UNet, CorrDiff, etc.)
against high-resolution reference data using various metrics and analyses.
"""

import logging
import os
from pathlib import Path

import iris
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr
from properscoring import crps_ensemble

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))
UNITS= {} #store units of variables

def _get_provenance_record(cfg, plot_file, caption, plot_types = ["map"], statistics= ["other"]):
    """Create a provenance record describing the diagnostic data and plot."""
    ancestor_files = []
    for dataset in cfg["input_data"].values():
        ancestor_files.append(dataset["filename"])
    
    record = {
        "caption": caption,
        "statistics": statistics,
        "domains": ["reg"],
        "plot_types": plot_types,
        "authors": ["debeire_kevin"],
        "references": [],
        "plot_file": plot_file,
        "ancestors": ancestor_files,
    }
    
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, record)


def load_ensemble_data(datasets_list):
    """Load ensemble members for a given method/dataset.
    
    Parameters
    ----------
    datasets_list : list of dict
        List of dataset metadata dictionaries
    
    Returns
    -------
    numpy.ndarray
        Array with shape (n_ensemble, time, lon, lat)
    """
    cubes = []
    for dataset in datasets_list:
        cube = iris.load_cube(dataset["filename"])
        # Remove height dimension if present (for uas, vas, tas, huss)
        if cube.ndim == 4:  # (time, lon, lat, height)
            cube = cube[:, :, :, 0] # Select first height level
            cube.transpose([0,2,1])  
        elif cube.ndim == 3:  # (time, lon, lat, height)
            cube.transpose([0,2,1])
        cubes.append(cube.data)
    
    # Stack along ensemble dimension
    ensemble_data = np.stack(cubes, axis=0)
    return ensemble_data

def update_variable_units(datasets: list[dict]) -> bool:
    """    
    Parameters
    ----------
    datasets : list[dict]
        Datasets
    Returns
    -------
    str
        Units list per variable
    """

    for dataset in datasets:
        units = dataset["units"]
        short_name = dataset["short_name"]
        if short_name in UNITS:
            pass
        else:
            UNITS.setdefault(short_name, units)

    return UNITS

def calculate_spatial_metrics(truth, pred_ens, var_name):
    """Calculate spatial metrics for evaluation.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lon, lat)
    var_name : str
        Variable name
    
    Returns
    -------
    dict
        Dictionary containing all calculated metrics
    """
    n_ens, time, lon, lat  = pred_ens.shape
    ens_mean = np.mean(pred_ens, axis=0)
    
    # Bias and Relative Bias
    bias = np.mean(ens_mean - truth, axis=0)
    if var_name in ["pr"]:
        relbias = bias / (np.mean(truth, axis=0) + 1e-6)
    else:
        relbias = None
    
    # CRPS
    crps_map = np.zeros((lat, lon))
    for i in range(lat):
        for j in range(lon):
            crps_map[i, j] = np.mean(
                crps_ensemble(truth[:, i, j], pred_ens[:, :, i, j].T)
            )
    
    # Correlation
    corr_map = np.zeros((lat, lon))
    for i in range(lat):
        for j in range(lon):
            if np.std(ens_mean[:, i, j]) > 0 and np.std(truth[:, i, j]) > 0:
                corr_map[i, j], _ = pearsonr(ens_mean[:, i, j], truth[:, i, j])
            else:
                corr_map[i, j] = np.nan
    
    # Variance ratio
    var_truth = np.var(truth, axis=0)
    variance_ratio_ens = []
    for n in range(n_ens):
        var_pred = np.var(pred_ens[n], axis=0)
        var_ratio = var_pred / (var_truth + 1e-8)
        variance_ratio_ens.append(var_ratio)
    variance_ratio = np.mean(variance_ratio_ens, axis=0)
    
    # Spread-Skill Ratio (SSR)
    ensemble_mean = np.mean(pred_ens, axis=0)
    ensemble_std = np.std(pred_ens, axis=0)
    mean_spread = np.mean(ensemble_std)
    rmse = np.sqrt(np.mean((ensemble_mean - truth)**2))
    ssr = mean_spread / rmse if rmse > 0 else np.nan
    
    # MAE
    mae = np.mean(np.abs(ens_mean - truth), axis=0)
    mae_mean = np.mean(mae)
    
    # Summary statistics
    bias_mean = np.mean(np.abs(bias))
    relbias_mean = np.mean(np.abs(relbias)) if relbias is not None else None
    crps_mean = np.mean(crps_map)
    corr_mean = np.nanmean(corr_map)
    
    return {
        "bias": bias,
        "relbias": relbias,
        "crps": crps_map,
        "corr": corr_map,
        "varratio": variance_ratio,
        "bias_mean": bias_mean,
        "relbias_mean": relbias_mean,
        "crps_mean": crps_mean,
        "corr_mean": corr_mean,
        "mae_mean": mae_mean,
        "ssr": ssr,
    }

def plot_panel_metrics(metrics_all_methods, var_name, method_names, extent, cfg):
    """Plot panel of spatial metrics for all methods.
    
    Parameters
    ----------
    metrics_all_methods : list of dict
        List of metrics dictionaries for each method
    var_name : str
        Variable name
    method_names : list of str
        List of method names
    cfg : dict
        Configuration dictionary
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    
    n_methods = len(method_names)
    
    # Get variable units
    var_units = UNITS[var_name]
    
    # Define metrics to plot based on variable
    if var_name in ["pr"]:
        metrics_to_plot = ["bias", "relbias", "crps", "corr", "varratio"]
        vmins = [-0.1, -0.4, 0., 0.6, 0.5]
        vmaxs = [0.1, 0.4, 0.2, 1., 1.5]
        titles = [
            f"Bias ({var_units})", 
            "Relative Bias (%)", 
            f"CRPS ({var_units})", 
            "Correlation", 
            "Variance Ratio (%)"
        ]
        cmaps = ["RdBu", "BrBG", "YlGnBu", "gist_ncar", "PiYG"]
        # Whether to show percentage in the mean text
        show_pct = [False, True, False, False, True]
    else:
        metrics_to_plot = ["bias", "crps", "corr", "varratio"]
        vmins = [-1, 0, 0.5, 0.5]
        vmaxs = [1, 2, 1, 1.5]
        titles = [
            f"Bias ({var_units})", 
            f"CRPS ({var_units})", 
            "Correlation", 
            "Variance Ratio (%)"
        ]
        cmaps = ["BrBG", "YlGnBu", "gist_ncar", "PiYG"]
        show_pct = [False, False, False, True]
    
    # Create figure with cartopy projection
    projection = ccrs.PlateCarree()
    
    fig = plt.figure(figsize=(3.5 * len(metrics_to_plot), 3.2 * n_methods))
    
    # Create grid for subplots
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(n_methods, len(metrics_to_plot), figure=fig, 
                  hspace=0.25, wspace=0.15, 
                  left=0.08, right=0.95, top=0.95, bottom=0.08)
    
    axes = []
    for row in range(n_methods):
        row_axes = []
        for col in range(len(metrics_to_plot)):
            ax = fig.add_subplot(gs[row, col], projection=projection)
            row_axes.append(ax)
        axes.append(row_axes)
    
    for row, method in enumerate(method_names):
        for col, metric in enumerate(metrics_to_plot):
            ax = axes[row][col]
            
            metric_map = metrics_all_methods[row][metric]
            
            if metric_map is None:
                ax.text(0.5, 0.5, "N/A", ha='center', va='center', fontsize=14,
                       transform=ax.transAxes)
                ax.axis('off')
                continue
            
            # Plot the metric map
            im = ax.imshow(metric_map, cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin='lower', aspect='auto',
                transform=projection,
                extent=extent  # Adjust based on your domain
            )
            
            # Add coastlines and features
            ax.coastlines(resolution='50m', linewidth=0.8, color='black', alpha=0.6)
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black', alpha=0.3)
            
            # Set extent to your region (adjust these values based on your data)
            # Example for a specific region - you may need to extract this from data
            # ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=projection)
            
            # Add gridlines
            gl = ax.gridlines(draw_labels=False, linewidth=0.5, 
                            color='gray', alpha=0.3, linestyle='--')
            
            # Add column title only on first row
            if row == 0:
                ax.set_title(titles[col], fontsize=11, fontweight='bold', pad=10)
            
            # Add method name on the left (outside the axis)
            if col == 0:
                # Position the text to the left of the axis
                ax.text(
                    -0.15, 0.5, method,
                    transform=ax.transAxes,
                    fontsize=12, fontweight='bold',
                    va='center', ha='right',
                    rotation=90
                )
            
            # Show mean value
            mean_val = np.nanmean(np.abs(metric_map))
            
            # Convert to percentage if needed (for display in label, not text box)
            if show_pct[col]:
                display_val = mean_val * 100
            else:
                display_val = mean_val
            
            ax.text(
                0.05, 0.95, f"Mean: {display_val:.3f}",
                ha='left', va='top', transform=ax.transAxes,
                fontsize=9, bbox=dict(
                    facecolor='white', edgecolor='black',
                    boxstyle='round,pad=0.5', alpha=0.9
                ),
                zorder=10
            )
    
    # Add colorbars at the bottom of each column
    for col in range(len(metrics_to_plot)):
        # Get the position of the last axis in this column
        ax_for_cbar = axes[-1][col]
        
        # Create colorbar
        norm = plt.cm.ScalarMappable(
            cmap=cmaps[col],
            norm=plt.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        )
        norm.set_array([])
        
        # Get axis position
        pos = ax_for_cbar.get_position()
        cbar_ax = fig.add_axes([pos.x0, pos.y0 - 0.04, pos.width, 0.02])
        
        cbar = fig.colorbar(
            norm, cax=cbar_ax,
            orientation='horizontal'
        )
        cbar.ax.tick_params(labelsize=9)
    
    # Save figure
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"spatial_metrics_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    caption = f"Spatial metrics for {var_name} across different ML methods"
    _get_provenance_record(cfg, plot_file, caption, ["metrics", "map"],["mean", "other", "corr"])
    
    logger.info("Saved spatial metrics plot: %s", plot_file)


def plot_energy_spectrum(truth, preds_methods, method_names, var_name, cfg):
    """Plot energy spectrum comparison.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    preds_methods : list of numpy.ndarray
        List of prediction ensembles for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name
    cfg : dict
        Configuration dictionary
    """
    def radial_average(psd2d):
        y, x = np.indices(psd2d.shape)
        center = np.array([
            (x.max() - x.min()) / 2.0,
            (y.max() - y.min()) / 2.0
        ])
        r = np.sqrt((x - center[0])**2 + (y - center[1])**2).astype(int)
        tbin = np.bincount(r.ravel(), psd2d.ravel())
        nr = np.bincount(r.ravel())
        return tbin / np.maximum(nr, 1)
    
    def compute_spectrum(data):
        spectra = []
        for t in range(data.shape[0]):
            field = data[t] - np.mean(data[t])
            fft2 = np.fft.fft2(field)
            psd2d = np.abs(fft2)**2
            psd2d = np.fft.fftshift(psd2d)
            spectra.append(radial_average(psd2d))
        return np.mean(spectra, axis=0)
    
    dx_km = cfg.get("spatial_scale", 6.3)
    # Compute truth spectrum
    truth_spectrum = compute_spectrum(truth)
    truth_spectrum /= truth_spectrum.max()
    
    # Compute prediction spectra
    pred_spectra = []
    for pred_ens in preds_methods:
        ens_spectra = [
            compute_spectrum(pred_ens[i])
            for i in range(pred_ens.shape[0])
        ]
        mean_spectrum = np.mean(ens_spectra, axis=0)
        mean_spectrum /= mean_spectrum.max()
        pred_spectra.append(mean_spectrum)
    
    # Compute wavelengths
    n = truth.shape[1]
    freqs = np.fft.fftfreq(n, d=dx_km)[:n//2]
    wavelengths = 1 / freqs[1:]
    
    truth_spectrum = truth_spectrum[1:n//2]
    pred_spectra = [s[1:n//2] for s in pred_spectra]
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, truth_spectrum, label="Reference", lw=2.5, color='black')
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(method_names)))
    for method, spectrum, color in zip(method_names, pred_spectra, colors):
        plt.plot(wavelengths, spectrum, label=method, linestyle="--", lw=2, color=color)
    
    plt.gca().invert_xaxis()
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Length scale (km)", fontsize=12)
    plt.ylabel("Normalized Energy", fontsize=12)
    plt.title(f"Energy Spectrum: {var_name}", fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    
    # Save
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"energy_spectrum_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    caption = f"Energy spectrum comparison for {var_name}"
    _get_provenance_record(cfg, plot_file, caption, ["line"],["spectrum"])
    
    logger.info("Saved energy spectrum plot: %s", plot_file)


def plot_log_pdf(truth, preds_methods, method_names, var_name, cfg):
    """Plot log PDF comparison.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    preds_methods : list of numpy.ndarray
        List of prediction ensembles for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name
    cfg : dict
        Configuration dictionary
    """
    plt.figure(figsize=(10, 6))
    bins = 200
    
    # Get variable units
    var_units = UNITS[var_name]
    
    # Truth histogram
    truth_flat = truth.flatten()
    truth_flat = truth_flat[np.isfinite(truth_flat)]
    hist_range = (np.min(truth_flat), np.max(truth_flat))
    truth_hist, edges = np.histogram(
        truth_flat, bins=bins, range=hist_range, density=True
    )
    centers = (edges[:-1] + edges[1:]) / 2
    epsilon = 1e-10
    truth_hist = np.maximum(truth_hist, epsilon)
    plt.plot(centers, np.log10(truth_hist), label="Reference", lw=2.5, color='black')
    
    # Prediction histograms
    colors = plt.cm.tab10(np.linspace(0, 1, len(method_names)))
    for method, pred_ens, color in zip(method_names, preds_methods, colors):
        ensemble_pdfs = []
        for i in range(pred_ens.shape[0]):
            pred_flat = pred_ens[i].flatten()
            pred_flat = pred_flat[np.isfinite(pred_flat)]
            pred_hist, _ = np.histogram(pred_flat, bins=edges, density=True)
            ensemble_pdfs.append(pred_hist)
        mean_pred_pdf = np.mean(ensemble_pdfs, axis=0)
        mean_pred_pdf = np.maximum(mean_pred_pdf, epsilon)
        plt.plot(
            centers, np.log10(mean_pred_pdf),
            label=method, linestyle='--', lw=2, color=color
        )
    
    # Add units to x-axis label
    xlabel = f"{var_name}"
    if var_units:
        xlabel += f" ({var_units})"
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel("log₁₀(PDF)", fontsize=12)
    plt.title(f"Distribution Comparison: {var_name}", fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    
    # Save
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"log_pdf_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    caption = f"Log PDF comparison for {var_name}"
    _get_provenance_record(cfg, plot_file, caption, ["probability"], ["pdf"])
    
    logger.info("Saved log PDF plot: %s", plot_file)


def create_metrics_table(metrics_all_methods, method_names, variables_list, cfg):
    """Create summary table of metrics.
    
    Parameters
    ----------
    metrics_all_methods : dict
        Dictionary of metrics for all methods and variables
    method_names : list of str
        List of method names
    variables_list : list of str
        List of variable names
    cfg : dict
        Configuration dictionary
    """
    rows = []
    for var in variables_list:
        for metric_name in ['mae_mean', 'crps_mean', 'ssr', 'corr_mean']:
            row = {'variable': var, 'metric': metric_name}
            for method, metrics in zip(method_names, metrics_all_methods[var]):
                row[method] = f"{metrics[metric_name]:.4f}"
            rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Save as CSV
    table_file = os.path.join(cfg["work_dir"], "summary_metrics_table.csv")
    df.to_csv(table_file, index=False)
    
    # Save as formatted text
    txt_file = os.path.join(cfg["work_dir"], "summary_metrics_table.txt")
    with open(txt_file, 'w') as f:
        f.write(df.to_string(index=False))
    
    logger.info("Saved metrics table: %s", table_file)
    
    return df

def calculate_acf_error(truth, pred_ens, var_name):
    """
    Calculate the error in lag-1 autocorrelation (ACF) for each variable.
    For ensemble methods, average ACF at the end.

    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lon, lat)
    var_name : str
        Variable name

    Returns
    -------
    float
        Average signed difference of samples' ACF minus target ACF over all locations
    """
    n_ens, time, lon, lat = pred_ens.shape
    acf_errors = []

    # For each location, calculate lag-1 ACF for truth and each ensemble member
    for i in range(lon):
        for j in range(lat):
            # Truth ACF
            truth_ts = truth[:, i, j]
            if np.std(truth_ts) < 1e-6:
                continue  # skip if no variance
            truth_acf = np.corrcoef(truth_ts[:-1], truth_ts[1:])[0, 1]

            # Ensemble ACF
            ens_acfs = []
            for n in range(n_ens):
                pred_ts = pred_ens[n, :, i, j]
                if np.std(pred_ts) < 1e-6:
                    continue
                pred_acf = np.corrcoef(pred_ts[:-1], pred_ts[1:])[0, 1]
                ens_acfs.append(pred_acf)
            ens_acf = np.mean(ens_acfs)

            # Signed difference
            acf_errors.append(ens_acf - truth_acf)

    return np.mean(acf_errors) if acf_errors else np.nan

def run_temporal_structure_analysis(truth, preds_methods, method_names, var_name, cfg):
    """
    Run temporal structure analysis (lag-1 ACF error) for a single variable.
    """
    acf_errors = []
    for pred_ens in preds_methods:
        acf_error = calculate_acf_error(truth, pred_ens, var_name)
        acf_errors.append(acf_error)

    return dict(zip(method_names, acf_errors))

def create_temporal_structure_table(results, method_names, variables_list, cfg):
    """
    Create a table of lag-1 ACF errors for each method and variable.
    """
    rows = []
    for var in variables_list:
        if var not in results:
            continue
        row = {'variable': var}
        for method in method_names:
            if method in results[var]:
                row[method] = f"{results[var][method]:.4f}"
            else:
                row[method] = "N/A"
        rows.append(row)

    df = pd.DataFrame(rows)

    # Save as CSV
    table_file = os.path.join(cfg["work_dir"], "temporal_structure_table.csv")
    df.to_csv(table_file, index=False)

    # Save as formatted text
    txt_file = os.path.join(cfg["work_dir"], "temporal_structure_table.txt")
    with open(txt_file, 'w') as f:
        f.write(df.to_string(index=False))

    logger.info("Saved temporal structure table: %s", table_file)

    return df

def create_animation(truth, truth_dates, preds_methods, method_names, var_name, extent, cfg):
    """Create animation comparing reference and downscaled fields over time.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    truth_dates : list
        List of datetime objects for each time step
    preds_methods : list of numpy.ndarray
        List of prediction ensembles for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name
    extent : list
        Spatial extent [lon_min, lon_max, lat_min, lat_max]
    cfg : dict
        Configuration dictionary
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.colors import LogNorm, Normalize

    # Get animation configuration
    start_time = cfg.get("animation_start_time", 0)
    end_time = cfg.get("animation_end_time", start_time+10)
    frame_duration = cfg.get("animation_frame_duration", 1)  # seconds per frame
    
    # Ensure valid time range
    end_time = min(end_time, truth.shape[0])
    if start_time >= end_time:
        logger.warning("Invalid time range for animation: start=%d, end=%d", 
                      start_time, end_time)
        return
    
    # Get variable units
    var_units = UNITS.get(var_name, "")
    
    # Determine colormap and value range
    if var_name in ["pr"]:
        cmap = "YlGnBu"
        vmin = 0
        vmax = np.ceil(np.percentile(truth[start_time:end_time], 99.5))
        norm = Normalize(vmin=vmin, vmax=vmax) #LogNorm(vmin=vmin, vmax=vmax) #for log scale
    elif var_name in ["tas", "tasmax", "tasmin"]:
        cmap = "RdYlBu_r"
        vmin = np.min(truth[start_time:end_time])
        vmax = np.max(truth[start_time:end_time])
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        cmap = "viridis"
        vmin = np.min(truth[start_time:end_time])
        vmax = np.max(truth[start_time:end_time])
        norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Setup figure
    projection = ccrs.PlateCarree()
    n_cols = len(method_names) + 1  # +1 for reference
    fig = plt.figure(figsize=(4 * n_cols, 5.5))
    
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(1, n_cols, figure=fig, 
                  hspace=0.1, wspace=0.25,
                  left=0.05, right=0.95, top=0.85, bottom=0.15)
    
    # Create axes
    axes = []
    for col in range(n_cols):
        ax = fig.add_subplot(gs[0, col], projection=projection)
        axes.append(ax)
    
    # Title for variable and units (will be updated with date)
    title_text = f"{var_name}"
    if var_units:
        title_text += f" ({var_units})"
    fig_title = fig.suptitle(title_text, fontsize=14, fontweight='bold', y=0.99)
    
    # Initialize plots
    images = []
    for col, ax in enumerate(axes):
        # Setup map features
        ax.coastlines(resolution='50m', linewidth=0.8, color='black', alpha=0.6)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black', alpha=0.3)
        ax.gridlines(draw_labels=False, linewidth=0.5, 
                    color='gray', alpha=0.3, linestyle='--')
        
        # Set extent
        if extent:
            ax.set_extent(extent, crs=projection)
        
        # Add method labels
        if col == 0:
            ax.set_title("Reference", fontsize=12, fontweight='bold', pad=10)
        else:
            ax.set_title(method_names[col-1], fontsize=12, fontweight='bold', pad=10)

        # Initialize with first frame
        im = ax.imshow(truth[start_time], cmap=cmap, norm=norm,
            origin='lower', aspect='auto',
            transform=projection,
            extent=extent if extent else None
        )
        images.append(im)
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.15, 0.08, 0.7, 0.03])
    cbar = fig.colorbar(images[0], cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=10)
    if var_units:
        cbar.set_label(var_units, fontsize=11)
    
    # Animation update function
    def update(frame_idx):
        time_idx = start_time + frame_idx
        
        # Update date in title
        if truth_dates and time_idx < len(truth_dates):
            date_str = truth_dates[time_idx].strftime("%Y-%m-%d %H:%M")
            title_with_date = f"{title_text}\n{date_str}"
        else:
            title_with_date = f"{title_text}\nTime step: {time_idx}"
        fig_title.set_text(title_with_date)
        
        # Update reference
        images[0].set_data(truth[time_idx])
        
        # Update predictions (use first ensemble member for each method)
        for col in range(1, n_cols):
            method_idx = col - 1
            if method_idx < len(preds_methods):
                pred_data = preds_methods[method_idx][0, time_idx]  # First ensemble member
                images[col].set_data(pred_data)
        
        return images + [fig_title]
    
    # Create animation
    n_frames = end_time - start_time
    anim = FuncAnimation(
        fig, update, frames=n_frames,
        interval=frame_duration * 1000,  # Convert to milliseconds
        blit=False
    )
    
    # Save as GIF
    output_file = os.path.join(
        cfg["plot_dir"],
        f"animation_{var_name}.gif"
    )
    
    writer = PillowWriter(fps=1.0/frame_duration)
    anim.save(output_file, writer=writer, dpi=100)
    plt.close()
    
    caption = f"Animation of {var_name} from time {start_time} to {end_time}"
    _get_provenance_record(cfg, output_file, caption, ["map"], ["mean"])
    
    logger.info("Saved animation: %s", output_file)

def calculate_wbgt(tas, huss, ps):
    """Calculate Wet-Bulb Globe Temperature (WBGT).
    
    Parameters
    ----------
    tas : numpy.ndarray
        Air temperature in Kelvin (time, lon, lat)
    huss : numpy.ndarray
        Specific humidity in kg/kg (time, lon, lat)
    ps : numpy.ndarray
        Surface pressure in Pa (time, lon, lat)
    
    Returns
    -------
    numpy.ndarray
        WBGT in degrees Celsius (time, lon, lat)
    """
    # Convert huss from g/kg to g/kg
    huss = huss/1000
    T = tas
    # Calculate vapor pressure in hPa
    epsilon = 0.622
    e = (ps * huss) / (epsilon + (1 - epsilon) * huss)
    
    # Calculate WBGT
    wbgt = 0.567 * T + 0.393 * e + 3.94
    
    return wbgt


def calculate_wind_speed(uas, vas):
    """Calculate wind speed from u and v components.
    
    Parameters
    ----------
    uas : numpy.ndarray
        Eastward wind component (time, lon, lat)
    vas : numpy.ndarray
        Northward wind component (time, lon, lat)
    
    Returns
    -------
    numpy.ndarray
        Wind speed (time, lon, lat)
    """
    return np.sqrt(uas**2 + vas**2)


def calculate_quantile_mae(truth, pred_ens, quantiles=None):
    """Calculate MAE of quantiles between truth and predictions.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lon, lat)
    quantiles : numpy.ndarray, optional
        Array of quantile values to compute. Default is np.linspace(0, 1, 101)
    
    Returns
    -------
    tuple
        (quantiles, mae_values) where mae_values is the MAE at each quantile
    """
    if quantiles is None:
        quantiles = np.linspace(0, 1, 101)
    
    # Flatten spatial dimensions
    truth_flat = truth.reshape(truth.shape[0], -1)  # (time, space)
    
    # Pool all ensemble members
    # pred_ens has shape (n_ens, time, lon, lat)
    n_ens = pred_ens.shape[0]
    pred_flat = pred_ens.reshape(n_ens, pred_ens.shape[1], -1)  # (n_ens, time, space)
    
    # Compute quantiles for truth
    truth_quantiles = np.quantile(truth_flat, quantiles, axis=1)  # (n_quantiles, space)
    truth_quantiles_mean = np.mean(truth_quantiles, axis=1)  # (n_quantiles,)
    
    # Compute quantiles for each ensemble member and average
    pred_quantiles_list = []
    for i in range(n_ens):
        pred_q = np.quantile(pred_flat[i], quantiles, axis=1)  # (n_quantiles, space)
        pred_quantiles_list.append(np.mean(pred_q, axis=1))  # (n_quantiles,)
    
    pred_quantiles_mean = np.mean(pred_quantiles_list, axis=0)  # (n_quantiles,)
    
    # Calculate MAE
    mae_values = np.abs(pred_quantiles_mean - truth_quantiles_mean)
    
    return quantiles, mae_values


def plot_quantile_mae(truth_data, preds_methods, method_names, var_name, cfg, range= [0,1]):
    """Plot quantile MAE comparison.
    
    Parameters
    ----------
    truth_data : numpy.ndarray
        Ground truth data (time, lon, lat)
    preds_methods : list of numpy.ndarray
        List of prediction ensembles for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name (e.g., 'wbgt', 'sfcWind', 'pr')
    cfg : dict
        Configuration dictionary
    """
    plt.figure(figsize=(10, 6))
    
    # Get variable units
    var_units = UNITS.get(var_name, "")
    
    # Define quantile range
    quantiles = np.linspace(range[0], range[1], 101)
    
    # Colors for different methods
    colors = plt.cm.tab10(np.linspace(0, 1, len(method_names)))
    
    # Calculate and plot quantile MAE for each method
    for method, pred_ens, color in zip(method_names, preds_methods, colors):
        q_vals, mae_vals = calculate_quantile_mae(truth_data, pred_ens, quantiles)
        plt.plot(q_vals, mae_vals, label=method, lw=2, color=color)
    
    # Formatting
    plt.xlabel("Quantile", fontsize=12)
    ylabel = "MAE"
    if var_units:
        ylabel += f" ({var_units})"
    plt.ylabel(ylabel, fontsize=12)
    
    # Create title with proper variable name
    var_display_name = var_name
    if var_name == "wbgt":
        var_display_name = "Summer WBGT"
    elif var_name == "sfcWind":
        var_display_name = "Fall Wind Speed"
    elif var_name == "pr":
        var_display_name = "Winter Precipitation"
    
    plt.title(f"Quantile MAE: {var_display_name}", fontsize=14, fontweight='bold')
    plt.legend(fontsize=11, loc='best')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.xlim(range[0], range[1])
    plt.tight_layout()
    
    # Save
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"quantile_mae_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    caption = f"Quantile MAE comparison for {var_display_name}"
    _get_provenance_record(cfg, plot_file, caption, ["line"], ["other"])
    
    logger.info("Saved quantile MAE plot: %s", plot_file)

def main(cfg):
    """Run ML downscaling evaluation diagnostic."""
    logger.setLevel(cfg["log_level"].upper())
    
    # Get configuration
    analysis_type = cfg.get("analysis_type", "spatial_metrics")
    ml_methods = cfg.get("ml_methods", [])
    reference_name = cfg.get("reference", "HIGHRES-REF")
    
    logger.info("Analysis type: %s", analysis_type)
    logger.info("ML methods: %s", ml_methods)
    logger.info("Reference: %s", reference_name)
    
    # Group input data
    input_data = cfg["input_data"].values()
    grouped_data = group_metadata(input_data, "short_name", sort="dataset")
    temporal_structure_results = {}
    all_metrics = {}
    # Process each variable
    for var_name in grouped_data:
        logger.info("Processing variable: %s", var_name)
        
        var_datasets = grouped_data[var_name]
        
        # Separate reference and ML methods
        reference_datasets = [
            d for d in var_datasets
            if d["dataset"] == reference_name
        ]
        UNITS = update_variable_units(reference_datasets)
    
        ml_datasets_by_method = {}
        for method in ml_methods:
            ml_datasets_by_method[method] = [
                d for d in var_datasets
                if d["dataset"] == method
            ]
        
        if not reference_datasets:
            logger.warning("No reference data found for %s", var_name)
            continue
        
        # Load reference data (assume single ensemble/version for reference)
        truth_cube = iris.load_cube(reference_datasets[0]["filename"])
        print(truth_cube)
        # Extracting dates
        truth_dates = None
        try:
            time_coord = truth_cube.coord('time')
            truth_dates = [time_coord.units.num2date(t) for t in time_coord.points]
        except (iris.exceptions.CoordinateNotFoundError, AttributeError):
            logger.warning("Could not extract dates from time coordinate")

        if truth_cube.ndim == 4:  # Remove height dimension if present
            truth_cube = truth_cube[:, :, :, 0]
            truth_cube.transpose([0,2,1])
        elif truth_cube.ndim == 3:
            truth_cube.transpose([0,2,1])

        truth = truth_cube.data
        
        #Get extent of regional domain
        # Extract lat/lon coordinates from the cube
        try:
            # Try to get latitude coordinate
            lat_coord = truth_cube.coord('latitude')
            lat_points = lat_coord.points
            
            # Try to get longitude coordinate
            lon_coord = truth_cube.coord('longitude')
            lon_points = lon_coord.points
            
            # Calculate extent
            # Handle both 1D and 2D coordinate arrays
            if lat_points.ndim == 1 and lon_points.ndim == 1:
                # 1D coordinates (regular grid)
                extent = [lon_points.min(), lon_points.max(), 
                        lat_points.min(), lat_points.max()]
            elif lat_points.ndim == 2 and lon_points.ndim == 2:
                # 2D coordinates (curvilinear grid)
                extent = [lon_points.min(), lon_points.max(), 
                        lat_points.min(), lat_points.max()]
            else:
                # Fallback to global extent
                extent = [-180, 180, -90, 90]
                logger.warning("Could not determine extent from coordinates, using global extent")
            
            logger.info(f"Data extent: lon=[{extent[0]:.2f}, {extent[1]:.2f}], "
                        f"lat=[{extent[2]:.2f}, {extent[3]:.2f}]")
            
        except iris.exceptions.CoordinateNotFoundError:
            # If coordinates not found, use global extent
            extent = [-180, 180, -90, 90]
            logger.warning("Latitude/longitude coordinates not found, using global extent")

        # Store extent in cfg for use in plotting functions
        if not hasattr(cfg, 'plot_extent'):
            cfg['plot_extent'] = extent

        # Load ML method ensembles
        preds_methods = []
        method_names_loaded = []
        for method in ml_methods:
            if method not in ml_datasets_by_method or not ml_datasets_by_method[method]:
                logger.warning("No data found for method %s", method)
                continue
            
            pred_ens = load_ensemble_data(ml_datasets_by_method[method])
            preds_methods.append(pred_ens)
            method_names_loaded.append(method)
        
        if not preds_methods:
            logger.warning("No ML method data loaded for %s", var_name)
            continue
        
        # Perform analysis based on type
        if analysis_type == "spatial_metrics":
            # Calculate metrics for all methods
            metrics_list = []
            for pred_ens in preds_methods:
                metrics = calculate_spatial_metrics(truth, pred_ens, var_name)
                metrics_list.append(metrics)
            all_metrics[var_name] = metrics_list
            # Plot spatial metrics
            plot_panel_metrics(metrics_list, var_name, method_names_loaded, extent, cfg)
            
            # # Create summary table (store for later aggregation)
            # if not hasattr(cfg, '_metrics_storage'):
            #     cfg['_metrics_storage'] = {}
            # cfg['_metrics_storage'][var_name] = (metrics_list, method_names_loaded)
        
        elif analysis_type == "energy_spectrum":
            plot_energy_spectrum(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
        
        elif analysis_type == "log_density":
            plot_log_pdf(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
                # Perform analysis based on type
        elif analysis_type == "temporal_structure":
            temporal_structure_results[var_name] = run_temporal_structure_analysis(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )        
        elif analysis_type == "create_animation":
            # Create animation
            create_animation(
                truth, truth_dates, preds_methods, 
                method_names_loaded, var_name, extent, cfg
            )
    
    # Handle quantile_MAE analysis separately (may need derived variables)
    if analysis_type == "quantile_MAE":
        # Check if we need to compute derived variables
        compute_wbgt = cfg.get("compute_wbgt", False)
        compute_wind = cfg.get("compute_wind", False)
        compute_pr = cfg.get("compute_pr", False)

        if compute_wbgt and all(v in grouped_data for v in ["tas", "huss", "ps"]):
            logger.info("Computing WBGT for reference and all methods")
            
            # Load reference data for all required variables
            ref_tas = iris.load_cube([d["filename"] for d in grouped_data["tas"] 
                                     if d["dataset"] == reference_name][0])
            ref_huss = iris.load_cube([d["filename"] for d in grouped_data["huss"] 
                                      if d["dataset"] == reference_name][0])
            ref_ps = iris.load_cube([d["filename"] for d in grouped_data["ps"] 
                                    if d["dataset"] == reference_name][0])
            
            # Remove height dimension if present
            if ref_tas.ndim == 4:
                ref_tas = ref_tas[:, :, :, 0]
                ref_tas.transpose([0,2,1])
            elif ref_tas.ndim == 3:
                ref_tas.transpose([0,2,1])
            if ref_huss.ndim == 4:
                ref_huss = ref_huss[:, :, :, 0]
                ref_huss.transpose([0,2,1])
            elif ref_huss.ndim == 3:
                ref_huss.transpose([0,2,1])
            if ref_ps.ndim == 4:
                ref_ps = ref_ps[:, :, :, 0]
                ref_ps.transpose([0,2,1])
            elif ref_ps.ndim == 3:
                ref_ps.transpose([0,2,1])
            
            # Calculate WBGT for reference
            truth_wbgt = calculate_wbgt(ref_tas.data, ref_huss.data, ref_ps.data)
            
            # Calculate WBGT for each method
            preds_wbgt = []
            method_names_loaded = []
            for method in ml_methods:
                # Load tas, huss, ps for this method
                method_tas_datasets = [d for d in grouped_data["tas"] if d["dataset"] == method]
                method_huss_datasets = [d for d in grouped_data["huss"] if d["dataset"] == method]
                method_ps_datasets = [d for d in grouped_data["ps"] if d["dataset"] == method]
                
                if not (method_tas_datasets and method_huss_datasets and method_ps_datasets):
                    logger.warning(f"Missing required variables for WBGT calculation for method {method}")
                    continue
                
                # Load ensembles
                tas_ens = load_ensemble_data(method_tas_datasets)
                huss_ens = load_ensemble_data(method_huss_datasets)
                ps_ens = load_ensemble_data(method_ps_datasets)
                
                # Calculate WBGT for each ensemble member
                wbgt_ens = np.zeros_like(tas_ens)
                for i in range(tas_ens.shape[0]):
                    wbgt_ens[i] = calculate_wbgt(tas_ens[i], huss_ens[i], ps_ens[i])
                
                preds_wbgt.append(wbgt_ens)
                method_names_loaded.append(method)
            
            # Set units for WBGT
            UNITS["wbgt"] = "°C"
            
            # Plot quantile MAE
            if preds_wbgt:
                plot_quantile_mae(truth_wbgt, preds_wbgt, method_names_loaded, "wbgt", cfg)

        if compute_wind and all(v in grouped_data for v in ["uas", "vas"]):
            logger.info("Computing wind speed for reference and all methods")
            
            # Load reference data
            ref_uas = iris.load_cube([d["filename"] for d in grouped_data["uas"] 
                                     if d["dataset"] == reference_name][0])
            ref_vas = iris.load_cube([d["filename"] for d in grouped_data["vas"] 
                                     if d["dataset"] == reference_name][0])
            
            # Remove height dimension if present
            if ref_uas.ndim == 4:
                ref_uas = ref_uas[:, :, :, 0]
                ref_uas.transpose([0,2,1])
            elif ref_uas.ndim == 3:
                ref_uas.transpose([0,2,1])
            if ref_vas.ndim == 4:
                ref_vas = ref_vas[:, :, :, 0]
                ref_vas.transpose([0,2,1])
            elif ref_vas.ndim == 3:
                ref_vas.transpose([0,2,1])
            
            # Calculate wind speed for reference
            truth_wind = calculate_wind_speed(ref_uas.data, ref_vas.data)
            
            # Calculate wind speed for each method
            preds_wind = []
            method_names_loaded = []
            for method in ml_methods:
                method_uas_datasets = [d for d in grouped_data["uas"] if d["dataset"] == method]
                method_vas_datasets = [d for d in grouped_data["vas"] if d["dataset"] == method]
                
                if not (method_uas_datasets and method_vas_datasets):
                    logger.warning(f"Missing required variables for wind speed calculation for method {method}")
                    continue
                
                # Load ensembles
                uas_ens = load_ensemble_data(method_uas_datasets)
                vas_ens = load_ensemble_data(method_vas_datasets)
                
                # Calculate wind speed for each ensemble member
                wind_ens = np.zeros_like(uas_ens)
                for i in range(uas_ens.shape[0]):
                    wind_ens[i] = calculate_wind_speed(uas_ens[i], vas_ens[i])
                
                preds_wind.append(wind_ens)
                method_names_loaded.append(method)
            
            # Set units for wind speed
            UNITS["sfcWind"] = "m s-1"
            
            # Plot quantile MAE
            if preds_wind:
                plot_quantile_mae(truth_wind, preds_wind, method_names_loaded, "sfcWind", cfg)
        
        # Process regular variables (e.g., precipitation)
        if compute_pr and all(v in grouped_data for v in ["pr"]):
            var_name= "pr"
            
            logger.info(f"Processing quantile MAE for variable: {var_name}")
            
            var_datasets = grouped_data[var_name]
            reference_datasets = [d for d in var_datasets if d["dataset"] == reference_name]
            
            if not reference_datasets:
                logger.warning(f"Missing Reference data for variable {var_name}")
            
            # Load reference
            truth_cube = iris.load_cube(reference_datasets[0]["filename"])
            if truth_cube.ndim == 4:
                truth_cube = truth_cube[:, :, :, 0]
                truth_cube.transpose([0,2,1])
            elif truth_cube.ndim == 3:
                truth_cube.transpose([0,2,1])
            truth_data = truth_cube.data
            
            # Load ML methods
            preds_methods = []
            method_names_loaded = []
            for method in ml_methods:
                method_datasets = [d for d in var_datasets if d["dataset"] == method]
                if not method_datasets:
                    continue
                pred_ens = load_ensemble_data(method_datasets)
                preds_methods.append(pred_ens)
                method_names_loaded.append(method)
            
            if preds_methods:
                plot_quantile_mae(truth_data, preds_methods, method_names_loaded, var_name, cfg, range = [0.8,1.0])

    # Create summary table if spatial_metrics was run
    if analysis_type == "spatial_metrics":
        if all_metrics:
            create_metrics_table(
                all_metrics,
                ml_methods,
                list(all_metrics.keys()),
                cfg
            )
        # Create summary table if temporal_structure was run
    elif analysis_type == "temporal_structure":
        create_temporal_structure_table(
            temporal_structure_results,
            ml_methods,
            list(temporal_structure_results.keys()),
            cfg
        )

if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)               