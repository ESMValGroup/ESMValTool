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
#color-blind friendly color palette
CB_COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
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


def avg_pool_2d_weighted(data, lat_points, kernel=10, stride=10):
    """Average pool 2D data with cos(lat) weighting.
    
    Parameters
    ----------
    data : numpy.ndarray
        Data array with shape (time, lat, lon)
    lat_points : numpy.ndarray
        Latitude values for each row (1D array of length lat)
    kernel : int
        Pooling kernel size
    stride : int
        Pooling stride
    
    Returns
    -------
    numpy.ndarray
        Pooled data with shape (time, lat_out, lon_out)
    numpy.ndarray
        Coarse latitude points (lat_out,)
    """
    t, h, w = data.shape
    h_out = (h - kernel) // stride + 1
    w_out = (w - kernel) // stride + 1
    
    pooled = np.zeros((t, h_out, w_out))
    coarse_lat_points = np.zeros(h_out)
    
    for i in range(h_out):
        lat_start = i * stride
        lat_end = lat_start + kernel
        
        # Get latitudes for this coarse cell
        cell_lats = lat_points[lat_start:lat_end]
        coarse_lat_points[i] = np.mean(cell_lats)
        
        # Compute cos(lat) weights for each row in the cell
        cos_weights = np.cos(np.deg2rad(cell_lats))
        # Normalize weights
        cos_weights = cos_weights / np.sum(cos_weights)
        
        for j in range(w_out):
            lon_start = j * stride
            lon_end = lon_start + kernel
            
            # Extract the block: (time, kernel, kernel)
            block = data[:, lat_start:lat_end, lon_start:lon_end]
            
            # Weight by cos(lat) along latitude dimension
            # First average over longitude, then weighted average over latitude
            lon_mean = np.mean(block, axis=2)  # (time, kernel_lat)
            weighted_mean = np.sum(lon_mean * cos_weights[np.newaxis, :], axis=1)  # (time,)
            pooled[:, i, j] = weighted_mean
    
    return pooled, coarse_lat_points


def calculate_conservation_error(truth, pred_ens, lat_points, pool_factor=10):
    """Calculate conservation error by comparing coarsened truth and predictions.
    
    The conservation error measures how well the downscaled predictions preserve
    the large-scale mean when aggregated back to a coarser resolution.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lat, lon)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lat, lon)
    lat_points : numpy.ndarray
        Latitude values for each row
    pool_factor : int
        Coarsening factor (default: 10)
    
    Returns
    -------
    dict
        Dictionary containing conservation error metrics
    """
    n_ens = pred_ens.shape[0]
    
    # Coarsen the ground truth
    truth_coarse, coarse_lats = avg_pool_2d_weighted(
        truth, lat_points, kernel=pool_factor, stride=pool_factor
    )
    
    # Coarsen ensemble mean prediction
    ens_mean = np.mean(pred_ens, axis=0)
    pred_coarse, _ = avg_pool_2d_weighted(
        ens_mean, lat_points, kernel=pool_factor, stride=pool_factor
    )
    
    # Conservation error: difference between coarsened prediction and coarsened truth
    # This represents how well the downscaling preserves large-scale averages
    conservation_error_map = np.mean(np.abs(pred_coarse - truth_coarse), axis=0)
    
    # Mean conservation error (spatial mean with cos(lat) weighting)
    cos_weights = np.cos(np.deg2rad(coarse_lats))
    cos_weights = cos_weights / np.sum(cos_weights)
    
    # Weighted spatial mean
    conservation_mean = np.sum(
        np.mean(conservation_error_map, axis=1) * cos_weights
    )
    
    # Also compute the mean absolute spatial difference 
    # (difference between domain-averaged coarse pred and truth at each timestep)
    domain_mean_truth = np.mean(truth_coarse, axis=(1, 2))
    domain_mean_pred = np.mean(pred_coarse, axis=(1, 2))
    conservation_domain_mean = np.mean(np.abs(domain_mean_pred - domain_mean_truth))
    
    return {
        "conservation": conservation_error_map,
        "conservation_mean": conservation_mean,
        "conservation_domain_mean": conservation_domain_mean,
        "coarse_lat_points": coarse_lats,
    }


def calculate_spatial_metrics(truth, pred_ens, var_name, lat_points=None, cfg=None):
    """Calculate spatial metrics for evaluation.
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lat, lon)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lat, lon)
    var_name : str
        Variable name
    lat_points : numpy.ndarray, optional
        Latitude values for cos(lat) weighting in conservation error
    cfg : dict, optional
        Configuration dictionary with options:
        - compute_relbias_pr: bool, whether to compute relative bias for pr
        - compute_conservation: list, variables for which to compute conservation error
        - pool_factor: int, coarsening factor for conservation error
    
    Returns
    -------
    dict
        Dictionary containing all calculated metrics
    """
    n_ens, time, lat, lon = pred_ens.shape
    ens_mean = np.mean(pred_ens, axis=0)
    
    # Get configuration options
    if cfg is None:
        cfg = {}
    compute_relbias_pr = cfg.get("compute_relbias_pr", True)
    compute_conservation_vars = cfg.get("compute_conservation", ["pr", "huss"])
    pool_factor = cfg.get("pool_factor", 10)
    
    # Bias
    bias = np.mean(ens_mean - truth, axis=0)
    
    # Relative Bias (optional for pr)
    if var_name in ["pr"] and compute_relbias_pr:
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
    
    result = {
        "bias": bias,
        "relbias": relbias,
        "crps": crps_map,
        "varratio": variance_ratio,
        "bias_mean": bias_mean,
        "relbias_mean": relbias_mean,
        "crps_mean": crps_mean,
        "mae_mean": mae_mean,
        "ssr": ssr,
    }
    
    # Conservation error (only for specified variables)
    if var_name in compute_conservation_vars and lat_points is not None:
        conservation_results = calculate_conservation_error(
            truth, pred_ens, lat_points, pool_factor
        )
        result.update(conservation_results)
    else:
        result["conservation"] = None
        result["conservation_mean"] = None
        result["conservation_domain_mean"] = None
    
    return result


def plot_panel_metrics(metrics_all_methods, var_name, method_names, extent, cfg, lat_points=None):
    """Plot panel of spatial metrics for all methods.
    
    Simplified version showing: Bias, Relative Bias (optional for pr), CRPS, 
    Variance Ratio, and Conservation Error (for pr and huss).
    
    Parameters
    ----------
    metrics_all_methods : list of dict
        List of metrics dictionaries for each method
    var_name : str
        Variable name
    method_names : list of str
        List of method names
    extent : list
        Spatial extent [lon_min, lon_max, lat_min, lat_max]
    cfg : dict
        Configuration dictionary
    lat_points : numpy.ndarray, optional
        Latitude points for conservation error plotting
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.gridspec import GridSpec
    
    n_methods = len(method_names)
    
    # Get variable units
    var_units = UNITS.get(var_name, "")
    
    # Get configuration options
    compute_relbias_pr = cfg.get("compute_relbias_pr", True)
    compute_conservation_vars = cfg.get("compute_conservation", ["pr", "huss"])
    
    # Build metrics to plot based on variable and config
    metrics_to_plot = ["bias"]
    vmins = [-0.1 if var_name == "pr" else -0.5]
    vmaxs = [0.1 if var_name == "pr" else 0.5]
    titles = [f"Bias ({var_units})"]
    cmaps = ["RdBu"]
    show_pct = [False]
    is_coarse = [False]  # Track which metrics are on coarse grid
    
    # Add relative bias for pr if enabled
    if var_name == "pr" and compute_relbias_pr:
        metrics_to_plot.append("relbias")
        vmins.append(-0.5)
        vmaxs.append(0.5)
        titles.append("Relative Bias (%)")
        cmaps.append("BrBG")
        show_pct.append(True)
        is_coarse.append(False)
    
    # Add CRPS
    metrics_to_plot.append("crps")
    vmins.append(0.)
    vmaxs.append(0.1 if var_name == "pr" else 0.5)
    titles.append(f"CRPS ({var_units})")
    cmaps.append("YlGnBu")
    show_pct.append(False)
    is_coarse.append(False)
    
    # Add variance ratio
    # metrics_to_plot.append("varratio")
    # vmins.append(0.5)
    # vmaxs.append(1.5)
    # titles.append("Variance Ratio")
    # cmaps.append("PiYG")
    # show_pct.append(True)
    # is_coarse.append(False)
    
    # Add conservation error for specified variables
    if var_name in compute_conservation_vars:
        metrics_to_plot.append("conservation")
        if var_name == "pr":
            vmins.append(0.)
            vmaxs.append(0.05)
        else:
            vmins.append(0.)
            vmaxs.append(0.5)
        titles.append(f"Conservation Error ({var_units})")
        cmaps.append("cool")
        show_pct.append(False)
        is_coarse.append(True)
    
    # Create figure with cartopy projection
    projection = ccrs.PlateCarree()
    
    fig = plt.figure(figsize=(3.5 * len(metrics_to_plot), 3.2 * n_methods))
    
    gs = GridSpec(n_methods, len(metrics_to_plot), figure=fig, 
                  hspace=0.05, wspace=0.15, 
                  left=0.05, right=0.95, top=0.95, bottom=0.08)
    
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
            
            metric_map = metrics_all_methods[row].get(metric)
            
            if metric_map is None:
                ax.text(0.5, 0.5, "N/A", ha='center', va='center', fontsize=14,
                       transform=ax.transAxes)
                ax.axis('off')
                continue
            
            # Determine extent for this metric (coarse grid has different extent)
            if is_coarse[col] and "coarse_lat_points" in metrics_all_methods[row]:
                # For coarse grid, we still use the same extent but the data is smaller
                plot_extent = extent
            else:
                plot_extent = extent
            
            # Plot the metric map
            im = ax.imshow(metric_map, cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin='lower', aspect='auto',
                transform=projection,
                extent=plot_extent
            )
            
            # Add coastlines and features
            ax.coastlines(resolution='50m', linewidth=1.0, color='black', alpha=0.8)
            ax.add_feature(cfeature.BORDERS, linewidth=0.7, edgecolor='black', alpha=0.5)
            
            # Add gridlines
            gl = ax.gridlines(draw_labels=False, linewidth=0.5, 
                            color='gray', alpha=0.3, linestyle='--')
            
            # Add column title only on first row
            if row == 0:
                ax.set_title(titles[col], fontsize=14, fontweight='bold', pad=10)
            
            # Add method name on the left (outside the axis)
            if col == 0:
                ax.text(
                    -0.15, 0.5, method,
                    transform=ax.transAxes,
                    fontsize=14, fontweight='bold',
                    va='center', ha='right',
                    rotation=90
                )
            
            # Show mean value
            mean_val = np.nanmean(np.abs(metric_map))
            
            # Convert to percentage if needed
            if show_pct[col]:
                display_val = mean_val * 100
            else:
                display_val = mean_val
            
            ax.text(
                0.05, 0.95, f"Mean: {display_val:.3f}",
                ha='left', va='top', transform=ax.transAxes,
                fontsize=10, bbox=dict(
                    facecolor='white', edgecolor='white',
                    boxstyle='round,pad=0.5', alpha=0.9
                ),
                zorder=10
            )
            
            # For conservation error, also show domain mean error
            # if metric == "conservation" and metrics_all_methods[row].get("conservation_domain_mean") is not None:
            #     domain_mean = metrics_all_methods[row]["conservation_domain_mean"]
            #     ax.text(
            #         0.05, 0.82, f"Domain: {domain_mean:.4f}",
            #         ha='left', va='top', transform=ax.transAxes,
            #         fontsize=10, bbox=dict(
            #             facecolor='white', edgecolor='white',
            #             boxstyle='round,pad=0.3', alpha=0.9
            #         ),
            #         zorder=10
            #     )
    
    # Add colorbars at the bottom of each column
    for col in range(len(metrics_to_plot)):
        ax_for_cbar = axes[-1][col]
        
        norm = plt.cm.ScalarMappable(
            cmap=cmaps[col],
            norm=plt.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        )
        norm.set_array([])
        
        pos = ax_for_cbar.get_position()
        cbar_ax = fig.add_axes([pos.x0, pos.y0 - 0.04, pos.width, 0.02])
        
        cbar = fig.colorbar(
            norm, cax=cbar_ax,
            orientation='horizontal'
        )
        cbar.set_ticks(np.linspace(vmins[col], vmaxs[col], 5))
        cbar.ax.tick_params(labelsize=10)
        
    # Save figure
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"spatial_metrics_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    caption = f"Spatial metrics for {var_name} across different ML methods"
    _get_provenance_record(cfg, plot_file, caption, ["metrics", "map"], ["mean", "other"])
    
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
    plt.plot(wavelengths, truth_spectrum, label="Reference", lw=3.5, color='black')
    
    colors = CB_COLORS[:len(method_names)]
    for method, spectrum, color in zip(method_names, pred_spectra, colors):
        plt.plot(wavelengths, spectrum, label=method, linestyle="--", lw=3.5, color=color)
    
    plt.gca().invert_xaxis()
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Length scale (km)", fontsize=22)
    plt.ylabel("Normalized Energy", fontsize=22)
    plt.title(f"Energy Spectrum", fontsize=32, fontweight='bold')
    plt.legend(fontsize=22)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tight_layout()
    
    # Save
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"energy_spectrum_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    caption = f"Energy spectrum comparison for {var_name}"
    _get_provenance_record(cfg, plot_file, caption, ["line"],["spectrum"])

    logger.info("Saved energy spectrum plot: %s", plot_file)


def calculate_ralsd(truth, pred_ens, cfg):
    """Calculate Relative Average Log Spectral Distance (RALSD).

    RALSD measures the difference between power spectral densities in log space.
    Following the formula from literature:
    RALSD = sqrt(mean_i((10 * log10(PSD_truth[i] / PSD_pred[i]))^2))

    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lon, lat)
    cfg : dict
        Configuration dictionary

    Returns
    -------
    float
        RALSD value (lower is better, 0 means perfect match)
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

    # Compute truth spectrum
    truth_spectrum = compute_spectrum(truth)

    # Compute prediction spectrum (average over ensemble)
    ens_spectra = [compute_spectrum(pred_ens[i]) for i in range(pred_ens.shape[0])]
    pred_spectrum = np.mean(ens_spectra, axis=0)

    # Use same wavelength range as in plotting
    n = truth.shape[1]
    truth_spectrum = truth_spectrum[1:n//2]
    pred_spectrum = pred_spectrum[1:n//2]

    # Avoid division by zero and log of zero
    epsilon = 1e-10
    truth_spectrum = np.maximum(truth_spectrum, epsilon)
    pred_spectrum = np.maximum(pred_spectrum, epsilon)

    # Calculate RALSD: sqrt(mean((10 * log10(PSD_truth / PSD_pred))^2))
    log_ratio = 10 * np.log10(truth_spectrum / pred_spectrum)
    ralsd = np.sqrt(np.mean(log_ratio**2))

    return ralsd


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
    plt.plot(centers, np.log10(truth_hist), label="Reference", lw=3.5, color='black')
    
    # Prediction histograms
    colors = CB_COLORS[:len(method_names)]
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
            label=method, linestyle='--', lw=3.5, color=color
        )
    
    # Add units to x-axis label
    xlabel = f"{var_name}"
    if var_units:
        xlabel += f" ({var_units})"
    plt.xlabel(xlabel, fontsize=22)
    plt.ylabel("log₁₀(PDF)", fontsize=22)
    plt.title(f"Probability Distribution", fontsize=32, fontweight='bold')
    plt.legend(fontsize=22)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tight_layout()
    
    # Save
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"log_pdf_{var_name}.png"
    )
    plt.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    caption = f"Log PDF comparison for {var_name}"
    _get_provenance_record(cfg, plot_file, caption, ["probability"], ["pdf"])

    logger.info("Saved log PDF plot: %s", plot_file)


def calculate_log_pdf_distance(truth, pred_ens, bins=200):
    """Calculate log PDF distance between truth and predictions.

    Computes the RMS distance between log-transformed probability density
    functions, similar to RALSD but for histograms.

    Formula: sqrt(mean((log10(PDF_truth) - log10(PDF_pred))^2))

    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lon, lat)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lon, lat)
    bins : int, optional
        Number of histogram bins (default: 200)

    Returns
    -------
    float
        Log PDF distance (lower is better, 0 means perfect match)
    """
    epsilon = 1e-10

    # Flatten and clean truth data
    truth_flat = truth.flatten()
    truth_flat = truth_flat[np.isfinite(truth_flat)]
    hist_range = (np.min(truth_flat), np.max(truth_flat))

    # Compute truth histogram
    truth_hist, edges = np.histogram(truth_flat, bins=bins, range=hist_range, density=True)
    truth_hist = np.maximum(truth_hist, epsilon)
    log_truth = np.log10(truth_hist)

    # Compute prediction histograms (average over ensemble)
    ensemble_pdfs = []
    for i in range(pred_ens.shape[0]):
        pred_flat = pred_ens[i].flatten()
        pred_flat = pred_flat[np.isfinite(pred_flat)]
        pred_hist, _ = np.histogram(pred_flat, bins=edges, density=True)
        ensemble_pdfs.append(pred_hist)

    mean_pred_pdf = np.mean(ensemble_pdfs, axis=0)
    mean_pred_pdf = np.maximum(mean_pred_pdf, epsilon)
    log_pred = np.log10(mean_pred_pdf)

    # Calculate RMS distance in log space
    log_pdf_distance = np.sqrt(np.mean((log_truth - log_pred)**2))

    return log_pdf_distance


def calculate_average_quantile_mae(truth, pred_ens, quantiles=None):
    """Calculate average MAE across all quantiles (area under the quantile MAE curve).

    This provides a single summary metric for quantile-based evaluation.

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
    float
        Average MAE across all quantiles
    """
    q_vals, mae_vals = calculate_quantile_mae(truth, pred_ens, quantiles)

    # Calculate area under the curve using trapezoidal rule, normalized by range
    # This gives the average MAE across the quantile range
    avg_mae = np.trapz(mae_vals, q_vals) / (q_vals[-1] - q_vals[0])

    return avg_mae


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
    # Get configuration for which metrics to include
    compute_conservation_vars = cfg.get("compute_conservation", ["pr", "huss"])
    compute_relbias_pr = cfg.get("compute_relbias_pr", True)
    
    rows = []
    for var in variables_list:
        # Define metrics to include based on variable
        metric_names = ['mae_mean', 'crps_mean', 'ssr']
        
        # Add relbias_mean for pr if enabled
        if var == "pr" and compute_relbias_pr:
            metric_names.insert(1, 'relbias_mean')
        
        # Add conservation_mean for specified variables
        if var in compute_conservation_vars:
            metric_names.append('conservation_mean')
        
        for metric_name in metric_names:
            row = {'variable': var, 'metric': metric_name}
            for method, metrics in zip(method_names, metrics_all_methods[var]):
                val = metrics.get(metric_name)
                row[method] = f"{val:.4f}" if val is not None else "N/A"
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

    Fully vectorized across all spatial locations for performance.

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
    n_ens, time_len, nx, ny = pred_ens.shape

    def _vectorized_lag1_acf(data):
        """Compute lag-1 ACF for all spatial locations at once.

        Parameters
        ----------
        data : numpy.ndarray
            Shape (time, nx, ny)

        Returns        -------
        numpy.ndarray
            Lag-1 ACF at each location, shape (nx, ny). NaN where variance is
            too low.
        """
        x = data[:-1]  # (T-1, nx, ny)
        y = data[1:]   # (T-1, nx, ny)

        # Mean and deviations
        x_mean = np.mean(x, axis=0)  # (nx, ny)
        y_mean = np.mean(y, axis=0)
        dx = x - x_mean[np.newaxis, :, :]
        dy = y - y_mean[np.newaxis, :, :]

        # Covariance and standard deviations
        n = x.shape[0]
        cov_xy = np.sum(dx * dy, axis=0) / n       # (nx, ny)
        std_x = np.sqrt(np.sum(dx**2, axis=0) / n)
        std_y = np.sqrt(np.sum(dy**2, axis=0) / n)
        denom = std_x * std_y

        # Mask low-variance locations
        acf = np.full((nx, ny), np.nan)
        valid = denom > 1e-6
        acf[valid] = cov_xy[valid] / denom[valid]
        return acf

    # Truth ACF: single vectorized call
    truth_acf = _vectorized_lag1_acf(truth)  # (nx, ny)

    # Prediction ACF: one vectorized call per ensemble member, then average
    pred_acf_sum = np.zeros((nx, ny))
    pred_acf_count = np.zeros((nx, ny))
    for n in range(n_ens):
        ens_acf = _vectorized_lag1_acf(pred_ens[n])  # (nx, ny)
        valid = np.isfinite(ens_acf)
        pred_acf_sum[valid] += ens_acf[valid]
        pred_acf_count[valid] += 1

    # Average ensemble ACF
    pred_acf_mean = np.full((nx, ny), np.nan)
    valid_ens = pred_acf_count > 0
    pred_acf_mean[valid_ens] = pred_acf_sum[valid_ens] / pred_acf_count[valid_ens]

    # Signed difference only at locations where both truth and pred are valid
    valid = np.isfinite(truth_acf) & np.isfinite(pred_acf_mean)
    if not np.any(valid):
        return np.nan

    acf_errors = pred_acf_mean[valid] - truth_acf[valid]
    return float(np.mean(acf_errors))


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


def calculate_rank_histogram(truth, pred_ens):
    """Calculate rank histogram for ensemble predictions.
    
    - Rank histogram computed by pooling ranks from all locations and times (for plotting)
    - MCB computed per-location, then averaged across all locations (EnScale method)
    
    For each location and time, compute the rank of the truth within the
    sorted ensemble members. A calibrated ensemble should produce a flat
    rank histogram (uniform distribution).
    
    The MCB (miscalibration) is computed as the sum of absolute deviations
    between the observed rank histogram and a uniform distribution.
    MCB = 0 indicates perfect calibration.
    
    Reference: Schillinger et al., "EnScale" (arXiv)
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lat, lon)
    pred_ens : numpy.ndarray
        Ensemble predictions (n_ens, time, lat, lon)
    
    Returns
    -------
    dict
        Dictionary containing:
        - rank_hist_global: Rank histogram pooled over all locations and times (for plotting)
        - mcb_mean: Mean MCB across all locations (EnScale metric)
        - n_bins: Number of bins in rank histogram
    """
    n_ens, n_time, n_lat, n_lon = pred_ens.shape
    n_bins = n_ens + 1  # Truth can fall in n_ens + 1 positions
    uniform = 1.0 / n_bins
    
    # =========================================================================
    # Compute ranks at all locations - VECTORIZED
    # =========================================================================
    # Sort ensemble members at each location and time
    # pred_ens shape: (n_ens, time, lat, lon)
    pred_ens_sorted = np.sort(pred_ens, axis=0)  # (n_ens, time, lat, lon)
    
    # Compute ranks vectorized using broadcasting
    # For each (time, lat, lon), count how many ensemble members are < truth
    # truth shape: (time, lat, lon), pred_ens_sorted shape: (n_ens, time, lat, lon)
    ranks = np.sum(pred_ens_sorted < truth[np.newaxis, :, :, :], axis=0)  # (time, lat, lon)
    
    # =========================================================================
    # 1. Global pooled rank histogram (for plotting)
    # =========================================================================
    # Pool all ranks from all locations and times into one histogram
    ranks_all = ranks.flatten()  # (time * lat * lon,)
    rank_hist_global = np.bincount(ranks_all, minlength=n_bins)[:n_bins]
    rank_hist_global_normalized = rank_hist_global / len(ranks_all)
    
    # =========================================================================
    # 2. Per-location MCB, then averaged (EnScale metric)
    # =========================================================================
    # Reshape for efficient bincount: flatten lat/lon, keep time separate
    ranks_flat = ranks.reshape(n_time, -1)  # (time, lat*lon)
    n_locations = n_lat * n_lon
    
    # Compute local histograms using vectorized bincount
    local_hists = np.zeros((n_locations, n_bins))
    for loc in range(n_locations):
        local_hists[loc] = np.bincount(ranks_flat[:, loc], minlength=n_bins)[:n_bins]
    
    # Normalize to get frequencies
    local_hists_normalized = local_hists / n_time  # (n_locations, n_bins)
    
    # Compute MCB at each location: sum of |hist - uniform|
    mcb_per_location = np.sum(np.abs(local_hists_normalized - uniform), axis=1)  # (n_locations,)
    
    # Average MCB across all locations (EnScale metric)
    mcb_mean = np.mean(mcb_per_location)
    
    return {
        "rank_hist_global": rank_hist_global_normalized,
        "mcb_mean": mcb_mean,
        "n_bins": n_bins,
    }


def plot_rank_histograms(calibration_results, method_names, var_name, cfg):
    """Plot rank histograms for all methods using global pooled ranks.
    
    Rank histograms pool all ranks from all locations and times into a single
    histogram for visualization. The MCB shown is the per-location average
    (EnScale metric), which is computed separately.
    
    Parameters
    ----------
    calibration_results : list of dict
        List of calibration results for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name
    cfg : dict
        Configuration dictionary
    """
    n_methods = len(method_names)
    n_bins = calibration_results[0]["n_bins"]
    
    # Create figure with 1 row and n_methods columns
    fig, axes = plt.subplots(1, n_methods, figsize=(3.5 * n_methods, 4), 
                              squeeze=False, sharey=True)
    
    # Uniform reference
    uniform = 1.0 / n_bins
    bins = np.arange(n_bins)
    
    colors = CB_COLORS[:n_methods]
    
    for col, (method, results, color) in enumerate(zip(method_names, calibration_results, colors)):
        ax = axes[0, col]
        
        rank_hist = results["rank_hist_global"]
        mcb_mean = results["mcb_mean"]  # Per-location average
        
        # Plot histogram bars
        ax.bar(bins, rank_hist, width=0.8, color=color, alpha=0.7, 
               edgecolor='black', linewidth=0.5)
        
        # Plot uniform reference line
        ax.axhline(y=uniform, color='red', linestyle='--', linewidth=2, 
                  label='Uniform')
        
        # Add MCB annotation (per-location average)
        ax.text(0.95, 0.95, f"MCB: {mcb_mean:.3f}", 
               transform=ax.transAxes, ha='right', va='top',
               fontsize=15, fontweight='bold',
               bbox=dict(facecolor='white', edgecolor='black', 
                        boxstyle='round,pad=0.3', alpha=0.9))
        
        # Labels and titles
        ax.set_title(method, fontsize=15, fontweight='bold')
        if col == 0:
            ax.set_ylabel("Frequency", fontsize=15)
        ax.set_xlabel("Rank", fontsize=15)
        
        # Set x-ticks
        ax.set_xticks(bins[::max(1, n_bins//10)])
        ax.set_xlim(-0.5, n_bins - 0.5)
        
        # Grid
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    
    # Add overall title
    fig.suptitle(f"Rank Histograms\n",
                 fontsize=21, fontweight='bold', y=1.01)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Save figure
    plot_file = os.path.join(cfg["plot_dir"], f"rank_histogram_{var_name}.png")
    plt.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close()
    
    caption = f"Rank histograms for {var_name} showing ensemble calibration"
    _get_provenance_record(cfg, plot_file, caption, ["histogram"], ["other"])
    
    logger.info("Saved rank histogram plot: %s", plot_file)


def create_calibration_table(calibration_all_vars, method_names, variables_list, cfg):
    """Create summary table of calibration metrics (MCB).
    
    Following EnScale: reports MCB averaged across all locations.
    
    Parameters
    ----------
    calibration_all_vars : dict
        Dictionary mapping var_name -> list of calibration results per method
    method_names : list of str
        List of method names
    variables_list : list of str
        List of variable names
    cfg : dict
        Configuration dictionary
    
    Returns
    -------
    pd.DataFrame
        Summary dataframe
    """
    rows = []
    
    for var in variables_list:
        if var not in calibration_all_vars:
            continue
        
        results_list = calibration_all_vars[var]
        
        # MCB mean (averaged over all locations) - main EnScale metric
        row = {'variable': var, 'metric': 'MCB'}
        for method, results in zip(method_names, results_list):
            row[method] = f"{results['mcb_mean']:.4f}"
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Save as CSV
    table_file = os.path.join(cfg["work_dir"], "calibration_mcb_table.csv")
    df.to_csv(table_file, index=False)
    
    # Save as formatted text with description
    txt_file = os.path.join(cfg["work_dir"], "calibration_mcb_table.txt")
    with open(txt_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CALIBRATION METRICS (MCB - Miscalibration)\n")
        f.write("=" * 80 + "\n\n")
        f.write("Methodology:\n")
        f.write("- Rank histograms computed at each grid point separately\n")
        f.write("- MCB = sum of absolute deviations from uniform rank histogram\n")
        f.write("- Final MCB = average across all spatial locations\n")
        f.write("- Visualization: global pooled histogram (all ranks pooled)\n\n")
        f.write("Interpretation:\n")
        f.write("- MCB = 0 indicates perfect calibration\n")
        f.write("- Higher MCB indicates worse calibration\n")
        f.write("- U-shaped histogram → overconfident (too little spread)\n")
        f.write("- Dome-shaped histogram → underconfident (too much spread)\n")
        f.write("-" * 80 + "\n\n")
        f.write(df.to_string(index=False))
    
    logger.info("Saved calibration table: %s", table_file)
    
    return df


def run_calibration_analysis(truth, preds_methods, method_names, var_name, extent, cfg):
    """Run calibration analysis for a single variable.
    
    - Computes global pooled rank histograms (for visualization)
    - Computes MCB per-location, then averages (EnScale calibration metric)
    
    Parameters
    ----------
    truth : numpy.ndarray
        Ground truth data (time, lat, lon)
    preds_methods : list of numpy.ndarray
        List of prediction ensembles for each method
    method_names : list of str
        List of method names
    var_name : str
        Variable name
    extent : list
        Spatial extent for plotting (unused, kept for API compatibility)
    cfg : dict
        Configuration dictionary
    
    Returns
    -------
    list of dict
        List of calibration results for each method
    """
    calibration_results = []
    
    for method, pred_ens in zip(method_names, preds_methods):
        logger.info("Computing rank histogram and MCB for %s (%s)", method, var_name)
        results = calculate_rank_histogram(truth, pred_ens)
        calibration_results.append(results)
        logger.info("MCB for %s (%s): %.4f", method, var_name, results['mcb_mean'])
    
    # Plot rank histograms (global pooled, with per-location MCB annotated)
    plot_rank_histograms(calibration_results, method_names, var_name, cfg)
    
    return calibration_results


def create_comprehensive_summary_table(all_results, method_names, cfg):
    """Create comprehensive summary table with all metrics across analysis types.

    This table aggregates metrics from:
    - Spatial metrics (bias, CRPS, MAE, SSR, conservation)
    - Energy spectrum (RALSD)
    - Log PDF distance
    - Temporal consistency (ACF error)
    - Quantile MAE (average)

    Parameters
    ----------
    all_results : dict
        Dictionary containing all collected metrics organized by analysis type
    method_names : list of str
        List of method names
    cfg : dict
        Configuration dictionary

    Returns
    -------
    pd.DataFrame
        Summary dataframe with all metrics
    """
    # Get configuration
    compute_conservation_vars = cfg.get("compute_conservation", ["pr", "huss"])
    compute_relbias_pr = cfg.get("compute_relbias_pr", True)
    
    rows = []

    # Collect all variables across all analysis types
    all_variables = set()
    for analysis_type in all_results:
        if all_results[analysis_type]:
            all_variables.update(all_results[analysis_type].keys())

    for var_name in sorted(all_variables):
        # Spatial metrics
        if 'spatial_metrics' in all_results and var_name in all_results['spatial_metrics']:
            metrics_list = all_results['spatial_metrics'][var_name]
            
            # Core metrics
            core_metrics = ['bias_mean', 'crps_mean', 'mae_mean', 'ssr']
            
            # Add relbias for pr if enabled
            if var_name == "pr" and compute_relbias_pr:
                core_metrics.insert(1, 'relbias_mean')
            
            # Add conservation for specified variables
            if var_name in compute_conservation_vars:
                core_metrics.append('conservation_mean')
            
            for metric_name in core_metrics:
                row = {'variable': var_name, 'metric': metric_name, 'analysis_type': 'spatial'}
                for i, method in enumerate(method_names):
                    if i < len(metrics_list) and metric_name in metrics_list[i]:
                        val = metrics_list[i][metric_name]
                        row[method] = f"{val:.4f}" if val is not None else "N/A"
                    else:
                        row[method] = "N/A"
                rows.append(row)

        # Energy spectrum (RALSD)
        if 'energy_spectrum' in all_results and var_name in all_results['energy_spectrum']:
            row = {'variable': var_name, 'metric': 'RALSD', 'analysis_type': 'spectrum'}
            ralsd_dict = all_results['energy_spectrum'][var_name]
            for method in method_names:
                if method in ralsd_dict:
                    row[method] = f"{ralsd_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)

        # Log PDF distance
        if 'log_density' in all_results and var_name in all_results['log_density']:
            row = {'variable': var_name, 'metric': 'log_pdf_distance', 'analysis_type': 'distribution'}
            lpd_dict = all_results['log_density'][var_name]
            for method in method_names:
                if method in lpd_dict:
                    row[method] = f"{lpd_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)

        # Temporal structure (ACF error)
        if 'temporal_structure' in all_results and var_name in all_results['temporal_structure']:
            row = {'variable': var_name, 'metric': 'acf_error', 'analysis_type': 'temporal'}
            acf_dict = all_results['temporal_structure'][var_name]
            for method in method_names:
                if method in acf_dict:
                    row[method] = f"{acf_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)

        # Quantile MAE average
        if 'quantile_mae' in all_results and var_name in all_results['quantile_mae']:
            row = {'variable': var_name, 'metric': 'avg_quantile_mae', 'analysis_type': 'quantile'}
            qmae_dict = all_results['quantile_mae'][var_name]
            for method in method_names:
                if method in qmae_dict:
                    row[method] = f"{qmae_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)

        # Calibration (MCB)
        if 'calibration' in all_results and var_name in all_results['calibration']:
            calib_list = all_results['calibration'][var_name]
            
            # MCB mean (averaged over locations) - main EnScale metric
            row = {'variable': var_name, 'metric': 'MCB', 'analysis_type': 'calibration'}
            for i, method in enumerate(method_names):
                if i < len(calib_list):
                    row[method] = f"{calib_list[i]['mcb_mean']:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)

    if not rows:
        logger.warning("No metrics collected for comprehensive summary table")
        return None

    df = pd.DataFrame(rows)

    # Reorder columns
    cols = ['variable', 'analysis_type', 'metric'] + method_names
    df = df[[c for c in cols if c in df.columns]]

    # Save as CSV
    table_file = os.path.join(cfg["work_dir"], "comprehensive_metrics_summary.csv")
    df.to_csv(table_file, index=False)

    # Save as formatted text
    txt_file = os.path.join(cfg["work_dir"], "comprehensive_metrics_summary.txt")
    with open(txt_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("COMPREHENSIVE ML DOWNSCALING EVALUATION METRICS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")
        f.write("-" * 80 + "\n")
        f.write("Metric descriptions:\n")
        f.write("-" * 80 + "\n")
        f.write("  bias_mean        : Mean absolute bias (lower is better)\n")
        f.write("  relbias_mean     : Mean absolute relative bias for pr (lower is better)\n")
        f.write("  crps_mean        : Continuous Ranked Probability Score (lower is better)\n")
        f.write("  mae_mean         : Mean Absolute Error (lower is better)\n")
        f.write("  ssr              : Spread-Skill Ratio (closer to 1 is better)\n")
        f.write("  conservation_mean: Conservation Error (lower is better)\n")
        f.write("  RALSD            : Relative Avg Log Spectral Distance (lower is better)\n")
        f.write("  log_pdf_distance : Log PDF Distance (lower is better)\n")
        f.write("  acf_error        : Lag-1 ACF Error (closer to 0 is better)\n")
        f.write("  avg_quantile_mae : Average Quantile MAE (lower is better)\n")
        f.write("  MCB              : Miscalibration - avg over locations (lower is better, 0=perfect)\n")

    logger.info("Saved comprehensive metrics summary: %s", table_file)

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
        norm = Normalize(vmin=vmin, vmax=vmax)
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
                  hspace=0.1, wspace=0.3,
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
        ax.coastlines(resolution='50m', linewidth=1.0, color='black', alpha=0.8)
        ax.add_feature(cfeature.BORDERS, linewidth=0.7, edgecolor='black', alpha=0.5)
        # ax.gridlines(draw_labels=False, linewidth=0.5, 
        #             color='gray', alpha=0.3, linestyle='--')
        
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
    colors = CB_COLORS[:len(method_names)]
    
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
    plt.savefig(plot_file, dpi=200, bbox_inches='tight')
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
    
    # Log new configuration options
    compute_relbias_pr = cfg.get("compute_relbias_pr", True)
    compute_conservation_vars = cfg.get("compute_conservation", ["pr", "huss"])
    pool_factor = cfg.get("pool_factor", 10)
    logger.info("Compute relative bias for pr: %s", compute_relbias_pr)
    logger.info("Compute conservation error for: %s", compute_conservation_vars)
    logger.info("Pool factor for conservation: %d", pool_factor)
    
    # Group input data
    input_data = cfg["input_data"].values()
    grouped_data = group_metadata(input_data, "short_name", sort="dataset")

    # Initialize storage for all metrics (for comprehensive summary table)
    all_results = {
        'spatial_metrics': {},
        'energy_spectrum': {},
        'log_density': {},
        'temporal_structure': {},
        'quantile_mae': {},
        'calibration': {}
    }
    temporal_structure_results = {}
    calibration_all_vars = {}
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
        
        # Get extent of regional domain and lat/lon coordinates
        lat_points = None
        try:
            # Try to get latitude coordinate
            lat_coord = truth_cube.coord('latitude')
            lat_points = lat_coord.points
            
            # Handle 2D coordinates
            if lat_points.ndim == 2:
                # Take mean along longitude for 1D lat array
                lat_points = np.mean(lat_points, axis=1)
            
            # Try to get longitude coordinate
            lon_coord = truth_cube.coord('longitude')
            lon_points = lon_coord.points
            
            # Calculate extent
            if lat_points.ndim == 1 and lon_points.ndim == 1:
                extent = [lon_points.min(), lon_points.max(), 
                        lat_points.min(), lat_points.max()]
            elif lon_points.ndim == 2:
                extent = [lon_points.min(), lon_points.max(), 
                        lat_points.min(), lat_points.max()]
            else:
                extent = [-180, 180, -90, 90]
                logger.warning("Could not determine extent from coordinates, using global extent")
            
            logger.info(f"Data extent: lon=[{extent[0]:.2f}, {extent[1]:.2f}], "
                        f"lat=[{extent[2]:.2f}, {extent[3]:.2f}]")
            
        except iris.exceptions.CoordinateNotFoundError:
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
                metrics = calculate_spatial_metrics(
                    truth, pred_ens, var_name, 
                    lat_points=lat_points, cfg=cfg
                )
                metrics_list.append(metrics)
            all_metrics[var_name] = metrics_list
            all_results['spatial_metrics'][var_name] = metrics_list
            # Plot spatial metrics
            plot_panel_metrics(
                metrics_list, var_name, method_names_loaded, extent, cfg,
                lat_points=lat_points
            )

        elif analysis_type == "energy_spectrum":
            # Plot energy spectrum
            plot_energy_spectrum(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
            # Calculate RALSD for each method
            ralsd_dict = {}
            for method, pred_ens in zip(method_names_loaded, preds_methods):
                ralsd = calculate_ralsd(truth, pred_ens, cfg)
                ralsd_dict[method] = ralsd
                logger.info("RALSD for %s (%s): %.4f", method, var_name, ralsd)
            all_results['energy_spectrum'][var_name] = ralsd_dict

        elif analysis_type == "log_density":
            # Plot log PDF
            plot_log_pdf(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
            # Calculate log PDF distance for each method
            lpd_dict = {}
            for method, pred_ens in zip(method_names_loaded, preds_methods):
                lpd = calculate_log_pdf_distance(truth, pred_ens)
                lpd_dict[method] = lpd
                logger.info("Log PDF distance for %s (%s): %.4f", method, var_name, lpd)
            all_results['log_density'][var_name] = lpd_dict
            
        elif analysis_type == "temporal_structure":
            temporal_structure_results[var_name] = run_temporal_structure_analysis(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
            all_results['temporal_structure'][var_name] = temporal_structure_results[var_name]
            
        elif analysis_type == "create_animation":
            # Create animation
            create_animation(
                truth, truth_dates, preds_methods, 
                method_names_loaded, var_name, extent, cfg
            )
        
        elif analysis_type == "calibration":
            # Run calibration analysis (rank histogram + MCB)
            calibration_results = run_calibration_analysis(
                truth, preds_methods, method_names_loaded, var_name, extent, cfg
            )
            all_results['calibration'][var_name] = calibration_results
    
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

            # Plot quantile MAE and calculate average
            if preds_wbgt:
                plot_quantile_mae(truth_wbgt, preds_wbgt, method_names_loaded, "wbgt", cfg)
                # Calculate average quantile MAE for each method
                qmae_dict = {}
                for method, pred_ens in zip(method_names_loaded, preds_wbgt):
                    avg_qmae = calculate_average_quantile_mae(truth_wbgt, pred_ens)
                    qmae_dict[method] = avg_qmae
                    logger.info("Avg Quantile MAE for %s (wbgt): %.4f", method, avg_qmae)
                all_results['quantile_mae']['wbgt'] = qmae_dict

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

            # Plot quantile MAE and calculate average
            if preds_wind:
                plot_quantile_mae(truth_wind, preds_wind, method_names_loaded, "sfcWind", cfg)
                # Calculate average quantile MAE for each method
                qmae_dict = {}
                for method, pred_ens in zip(method_names_loaded, preds_wind):
                    avg_qmae = calculate_average_quantile_mae(truth_wind, pred_ens)
                    qmae_dict[method] = avg_qmae
                    logger.info("Avg Quantile MAE for %s (sfcWind): %.4f", method, avg_qmae)
                all_results['quantile_mae']['sfcWind'] = qmae_dict
        
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
                plot_quantile_mae(truth_data, preds_methods, method_names_loaded, var_name, cfg, range = [0.95,1.0])
                # Calculate average quantile MAE for each method (using same range)
                qmae_dict = {}
                quantiles = np.linspace(0.8, 1.0, 101)
                for method, pred_ens in zip(method_names_loaded, preds_methods):
                    avg_qmae = calculate_average_quantile_mae(truth_data, pred_ens, quantiles)
                    qmae_dict[method] = avg_qmae
                    logger.info("Avg Quantile MAE for %s (pr): %.4f", method, avg_qmae)
                all_results['quantile_mae']['pr'] = qmae_dict

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

    # Create energy spectrum summary table with RALSD
    if analysis_type == "energy_spectrum" and all_results['energy_spectrum']:
        rows = []
        for var_name, ralsd_dict in all_results['energy_spectrum'].items():
            row = {'variable': var_name, 'metric': 'RALSD'}
            for method in ml_methods:
                if method in ralsd_dict:
                    row[method] = f"{ralsd_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)
        df = pd.DataFrame(rows)
        table_file = os.path.join(cfg["work_dir"], "energy_spectrum_ralsd_table.csv")
        df.to_csv(table_file, index=False)
        logger.info("Saved energy spectrum RALSD table: %s", table_file)

    # Create log density summary table
    if analysis_type == "log_density" and all_results['log_density']:
        rows = []
        for var_name, lpd_dict in all_results['log_density'].items():
            row = {'variable': var_name, 'metric': 'log_pdf_distance'}
            for method in ml_methods:
                if method in lpd_dict:
                    row[method] = f"{lpd_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)
        df = pd.DataFrame(rows)
        table_file = os.path.join(cfg["work_dir"], "log_density_distance_table.csv")
        df.to_csv(table_file, index=False)
        logger.info("Saved log density distance table: %s", table_file)

    # Create quantile MAE summary table
    if analysis_type == "quantile_MAE" and all_results['quantile_mae']:
        rows = []
        for var_name, qmae_dict in all_results['quantile_mae'].items():
            row = {'variable': var_name, 'metric': 'avg_quantile_mae'}
            for method in ml_methods:
                if method in qmae_dict:
                    row[method] = f"{qmae_dict[method]:.4f}"
                else:
                    row[method] = "N/A"
            rows.append(row)
        df = pd.DataFrame(rows)
        table_file = os.path.join(cfg["work_dir"], "quantile_mae_table.csv")
        df.to_csv(table_file, index=False)
        logger.info("Saved quantile MAE table: %s", table_file)

    # Create calibration (MCB) summary table
    if analysis_type == "calibration" and all_results['calibration']:
        create_calibration_table(
            all_results['calibration'],
            ml_methods,
            list(all_results['calibration'].keys()),
            cfg
        )

    # Create comprehensive summary table if any metrics were collected
    has_any_results = any(all_results[key] for key in all_results)
    if has_any_results:
        create_comprehensive_summary_table(all_results, ml_methods, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)