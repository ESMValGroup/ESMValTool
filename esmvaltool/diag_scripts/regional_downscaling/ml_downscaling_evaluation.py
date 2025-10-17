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
            cube = cube[:, :, :, 0]  # Select first height level
        cubes.append(cube.data)
    
    # Stack along ensemble dimension
    ensemble_data = np.stack(cubes, axis=0)
    return ensemble_data


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


def plot_panel_metrics(metrics_all_methods, var_name, method_names, cfg):
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
    n_methods = len(method_names)
    
    # Define metrics to plot based on variable
    if var_name in ["pr"]:
        metrics_to_plot = ["bias", "relbias", "crps", "corr", "varratio"]
        vmins = [-1, -0.4, 0., 0.6, 0.5]
        vmaxs = [1, 0.4, 0.5, 1., 1.5]
        titles = ["Bias", "Relative Bias", "CRPS", "Correlation", "Variance Ratio"]
        cmaps = ["RdBu", "BrBG", "YlGnBu", "gist_ncar", "PiYG"]
    else:
        metrics_to_plot = ["bias", "crps", "corr", "varratio"]
        vmins = [-1, 0, 0.5, 0.5]
        vmaxs = [1, 2, 1, 1.5]
        titles = ["Bias", "CRPS", "Correlation", "Variance Ratio"]
        cmaps = ["BrBG", "YlGnBu", "gist_ncar", "PiYG"]
    
    fig, axes = plt.subplots(
        n_methods, len(metrics_to_plot),
        figsize=(3.5 * len(metrics_to_plot), 3 * n_methods),
        squeeze=False, constrained_layout=True
    )
    
    for row, method in enumerate(method_names):
        for col, metric in enumerate(metrics_to_plot):
            ax = axes[row, col]
            
            metric_map = metrics_all_methods[row][metric]
            
            if metric_map is None:
                ax.text(0.5, 0.5, "N/A", ha='center', va='center', fontsize=14)
                ax.axis('off')
                continue
            
            im = ax.imshow(
                metric_map, cmap=cmaps[col],
                vmin=vmins[col], vmax=vmaxs[col],
                origin='lower'
            )
            
            if row == 0:
                ax.set_title(f"{titles[col]}", fontsize=11, fontweight='bold')
            
            # Add method name on left
            if col == 0:
                ax.set_ylabel(method, fontsize=11, fontweight='bold')
            
            # Show mean value
            mean_val = np.nanmean(np.abs(metric_map))
            unit = ""
            if metric in ["relbias", "varratio"]:
                unit = "%"
                mean_val *= 100
            
            ax.text(
                0.05, 0.95, f"Mean: {mean_val:.2f}{unit}",
                ha='left', va='top', transform=ax.transAxes,
                fontsize=9, bbox=dict(
                    facecolor='white', edgecolor='black',
                    boxstyle='round,pad=0.5', alpha=0.8
                )
            )
            
            ax.axis("off")
    
    # Add colorbars
    for col in range(len(metrics_to_plot)):
        ax_for_cbar = axes[-1, col]
        norm = plt.cm.ScalarMappable(
            cmap=cmaps[col],
            norm=plt.Normalize(vmin=vmins[col], vmax=vmaxs[col])
        )
        norm.set_array([])
        cbar = fig.colorbar(
            norm, ax=ax_for_cbar,
            orientation='horizontal', pad=0.05, fraction=0.05
        )
        cbar.set_label(titles[col], fontsize=10)
    
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
    dx_km : float
        Grid spacing in km
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
    
    # Truth histogram
    truth_flat = truth.flatten()
    truth_flat = truth_flat[np.isfinite(truth_flat)]
    hist_range = (np.percentile(truth_flat, 0.1), np.percentile(truth_flat, 99.9))
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
    
    plt.xlabel(f"{var_name} value", fontsize=12)
    plt.ylabel("log10(PDF)", fontsize=12)
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
    
    # Process each variable
    for var_name in grouped_data:
        logger.info("Processing variable: %s", var_name)
        
        var_datasets = grouped_data[var_name]
        
        # Separate reference and ML methods
        reference_datasets = [
            d for d in var_datasets
            if d["dataset"] == reference_name
        ]
        
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
        if truth_cube.ndim == 4:  # Remove height dimension if present
            truth_cube = truth_cube[:, 0, :, :]
        truth = truth_cube.data
        
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
            
            # Plot spatial metrics
            plot_panel_metrics(metrics_list, var_name, method_names_loaded, cfg)
            
            # Create summary table (store for later aggregation)
            if not hasattr(cfg, '_metrics_storage'):
                cfg['_metrics_storage'] = {}
            cfg['_metrics_storage'][var_name] = (metrics_list, method_names_loaded)
        
        elif analysis_type == "energy_spectrum":
            plot_energy_spectrum(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
        
        elif analysis_type == "log_density":
            plot_log_pdf(
                truth, preds_methods, method_names_loaded, var_name, cfg
            )
    
    # Create summary table if spatial_metrics was run
    if analysis_type == "spatial_metrics" and hasattr(cfg, '_metrics_storage'):
        # Aggregate all variables
        all_metrics = {}
        all_methods = None
        for var_name, (metrics_list, methods) in cfg['_metrics_storage'].items():
            all_metrics[var_name] = metrics_list
            if all_methods is None:
                all_methods = methods
        
        if all_metrics:
            create_metrics_table(
                all_metrics,
                all_methods,
                list(all_metrics.keys()),
                cfg
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)