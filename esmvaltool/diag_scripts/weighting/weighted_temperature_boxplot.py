"""Implementation of the weighting scheme.

Boxplot
"""
import logging
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import xarray as xr
from climwip.io_functions import log_provenance, read_model_data
from climwip.core_functions import area_weighted_mean
from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)
from esmvaltool.diag_scripts.weighting.plot_utilities import (
    calculate_percentiles,
    read_metadata,
    read_weights,
)
import scienceplots
plt.style.use(['science','nature'])

logger = logging.getLogger(os.path.basename(__file__))

def visualize_and_save_boxplot(data_df1: pd.DataFrame, data_df2: pd.DataFrame, cfg: dict, ancestors: list):
    #Create figure and axes
    fig, ax = plt.subplots(dpi=300)
    data_list= [data_df1,data_df2]

    for idx in range(2):
        # Create darkgrey box for unweighted projections
        pos_unweighted = idx * 1.5 + 0.75
        pos_weighted = idx * 1.5 + 1
        pos_weighted_era5 = idx * 1.5 + 1.25
        bp_unweighted = ax.boxplot(
            [data_list[idx].loc['Unweighted', ['min', 'Q1', 'mean', 'Q3', 'max']]],
            positions=[pos_unweighted],
            widths=0.2,
            patch_artist=True,
            boxprops=dict(facecolor='darkgrey', color='darkgrey'),
            medianprops=dict(color='black'),
            whiskerprops=dict(color='darkgrey'),
            capprops=dict(color='darkgrey')
        )
        # Create tab:green box for weighted projections
        bp_weighted = ax.boxplot(
            [data_list[idx].loc['Weighted', ['min', 'Q1', 'mean', 'Q3', 'max']]],
            positions=[pos_weighted],
            widths=0.2,
            patch_artist=True,
            boxprops=dict(facecolor='tab:green', color='tab:green'),
            medianprops=dict(color='black'),
            whiskerprops=dict(color='tab:green'),
            capprops=dict(color='tab:green')
        )
        # Create tab:blue box for era5 weighted projections
        bp_weighted_era5 = ax.boxplot(
            [data_list[idx].loc['Weighted_era5', ['min', 'Q1', 'mean', 'Q3', 'max']]],
            positions=[pos_weighted_era5],
            widths=0.2,
            patch_artist=True,
            boxprops=dict(facecolor='tab:blue', color='tab:blue'),
            medianprops=dict(color='black'),
            whiskerprops=dict(color='tab:blue'),
            capprops=dict(color='tab:blue')
        )

    # Customize the plot
    ax.set_xticks([1,2.5])
    caption = cfg['title']
    ax.set_title(caption)
    ax.set_ylabel(cfg['ylabel'])
    ax.set_xticklabels(cfg['xlabel'])
    # Add darkgrey grid to the plot
    ax.grid(True, which='major', color='darkgrey', linestyle='--', linewidth=0.7,axis="y")
    #Switch off minor xticks
    ax.xaxis.set_minor_locator(plt.NullLocator())
    # Add proxy artists for the legend
    legend_elements = [ #Line2D([0], [0], color='white', label=r'Mean, 66$\%$, 90$\%$ boxplot'),
                    Patch(facecolor='darkgrey', edgecolor='darkgrey', label='Unweighted'),
                    Patch(facecolor='tab:green', edgecolor='tab:green', label='Weighted (NCEP/NCAR)'),
                    Patch(facecolor='tab:blue', edgecolor='tab:blue', label='Weighted (ERA5)')
                    ]
    # Add the legend to the plot
    ax.legend(handles=legend_elements, loc='upper left',title=r'Mean, 66$\%$, 90$\%$ boxplot',alignment='left')

    filename_plot = get_plot_filename('boxplot', cfg)
    fig.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(fig)

    filename_data1 = get_diagnostic_filename('boxplot_period_%s'%cfg['xlabel'][0],
                                             cfg,
                                             extension='csv')
    filename_data2 = get_diagnostic_filename('boxplot_period_%s'%cfg['xlabel'][1],
                                             cfg,
                                             extension='csv')
    data_df1.to_csv(filename_data1, index=False)
    data_df2.to_csv(filename_data2, index=False)
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data1, cfg, ancestors)
    log_provenance(caption, filename_data2, cfg, ancestors)

def get_central_and_range(model_data,central_estimate_var,percentiles,weights=None):
    if isinstance(central_estimate_var, (int, float)):
        if weights is None:
            central_estimate = calculate_percentiles(model_data,np.array([central_estimate_var]),)
        else:
            central_estimate = calculate_percentiles(model_data,np.array([central_estimate_var]),weights=weights)
    elif central_estimate_var == 'mean':
        if weights is None:
            central_estimate = model_data.mean('model_ensemble')
        else:
            central_estimate = model_data.weighted(weights).mean('model_ensemble')

    if weights is None:
        uncertainty_range = calculate_percentiles(
            model_data,percentiles,)
    else: 
        uncertainty_range = calculate_percentiles(
            model_data,percentiles,weights=weights)
    return central_estimate,uncertainty_range

def main(cfg):
    """Plot weighted and unweighted boxplot."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    #Weights for period 2041-2060
    weights_path = Path(input_files[0]) / filename
    weights_p1 = read_weights(weights_path)
    weights_path_era5 = Path(input_files[1]) / filename
    weights_era5_p1 = read_weights(weights_path_era5)
    #Weights for period 2081-2100
    weights_path = Path(input_files[2]) / filename
    weights_p2 = read_weights(weights_path)
    weights_path_era5 = Path(input_files[3]) / filename
    weights_era5_p2 = read_weights(weights_path_era5)

    metadata = read_metadata(cfg)
    models1 = metadata['pr_CLIM_future1']
    model_data1, model_data_files = read_model_data(models1)
    models2 = metadata['pr_CLIM_future2']
    model_data2, model_data_files = read_model_data(models2)
    # if a historical period is given calculate the change
    if 'pr_CLIM_reference' in metadata:
        models_hist = metadata['pr_CLIM_reference']
        model_data_hist, model_data_files_hist = read_model_data(models_hist)
        model_data_files += model_data_files_hist
        model_data1 = model_data1 - model_data_hist
        model_data2 = model_data2 - model_data_hist

    model_data1 = area_weighted_mean(model_data1)
    model_data2 = area_weighted_mean(model_data2)

    settings = cfg['settings']
    central_estimate_var = settings.get('central_estimate', 50)
    percentiles = np.array(
        [settings.get('min_bound', 5),settings.get('q1_bound', 17),
         settings.get('q3_bound', 83),settings.get('max_bound', 95)])

    #Central estimate, ranges for period 1 (2041-2060) with associated weights
    central_estimate1,uncertainty_range1 = get_central_and_range(model_data1,central_estimate_var,percentiles,weights=None)
    central_estimate_weighted1,uncertainty_range_weighted1 = get_central_and_range(model_data1,central_estimate_var,
                                                                                   percentiles,weights=weights_p1)
    central_estimate_weighted_era5,uncertainty_range_weighted_era5 = get_central_and_range(model_data1,central_estimate_var,
                                                                                   percentiles,weights=weights_era5_p1)
    
    #Central estimate, ranges for period 1 (2041-2060) with associated weights
    central_estimate2,uncertainty_range2 = get_central_and_range(model_data2,central_estimate_var,percentiles,weights=None)
    central_estimate_weighted2,uncertainty_range_weighted2 = get_central_and_range(model_data2,central_estimate_var,
                                                                                   percentiles,weights=weights_p2)
    central_estimate_weighted2_era5,uncertainty_range_weighted2_era5 = get_central_and_range(model_data2,central_estimate_var,
                                                                                   percentiles,weights=weights_era5_p2)

    data1 = { #data for period1: 2041-2060
    'mean': [central_estimate1.values, central_estimate_weighted1.values, central_estimate_weighted_era5.values],
    'min': [uncertainty_range1[0].values, uncertainty_range_weighted1[0].values, uncertainty_range_weighted_era5[0].values],
    'Q1': [uncertainty_range1[1].values, uncertainty_range_weighted1[1].values, uncertainty_range_weighted_era5[1].values],
    'Q3': [uncertainty_range1[2].values, uncertainty_range_weighted1[2].values, uncertainty_range_weighted_era5[2].values],
    'max': [uncertainty_range1[3].values, uncertainty_range_weighted1[3].values, uncertainty_range_weighted_era5[3].values],
    }
    data2 = { #data for period2: 2081-2100
    'mean': [central_estimate2.values, central_estimate_weighted2.values, central_estimate_weighted2_era5.values],
    'min': [uncertainty_range2[0].values, uncertainty_range_weighted2[0].values, uncertainty_range_weighted2_era5[0].values],
    'Q1': [uncertainty_range2[1].values, uncertainty_range_weighted2[1].values, uncertainty_range_weighted2_era5[1].values],
    'Q3': [uncertainty_range2[2].values, uncertainty_range_weighted2[2].values, uncertainty_range_weighted2_era5[2].values],
    'max': [uncertainty_range2[3].values, uncertainty_range_weighted2[3].values, uncertainty_range_weighted2_era5[3].values],
    }

    data_df1 = pd.DataFrame(data1, index=['Unweighted', 'Weighted', 'Weighted_era5'])
    data_df2 = pd.DataFrame(data2, index=['Unweighted', 'Weighted', 'Weighted_era5'])
    print(data_df1)
    visualize_and_save_boxplot(
        data_df1,data_df2,
        cfg,
        model_data_files,
    )

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
