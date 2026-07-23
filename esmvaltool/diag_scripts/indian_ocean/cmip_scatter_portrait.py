import logging
from pathlib import Path
from pprint import pformat
import os
import iris # type: ignore
from iris.cube import Cube  # type: ignore
from iris.coords import AuxCoord, DimCoord  # type: ignore
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from basic_functions import (compute_cube_skew, get_provenance_record, iso_depth_3d, 
                             load_data, 
                             compute_cube_diff, 
                             load_and_update_dict,
                             get_prefix,
                             to_set)

from esmvaltool.diag_scripts.shared import (   # type: ignore
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
    Datasets,
    Variables,
)
import esmvaltool.diag_scripts.shared as e   # type: ignore

logger = logging.getLogger(Path(__file__).stem)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

OBS_CANDIDATES = {
    'sst': ['HadISST'],
    'wind': ['NCEP', 'ERA5'],
    'ocean': ['EN4', 'SODA', 'ORAS5'],
    'dmi': ['HadISST'],
}


def as_scalar(cube, method='mean'):
    """Return a scalar float from a cube using a selected reduction method."""
    if method == 'mean':
        if cube.ndim == 0:
            return float(cube.data)
        return float(cube.collapsed(cube.dim_coords, iris.analysis.MEAN).data)
    if method == 'std':
        return float(cube.collapsed('time', iris.analysis.STD_DEV).data)
    if method == 'raw':
        return float(cube.data)
    raise ValueError(f"Unsupported scalar reduction method: {method}")


def get_obs_dataset_name(data_dict, candidates):
    """Return the first observational dataset name available in the dictionary."""
    if not data_dict:
        return None
    for candidate in candidates:
        if candidate in data_dict:
            return candidate
    return None


def save_scalar_bias(cfg, metric_key, dataset, value, ancestors):
    """Store a scalar bias value as a cube using save_data."""
    cube = Cube(np.array(value), long_name=f"{metric_key} bias", units='1')
    save_name = f"bias_{metric_key}_{dataset}"
    provenance_record = get_provenance_record(save_name, sorted(list(ancestors)))
    save_data(save_name, provenance_record, cfg, cube)


def build_bias_table(cfg, diagnostics):
    """Build and save model-by-diagnostic bias matrix and plot a heatmap."""
    metric_names = [
        'CEIO zonal wind bias',
        'EEIO mean SST bias',
        'WEIO mean SST bias',
        'WEIO-EEIO SST gradient bias',
        'DMI standard deviation bias',
        'DMI skewness bias',
        'equatorial D20 tilt bias',
        'SCTR D20 bias',
    ]

    model_names = sorted(
        diagnostics.keys(),
        key=lambda model: (
            diagnostics[model].get('CEIO zonal wind bias', np.inf),
            diagnostics[model].get('EEIO mean SST bias', np.inf),
            model,
        ),
    )
    if not model_names:
        logger.warning("No model diagnostics available for bias portrait plot.")
        return

    matrix = np.full((len(model_names), len(metric_names)), np.nan)
    ancestors = set()

    for i, model in enumerate(model_names):
        model_vals = diagnostics[model]
        ancestors.update(model_vals.get('ancestors', set()))
        for j, metric in enumerate(metric_names):
            matrix[i, j] = model_vals.get(metric, np.nan)

    # Standardize each diagnostic across models (z-score) so all columns share a
    # comparable, unitless scale for plotting and table export.
    standardized_matrix = np.full_like(matrix, np.nan, dtype=float)
    for j in range(matrix.shape[1]):
        col = matrix[:, j]
        valid = np.isfinite(col)
        if np.count_nonzero(valid) < 2:
            standardized_matrix[valid, j] = 0.0
            continue
        col_mean = np.mean(col[valid])
        col_std = np.std(col[valid])
        if col_std == 0:
            standardized_matrix[valid, j] = 0.0
        else:
            standardized_matrix[valid, j] = (col[valid] - col_mean) / col_std

    # Save a plain-text table for easy downstream use.
    table_path = Path(cfg['work_dir']) / 'iod_bias_portrait_table.csv'
    with table_path.open('w', encoding='utf-8') as handle:
        handle.write('model,' + ','.join(metric_names) + '\n')
        for i, model in enumerate(model_names):
            values = [
                f"{standardized_matrix[i, j]:.6g}"
                if np.isfinite(standardized_matrix[i, j])
                else 'nan'
                for j in range(standardized_matrix.shape[1])
            ]
            handle.write(model + ',' + ','.join(values) + '\n')

    # Save standardized bias matrix as a cube.
    bias_cube = Cube(
        standardized_matrix,
        long_name='IOD standardized model-by-diagnostic bias matrix',
        units='1',
        dim_coords_and_dims=[
            (DimCoord(np.arange(len(model_names)), long_name='model_index'), 0),
            (DimCoord(np.arange(len(metric_names)), long_name='diagnostic_index'), 1),
        ],
    )
    bias_cube.add_aux_coord(AuxCoord(np.array(model_names), long_name='model_name'), 0)
    bias_cube.add_aux_coord(AuxCoord(np.array(metric_names), long_name='diagnostic_name'), 1)
    save_name = 'iod_bias_portrait_matrix'
    provenance_record = get_provenance_record(save_name, sorted(list(ancestors)))
    save_data(save_name, provenance_record, cfg, bias_cube)

    plot_bias_heatmap(cfg, standardized_matrix, model_names, metric_names, ancestors)


def plot_bias_heatmap(cfg, matrix, model_names, metric_names, ancestors):
    """Plot a model-by-diagnostic bias heatmap with fixed model ordering."""
    ordered_matrix = matrix.copy()
    ordered_models = list(model_names)

    fig_width = max(8, len(metric_names) * 1.1)
    fig_height = max(8, len(ordered_models) * 0.45)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    finite_vals = ordered_matrix[np.isfinite(ordered_matrix)]
    vmax = np.nanmax(np.abs(finite_vals)) if finite_vals.size else 1.0
    vmax = 1.0 if vmax == 0 else vmax
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    img = ax.imshow(ordered_matrix, aspect='equal', cmap='BrBG', norm=norm)
    ax.set_xticks(np.arange(len(metric_names)))
    ax.set_xticklabels(metric_names, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(np.arange(len(ordered_models)))
    ax.set_yticklabels(ordered_models, fontsize=9)
    ax.set_title('Model-by-diagnostic bias portrait', fontsize=14)
    ax.set_xlabel('Diagnostics')
    ax.set_ylabel('Models')

    cbar = fig.colorbar(img, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label('Standardized bias (z-score)')

    fig.tight_layout()
    save_name = 'iod_bias_portrait_heatmap'
    provenance_record = get_provenance_record(save_name, sorted(list(ancestors)))
    save_figure(save_name, provenance_record, cfg, bbox_inches='tight')
    plt.close(fig)

def get_iso_data(cfg, group_md, variable):
    """
    Compute T20D from 3D temperature fields and return dataset-mean cube dictionary.
    """
    logger.info(f"Processing iso data for variable: {variable}")
    iso_results = {}

    ISO_LEVEL = 20  # Define the isotherm level for T20D
    FILL_VALUE = 1e+20  # Define the fill value for masking

    # Load temperatures using load_data
    temp_data = load_data(group_md, variable, get_filenames=True)

    if not temp_data or not list(temp_data.keys()):
        logger.warning("No temperature data found. Skipping.")
        return None  

    (dataset, info), = temp_data.items()
    cube = info['cube']
    file = info['filename']

    t20d_cube = iso_depth_3d(cube, ISO_LEVEL)

    # Apply masking
    masked_data = np.ma.masked_values(t20d_cube.data, FILL_VALUE)
    #masked_mask = np.ma.getmask(masked_data)

    # Assign masked data back to cube
    t20d_cube.data = masked_data

    # Collapse the cube
    t20d_mean = t20d_cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)
    print(f"Processed T20D for dataset: {dataset}, mean value: {t20d_mean.data}")  # Debug statement
    iso_results[dataset] = {
                'cube': t20d_mean,
                'filename': file
            }

    # Save output
    save_prefix = get_prefix(file)
    output_basename = f"{save_prefix}_{dataset}_{variable}_t20d"
    provenance_record = get_provenance_record(output_basename, [file])
    save_data(output_basename, provenance_record, cfg, t20d_cube)  ## or t20d_mean?

    logger.info(f"T20D processed and saved for {dataset}")
    return iso_results
    
def scat_mean_vs_std(cfg, xy_dict, x_label, y_label, title, output_basename):
    """
    Create mean vs standard deviation for multiple datasets.
    """
    logger.info("Creating mean vs std scatter plot...")
    # Define plot
    fig = plt.figure(figsize=(10, 6))
    colors = plt.cm.tab20(np.linspace(0, 1, len(xy_dict)))  # Generate distinct colors
    input_filenames = set()

    # Lists to collect all mean and std values for correlation
    all_means = []
    all_stds = []

    for i, (dataset, xy_info) in enumerate(xy_dict.items()):
        file = xy_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])
        cube = xy_info['cube']
        print(cube)  # Debug statement
        mean_val = cube.collapsed('time', iris.analysis.MEAN).data
        std_val = cube.collapsed('time', iris.analysis.STD_DEV).data
        print(f"Processing dataset: {dataset}, file: {file}, mean: {mean_val}, std: {std_val}")  # Debug statement
        color = 'k' if 'HadISST' in dataset else colors[i]
        label = 'Obs' if 'HadISST' in dataset else dataset
        plt.scatter(mean_val, std_val, color=color, label=label, alpha=0.7,s=80)

        # Append to lists for correlation calculation
        all_means.append(mean_val)
        all_stds.append(std_val)

    print(f"All means: {all_means}")  # Debug statement
    print(f"All stds: {all_stds}")  # Debug statement
    # Calculate and log correlation
    if all_means and all_stds:
        correlation = np.corrcoef(all_means, all_stds)[0, 1]
        logger.info(f"Correlation between mean and std: {correlation:.2f}")

    # Set plot labels and limits
    plt.xlabel(x_label, fontsize=18)
    plt.ylabel(y_label, fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title(title, fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=10)
    plt.grid(True)
     # Add correlation text in top right corner
    plt.text(0.96, 1.07, f'r = {correlation:.3f}', 
             transform=plt.gca().transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    fig.tight_layout()

    # Save figure
    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info(f"Plot saved: {output_basename}")

def find_matching_y_data(dataset, x_prefix, y_dict, colors, i):
    """Match x data with corresponding y data, handling OBS specially."""    

    OBS_PREFIX = 'OBS'
    OBS_COLOUR = 'k'
    OBS_LABEL = 'OBS'
    print(f"Finding match for dataset: {dataset} with prefix: {x_prefix}")  # Debug statement
    
    if x_prefix.startswith(OBS_PREFIX):
        for y_dataset, y_candidate in y_dict.items():
            y_file = y_candidate['filename']
            y_file_base = os.path.basename(y_file) if isinstance(y_file, str) else os.path.basename(y_file[0])
            if y_file_base.startswith(x_prefix):
                print(f"Found matching OBS entry of {y_candidate}")  # Debug statement
                return y_candidate, OBS_LABEL, OBS_COLOUR
        else:
            # No matching OBS entry found
            logger.warning(f"No matching OBS entry found in y_dict for {dataset}.")
            return None, OBS_LABEL, OBS_COLOUR

    elif dataset in y_dict:
        print(f"Found matching dataset in {y_dict[dataset]}")  # Debug statement
        return y_dict[dataset], dataset, colors[i]
    
    else:
        logger.warning(f"{dataset} not found in y_dict. Skipping...")
        return None, None, None
    

def scat_plot(cfg, x_dict, y_dict, x_label, y_label, title, output_basename):
    """
    Create a scatter plot for multiple datasets.
    """
    logger.info("Creating general scatter plot...")
    # Define plot
    fig = plt.figure(figsize=(10, 6))
    colors = plt.cm.tab20(np.linspace(0, 1, len(x_dict)))  # Generate distinct colors

    input_filenames = set()

    # Lists to collect all mean and std values for correlation
    all_x = []
    all_y = []

    for i, (dataset, x_info) in enumerate(x_dict.items()):
        logger.info(f'Plotting {dataset}')
        x_cube = x_info['cube']
        x_file = x_info['filename']

        input_filenames.update(x_file if isinstance(x_file, list) else [x_file])

        x_prefix = get_prefix(x_file)


        y_info, label, color = find_matching_y_data(dataset, x_prefix, y_dict, colors, i)

        if y_info:
            y_cube = y_info['cube']
            y_file = y_info['filename']
            print(np.shape(x_cube.data), np.shape(y_cube.data))  # Debug statement
            plt.scatter(x_cube.data, y_cube.data, color=color,
                        label=label, alpha=0.7,s=80)
            
            input_filenames.update(y_file if isinstance(y_file, list) else [y_file])

        # Store values for correlation calculation
        all_x.append(x_cube.data)
        all_y.append(y_cube.data)

    # Set plot labels and limits
    plt.xlabel(x_label, fontsize=18)
    plt.ylabel(y_label, fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title(title, fontsize=18)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=10)
    plt.grid(True)
    fig.tight_layout()

    # Calculate and log correlation
    if all_x and all_y:
        correlation = np.corrcoef(all_x, all_y)[0, 1]
        logger.info(f"Correlation between x and y: {correlation:.2f}")
        # Add correlation text in top right corner
        plt.text(0.96, 1.07, f'r = {correlation:.3f}', 
                 transform=plt.gca().transAxes,
                 fontsize=14, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Save figure
    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info(f"Scatter plot saved: {output_basename}")
    plt.close()

def main(cfg):
    """
     Main function to compute and plot scatter diagnostics.
    """
    logger.info("Starting main diagnostic process.")
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')

    therm_tilt, sst_grad, eq_winds, skew_dmi, east_sst_son, west_sst_son, sctr_t20d, nino_ssts, dmi, east_sst_son_ts, west_sst_son_ts = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    for group_name, group_md in grouped_data.items():
        logger.info(f"Processing group: {group_name}")
        # If checks are necessary as not all models/obs (group names) have all variables.

        # Load west and east SST using load_data and calculate gradient
        west_sst_son_item = load_data(group_md, 'west_sst_son', get_filenames=True)
        east_sst_son_item = load_data(group_md, 'east_sst_son', get_filenames=True)
        if west_sst_son_item:
            west_sst_son.update(west_sst_son_item)
        if east_sst_son_item:
            east_sst_son.update(east_sst_son_item)
        if west_sst_son_item and east_sst_son_item:
            sst_grad_item = compute_cube_diff(cfg, west_sst_son_item, east_sst_son_item, 'sst_grad')
            if sst_grad_item:
                sst_grad.update(sst_grad_item)
        
        # Load 3d ocean temps, calculate and find area average t20d
        # Calculate gradient for thermocline tilt
        west_t20d = get_iso_data(cfg, group_md, 'west_temps')
        east_t20d = get_iso_data(cfg, group_md, 'east_temps')
        
        if west_t20d and east_t20d:
            t20d_tilt = compute_cube_diff(cfg, west_t20d, east_t20d, 'therm_tilt')
            if t20d_tilt:
                therm_tilt.update(t20d_tilt)

        load_and_update_dict(group_md, 'east_sst_son_ts', east_sst_son_ts)
        load_and_update_dict(group_md, 'west_sst_son_ts', west_sst_son_ts)
        load_and_update_dict(group_md, 'eq_winds', eq_winds)
        load_and_update_dict(group_md, 'nino_ssts', nino_ssts)
        sctr_t20d_item = get_iso_data(cfg, group_md, 'sctr_temps')
        if sctr_t20d_item:
            sctr_t20d.update(sctr_t20d_item) 
              

        # Load SST anomalies, calculate diff for DMI and find skewness
        west_anoms = load_data(group_md, 'west_anoms_ts', get_filenames=True)
        east_anoms = load_data(group_md, 'east_anoms_ts', get_filenames=True)

        if west_anoms and east_anoms:
            dmi_item = compute_cube_diff(cfg, west_anoms, east_anoms, 'dmi')
            if dmi_item:
                dmi.update(dmi_item)
                skew_dmi_item = compute_cube_skew(cfg, dmi_item, 'skew_dmi')
                skew_dmi.update(skew_dmi_item)

    # Build scalar diagnostics and biases (model - observation).
    obs_sst = get_obs_dataset_name(east_sst_son, OBS_CANDIDATES['sst'])
    obs_wind = get_obs_dataset_name(eq_winds, OBS_CANDIDATES['wind'])
    obs_ocean = get_obs_dataset_name(therm_tilt, OBS_CANDIDATES['ocean'])
    obs_dmi = get_obs_dataset_name(dmi, OBS_CANDIDATES['dmi'])

    obs_values = {}
    if obs_sst and obs_sst in east_sst_son and obs_sst in west_sst_son and obs_sst in sst_grad:
        obs_values['EEIO mean SST bias'] = as_scalar(east_sst_son[obs_sst]['cube'], 'mean')
        obs_values['WEIO mean SST bias'] = as_scalar(west_sst_son[obs_sst]['cube'], 'mean')
        obs_values['WEIO-EEIO SST gradient bias'] = as_scalar(sst_grad[obs_sst]['cube'], 'mean')
    if obs_wind and obs_wind in eq_winds:
        obs_values['CEIO zonal wind bias'] = as_scalar(eq_winds[obs_wind]['cube'], 'mean')
    if obs_ocean and obs_ocean in therm_tilt and obs_ocean in sctr_t20d:
        obs_values['equatorial D20 tilt bias'] = as_scalar(therm_tilt[obs_ocean]['cube'], 'mean')
        obs_values['SCTR D20 bias'] = as_scalar(sctr_t20d[obs_ocean]['cube'], 'mean')
    if obs_dmi and obs_dmi in dmi and obs_dmi in skew_dmi:
        obs_values['DMI standard deviation bias'] = as_scalar(dmi[obs_dmi]['cube'], 'std')
        obs_values['DMI skewness bias'] = as_scalar(skew_dmi[obs_dmi]['cube'], 'raw')

    diagnostics_bias = {}
    candidate_models = sorted(
        set(east_sst_son.keys()) |
        set(west_sst_son.keys()) |
        set(sst_grad.keys()) |
        set(eq_winds.keys()) |
        set(therm_tilt.keys()) |
        set(sctr_t20d.keys()) |
        set(dmi.keys()) |
        set(skew_dmi.keys())
    )
    excluded = set([obs_sst, obs_wind, obs_ocean, obs_dmi])

    for model in candidate_models:
        if model in excluded:
            continue

        model_metrics = {}
        model_ancestors = set()

        if model in east_sst_son and 'EEIO mean SST bias' in obs_values:
            val = as_scalar(east_sst_son[model]['cube'], 'mean') - obs_values['EEIO mean SST bias']
            model_metrics['EEIO mean SST bias'] = val
            model_ancestors.update(to_set(east_sst_son[model]['filename']))
            save_scalar_bias(cfg, 'eeio_mean_sst', model, val, model_ancestors)

        if model in west_sst_son and 'WEIO mean SST bias' in obs_values:
            val = as_scalar(west_sst_son[model]['cube'], 'mean') - obs_values['WEIO mean SST bias']
            model_metrics['WEIO mean SST bias'] = val
            model_ancestors.update(to_set(west_sst_son[model]['filename']))
            save_scalar_bias(cfg, 'weio_mean_sst', model, val, model_ancestors)

        if model in sst_grad and 'WEIO-EEIO SST gradient bias' in obs_values:
            val = as_scalar(sst_grad[model]['cube'], 'mean') - obs_values['WEIO-EEIO SST gradient bias']
            model_metrics['WEIO-EEIO SST gradient bias'] = val
            model_ancestors.update(to_set(sst_grad[model]['filename']))
            save_scalar_bias(cfg, 'weio_minus_eeio_sst_grad', model, val, model_ancestors)

        if model in eq_winds and 'CEIO zonal wind bias' in obs_values:
            val = as_scalar(eq_winds[model]['cube'], 'mean') - obs_values['CEIO zonal wind bias']
            model_metrics['CEIO zonal wind bias'] = val
            model_ancestors.update(to_set(eq_winds[model]['filename']))
            save_scalar_bias(cfg, 'ceio_zonal_wind', model, val, model_ancestors)

        if model in therm_tilt and 'equatorial D20 tilt bias' in obs_values:
            val = as_scalar(therm_tilt[model]['cube'], 'mean') - obs_values['equatorial D20 tilt bias']
            model_metrics['equatorial D20 tilt bias'] = val
            model_ancestors.update(to_set(therm_tilt[model]['filename']))
            save_scalar_bias(cfg, 'equatorial_d20_tilt', model, val, model_ancestors)

        if model in sctr_t20d and 'SCTR D20 bias' in obs_values:
            val = as_scalar(sctr_t20d[model]['cube'], 'mean') - obs_values['SCTR D20 bias']
            model_metrics['SCTR D20 bias'] = val
            model_ancestors.update(to_set(sctr_t20d[model]['filename']))
            save_scalar_bias(cfg, 'sctr_d20', model, val, model_ancestors)

        if model in dmi and 'DMI standard deviation bias' in obs_values:
            val = as_scalar(dmi[model]['cube'], 'std') - obs_values['DMI standard deviation bias']
            model_metrics['DMI standard deviation bias'] = val
            model_ancestors.update(to_set(dmi[model]['filename']))
            save_scalar_bias(cfg, 'dmi_std', model, val, model_ancestors)

        if model in skew_dmi and 'DMI skewness bias' in obs_values:
            val = as_scalar(skew_dmi[model]['cube'], 'raw') - obs_values['DMI skewness bias']
            model_metrics['DMI skewness bias'] = val
            model_ancestors.update(to_set(skew_dmi[model]['filename']))
            save_scalar_bias(cfg, 'dmi_skewness', model, val, model_ancestors)

        if model_metrics:
            model_metrics['ancestors'] = model_ancestors
            diagnostics_bias[model] = model_metrics

    build_bias_table(cfg, diagnostics_bias)

    print('therm_tilt', therm_tilt)
    scat_plot(cfg, sst_grad, therm_tilt, 'SST gradient / $^\\circ$C', 'Thermocline tilt / m', 'SON', 'sst_vs_tilt')
    scat_plot(cfg, sst_grad, eq_winds, 'SST gradient / $^\\circ$C', 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'SON', 'sst_vs_winds')
    scat_plot(cfg, eq_winds, therm_tilt, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'Thermocline tilt / m', 'SON', 'winds_vs_tilt')
    scat_plot(cfg, eq_winds, skew_dmi, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'Skewness of DMI', 'SON', 'winds_vs_skew')
    scat_plot(cfg, eq_winds, sctr_t20d, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'SCTR 20$^\\circ$C isotherm depth / m', 'SON', 'winds_vs_sctr')
    scat_plot(cfg, sst_grad, nino_ssts, 'SST gradient / $^\\circ$C', 'Nino 3.4 SST / $^\\circ$C', 'SST gradient-SON, Nino-DJF', 'dmi_vs_nino')
    scat_plot(cfg, east_sst_son, west_sst_son, 'EEIO SST / $^\\circ$C', 'WEIO SST / $^\\circ$C', 'SON', 'east_sst_vs_west_sst')
    scat_plot(cfg, eq_winds, east_sst_son, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'EEIO SST / $^\\circ$C', 'SON', 'winds_vs_east_sst')
    scat_plot(cfg, eq_winds, west_sst_son, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'WEIO SST / $^\\circ$C', 'SON', 'winds_vs_west_sst')
    scat_mean_vs_std(cfg, east_sst_son_ts, 'SST mean in EEIO / $^\\circ$C', 'STD of EEIO SST / $^\\circ$C', 'SON', 'east_sst_mean_vs_std')
    scat_mean_vs_std(cfg, west_sst_son_ts, 'SST mean in WEIO / $^\\circ$C', 'STD of WEIO SST / $^\\circ$C', 'SON', 'west_sst_mean_vs_std')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)

