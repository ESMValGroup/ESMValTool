import logging
from pathlib import Path
from pprint import pformat
import os
import iris # type: ignore
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from basic_functions import (compute_cube_skew, get_provenance_record, iso_depth_3d, 
                             load_data, 
                             compute_cube_diff, 
                             get_prefix,
                             load_and_update_dict)

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
        mean_val = cube.collapsed('time', iris.analysis.MEAN).data
        std_val = cube.collapsed('time', iris.analysis.STD_DEV).data

        color = 'k' if 'HadISST' in dataset else colors[i]
        label = 'Obs' if 'HadISST' in dataset else dataset
        plt.scatter(mean_val, std_val, color=color, label=label, alpha=0.7,s=80)

        # Append to lists for correlation calculation
        all_means.append(mean_val)
        all_stds.append(std_val)

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
            if dataset == 'NCEP':
                logger.debug(f"{y_file}")
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

    therm_tilt, sst_grad, eq_winds, east_ts, skew_dmi, sctr_t20d, nino_anoms, dmi = {}, {}, {}, {}, {}, {}, {}, {}
    for group_name, group_md in grouped_data.items():
        logger.info(f"Processing group: {group_name}")
        # If checks are necessary as not all models/obs (group names) have all variables.

        # Load west and east SST using load_data and calculate gradient
        west_ssts = load_data(group_md, 'west_ssts', get_filenames=True)
        east_ssts = load_data(group_md, 'east_ssts', get_filenames=True)
        if west_ssts and east_ssts:
            sst_grad_item = compute_cube_diff(cfg, west_ssts, east_ssts, 'sst_grad')
            sst_grad.update(sst_grad_item)
        
        # Load 3d ocean temps, calculate and find area average t20d
        # Calculate gradient for thermocline tilt
        west_t20d = get_iso_data(cfg, group_md, 'west_temps')
        east_t20d = get_iso_data(cfg, group_md, 'east_temps')

        if west_t20d and east_t20d:
            t20d_tilt = compute_cube_diff(cfg, west_t20d, east_t20d, 'therm_tilt')
            therm_tilt.update(t20d_tilt)

        load_and_update_dict(group_md, 'eq_winds', eq_winds)
        load_and_update_dict(group_md, 'east_ssts_ts', east_ts)
        load_and_update_dict(group_md, 'nino_anoms', nino_anoms)
        sctr_t20d_item = get_iso_data(cfg, group_md, 'sctr_temps')
        if sctr_t20d_item:
            sctr_t20d.update(sctr_t20d_item) 
              

        # Load SST anomalies, calculate diff for DMI and find skewness
        west_anoms = load_data(group_md, 'west_anoms', get_filenames=True)
        east_anoms = load_data(group_md, 'east_anoms', get_filenames=True)

        if west_anoms and east_anoms:
            dmi_item = compute_cube_diff(cfg, west_anoms, east_anoms, 'dmi')
            dmi.update(dmi_item)
            skew_dmi_item = compute_cube_skew(cfg, dmi_item, 'skew_dmi')
            skew_dmi.update(skew_dmi_item)

    
    scat_plot(cfg, sst_grad, therm_tilt, 'SST gradient / $^\\circ$C', 'Thermocline tilt / m', 'SON', 'sst_vs_tilt')
    scat_plot(cfg, sst_grad, eq_winds, 'SST gradient / $^\\circ$C', 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'SON', 'sst_vs_winds')
    scat_plot(cfg, eq_winds, therm_tilt, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'Thermocline tilt / m', 'SON', 'winds_vs_tilt')
    scat_plot(cfg, eq_winds, skew_dmi, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'Skewness of DMI', 'SON', 'winds_vs_skew')
    scat_plot(cfg, eq_winds, sctr_t20d, 'Zonal wind speed in CEIO / m $\\mathregular{s^{-1}}$', 'SCTR 20$^\\circ$C isotherm depth / m', 'SON', 'winds_vs_sctr')
    scat_plot(cfg, sst_grad, nino_anoms, 'SST gradient / $^\\circ$C', 'Nino 3.4 SST / $^\\circ$C', 'SST gradient-SON, Nino-DJF', 'dmi_vs_nino')


    scat_mean_vs_std(cfg, east_ts, 'SST mean in EEIO / $^\\circ$C', 'STD of EEIO SST / $^\\circ$C', 'SON', 'east_sst_mean_vs_std')

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)

