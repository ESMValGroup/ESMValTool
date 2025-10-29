import logging
from pathlib import Path
from pprint import pformat
import os
import iris # type: ignore
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from basic_functions import (get_provenance_record, 
                             load_data, 
                             compute_cube_diff, 
                             iso_depth_3d,
                             get_prefix)

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

    # Load temperatures using load_data
    temp_data = load_data(group_md, variable, get_filenames=True)

    if not temp_data or not list(temp_data.keys()):
        logger.warning("No temperature data found. Skipping.")
        return None  

    (dataset, info), = temp_data.items()
    cube = info['cube']
    file = info['filename']
    logger.debug(f"Selected dataset: {dataset}")

    t20d_cube = iso_depth_3d(cube, 20)

    # Check for invalid values before masking
    has_nan = np.isnan(t20d_cube.data).any()
    has_inf = np.isinf(t20d_cube.data).any()
    mask = np.ma.getmask(t20d_cube.data)
    logger.warning(f"Before masking: NaN={has_nan}, Inf={has_inf}, Writeable={t20d_cube.data.flags.writeable}")
    print("Any masked values before masking:", np.any(mask))
    print("Mask shape matches data before masking:", mask.shape == t20d_cube.data.shape)

    # Apply masking
    #masked_data = np.ma.masked_invalid(t20d_cube.data)
    masked_data = np.ma.masked_values(t20d_cube.data, 1e+20)
    masked_mask = np.ma.getmask(masked_data)
    print("Any masked values after masking:", np.any(masked_mask))
    print("Mask shape matches data after masking:", masked_mask.shape == masked_data.data.shape)
    logger.warning(f"After masking: Writeable={masked_data.flags.writeable}, IsMasked={np.ma.isMaskedArray(masked_data)}")

    # Assign masked data back to cube
    t20d_cube.data = masked_data

    # Collapse the cube
    t20d_mean = t20d_cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)

    # Check writeability after collapse
    logger.warning(f"After collapse: Writeable={t20d_mean.data.flags.writeable}")


    # 
    # t20d_cube.data = np.ma.masked_invalid(t20d_cube.data)
    # t20d_mean = t20d_cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)
    # if not t20d_mean.data.flags.writeable:
    #     t20d_mean = t20d_mean.copy()    
    # print(file, cube.data.flags.writeable, t20d_mean.data.flags.writeable)
    iso_results[dataset] = {
                'cube': t20d_mean,
                'filename': file
            }

    # Save output
    save_prefix = get_prefix(file)
    output_basename = f"{save_prefix}_{dataset}_{variable}_t20d"
    provenance_record = get_provenance_record(output_basename, [file])
    save_data(output_basename, provenance_record, cfg, t20d_cube)

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

    for i, (dataset, xy_info) in enumerate(xy_dict.items()):
        file = xy_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])
        cube = xy_info['cube']
        mean_val = cube.collapsed('time', iris.analysis.MEAN).data
        std_val = cube.collapsed('time', iris.analysis.STD_DEV).data

        color = 'k' if 'HadISST' in dataset else colors[i]
        label = 'Obs' if 'HadISST' in dataset else dataset
        plt.scatter(mean_val, std_val, color=color, label=label, alpha=0.7,s=80)

    # Set plot labels and limits
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title(title, fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=10)
    plt.grid(True)
    fig.tight_layout()

    # Save figure
    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info(f"Plot saved: {output_basename}")


def scat_plot(cfg, x_dict, y_dict, x_label, y_label, title, output_basename):
    """
    Create a scatter plot for multiple datasets.
    """
    logger.info("Creating general scatter plot...")
    # Define plot
    fig = plt.figure(figsize=(10, 6))
    colors = plt.cm.tab20(np.linspace(0, 1, len(x_dict)))  # Generate distinct colors

    input_filenames = set()

    for i, (dataset, x_info) in enumerate(x_dict.items()):
        logger.info(f'Plotting {dataset}')
        x_cube = x_info['cube']
        x_file = x_info['filename']

        input_filenames.update(x_file if isinstance(x_file, list) else [x_file])

        x_prefix = get_prefix(x_file)
        y_info = None
            
        if x_prefix.startswith('OBS'):
            for y_dataset, y_candidate in y_dict.items():
                y_file = y_candidate['filename']
                y_file_base = os.path.basename(y_file) if isinstance(y_file, str) else os.path.basename(y_file[0])
                if y_file_base.startswith(x_prefix):
                    y_info = y_candidate
                    label = 'Obs'
                    color = 'k'
                    break
            else:
                logger.warning(f"No matching OBS entry found in y_dict for {dataset}. Skipping...")
                plt.scatter(x_cube.data, 0, color='k', label='Obs', alpha=0.7,s=80)
                #plt.scatter(x_cube.data, -8.43, color='k', label='Obs', alpha=0.7,s=80)
                continue  # Skip if no OBS match found
        elif dataset in y_dict:
            y_info = y_dict[dataset]
            label = dataset
            color = colors[i]


        else:
            logger.warning(f"{dataset} not found in y_dict. Skipping...")
            continue

        if y_info:
            y_cube = y_info['cube']
            y_file = y_info['filename']
            plt.scatter(x_cube.data, y_cube.data, color=color,
                        label=label, alpha=0.7,s=80)
            
            input_filenames.update(y_file if isinstance(y_file, list) else [y_file])

    # Set plot labels and limits
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title(title, fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',fontsize=10)
    plt.grid(True)
    fig.tight_layout()

    # Save figure
    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg, bbox_inches='tight')
    logger.info(f"Scatter plot saved: {output_basename}")
    plt.close()

def check_item_load(metadata, variable, dict_name):
    if (item := load_data(metadata, variable, get_filenames=True)):
            dict_name.update(item)

def main(cfg):
    """
     Main function to compute and plot scatter diagnostics.
    """
    logger.info("Starting main diagnostic process.")
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')

    therm_tilt, sst_grad, eq_winds, east_sst, skew_dmi, sctr_t20d = {}, {}, {}, {}, {}, {}
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

        check_item_load(group_md, 'eq_winds',eq_winds)
        check_item_load(group_md, 'east_ssts_ts', east_sst)

        # eq_winds_item = load_data(group_md, 'eq_winds', get_filenames=True)
        # if eq_winds_item:
        #     eq_winds.update(eq_winds_item)
        # east_sst.update(load_data(group_md, 'east_ssts_ts', get_filenames=True))

        sctr_t20d_item = get_iso_data(cfg, group_md, 'sctr_temps')
        if sctr_t20d_item:
            sctr_t20d.update(sctr_t20d_item) 
              

        # Load SST anomalies, calculate diff for DMI and find skewness
        west_anoms = load_data(group_md, 'west_anoms', get_filenames=True)
        east_anoms = load_data(group_md, 'east_anoms', get_filenames=True)

        if west_anoms and east_anoms:
            dmi = compute_cube_diff(cfg, west_anoms, east_anoms, 'skew_dmi', skew=True)
            skew_dmi.update(dmi)

    
    scat_plot(cfg, sst_grad, therm_tilt, 'SST gradient / $^\circ$C', 'Thermocline tilt / m', 'SON', 'sst_vs_tilt')
    scat_plot(cfg, sst_grad, eq_winds, 'SST gradient / $^\circ$C', 'Zonal wind speed in CEIO / m $\mathregular{s^{-1}}$', 'SON', 'sst_vs_winds')
    scat_plot(cfg, eq_winds, therm_tilt, 'Zonal wind speed in CEIO / m $\mathregular{s^{-1}}$', 'Thermocline tilt / m', 'SON', 'winds_vs_tilt')
    scat_plot(cfg, eq_winds, skew_dmi, 'Zonal wind speed in CEIO / m $\mathregular{s^{-1}}$', 'Skewness of DMI', 'SON', 'winds_vs_skew')
    scat_plot(cfg, eq_winds, sctr_t20d, 'Zonal wind speed in CEIO / m $\mathregular{s^{-1}}$', 'SCTR 20$^\circ$C isotherm depth / m', 'SON', 'winds_vs_sctr')


    scat_mean_vs_std(cfg, east_sst, 'SST mean in EEIO / $^\circ$C', 'STD of EEIO SST / $^\circ$C', 'SON', 'east_sst_mean_vs_std')

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)

