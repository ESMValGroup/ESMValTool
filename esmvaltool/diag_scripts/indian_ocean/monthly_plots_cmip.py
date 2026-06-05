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


def is_mohc_dataset(dataset):
    """Return True if dataset likely belongs to MOHC."""
    dataset_upper = dataset.upper()
    mohc_markers = ("MOHC", "HADGEM", "HADCM", "UKESM")
    return any(marker in dataset_upper for marker in mohc_markers)

def plot_ts(cfg, plot_dict, title, y_title, output_basename):
    """
    Plot monthly values for all datasets in a single figure.
    """

    logger.info("Plotting monthly timeseries")
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    plt.figure(figsize=(10, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, len(plot_dict)))
    input_filenames = set()
    max_months = 12
    x_labels = month_labels

    for i, (dataset, dict_info) in enumerate(plot_dict.items()):
        cube = dict_info['cube']
        file = dict_info['filename']
        input_filenames.update(file if isinstance(file, list) else [file])

        values = np.asarray(cube.data, dtype=float).squeeze()
        if values.ndim != 1:
            logger.warning(f"Skipping {dataset}: expected 1D monthly data, got shape {values.shape}.")
            continue

        n_months = values.shape[0]
        x = np.arange(n_months)

        if n_months == 12:
            x_labels = month_labels
        else:
            x_labels = [str(i + 1) for i in range(n_months)]
        max_months = max(max_months, n_months)

        if dataset == 'ERA5' or dataset == 'NCEP':
            color = 'black'
            linewidth = 2.2
            alpha = 1.0
        elif is_mohc_dataset(dataset):
            color = colors[i]
            linewidth = 1.4
            alpha = 0.5
        else:
            color = colors[i]
            linewidth = 1.2
            alpha = 0.35

        plt.plot(x, values, label=dataset, color=color, linewidth=linewidth, alpha=alpha)

    plt.xlim(-0.5, max_months - 0.5)
    if max_months == 12:
        tick_labels = month_labels
    else:
        tick_labels = [str(i + 1) for i in range(max_months)]
    plt.xticks(np.arange(max_months), tick_labels)
    plt.xlabel("Month", fontsize=12)
    plt.ylabel(y_title, fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()

    provenance_record = get_provenance_record(output_basename, list(input_filenames))
    save_figure(output_basename, provenance_record, cfg)
    logger.info(f"Monthly plot saved: {output_basename}")
    plt.close()


def main(cfg):
    """Plot monthly CEIO windspeeds for multiple datasets and observations."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'dataset')
    ceio_wind_monthly, east_sst_monthly = {}, {}
    for group_name, group_md in grouped_data.items():
        load_and_update_dict(group_md, 'ceio_wind_monthly', ceio_wind_monthly)
        load_and_update_dict(group_md, 'east_sst_monthly', east_sst_monthly)
    
    logger.info("Thermocline calculated, now plotting.")
    # Plot results for all datasets
    plot_ts(
        cfg,
        ceio_wind_monthly,
        'Monthly CEIO zonal wind speed',
        'Zonal wind speed / m $\\mathregular{s^{-1}}$',
        'ceio_winds_monthly',
    )
    plot_ts(
        cfg,
        east_sst_monthly,
        'Monthly SST in the EEIO',
        'SST / $^\\circ$C',
        'east_sst_monthly',
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

