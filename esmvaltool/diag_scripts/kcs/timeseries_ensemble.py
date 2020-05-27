"""
Visualize temperature accross datasets
"""

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr


from esmvaltool.diag_scripts.shared import (run_diagnostic, get_plot_filename,
                                            get_diagnostic_filename)

def plot_datasets(cfg):

    # a list of dictionaries describing all datasets passed on to the recipe
    metadata = cfg['input_data'].values()

    # Create a figure and a list that we'll use to store output data to csv
    fig, ax = plt.subplots()
    statistics = []

    # Loop over all datasets
    for info_dict in metadata:

        # Open file and add timeseries to figure
        ds = xr.open_dataset(info_dict['filename'])
        ds.tas.plot(ax=ax, label=info_dict['alias'])

        # Add statistics to list for later saving
        if 'MultiModel' in info_dict['alias']:
            s = ds.tas.to_series()
            s.name = info_dict['alias']
            statistics.append(s)

    # Save figure
    ax.legend()
    filename = get_plot_filename('temperature_change_pdf', cfg)
    fig.savefig(filename, bbox_inches='tight', dpi=300)

    # Save csv file
    filename = get_diagnostic_filename('tas_annual_change', cfg, extension='csv')
    df = pd.concat(statistics, axis=1)
    df.to_csv(filename)

    return 'Ready'


if __name__ == '__main__':
    with run_diagnostic() as config:
        plot_datasets(config)
