import logging
from pathlib import Path
import iris
import matplotlib.pyplot as plt
from esmvalcore.preprocessor import (
    extract_levels,
    anomalies,
)

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)

def plotting_support(cube, key, linestyle='-', marker=None, color=None):
    """Helper function for pretty plotting."""
    print(f"Plotting {key} - Cube data: {cube.data}")

    if cube.coords('time', dim_coords=True):
        ih.unify_time_coord(cube)

    plt.plot(cube.coord('time').points, cube.data, label=key, linestyle=linestyle, marker=marker, color=color)

def plot_anomalies(input_file, name, variable_group, color=None):
    """
    Plot anomalies over time using the provided plotting support function.
    """
    try:
        # Load data
        cube = iris.load_cube(input_file)

        # Extract time and anomalies
        if cube.coords('time', dim_coords=True):
            ih.unify_time_coord(cube)

        # Plotting function by using the plotting_support function
        plotting_support(cube, key=name, linestyle='-', marker=None, color=color)

        return {'cube': cube, 'name': name}

    except Exception as e:
        logger.error(f"Error plotting: {e}")
        return None

def create_subplot(ax, data_list, variable_group):
    """Plot subplot over time considering minimum and maximum time across dataset."""
    for i, data in enumerate(data_list):
        color = plt.cm.tab10(i % 10)
        plotting_support(data['cube'], key=data['name'], linestyle='-', marker=None, color=color)
   
    ax.set_ylabel(f"{data_list[0]['cube'].var_name.upper()} Anomalies, " + str(data_list[0]['cube'].units))
    ax.set_title(f"Time series of monthly mean {data_list[0]['cube'].var_name.upper()} anomalies - {variable_group}")
    ax.legend([data['name'] for data in data_list], loc="center", bbox_to_anchor=(1, 0.5), borderaxespad=0.)
    ax.set_xlim('2010', '2014')
    ax.set_xlabel('Time')
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(plt.matplotlib.dates.AutoDateLocator())

def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input
    input_data = cfg['input_data'].values()

    # Dictionary to store data for each variable_group
    data_by_variable_group = {}

    # loop for variables/dataset
    for i, dataset in enumerate(input_data):
        
        # Load data
        input_file = dataset['filename']
        name = dataset['dataset']
        variable_group = dataset['variable_group']

        # Check whether variable_group is alreadypresent
        if variable_group not in data_by_variable_group:
            data_by_variable_group[variable_group] = []

        # color for the dataset 
        color = plt.cm.tab10(i % 10)

        # Plotting function for each variable and appending data to the dictionary
        data = plot_anomalies(input_file=input_file, name=name, variable_group=variable_group, color=color)

        if data:
            data_by_variable_group[variable_group].append(data)

    # Check the data_list
    for variable_group, data_list in data_by_variable_group.items():
        if not data_list:
            continue

        plt.figure(figsize=(8, 6))

        create_subplot(plt.gca(), data_list, variable_group)
         
        # Save separate figures for each variable
        plt.savefig(f'plot_{variable_group}.png')

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)



