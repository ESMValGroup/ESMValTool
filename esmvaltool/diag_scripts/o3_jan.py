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

    plt.plot(cube.coord('time').points, cube.data, label=key,
 linestyle=linestyle, marker=marker, color=color)

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

        # Plotting using the provided plotting_support function
        plotting_support(cube, key=name, linestyle='-',
 marker=None, color=color)

        return {'cube': cube, 'name': name}

    except Exception as e:
        logger.error(f"Error plotting: {e}")
        return None

def create_subplot(ax, data_list, variable_group, min_time, max_time):
    """Plot subplot by considering min & max time."""
    for i, data in enumerate(data_list):
        color = plt.cm.tab10(i % 10)  
        plotting_support(data['cube'], key=data['name'], linestyle='-', marker=None, color=color)

    ax.set_ylabel(f"{data_list[0]['cube'].var_name.upper()} Anomalies, " + str(data_list[0]['cube'].units))
    ax.set_title(f"Time series of monthly mean {data_list[0]['cube'].var_name.upper()} anomalies - {variable_group}")
    ax.legend([data['name'] for data in data_list], loc="center", bbox_to_anchor=(1, 0.5), borderaxespad=0.)
    ax.set_xlim(min_time, max_time)
    ax.set_xlabel('Time')
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(plt.matplotlib.dates.AutoDateLocator())


def main(cfg, start_year, end_year):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input
    input_data = cfg['input_data'].values()

    # Dictionary to store data for each variable_group
    data_by_variable_group = {}

    # Loop for variables/datasets
    for i, dataset in enumerate(input_data):
        # Load the data
        input_file = dataset['filename']
        name = dataset['dataset']
        variable_group = dataset['variable_group']

        # Check whether variable_group is already in the dictionary
        if variable_group not in data_by_variable_group:
            data_by_variable_group[variable_group] = []

        # Select color for the dataset
        color = plt.cm.tab10(i % 10)
        # Plotting function for each variable and append data to the dictionary
        data = plot_anomalies(input_file=input_file, name=name, variable_group=variable_group, color=color)

        if data:
            data_by_variable_group[variable_group].append(data)

    # Minimum and maximum time values across all datasets
    min_time = min(data['cube'].coord('time').points.min()
 for data_list in data_by_variable_group.values() for data in data_list)
    max_time = max(data['cube'].coord('time').points.max()
 for data_list in data_by_variable_group.values() for data in data_list)


    for variable_group, data_list in data_by_variable_group.items():
        if not data_list:  
            continue

        plt.figure(figsize=(8, 6))

        create_subplot(plt.gca(), data_list,
 variable_group, min_time, max_time)

        # Save the figures for each variable
        plt.savefig(f'plot_{variable_group}.png')


if __name__ == '__main__':
    with run_diagnostic() as config:
        start_year = 2010  # Specifying the start year
        end_year = 2014    # Specifying the end year
        main(config, start_year, end_year)
