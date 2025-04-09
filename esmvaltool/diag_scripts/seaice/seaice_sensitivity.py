import logging
from pathlib import Path
from pprint import pformat
import iris
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import iris.coord_categorisation as cc  # we could avoid using this
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
    select_metadata,
)


logger = logging.getLogger(Path(__file__).stem)


# This is stolen directly from the example in the diag_scripts
# TODO: Work out what should be written here
def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    # Associated recipe uses contains a caption string with placeholders
    # like {long_name} that are now populated from attributes dictionary.
    # Note that for simple recipes, caption can be set here as a simple string
    # caption = attributes["caption"].format(**attributes)

    record = {
        "caption": 'I CHANGED THIS TO A STRING caption',
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["zonal"],
        "authors": [
            "andela_bouwe",
            "righi_mattia",
        ],
        "references": [
            "acknow_project",
        ],
        "ancestors": ancestor_files,
    }
    return record


# Ed Blockley's observations dictionary
observations = {
    '1979-2014': {
        # dSIA/dT values in million sq km per Kelvin
        'mean': -4.01,
        'std_dev': 0.32,
        'plausible': 1.28,
    }
}


def list_datasets(data):
    """ Actually returns a set of datatsets, to avoid duplication. """
    logger.debug('listing datasets')
    datasets = set()
    for element in data:
        datasets.add(element['dataset'])
    return datasets


def extract_cube(data, variable_group):
    """
    Selects records from data based on the variable_group,
    asserts that there is only one matching record,
    then returns an iris cube for that record
    """
    logger.debug(f'extracting {variable_group} cube from {data}')

    # Load any data in the variable_group
    selection = select_metadata(data, variable_group=variable_group)

    # Ensure there is only one file in the list
    assert len(selection) == 1, f'Too many matching files found for {variable_group}'

    # Load the cube, [0] is because selection returns a list
    cube = iris.load_cube(selection[0]['filename'])

    return cube


# Stolen from Ed Blockley - do we ever need to return more than the slope?
def calculate_trend(times, timeseries, slope_only=True):
    """
    Calculate linear trend for a timeseries
    Returns a numpy array containing the linear fit for that trend, which can
    be subtracted from another run to de-trend it.
    """
    logger.debug(f'calculating trend for {timeseries}')

    # Use SciPy stats to calculate the slope
    slope, intercept, r, p, stderr = stats.linregress(times, timeseries)

    # Return either the slope or the array [0, 1*slope, 2*slope, 3*slope, ...]
    if slope_only:
        return slope
    else:
        return slope * np.arange(len(times))


def calculate_sensitivity(data):
    """
    Calculates change in sea ice area per year
    divided by change in global mean temperature per year
    """
    logger.debug('calculating sensitivity')

    # Load the preprocessed cubes
    si_cube = extract_cube(data, 'arctic_sea_ice')
    tas_cube = extract_cube(data, 'avg_ann_global_temp')

    # Changing the sea ice cube's units to match literature
    si_cube.convert_units('1e6 km2')

    # Add years to si cube (tas cube already has years)
    cc.add_year(si_cube, 'time', name='year')
    # Check that the years match between the two cubes
    assert (si_cube.coord('year').points == tas_cube.coord('year').points).all()
    # Name the years array for later use
    years = tas_cube.coord('year').points

    # Calculate trends and sensitivity, both via time as in Ed Blockley's code
    si_trend = calculate_trend(years, si_cube.data)
    tas_trend = calculate_trend(years, tas_cube.data)
    sensitivity = si_trend / tas_trend

    return sensitivity


def plot_from_dict(dictionary, filename, cfg):
    """
    Saves a plot of sensitivities and observations to the given filename
    """
    logger.debug(f'creating plot {filename}')

    # Reading from Ed Blockley's dictionary
    obs_years = list(observations)[0]
    obs_mean = observations[obs_years]['mean']
    obs_std_dev = observations[obs_years]['std_dev']
    obs_plausible = observations[obs_years]['plausible']

    # Set up the figure
    fig, ax = plt.subplots(figsize=(4, 6), layout='constrained')
    fig.suptitle('September Arctic Sea Ice Sensitivity')
    ax.set_title(f'dSIA/dGMST {obs_years}')  # Ed's title

    # Iterate over the dictionary
    for dataset, sensitivity in dictionary.items():
        ax.plot(0.5, sensitivity, label=dataset, marker='_', markersize=20)

    # Add observations (style taken from Ed's code)
    ax.hlines(obs_mean, 0, 1, linestyle='--', color='black', linewidth=2)
    ax.hlines(obs_mean + obs_std_dev, 0, 1, linestyle=':', color='black', linewidth=1)
    ax.hlines(obs_mean - obs_std_dev, 0, 1, linestyle=':', color='black', linewidth=1)
    ax.hlines(obs_mean + obs_plausible, 0, 1, linestyle=':', color='0.5', linewidth=1)
    ax.hlines(obs_mean - obs_plausible, 0, 1, linestyle=':', color='0.5', linewidth=1)

    # Tidying the figure
    ax.set_xticks([])
    ax.set_ylabel('dSIA/dGMST [million km$^2$ K$^{-1}$]')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save the figure (also closes it)
    # TODO: work out what provedence file is about. Not the dictionary!
    provenance_record = get_provenance_record(dictionary, ancestor_files=[])
    save_figure(filename, provenance_record, cfg, figure=fig, close=True)


def main(cfg):
    # Get the data from the cfg
    logger.info('Getting data from the config')
    input_data = cfg['input_data'].values()

    # Initialize blank dictionary to send to plot later
    sensitivity_dict = {}

    # Get list of datasets from cfg
    logger.info('Listing datasets in the data')
    datasets = list_datasets(input_data)

    # Iterate over each dataset
    for dataset in datasets:

        # Select only data from that dataset
        logger.info(f'Selecting data from {dataset}')
        selection = select_metadata(input_data, dataset=dataset)

        # Calculate the sensitivity
        logger.info(f'Calculating sensitivity for {dataset}')
        sensitivity = calculate_sensitivity(selection)

        # Add the sensitivity to the dictionary
        sensitivity_dict[dataset] = sensitivity

    # Plot the sensitivities (and save and close the plot)
    logger.info(f'Creating plot')
    plot_from_dict(sensitivity_dict, 'Sea_ice_sensitivity', cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
