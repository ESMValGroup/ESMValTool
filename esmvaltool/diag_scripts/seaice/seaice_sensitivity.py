import logging
from pathlib import Path
from pprint import pformat
import iris
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
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
        "caption": 'I CHANGED THIS TO A STRING caption',  # TODO: remember edited
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


def extract_cube(data, variable_group):  # TODO: should I tweak the sea ice cube here? Or in preprocessing?
    """
    Selects records from data based on the variable_group,
    asserts that there is only one matching record,
    then returns an iris cube for that record
    """
    logger.debug(f'extracting {variable_group} cube from {data}')

    # Load any data in the variable_group
    selection = select_metadata(data, variable_group=variable_group)

    # Ensure there is only one file in the list
    assert len(selection) == 1, f'None or too many matching files found for {variable_group}'

    # Load the cube, [0] is because selection returns a list
    cube = iris.load_cube(selection[0]['filename'])

    # Amend units if it's a sea ice cube
    if cube.var_name == 'siconc':
        cube.convert_units('1e6 km2')  # to match literature

    # Add years if it hasn't been annually averaged
    if 'year' not in [coord.name() for coord in cube.coords()]:
        cc.add_year(cube, 'time', name='year')

    return cube


def calculate_slope(independent, dependent):
    """
    Use SciPy stats to calculate the least-squares regression
    """
    logger.debug(f'Calculating linear relationship between {dependent} and {independent}')

    # Use SciPy stats to calculate the regression
    slope, intercept, r, p, stderr = stats.linregress(independent, dependent)

    # Return ony the slope
    return slope


def calculate_rvalue(independent, dependent):
    """
    Use SciPy stats to calculate the Pearson correlation coefficient
    """
    logger.debug(f'Calculating rvalue between {dependent} and {independent}')

    # Use SciPy stats to calculate the regression
    slope, intercept, r, p, stderr = stats.linregress(independent, dependent)

    # Return ony the rvalue
    return r


def calculate_annual_trends(data):
    """
    Calculates ...
    """
    logger.debug('calculating annual trends')

    # Load the preprocessed cubes
    ann_si_cube = extract_cube(data, 'avg_ann_arctic_sea_ice')
    ann_tas_cube = extract_cube(data, 'avg_ann_global_temp')

    # Check that the years match between the two cubes
    assert (ann_si_cube.coord('year').points == ann_tas_cube.coord('year').points).all()
    years = ann_si_cube.coord('year').points

    # Calculate the trends over time
    si_trend = calculate_slope(years, ann_si_cube.data)
    tas_trend = calculate_slope(years, ann_tas_cube.data)

    # Calculate the direct correlation coefficient
    direct_r_val = calculate_rvalue(ann_tas_cube.data, ann_si_cube.data)

    dictionary = {
        'si_ann_trend': si_trend,
        'tas_ann_trend': tas_trend,
        'direct_r_val': direct_r_val
    }

    return dictionary


def calculate_sept_sensitivity(data):
    """
    Calculates change in sea ice area per year
    divided by change in global mean temperature per year
    """
    logger.debug('calculating sensitivity')

    # Load the preprocessed cubes
    si_cube = extract_cube(data, 'arctic_sept_sea_ice')
    tas_cube = extract_cube(data, 'avg_ann_global_temp')

    # Check that the years match between the two cubes
    assert (si_cube.coord('year').points == tas_cube.coord('year').points).all()

    # Calculate the sensitivity, NOT via time, so unlike Ed Blockley's code
    sensitivity = calculate_slope(tas_cube.data, si_cube.data)

    return sensitivity


def notz_style_plot_from_dict(dictionary, filename, cfg):
    """
    Saves a plot of sensitivities and observations to the given filename
    """
    # Read from Ed Blockley's dictionary
    obs_years = list(observations)[0]
    obs_mean = observations[obs_years]['mean']
    obs_std_dev = observations[obs_years]['std_dev']
    obs_plausible = observations[obs_years]['plausible']

    # Set up the figure
    fig, ax = plt.subplots(figsize=(4, 6), layout='constrained')
    fig.suptitle('September Arctic Sea Ice Sensitivity')
    ax.set_title(f'dSIA/dGMST {obs_years}')  # Ed's title

    # Iterate over the dictionary
    for dataset, inner_dict in dictionary.items():
        ax.plot(0.5, inner_dict['direct_sept_sensitivity'], label=dataset, marker='_', markersize=20)

    # Add observations (style taken from Ed's code)
    ax.hlines(obs_mean, 0, 1, linestyle='--', color='black', linewidth=2)
    ax.fill_between([0, 1], obs_mean - obs_std_dev, obs_mean + obs_std_dev, facecolor='k', alpha=0.15)
    ax.hlines(obs_mean + obs_plausible, 0, 1, linestyle=':', color='0.5', linewidth=1)
    ax.hlines(obs_mean - obs_plausible, 0, 1, linestyle=':', color='0.5', linewidth=1)

    # Tidy the figure
    ax.set_xticks([])
    ax.set_ylabel('dSIA/dGMST (million km$^2$ K$^{-1}$)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save the figure (also closes it)
    # TODO: work out what provenance file is about. Not the dictionary!
    provenance_record = get_provenance_record(dictionary, ancestor_files=[])
    save_figure(filename, provenance_record, cfg, figure=fig, close=True)


# TODO: the published plot uses a trend of K/decade, not per year, and also 1979-2018
def roach_style_plot_from_dict(dictionary, filename, cfg):
    """
    Saves a plot of trend in SIA against trend in GMST to the given filename
    """
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
    fig.suptitle('Trends in Annual Mean Temperature And Arctic Sea Ice')
    ax.set_title('Decades are scaled up from annual data')

    # Set up for colouring the points
    norm = Normalize(vmin=-1, vmax=1)
    cmap = plt.get_cmap('PiYG_r')

    # Set up the axes
    ax.set_xlim(-0.02, 0.505)
    ax.set_ylim(-1.1, 0.21)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_yticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2])
    ax.hlines(0, -0.02, 0.505, color='black', alpha=0.5)
    ax.vlines(0, -1.1, 0.21, color='black', alpha=0.5)
    ax.set_xlabel('Trend in GMST (K decade$^{-1}$)')
    ax.set_ylabel('Trend in SIA (million km$^2$ decade$^{-1}$)')

    # Iterate over the dictionary
    for dataset, inner_dict in dictionary.items():
        x = 10 * inner_dict['tas_trend']  # x10 to approximate decades
        y = 10 * inner_dict['siconc_trend']  # x10 to approximate decades
        r = inner_dict['direct_r_val']**2 * np.sign(inner_dict['direct_r_val'])  # TODO: check this is right
        plt.scatter(x, y, marker='o', s=100, c=[r], cmap=cmap, norm=norm)
        plt.annotate(dataset, xy=(x, y), xytext=(x+0.02, y+0.02))

    # Add a colour bar
    plt.colorbar(label='r value')

    # Save the figure (also closes it)
    # TODO: provenance queries as above
    provenance_record = get_provenance_record(dictionary, ancestor_files=[])
    save_figure(filename, provenance_record, cfg, figure=fig, close=True)


def main(cfg):
    # Get the data from the cfg
    logger.info('Getting data from the config')
    input_data = cfg['input_data'].values()

    # Initialize blank dictionary to send to plot later
    trends_dict = {}

    # Get list of datasets from cfg
    logger.info('Listing datasets in the data')
    datasets = list_datasets(input_data)

    # Iterate over each dataset
    for dataset in datasets:

        # Select only data from that dataset
        logger.info(f'Selecting data from {dataset}')
        selection = select_metadata(input_data, dataset=dataset)

        # Add the dataset to the dictionary with a blank inner dictionary
        trends_dict[dataset] = {}

        # Calculations for the Notz-style plot
        logger.info('Calculating data for Notz-style plot')
        sensitivity = calculate_sept_sensitivity(selection)
        # Add to dictionary
        trends_dict[dataset]['direct_sept_sensitivity'] = sensitivity

        # Calculations for the Roach-style plot
        logger.info('Calculating data for Roach-style plot')
        trends = calculate_annual_trends(input_data)
        # Add to dictionary
        trends_dict[dataset]['siconc_trend'] = trends['si_ann_trend']
        trends_dict[dataset]['tas_trend'] = trends['tas_ann_trend']
        trends_dict[dataset]['direct_r_val'] = trends['direct_r_val']

    # Plot the sensitivities (and save and close the plot)
    logger.info(f'Creating Notz-style plot')
    notz_style_plot_from_dict(trends_dict, 'Notz-style_sea_ice_sensitivity', cfg)
    logger.info(f'Creating Roach-style plot')
    roach_style_plot_from_dict(trends_dict, 'Roach-style_avg_annual_trends', cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
