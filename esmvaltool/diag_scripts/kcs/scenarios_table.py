"""
Visualize temperature accross datasets

- Make a plot of the global mean temperature change according to all datasets (defined above)
- Get the global mean temperature change for specified years and specified percentiles (Delta T). These define our scenarios.
- Select the 30-year period from the target model (all ensemble members) where they match the Delta T for each scenario.
"""

from itertools import combinations
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr


from esmvaltool.diag_scripts.shared import (run_diagnostic, get_plot_filename,
                                            get_diagnostic_filename)


def get_target_delta_T(cfg, year):
    """Compute target delta T for KNMI scenarios."""
    # Read the 10th percentile dataset from the preprocessed data
    # Return the value for the scenario year.
    return multimodelstats[year]


def get_target_model_time_bounds(target_dataset, delta_T):
    """Return 30-year time bounds.

    Get the 30-year time bounds for which the target model delta T
    matches the target delta T for a specific scenario.
    """
    # Compute rolling average delta T of the target dataset
    # Get the year for which this delta T most closely matches the target delta T
    # Get the start- and end year belonging to that 30-year period
    # Possibly, apply/report a (pattern) scaling factor in order to match the target delta T more closely.


    return start, end, target_dT


def make_plot():
    # Create a figure and a list that we'll use to store output data to csv
    fig, ax = plt.subplots()

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

    return



def main():

    # a list of dictionaries describing all datasets passed on to the recipe
    metadata = cfg['input_data'].values()

    # Define the different scenario's
    scenarios = []
    for year, percentile in combinations(cfg['scenarios_years'], cfg['scenarios_percentiles']):
        delta_T = get_target_delta_T
        scenario = {
            'year': year,
            'percentile': percentile,
            'CMIP_delta_T_global': delta_T
            }
        scenarios.append(scenario)

    print(scenarios)

    target_datasets = select_metadata(metadata, dataset.contains(cfg['target_model']))

    # find all ensemble members (datasets) of the target_model
    # average them
    # use the average temperature change to get the time bounds
    for scenario in scenarios:
        cmip_dT = scenario['CMIP_delta_T_global']
        period_bounds, target_dT = get_target_model_time_bounds(average, scenario['CMIP_delta_T_global'])
        scenario['period_bounds'] = period_bounds
        scenario['target_model_delta_T'] = target_dT
        scenario['pattern_scaling_factor'] = cmip_dT / target_dT


    # Make a nice plot
    if cfg['write_plots'] == True:
        make_plot()

    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
