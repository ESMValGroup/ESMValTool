"""Visualize temperature accross datasets.

- Make a plot of the global mean temperature change
according to all datasets (defined above),

- Get the global mean temperature change for specified years
and specified percentiles (Delta T). These define our scenarios.

- Select the 30-year period from the target model (all ensemble members)
where they match the Delta T for each scenario.
"""

from itertools import product
from datetime import timedelta
import matplotlib.pyplot as plt
import iris
import xarray as xr

from esmvaltool.diag_scripts.shared import (run_diagnostic,
                                            get_plot_filename,
                                            select_metadata)


def get_target_delta_t(metadata, year, percentile):
    """Compute target delta T for KNMI scenarios."""
    # Read the 10th percentile dataset from the preprocessed data
    # Return the value for the scenario year.
    attribute = f'MultiModel{percentile}'
    multimodelstat = select_metadata(metadata, alias=attribute)
    cube = iris.load_cube(multimodelstat[0]['filename'])
    time = iris.Constraint(time=lambda cell: cell.point.year == year)
    return cube.extract(time).data.data


def get_mean_ensemble(metadata, model):
    """Compute average of ensembles for target model."""
    target_dataset = select_metadata(metadata, dataset=model)
    sum_cube = 0
    for ensemble in target_dataset:
        cube = iris.load_cube(ensemble['filename'])
        sum_cube += cube
    return sum_cube / len(target_dataset)


def get_target_model_time_bounds(target_dataset, delta_t, span=30):
    """Return 30-year time bounds.

    Get the 30-year time bounds for which the target model delta T
    matches the target delta T for a specific scenario.
    """
    # Compute rolling average delta T of the target dataset?!
    # Possibly, apply/report a (pattern) scaling factor
    # in order to match the target delta T more closely.

    cube = target_dataset - delta_t
    min_cube = cube.collapsed('time', iris.analysis.MIN)
    coord = min_cube.coord('time')
    start = coord.cell(0).point - timedelta(span / 2 * 365.24)
    end = coord.cell(0).point + timedelta(span / 2 * 365.24)
    return [start.year, end.year], min_cube.data


def make_plot(metadata, scenario):
    """Make figure 3, left graph.

    Multimodel values as line, reference value in black square,
    steering variables in dark dots.
    """
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


def main(cfg):
    """Return scenarios tables."""
    # a list of dictionaries describing all datasets passed on to the recipe
    metadata = cfg['input_data'].values()
    print(metadata)
    # Define the different scenario's
    scenarios = []
    combinations = product(cfg['scenario_years'], cfg['scenario_percentiles'])
    for year, percentile in combinations:
        cmip_dt = get_target_delta_t(metadata, year, percentile)
        scenario = {
            'year': year,
            'percentile': percentile,
            'CMIP_delta_T_global': cmip_dt
        }
        scenarios.append(scenario)
    print(scenarios)
    # find all ensemble members (datasets) of the target_model
    # average them TODO check for ensembles in different rcp
    average = get_mean_ensemble(metadata, cfg['target_model'])

    # use the average temperature change to get the time bounds
    for scenario in scenarios:
        cmip_dt = scenario['CMIP_delta_T_global']
        period_bounds, target_dt = get_target_model_time_bounds(
            average, scenario['CMIP_delta_T_global']
        )
        scenario['period_bounds'] = period_bounds
        scenario['target_model_delta_T'] = target_dt
        scenario['pattern_scaling_factor'] = cmip_dt / target_dt
    print(scenarios)

    # TODO Make a nice plot
    # if cfg['write_plots']:
    #     make_plot()
    # TODO add provenance


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
