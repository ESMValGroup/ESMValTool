"""Resampling.

In the second step (see sections 4.2 and 4.3) we resample the results of
ENS-EC for the selected time periods (from the previous step). This is done
by recombining 5 yr periods of the eight members of ENSEC into new resampled
climates, and selecting combinations that match with the spread in CMIP5.
This provides eight resampled EC-Earth time series for each of the scenarios.

"""
import logging
from pathlib import Path

from itertools import product

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np

from esmvaltool.diag_scripts.shared import (run_diagnostic,
                                            select_metadata,
                                            get_diagnostic_filename,
                                            ProvenanceLogger,
                                            get_plot_filename)

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption':
        "Resampling of local climate model.",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            'alidoost_sarah',
            # 'rol_evert',
        ],
        'ancestors': ancestor_files,
    }
    return record


def get_input_data(cfg):
    """Load ensembles for each variable and merge."""
    dataset_dicts = cfg['input_data'].values()
    target_model_metadata = select_metadata(
        dataset_dicts, dataset=cfg['target_model']
    )
    dataset = []
    ancestor_files = []
    for short_name in "pr", "tas":
        group = f'{short_name}_target'
        var = select_metadata(target_model_metadata, variable_group=group)
        files = [metadata['filename'] for metadata in var]
        dataset.append(
            xr.open_mfdataset(files, concat_dim='ensemble_member',
                              combine='nested')
        )
        ancestor_files.extend(files)
    provenance = create_provenance_record(ancestor_files)
    return xr.merge(dataset), provenance


def get_season_means_per_segment(dataset, period, step=5):
    """Segment (part of) a dataset into x blocks of 5 years each.

    Returns a new dataset with additional coordinate 'segment',
    and values as season means.
    """
    segments = []
    for year in range(*period, step):
        segments.append(
            dataset.sel(time=slice(str(year), str(year + step - 1)))
        )
    segmented_dataset = xr.concat(segments, dim='segment')
    season_means = segmented_dataset.groupby('time.season').mean()
    return season_means


def all_possible_combinations(n_ensemble_members, n_segments):
    """Generate indexer.

    Returns indices for all possible combinations of
    segments and ensembles.
    """
    # Create a DataArray once...
    indices = xr.DataArray(
        data=np.zeros(n_segments),
        dims=['segment'],
        coords={'segment': np.arange(n_segments)})

    # ...and update its values for all possible combinations
    for combination in product(range(n_ensemble_members), repeat=n_segments):
        indices.values = list(combination)
        yield indices


def selected_combinations(combinations):
    """Generate indexer for all selected combinations of segmented dataset.

    combinations: a list of combinations
    """
    n_segments = len(combinations[0])

    # Create a DataArray once...
    indices = xr.DataArray(
        data=np.zeros(n_segments),
        dims=['segment'],
        coords={'segment': np.arange(n_segments)})

    # ...and update its values for each selected combinations
    for combination in combinations:
        indices.values = combination
        yield indices


def most_promising_combinations(
        segment_means,
        target,
        n_ensemble_members,
        n_combination=1000):
    """Select n_combination that are closest to the target.

    First, get all possible combinations, then select n_combination.
    """
    promising_combinations = []
    for combination in all_possible_combinations(
            n_ensemble_members, n_segments=6):
        recombined_segment_means = segment_means.sel(
            ensemble_member=combination)
        new_overall_mean = recombined_segment_means.mean('segment')
        distance_to_target = abs(new_overall_mean - target)
        promising_combinations.append([combination.values, distance_to_target])
    top_all = pd.DataFrame(promising_combinations,
                           columns=['combination', 'distance_to_target'])
    return top_all.sort_values('distance_to_target').head(n_combination)


def seasonal_mean_characteristics(top1000, segment_means):
    """Compute summer pr,and summer and winter tas.

    Computation is done for each set of 1000 samples from 1a-b.
    """
    filter2_variables = []
    columns = ['combination', 'pr_summer', 'tas_winter', 'tas_summer']
    for combination in selected_combinations(top1000):
        recombined_segment_means = segment_means.sel(
            ensemble_member=combination)
        new_overall_season_means = recombined_segment_means.mean('segment')

        filter2_variables.append([
            combination.values,
            new_overall_season_means.pr.sel(season='JJA'),
            new_overall_season_means.tas.sel(season='DJF'),
            new_overall_season_means.tas.sel(season='JJA')
        ])
    return pd.DataFrame(filter2_variables, columns=columns)


def within_bounds(values, bounds):
    """Return true if value is within percentile bounds."""
    low, high = np.percentile(values, bounds)
    return values.between(low, high)


def select_percentile_bounds(top1000, info, period):
    """Select samples based on the percentile bounds.

    Select samples for which summer pr, and summer and winter tas
    are within the percentile bounds specified in the recipe.
    """
    pr_summer = within_bounds(
        top1000['pr_summer'], info[f'pr_summer_{period}'])
    tas_winter = within_bounds(
        top1000['tas_winter'], info[f'tas_winter_{period}'])
    tas_summer = within_bounds(
        top1000['tas_summer'], info[f'tas_summer_{period}'])
    subset = top1000[
        np.all([pr_summer, tas_winter, tas_summer], axis=0)
    ]
    return subset


def determine_penalties(overlap):
    """Determine penalties dependent on the number of overlaps."""
    return np.piecewise(
        overlap,
        condlist=[overlap < 3, overlap == 3, overlap == 4, overlap > 4],
        funclist=[0, 1, 5, 100])


def select_final_subset(combinations, n_sample=8):
    """Find n samples with minimal reuse of ensemble members per segment.

    combinations: a pandas series with the remaining candidates
    n: the final number of samples drawn from the remaining set.
    """
    # Convert series of tuples to 2d numpy array (much faster!)
    combinations = np.array(
        [list(combination) for combination in combinations]
    )

    # Store the indices in a nice dataframe
    _, n_segments = combinations.shape
    best_combination = pd.DataFrame(
        data=None,
        columns=[f'Segment {x}' for x in range(n_segments)],
        index=[f'Combination {x}' for x in range(n_sample)]
    )

    # Random number generator
    rng = np.random.default_rng()

    lowest_penalty = 500  # just a random high value
    for _ in range(10000):
        sample = rng.choice(combinations, size=n_sample)
        penalty = 0
        for segment in sample.T:
            _, counts = np.unique(segment, return_counts=True)
            penalty += determine_penalties(counts).sum()
        if penalty < lowest_penalty:
            lowest_penalty = penalty
            best_combination.loc[:, :] = sample

    return best_combination


def make_scenario_tables(dataset):
    """Make a table of final set of eight samples.

    Selection is done based on minimal reuse of the same
    ensemble member for the same period.
    """
    scenario_output_tables = []
    for period in ['control', 'future']:
        remaining_combinations = dataset[period].combination
        result = select_final_subset(remaining_combinations)
        scenario_output_tables.append(result)

    scenario_output_tables = pd.concat(
        scenario_output_tables, axis=1, keys=['control', 'future'])
    return scenario_output_tables


def save_scenario_tables(scenario, tables, cfg, provenance):
    """Save the scenario tables as a csv file."""
    basename = f'indices_{scenario}'
    output_file = get_diagnostic_filename(basename, cfg, extension='csv')
    tables.to_csv(output_file)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)


def make_plot(cfg):
    """Reproduce figure 5 from the paper."""
    # note that quantile is applied twice!
    # Once to get the pdf's of seasonal tas/pr
    # and once to get the multimodel pdf of the quantile changes
    # TODO: Add variable group tas_cmip to recipe:
    # variables:
    #   tas_cmip:
    #     short_name: tas
    #     preprocessor: preproc2 (extract NL point)
    #     additional_datasets: *cmip5_datasets

    metadata = cfg['input_data'].values()
    tas_cmip = select_metadata(metadata, variable_group='tas_cmip')
    envelope = []
    for data_dict in tas_cmip:
        dataset = xr.open_dataset(data_dict['filename'])
        grouped_control = (dataset
                           .sel(time=slice('1981', '2010'))
                           .groupby('time.season'))
        grouped_future = (dataset
                          .sel(time='2050')
                          .groupby('time.season'))

        quantiles = [.05, .1, .25, .5, .75, .90, .95]
        qcontrol = grouped_control.quantile(quantiles)
        qfuture = grouped_future.quantile(quantiles)
        qchange = qfuture - qcontrol
        envelope.append(qchange)
    cmip5 = xr.concat(envelope, dim='multimodel')

    fig, ax = plt.subplots()
    for high, low in [[0.9, 0.1], [0.75, 0.25]]:
        ax.fill_between(quantiles,
                        cmip5.quantile(high, dim='multimodel'),
                        cmip5.quantile(low, dim='multimodel'),
                        color='k',
                        alpha=0.5)
    filename = get_plot_filename('envelope_figure', cfg)
    fig.savefig(filename, bbox_inches='tight', dpi=300)
    return filename


def main(cfg):
    """Resample the model of interest.

    Step 0: Read the data, extract segmented subsets for both the control
    and future periods, and precompute seasonal means for each segment.

    Step 1a: Get 1000 combinations for the control period.
    These samples should have the same mean winter precipitation as the
    overall mean of the x ensemble members.

    Step 1b: Get 1000 combinations for the future period.
    Now, for future period, the target value is a relative change
    with respect to the overall mean of the control period.

    Step 2a: For each set of 1000 samples from 1a-b, compute
    summer precipitation, and summer and winter temperature.

    Step 2b: For each scenario, select samples for which
    summer precipitation, and summer and winter temperature are
    within the percentile bounds specified in the recipe.

    Step 3: Select final set of eight samples with minimal reuse of the same
    ensemble member for the same period.
    From 10.000 randomly selected sets of 8 samples, count
    and penalize re-used segments (1 for 3*reuse, 5 for 4*reuse).
    Chose the set with the lowest penalty.
    """
    # Step 0:
    # Read the data
    dataset, provenance = get_input_data(cfg)

    # Precompute the segment season means for the control period ...
    segments_season_means = {}
    segments_season_means['control'] = get_season_means_per_segment(
        dataset, cfg['control_period'], step=5
    )
    # ... and for all the future scenarios
    for scenario, info in cfg['scenarios'].items():
        segments_season_means[scenario] = get_season_means_per_segment(
            dataset, info['resampling_period'], step=5
        )

    # Step 1:
    # Create a dictionary to store the selected indices later on.
    selected_indices = {}

    # Step 1a: Get 1000 combinations for the control period
    # Get mean winter precipitation per segment
    winter_mean_pr_segments = segments_season_means['control'].pr.sel(
        season='DJF')
    # Get overall mean of segments
    control_period_winter_mean_pr = winter_mean_pr_segments.mean()
    # Select top1000 values
    selected_indices['control'] = most_promising_combinations(
        winter_mean_pr_segments, control_period_winter_mean_pr,
        len(dataset.ensemble_member)
    )

    # Step 1b: Get 1000 combinations for the future period
    for scenario, info in cfg['scenarios'].items():
        # The target value is a relative change with respect to
        # the overall mean of the control period
        target = control_period_winter_mean_pr * (1 + info['dpr_winter'] / 100)
        winter_mean_pr_segments = segments_season_means[scenario].pr.sel(
            season='DJF')
        selected_indices[scenario] = most_promising_combinations(
            winter_mean_pr_segments, target, len(dataset.ensemble_member)
        )

    # Step 2a: For each set of 1000 samples from 1a-b, compute summer pr,
    # and summer and winter tas.
    for name, dataframe in selected_indices.items():
        segment_means = segments_season_means[name]
        selected_indices[name] = seasonal_mean_characteristics(
            dataframe.combination, segment_means
        )

    # Step 2b: For each scenario, select samples for which summer pr, and
    # summer and winter tas are within the percentile bounds specified
    # in the recipe
    for scenario, info in cfg['scenarios'].items():
        selected_indices[scenario] = {
            # Get relatively high/low values for the control period ...
            'control': select_percentile_bounds(
                selected_indices['control'], info, 'control'
            ),
            # Get relatively high/low values for the future period
            'future': select_percentile_bounds(
                selected_indices[scenario], info, 'future'
            )
        }

        del selected_indices['control']  # No longer needed

    # Step 3: Select final set of eight samples with minimal reuse of the same
    # ensemble member for the same period.
    for scenario, dataframes in selected_indices.items():
        scenario_output_tables = make_scenario_tables(dataframes)
        print(f"Selected recombinations for scenario {scenario}:")
        print(scenario_output_tables)
        save_scenario_tables(scenario, scenario_output_tables, cfg, provenance)

    # Step 4: Plot the results
    if cfg['write_plots']:
        make_plot(cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
