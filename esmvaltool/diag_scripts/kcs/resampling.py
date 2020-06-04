
"""
In the second step (see sections 4.2 and 4.3) we resample the results of
ENS-EC for the selected time periods (from the previous step). This is done
by recombining 5 yr periods of the eight members of ENSEC into new resampled
climates, and selecting combinations that match with the spread in CMIP5.
This provides eight resampled EC-Earth time series for each of the scenarios.
"""

from itertools import combinations

import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.shared import (get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata)


def segment(dataset, start_year, end_year, block_size):
    """Segment (part of) a dataset into x blocks of n years each.

    Returns a new dataset with additional coordinate 'segment'.
    """
    segmented_dataset = xr.concat([
        ds.sel(time = slice(str(year), str(year+block_size)))
        for year in range(start_year, end_year, block_size)
    ], dim='segment')
    return segmented_dataset


def all_possible_combinations(n_ensemble_members, n_segments):
    """Generate indexer for all possible combinations of segmented dataset."""

    # Create a DataArray once...
    indices = xr.DataArray(
        data = np.zeros(n_segments),
        dims = ['segment'],
        coords = {'segment': np.arange(n_segments)})

    # ...and update its values for all possible combinations
    for combination in product(range(n_ensemble_members), repeat=n_segments):
        indices.values = list(update)
        yield indices


def selected_combinations(combinations):
    """Generate indexer for all selected combinations of segmented dataset.

    combinations: a list of combinations
    """
    n_segments = len(combinations[0])

    # Create a DataArray once...
    indices = xr.DataArray(
        data = np.zeros(n_segments),
        dims = ['segment'],
        coords = {'segment': np.arange(n_segments)})

    # ...and update its values for each selected combinations
    for combination in combinations:
        indices.values = combination
        yield indices


def most_promising_combinations(segment_means, target, n=1000):
    """Get n out of all possible combinations that are closest to the target."""
    promising_combinations = []
    for combination in all_possible_combinations(n_ensemble_members, n_segments=6):
        recombined_segment_means = segment_means.sel(ensemble_member=combination)
        new_overall_mean = recombined_segment_means.mean('segment')
        distance_to_target = (new_overall_mean - target).abs()
        promising_combinations.append([combination.values, distance_to_target])
    top_all = pd.DataFrame(promising_combinations, columns=['combination', 'distance_to_target'])
    return top_all.sort_by('distance_to_target').head(n)


def within_bounds(values, bounds):
        """Return true if value is within percentile bounds."""
        low, high = np.percentile(values, bounds)
        return (low <= values) & (values <= high)


def determine_penalties(overlap):
    """Determine penalties dependent on the number of overlaps."""
    return np.piecwise(x = overlap,
        condlist = [x<3, x==3, x==4, x>4],
        funclist = [0, 1, 5, 100])


def select_final_subset(combinations, n=8):
    """Find n samples with minimal reuse of ensemble members per segment.

    combinations: a pandas series with the remaining candidates
    n: the final number of samples drawn from the remaining set.
    """
    # Convert series of tuples to 2d numpy array (much faster!)
    combinations = np.array(
        [list(combination) for combination in combinations]
    )

    # Random number generator
    rng = np.random.default_rng()

    lowest_penalty = 500  # just a random high value
    for i in range(10000):
        sample = rng.choice(combinations, size=n)
        penalty = 0
        for segment in sample.T:
            _, counts = np.unique(segment, return_counts=True)
            penalty += determine_penalties(counts).sum()
        if penalty < lowest_penalty:
            lowest_penalty = penalty
            best_combination = sample

    # Store the indices in a nice dataframe
    _, n_segments = best_combination.shape
    best_combination = pd.DataFrame(
        data = best_combination,
        columns = [f'Segment {x}' for x in range(n_segments)],
        index = [f'Combination {x}' for x in range(n)]
    )

    return best_combination


def main(cfg):

    # Step 0: Read the data, extract segmented subsets for both the control
    # and future periods, and precompute seasonal means for each segment.
    dataset_dicts = cfg['input_data'].values()
    target_model_metadata = select_metadata(dataset_dicts, dataset=cfg['target_model'])
    files = [metadata['filename'] for metadata in target_model_metadata]
    dataset = xr.open_mfdataset(files, concat_dim='ensemble_member', combine='nested')[['pr', 'tas']]
    segments_season_means = {}

    # Precompute the segment season means for the control period ...
    segments = segment(dataset, *cfg['control_period'], step=5)
    season_means = segments.groupby('time.season').mean()
    segments_season_means['control'] = season_means

    # ... and for all the future scenarios
    for name, info, in cfg['scenarios']:
        segments = segment(dataset, *info['resampling_period'], step=5)
        season_means = segments.groupby('time.season').mean()
        segments_season_means[name] = season_means

    # Create a dictionary to store the selected indices later on.
    selected_indices = {}


    # Step 1a: Get 1000 combinations for the control period
    # These samples should have the same mean winter precipitation as the
    # overall mean of the x ensemble members
    winter_mean_pr_segments = segments_season_means['control'].pr.sel(season='DJF')
    target = winter_mean_pr_segments.mean()
    top1000 = most_promising_combinations(winter_mean_pr_segments, target)
    selected_indices['control'] = top1000

    # Step 1b: Get 1000 combinations for the future period
    # The target value is a relative change wrt the overall mean of the control period
    control_period_winter_mean_pr = target

    for scenario, info, in cfg['scenarios']:
        target = control_period_winter_mean_pr * (1 + info['dpr_winter']/100)
        winter_mean_pr_segments = segments_season_means[scenario].pr.sel(season='DJF')
        top1000 = most_promising_combinations(winter_mean_pr_segments, target)
        selected_indices[scenario] = top1000


    # Step 2a: For each set of 1000 samples from 1a-b, compute summer pr,
    # and summer and winter tas.
    for name, dataframe in selected_indices.items():
        top1000 = dataframe.combination
        segment_means = segments_season_means[name]

        filter2_variables = []
        columns = ['combination', 'pr_summer', 'tas_winter', 'tas_summer']
        for combination in selected_combinations(top1000):
            recombined_segment_means = segment_means.sel(ensemble_member=combination)
            new_overall_season_means = recombined_segment_means.mean('segment')

            filter2_variables.append([
                combination.values,
                new_overall_season_means.pr.sel(season='JJA'),
                new_overall_season_means.tas.sel(season='DJF'),
                new_overall_season_means.tas.sel(season='JJA')
                ])

        top1000 = pd.DataFrame(filter2_variables, columns=columns)
        selected_indices[scenario] = top1000

    # Step 2b: For each scenario, select samples for which summer pr,and
    # summer and winter tas are within the percentile bounds specified
    # in the recipe
    top1000_control = selected_indices['control']
    for scenario, info, in cfg['scenarios']:
        top1000 = selected_indices[scenario]

        # Get relatively high/low values for the control period ...
        pr_summer = within_bounds(top1000_control['pr_summer'], info['pr_summer_control'])
        tas_winter = within_bounds(top1000_control['tas_winter'], info['tas_winter_control'])
        tas_summer = within_bounds(top1000_control['tas_summer'], info['tas_summer_control'])
        subset_control = top1000_control[np.all(pr_summer, tas_winter, tas_summer)]

        # ... combined with relatively high/low values for the future period
        pr_summer = within_bounds(top1000['pr_summer'], info['pr_summer_future'])
        tas_winter = within_bounds(top1000['tas_winter'], info['tas_winter_future'])
        tas_summer = within_bounds(top1000['tas_summer'], info['tas_summer_future'])
        subset_future = top1000[np.all(pr_summer, tas_winter, tas_summer)]

        selected_indices[scenario] = {
            'control': subset_control,
            'future': subset_future
        }

        del selected_indices['control']  # No longer needed


    # Step 3:
    # Select final set of eight samples with minimal reuse of the same
    # ensemble member for the same period.
    #
    # From 10.000 randomly selected sets of 8 samples, count
    # and penalize re-used segments (1 for 3*reuse, 5 for 4*reuse).
    # Chose the set with the lowest penalty.
    for scenario, dataframes in selected_indices.items():
        scenario_output_tables = []
        for period in ['control', 'future']:
            remaining_combinations = dataframes[period].combination
            result = select_final_subset(remaining_combinations)
            scenario_output_tables.append(result)

        scenario_output_tables = pd.concat(scenario_output_tables, axis=1, keys=['control', 'future'])
        print(f"Selected recombinations for scenario {scenario}:")
        print(scenario_output_tables)
        filename = f'indices_{scenario}'
        scenario_output_tables.to_csv(filename)

    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
