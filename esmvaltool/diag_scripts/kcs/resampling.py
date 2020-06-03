
"""
- Divide the 30-year dataset into 5-year blocks
- Create all possible combinations out of these 30/5 = 6 periods and x ensemble members (may contain the same block multiple times, but less is better and maximum three times (see third item below))
- Determine the 1000 best ...
- Determine the final best

From the paper:
1. We first select 1000 re-samples (out of the total 86) that have a change in mean winter precipitation of 4% (for the 'L' scenarios) and 8% (for the 'H' scenarios) per degree rise in global mean temperature. These two rates of the precipitation change are primarily based results with EC-Earth shown in figure 1. They are also consistent with expectations related to the increased moisture of the atmosphere following from the Clausius–Clapeyron relation (Trenberth et al 2003). This would imply an increase of 6–7% per degree temperature rise, considering that the local temperature rise is approximately equal to global temperature rise and that changes in relative humidity are small (see also Lenderink and Attema, this issue).
2. From the subset of 1000 re-samples, we select approximately 50 re-samples based on changes in summer precipitation, and changes in summer and winter temperature (i.e. targets 3–5 in table 2). This is not done by imposing an absolute constraint (like in step 1), but by selecting a percentage range from the remaining re-samples (see supplementary material). These percentile ranges are chosen in an iterative way doing the evaluation with the CMIP5 range (as discussed in the next section) a number of times.
3. From approximately 50 re-samples obtained above, we finally sub-select eight re-samples that have a minimal re-use of the same model data. It is e.g. not allowed that a 5 yr period from an ensemble member is re-used more than three times in the set of eight re-samples that form a scenario.
"""

from itertools import combinations
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.shared import (run_diagnostic, get_plot_filename,
                                            get_diagnostic_filename, select_metadata)


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


def filter3(df, n=8):
    """Sort by the number of unique elements in 'combination' and return the top n."""
    # This would take the lowest penalty directly for each resampled version, but
    # instead we need to select the (random) combination of 8 samples that has the
    # overall lowest re-use of individual members.
    df['n_datasets'] = df.combination.map(lambda x: len(set(x)))
    return df.sort('n_datasets', descending=True).head(n)



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
                new_overall_season_means.pr.sel(season='JJA')
                new_overall_season_means.tas.sel(season='DJF')
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

        del selected_indices['control']  # No longer needed; we now have a control subset for each scenario


    # Step 3:
   """Calculate the subset S3: find a subset with the least re-use of
    segments, using random sampling"""

    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
