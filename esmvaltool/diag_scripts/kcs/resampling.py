
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


def filter1(pr_diffs_winter, target_dpr):
    """Get the 1000 combinations that are closest to the target.

    pr_diffs_winter: a pandas dataframe with two columns:
        1. the code for the combination
        2. the winter precipitation differences

    targed_dpr: target precipitation change (dependent on the scenario)
    """

    top1000 (pr_diffs_winter - target_dpr).sort('pr_diff_winter').head(1000)
    return top1000


def filter2():
    # These are needed only for the 1000 best samples from the step above, so may do it later, outside this superloop
    # pr_diff_summer = get_precipitation_change(resampled_future, resample_control, variable='tas', season='JJA') # Maybe use esmvaltool's seasonal statistics preprocessor
    # tas_diff_winter = ..
    # tas_diff_summer = ..
    return ...


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

    # First precompute the segment season means for the control period ...
    segments = segment(dataset, *cfg['control_period'], step=5)
    season_means = segments.groupby('time.season').mean()
    segments_season_means['control'] = season_means

    # ... then for all the future scenarios
    for name, info, in cfg['scenarios']:
        segments = segment(dataset, *info['resampling_period'], step=5)
        season_means = segments.groupby('time.season').mean()
        segments_season_means[name] = season_means

    # Step 1a: Get 1000 combinations for the control period
    # These samples should have the same mean winter precipitation as the mean of the x ensemble members
    winter_mean_pr_segments = segments_season_means['control'].pr.sel(season='DJF')
    target_winter_mean_pr = winter_mean_pr_segments.mean()
    top1000 = most_promising_combinations(winter_mean_pr_segments, target_winter_mean_pr)

    # Step 1b: Get 1000 combinations for the future period
    # The target value is a relative change wrt the overall mean of the control period
    control_period_winter_mean_pr = target_winter_mean_pr
    for scenario, info, in cfg['scenarios']:
        target_winter_mean_pr = control_period_winter_mean_pr * (1 + info['dpr_winter']/100)
        winter_mean_pr_segments = segments_season_means[scenario].pr.sel(season='DJF')
        top1000 = most_promising_combinations(winter_mean_pr_segments, target_winter_mean_pr)

    # Step 2a: From the 1000 samples in 1a, select those for which summer pr,
    # and summer and winter tas are within the percentile bounds specified
    # in the recipe for each scenario for the control period


    # Step 2b: From the 1000 samples in 1a, select those for which summer pr,
    # and summer and winter tas are within the percentile bounds specified
    # in the recipe for each scenario for the future period


    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)



# Speed test for generating xarray indexer arrays over and over, or updating them on the fly.
###############################################################################################3
def _get_helper_array(combination):
    # use a helper 'indices' dataarray to select one single ensemble member for each segment.
    indices = xr.DataArray(
        data = list(combination),
        dims = ['segment'],
        coords = {'segment': np.arange(6)})
    return indices

# %timeit [_get_helper_array(x) for x in product(range(8), repeat=6)]
# 35 s ± 4.45 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
###############################################################################################3
indices = xr.DataArray(
    data = np.zeros(6),
    dims = ['segment'],
    coords = {'segment': np.arange(6)})

def update(indices, values):
    indices.values = list(values)
    return indices

# %timeit [update(indices, x) for x in product(range(8), repeat=6)]
# 1.23 s ± 35 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)


###############################################################################################3
def all_possible_combinations(n_ensemble_members, n_segments):
    """Generate indexer for all possible combinations of segmented dataset."""

    # Create a DataArray once...
    indices = xr.DataArray(
        data = np.zeros(n_segments),
        dims = ['segment'],
        coords = {'segment': np.arange(n_segments)})

    # ...and update its values for all possible combinations
    for x in product(range(n_ensemble_members), repeat=n_segments):
        indices.values = list(x)
        yield indices

# %timeit [x for x in all_possible_combinations(8, 6)]
# 1.21 s ± 11.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
