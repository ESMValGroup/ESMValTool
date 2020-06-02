
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



def segment_datasets(datasetlist, start=1981, end=2010, step=5):
    """Return a nested list with [[x segments of n-year periods] for each input dataset]."""
    segmented_datasets = []
    for dataset in datasetlist:
        filename = dataset['input_file']
        cube = iris.load_cube('filename')
        segments = [
            cube.extract(iris.Constraint(time=lambda cell: year <= cell.point.year < year+step))
            for year in range(start, end, step)
        ]
        segmented_datasets.append(segments)

    return segmented_datasets


def recombine_cubes(segments, combination):
    """

    segments: a nested list with [6 segments of 5-year periods] for each ensemble member.
    combination: a specific order to reconstruct a new climate, e.g.
        [1, 3, 4, 1, 6, 8] means:
        - take the first 5 years of the 1st ensemble member,
        - the second 5 years of the 3rd ensemble member,
        - etc.
    """
    cubelist = iris.cube.CubeList(controls[i][j] for i, j in enumerate(combination))
    return cubelist.concatenate()


def get_change(control, future, season):
    """Return the difference in mean between future and control period.

    e.g. control = precipitation cubes for 1981-2010
         future = precipitation cubes for 2030-2059
         season = 'DJF'
         returns change in mean winter precipitation between future and control
    """

    dpr = control.extract_season(season).mean() - future.extract_season(season).mean()
    return dpr

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


def get_filter2_variables(combination):
    resampled_control_tas = recombine_cubes(tas_controls, combination)
    resampled_future_tas = recombine_cubes(tas_futures, combination)

    resampled_control_tas = recombine_cubes(tas_controls, combination)
    resampled_future_tas = recombine_cubes(tas_futures, combination)

    pr_diff_summer = get_change(resampled_future_pr, resampled_control_pr, season='JJA')
    tas_diff_summer = get_change(resampled_future_tas, resampled_control_tas, season='JJA')
    tas_diff_winter = get_change(resampled_future_tas, resampled_control_tas, season='DJF')

    return pr_diff_summer, tas_diff_summer, tas_diff_winter


def filter3(df, n=8):
    """Sort by the number of unique elements in 'combination' and return the top n."""
    df['n_datasets'] = df.combination.map(lambda x: len(set(x)))
    return df.sort('n_datasets', descending=True).head(n)


def main():

    # a list of dictionaries describing all datasets passed on to the recipe
    dataset_dicts = cfg['input_data'].values()

    # Only keep datasets from the target model:
    pr = select_metadata(dataset_dicts, dataset=cfg['target_model'], short_name='pr')
    tas = select_metadata(dataset_dicts, dataset=cfg['target_model'], short_name='tas')

    pr_controls = segment_datasets(pr, *cfg['control_period'])
    tas_controls = segment_datasets(tas, *cfg['control_period'])

    for scenario_name, info,  in cfg['scenarios']:
        print('Working on scenario:', scenario_name)

        pr_futures = segment_datasets(target_pr_datasets, *info['resampling_period'])
        tas_futures = segment_datasets(target_tas_datasets, *info['resampling_period'])

    # Make a lot of different re-combinations of the datasets
    n_ensemble_members = len(datasetlist)
    pr_diffs_winter = pd.DataFrame(columns=['combination', 'pr_diff_winter'])
    for combination in itertools.product(n_ensemble_members, repeat=6):
        resampled_control = recombine_cubes(pr_controls, combination)
        resampled_future = recombine_cubes(pr_futures, combination)

        pr_diff_winter = get_change(resampled_future, resampled_control, season='DJF')
        pr_diffs_winter.append([combination, pr_diff_winter])


    # Filter 1
    top1000 = filter1(pr_diffs_winter, target_dpr)

    # Filter 2
    # Create additional columns in the dataframe by applying the following steps:
    for combi in top1000.combinations:
        pr_diff_summer, tas_diff_summer, tas_diff_winter = get_filter2_variables(combination)
        # somehow append to dataframe: ['pr_diff_summer, tas_diff_summer, tas_diff_winter]

    # I think this function can also be applied using (avoiding the loop):
    # top1000['pr_diff_summer, tas_diff_summer, tas_diff_winter] = top1000.apply(get_filter2_variables)
    top1000.filter(tas_diff_summer - target_percentile_tas_summer < ...)
    top1000.filter(tas_diff_winter - target_percentile_tas_winter < ...)
    top1000.filter(pr_diff_summer - target_percentile_pr_summer < ...)



    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
