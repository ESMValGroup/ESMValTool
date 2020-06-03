
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
    segmented_dataset = xr.concat([
        ds.sel(time = slice(str(year), str(year+block_size)))
        for year in range(start_year, end_year, block_size)
    ], dim='segment')
    return segmented_dataset


def _recombine(segments, combination):
    # use a helper 'indices' dataarray to select one single ensemble member for each segment.
    indices = xr.DataArray(
        data = combination,
        dims = ['segment'],
        coords = {'segment': segments.segment})
    recombined = segments.sel(ensemble_member=indices)
    return recombined


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



def main():

    # a list of dictionaries describing all datasets passed on to the recipe
    dataset_dicts = cfg['input_data'].values()

    # Get all datasets for the target model; load all at once; concatenate along
    # a new dimension 'ensemble_member'; and only keep variables pr and tas
    target_model_metadata = select_metadata(dataset_dicts, dataset=cfg['target_model'])
    files = [metadata['filename'] for metadata in target_model_metadata]
    datasets = xr.open_mfdataset(files, concat_dim='ensemble_member', combine='nested')[['pr', 'tas']]

    segments_control = segment(dataset, *cfg['control_period'], step=5)

    winter_pr_mean_segments = ds_segmented.pr.groupby('time.season').mean().sel(season='DJF')

    # Find combinations of the control segments where the average is very close to the overall mean of all ensemble members in the control period
    controlmean = segments_control.pr.mean()
    recombined_climates = pd.DataFrame(columns=['combination', 'winter_mean_pr'])

    for combination in itertools.product(n_ensemble_members, repeat=6):

        # For now, only look at winter mean precipitation:
        recombined_winter_mean_pr_segments = _recombine(winter_pr_mean_segments, combination)
        recombined_winter_mean_pr_overal = resampled.mean('segment')

        # store them in a dataframe
        # keep the 1000 that are closest to controlmean.
        recombined_climates.append([combination, recombined_winter_mean)

    # Only keep the 1000 that are closest to controlmean
    recombined_climates.apply(lambda row: np.abs(row.winter_mean_pr - controlmean)).sort(ascending=True).head(1000)







    for scenario_name, info, in cfg['scenarios']:
        print('Working on scenario:', scenario_name)
        segments_future = segment(dataset, *info['resampling'], step=5)


        segmented_seasonal_means = ds_segmented.groupby('time.season').mean()
        segmented_seasonal_means.assign_coords(period=period, variable=variable)
        lookup_table.append(segmented_seasonal_means)
    lookup_table = xr.concatenate(lookup_table, dim=period)




    return


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

