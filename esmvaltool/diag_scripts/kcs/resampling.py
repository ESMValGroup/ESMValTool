
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





def main():


    # a list of dictionaries describing all datasets passed on to the recipe
    dataset_dicts = cfg['input_data'].values()

    # Only keep datasets from the target model:
    datasetlist = select_metadata(dataset_dicts, dataset=target_model, short_name='pr')

    controls = []
    futures = []
    for dataset in datasetlist():
        filename = dataset['input_file']
        cube = iris.read_cube('filename')
        control = cube.select(....)
        future =

        control_segments = [control.select(...), control.select(...)]
        future_segments = [control.select(...), control.select(...)]
        controls.append(cube)
        futures.append(cube)

    # Make a lot of different combinations of the datasets
    n_ensemble_members = len(datasetlist)
    pr_diffs_winter = pd.DataFrame(columns=['combination', 'pr_diff_winter'] )
    for combination in itertools.product(n_ensamble_members, repeat=6):
        resampled_control = iris.cube.CubeList(controls[i][j] for i, j in enumerate(combination)).concatenate()
        resampled_future = iris.cube.CubeList(futures[i] for i, j in enumerate(combination)).concatenate()

        pr_diff_winter = get_precipitation_change(resampled_future, resample_control, variable='pr', season='DJF')
        pr_diffs_winter.append([combination, pr_diff_winter])

        # These are needed only for the 1000 best samples from the step above, so may do it later, outside this superloop
        pr_diff_summer = get_precipitation_change(resampled_future, resample_control, variable='tas', season='JJA') # Maybe use esmvaltool's seasonal statistics preprocessor
        tas_diff_winter = ..
        tas_diff_summer = ..


    # Filter 1
    filter1 = pr_diffs_winter.sort('pr_diff_winter').head(1000)

    # Filter 2
    filter2 = filter1.#dostuff




    return

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
