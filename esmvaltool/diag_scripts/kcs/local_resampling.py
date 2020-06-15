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

from esmvaltool.diag_scripts.shared import (run_diagnostic, select_metadata,
                                            get_diagnostic_filename,
                                            ProvenanceLogger,
                                            get_plot_filename)

from tqdm import tqdm

LOGGER = logging.getLogger(Path(__file__).name)


def create_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption': "Resampling of local climate model.",
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
    LOGGER.info("Reading input data for target model")
    dataset_dicts = cfg['input_data'].values()
    target_model_metadata = select_metadata(dataset_dicts,
                                            dataset=cfg['target_model'])
    dataset = []
    ancestor_files = []
    for short_name in "pr", "tas":
        group = f'{short_name}_target'
        var = select_metadata(target_model_metadata, variable_group=group)
        files = [metadata['filename'] for metadata in var]
        dataset.append(
            xr.open_mfdataset(files,
                              concat_dim='ensemble_member',
                              combine='nested'))
        ancestor_files.extend(files)
    provenance = create_provenance_record(ancestor_files)
    return xr.merge(dataset), provenance


def compute_segment_season_means(dataset, period, step, filename):
    """Compute season means for each n-year segment of input dataset.

    Saves the seasonal means to filename
    """
    LOGGER.info("Computing seasonal means for segmented dataset")
    segments = []
    for year in range(*period, step):
        segments.append(
            dataset.sel(time=slice(str(year), str(year + step - 1))))
    segmented_dataset = xr.concat(segments, dim='segment')
    season_means = segmented_dataset.groupby('time.season').mean()
    season_means.to_netcdf(filename)
    LOGGER.info("Intermediate results stored as %s.", filename)
    # TODO store provenance?


def all_possible_combinations(n_ensemble_members, n_segments):
    """Generate indexer.

    Returns indices for all possible combinations of
    segments and ensembles.
    """
    # Create a DataArray once...
    indices = xr.DataArray(data=np.zeros(n_segments),
                           dims=['segment'],
                           coords={'segment': np.arange(n_segments)})

    # ...and update its values for all possible combinations
    for combination in product(range(n_ensemble_members), repeat=n_segments):
        indices.values = list(combination)
        name = ''.join([str(x) for x in combination])
        yield name, indices


def find_and_save_top1000(segment_means, target, filename):
    """Select n_combinations that are closest to the target.

    First, get all possible combinations, then select n_combinations.
    Store the result under 'filename'
    """
    n_segments = len(segment_means.segment)
    n_members = len(segment_means.ensemble_member)
    total = n_members**n_segments
    all_combinations = all_possible_combinations(n_members, n_segments)

    results = []
    for combination, idx in tqdm(all_combinations, total=total):
        results.append([
            combination,
            abs(
                segment_means.sel(ensemble_member=idx).mean('segment') -
                target).values.tolist()
        ])

    # Create a pandas dataframe with the combinations and distance to target
    top1000 = pd.DataFrame(results,
                           columns=['combination', 'distance_to_target'
                                    ]).sort_values('distance_to_target').head(
                                        1000).set_index('combination')
    top1000.to_csv(filename)
    LOGGER.info("Intermediate results stored as %s.", filename)
    # TODO: store provenance?


def selected_combinations(combinations):
    """Generate indexer for all selected combinations of segmented dataset.

    combinations: a list of combinations (dtype: str !)
    """
    n_segments = len(combinations[0])

    # Create a DataArray once...
    indices = xr.DataArray(data=np.zeros(n_segments),
                           dims=['segment'],
                           coords={'segment': np.arange(n_segments)})

    # ...and update its values for each selected combination
    for combination in combinations:
        indices.values = [int(x) for x in combination]
        yield indices


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
            new_overall_season_means.pr.sel(season='JJA').values,
            new_overall_season_means.tas.sel(season='DJF').values,
            new_overall_season_means.tas.sel(season='JJA').values
        ])
    return pd.DataFrame(filter2_variables, columns=columns)


def within_bounds(values, bounds):
    """Return true if value is within percentile bounds."""
    low, high = np.percentile(values, bounds)
    return values.between(low, high)


def percentile_based_subset(top1000, info, period):
    """Select samples based on the percentile bounds.

    Select samples for which summer pr, and summer and winter tas
    are within the percentile bounds specified in the recipe.
    """
    pr_summer = within_bounds(top1000['pr_summer'],
                              info[f'pr_summer_{period}'])
    tas_winter = within_bounds(top1000['tas_winter'],
                               info[f'tas_winter_{period}'])
    tas_summer = within_bounds(top1000['tas_summer'],
                               info[f'tas_summer_{period}'])
    subset = top1000[np.all([pr_summer, tas_winter, tas_summer], axis=0)]
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
        [list(combination) for combination in combinations])

    # Store the indices in a nice dataframe
    _, n_segments = combinations.shape
    best_combination = pd.DataFrame(
        data=None,
        columns=[f'Segment {x}' for x in range(n_segments)],
        index=[f'Combination {x}' for x in range(n_sample)])

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

    scenario_output_tables = pd.concat(scenario_output_tables,
                                       axis=1,
                                       keys=['control', 'future'])
    return scenario_output_tables


def save_scenario_tables(scenario, tables, cfg, provenance):
    """Save the scenario tables as a csv file."""
    output_file = get_diagnostic_filename(f'final_indices_{scenario}',
                                          cfg,
                                          extension='csv')
    tables.to_csv(output_file)
    LOGGER.info('Final output for scenario %s stored as %s', scenario,
                output_file)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(output_file, provenance)


def make_plot(cfg):
    """Reproduce figure 5 from the paper."""
    # note that quantile is applied twice!
    # Once to get the pdf's of seasonal tas/pr
    # and once to get the multimodel pdf of the quantile changes
    metadata = cfg['input_data'].values()
    tas_cmip = select_metadata(metadata, variable_group='tas_cmip')
    envelope = []
    for data_dict in tas_cmip:
        dataset = xr.open_dataset(data_dict['filename']).tas
        grouped_control = dataset.sel(
            time=slice('1981', '2010')).groupby('time.season')
        grouped_future = dataset.sel(time='2050').groupby('time.season')

        quantiles = [.05, .1, .25, .5, .75, .90, .95]
        qcontrol = grouped_control.quantile(quantiles)
        qfuture = grouped_future.quantile(quantiles)
        qchange = qfuture - qcontrol
        envelope.append(qchange)
    cmip5 = xr.concat(envelope, dim='multimodel')
    # prevent confusion between dimension 'quantile' and method 'quantile'
    cmip5 = cmip5.rename({'quantile': 'percentile'})

    fig, ax = plt.subplots()
    x_locs = np.arange(len(quantiles))
    for high, low in [[0.9, 0.1], [0.75, 0.25]]:
        ax.fill_between(x_locs,
                        cmip5.sel(season='DJF').quantile(high,
                                                         dim='multimodel'),
                        cmip5.sel(season='DJF').quantile(low,
                                                         dim='multimodel'),
                        color='k',
                        alpha=0.3)

    ax.set_xticks(x_locs)
    ax.set_xticklabels([f'P{100*x:02.0f}' for x in quantiles])
    filename = get_plot_filename('envelope_figure', cfg)
    fig.savefig(filename, bbox_inches='tight', dpi=300)
    LOGGER.info("Envelope figure stored as %s", filename)
    # TODO: Add multiple subplots for summer/winter and pr/tas
    # TODO: Add final selection of scenario to figure
    # TODO: Add provenance for figure?


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

    # Step 0: Read the data and compute segment season means
    dataset, provenance = get_input_data(cfg)
    segments_season_means = {}

    # First for control ...
    LOGGER.info("Precompute segment season means for the control period ...")
    filename = get_diagnostic_filename('season_means_control', cfg, 'nc')
    if Path(filename).exists():
        LOGGER.info("Found intermediate file %s", filename)
    else:
        compute_segment_season_means(dataset,
                                     cfg['control_period'],
                                     step=5,
                                     filename=filename)
    segments_season_means['control'] = xr.open_dataset(filename)

    # ... then for each scenario
    for scenario, info in cfg['scenarios'].items():
        LOGGER.info("Precompute segment season means for the scenario %s",
                    scenario)
        filename = get_diagnostic_filename(f'season_means_{scenario}', cfg,
                                           'nc')
        if Path(filename).exists():
            LOGGER.info(f"Found intermediate file %s", filename)
        else:
            compute_segment_season_means(dataset,
                                         info['resampling_period'],
                                         step=5,
                                         filename=filename)
        segments_season_means[scenario] = xr.open_dataset(filename)

    # Step 1: Get 1000 recombinations and compute their winter mean pr
    selected_indices = {}

    # First for the control ...
    LOGGER.info('Get 1000 recombinations for the control'
                ' period based on mean winter precipitation')
    djf_pr_control = segments_season_means['control'].pr.sel(season='DJF')
    djf_pr_control_mean = djf_pr_control.mean().values  # the target value
    filename = get_diagnostic_filename("top1000_control", cfg, 'csv')
    if Path(filename).exists():
        LOGGER.info("Found intermediate file %s", filename)
    else:
        find_and_save_top1000(djf_pr_control, djf_pr_control_mean, filename)
    selected_indices['control'] = pd.read_csv(filename,
                                              dtype={'combination': 'str'})

    # ... then for each scenario
    for scenario, info in cfg['scenarios'].items():
        LOGGER.info(
            'Get 1000 recombinations for the'
            ' resampling period of scenario %s', scenario)
        filename = get_diagnostic_filename(f"top1000_{scenario}", cfg, 'csv')
        if Path(filename).exists():
            LOGGER.info("Found intermediate file %s", filename)
        else:
            # Target is a relative change with respect to the control period
            target = djf_pr_control_mean * (1 + info['dpr_winter'] / 100)
            djf_pr_segments = segments_season_means[scenario].pr.sel(
                season='DJF')
            find_and_save_top1000(djf_pr_segments, target, filename)
        selected_indices[scenario] = pd.read_csv(filename,
                                                 dtype={'combination': 'str'})

    # Step 2a: For each set of 1000 samples from 1a-b
    for name, dataframe in selected_indices.items():
        LOGGER.info("Compute summer mean pr and summer and winter "
                    "mean tas for 1000 selected combinations for %s" % name)
        segment_means = segments_season_means[name]
        selected_indices[name] = seasonal_mean_characteristics(
            dataframe.combination, segment_means)

    # Step 2b: For each scenario, select samples for which summer pr, and
    # summer and winter tas are within the percentile bounds specified
    # in the recipe
    for scenario, info in cfg['scenarios'].items():
        LOGGER.info("Getting percentile-based subsets for scenario %s",
                    scenario)
        selected_indices[scenario] = {
            # Get relatively high/low values for the control period ...
            'control':
            percentile_based_subset(selected_indices['control'], info,
                                    'control'),
            # Get relatively high/low values for the future period
            'future':
            percentile_based_subset(selected_indices[scenario], info, 'future')
        }

        del selected_indices['control']  # No longer needed

    # Step 3: Select final set of eight samples with minimal reuse of the same
    # ensemble member for the same period.
    for scenario, dataframes in selected_indices.items():
        LOGGER.info("Selecting 8 final samples for scenario %s", scenario)
        scenario_output_tables = make_scenario_tables(dataframes)
        LOGGER.info("Selected recombinations for scenario %s: \n %s", scenario,
                    scenario_output_tables)
        save_scenario_tables(scenario, scenario_output_tables, cfg, provenance)

    # Step 4: Plot the results
    if cfg['write_plots']:
        make_plot(cfg)
        plt.show()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
