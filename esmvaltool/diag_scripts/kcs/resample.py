"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
import os
from pathlib import Path
from datetime import datetime
import functools
import logging
import warnings
from pprint import pformat
import json
import yaml
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvalcore.preprocessor import area_statistics
import kcsutils
from math import factorial, inf
import itertools
from datetime import datetime
import logging
import numpy as np
import pandas as pd
import iris
import h5py


REF_PERIOD = (1980, 2009)
RELATIVE = ('pr',)
ALLSEASONS = ['djf', 'mam', 'jja', 'son']
STATS = ['mean', '5', '10', '25', '50', '75', '90', '95']
NSTEP1 = 1000
NSTEP3 = 8
NSAMPLE = 10_000
NSECTIONS = 6

logger = logging.getLogger(__name__)   # pylint: disable=invalid-name


def read_steering_table(filename, scenarios, scenario_selection=None):
    """Read and amend the steering table

    Reads the steering table from CSV, adds the precip scenarios, and
    selects any given scenarios.

    """

    table = pd.read_csv(filename)
    table = table.rename(columns={'target year': 'epoch', 'resampling period': 'period',
                                  'name': 'scenario'})
    table['epoch'] = table['epoch'].astype(str)
    # Convert the "(start, stop)" string to a proper tuple of integer years.
    table['period'] = table['period'].apply(lambda x: tuple(map(int, x[1:-1].split(','))))

    # Add the precip scenarios.
    table2 = pd.DataFrame({'subscenario': list(scenarios.keys()),
                           'precip_change_per_t': list(scenarios.values())})
    table['key'] = table2['key'] = 0
    table = pd.merge(table, table2, on='key', how='outer').drop(columns='key')
    table['precip_change'] = table['model_delta_t'] * table['precip_change_per_t']

    # If there is any scenarion selection, remove the not-selected ones
    if scenario_selection:
        selection = np.zeros(len(table), dtype=np.bool)
        for scen, epoch, subscen in scenario_selection:
            selection = selection | ((table['scenario'] == scen) &
                                     (table['subscenario'] == subscen) &
                                     (table['epoch'] == epoch))
        table = table[selection]

    logger.debug("Scenario's and target precipiation for S1: %s",
                 table[['epoch', 'scenario', 'subscenario', 'period', 'model_delta_t',
                        'precip_change_per_t', 'precip_change']].sort_values(
                            ['subscenario', 'scenario', 'epoch'], ascending=[True, False, False]))

    return table


def save_indices_h5(filename, indices):
    """Save the (resampled) array of indices in a HDF5 file"""
    h5file = h5py.File(filename, 'a')
    for key, value in indices.items():
        name = "/".join(key)
        try:
            group = h5file[name]
        except KeyError:
            group = h5file.create_group(name)
        for k, v in value['meta'].items():  # pylint: disable=invalid-name
            # Transforms a dict to str so it can be saved in HDF 5
            group.attrs[k] = json.dumps(v) if isinstance(v, (dict, list)) else v
        if 'control' in group:
            del group['control']
        group.create_dataset('control', data=value['data']['control'])
        if 'future' in group:
            del group['future']
        group.create_dataset('future', data=value['data']['future'])


def save_resamples(filename, diffs, as_csv=True):
    """Save the resampled data, that is, the changes (differences) between
    epochs, to a HDF5 file

    With `as_csv` set to `True`, save each individual
    scenario-variable-season combination as a separate CSV file, named
    after the combination (and `resampled_` prepended). These CSV
    files can be used with kcs.change_perc.plot, with the --scenario
    option.

    Files are always overwritten if they exist.

    """

    h5file = h5py.File(filename, 'w')
    for key, value in diffs.items():
        keyname = "/".join(key)
        for var, value2 in value.items():
            for season, diff in value2.items():
                name = f"{keyname}/{var}/{season}"
                if name not in h5file:
                    h5file.create_group(name)
                group = h5file[name]

                # Remove existing datasets, to avoid problems
                # (we probably could overwrite them; this'll work just as easily)
                for k in {'diff', 'mean', 'std', 'keys'}:
                    if k in group:
                        del group[k]

                group.create_dataset('diff', data=diff.values)
                group.create_dataset('mean', data=diff.mean(axis=0))
                group.create_dataset('std', data=diff.std(axis=0))
                # pylint: disable=no-member
                dataset = group.create_dataset('keys', (len(STATS),), dtype=h5py.string_dtype())
                dataset[:] = STATS

                assert len(STATS) == len(diff.mean(axis=0))

                if as_csv:
                    csvfile = "_".join(key)
                    csvfile = f"resampled_{csvfile}_{var}_{season}.csv"
                    diff.to_csv(csvfile, index=False)



def segment_data(cubes, period, control_period, nsections, seasons=None):
    """Given a list of cubes (or CubeList), return a dict with periods and seasons extracted"""

    if seasons is None:
        seasons = ALLSEASONS

    data = {season: {} for season in seasons}
    # pragma pylint: disable=cell-var-from-loop
    for key, years in zip(['control', 'future'], [control_period, period]):
        # Chop each cube into n-year segments
        span = (years[1] - years[0] + 1) // nsections
        constraint = iris.Constraint(year=lambda point: years[0] <= point <= years[1])
        logger.debug("Extraction %s period %s for all datasets", key, years)
        excubes = [constraint.extract(cube) for cube in cubes]
        for season in seasons:
            logger.debug("Extraction season %s for all datasets", season)
            constraint = iris.Constraint(season=lambda point: point == season)
            season_cubes = [constraint.extract(cube) for cube in excubes]
            logger.debug("Extracting %d-year segments for all datasets", span)
            data[season][key] = np.array([
                [iris.Constraint(year=lambda point: year <= point < year+span).extract(cube)
                 for year in range(years[0], years[1], span)]
                for cube in season_cubes], dtype=np.object)
    # pragma pylint: enable=cell-var-from-loop

    return data


def create_indices(nsets, nsections):
    """Create an array of indices for all possible resamples

    This results in an array with a shape of (nsets^nsections,
    nsections). Each row is a separate resample, indexing the relevant
    segment from a dataset i.

    """

    indices = list([range(nsets) for _ in range(nsections)])
    indices = np.array(list(itertools.product(*indices)))
    assert len(indices) == nsets**nsections, "incorrect number of generated indices"

    return indices


class Calculation:
    """Class to calculate the difference between a calculated precipiation
    change, and a target value

    We use a class, so it can be used with multiprocessing: fixed
    arguments are passed to the constructor, while the variable
    argument (the subselection of the data) is passed to the call
    method.

    An alternative is the use of functools.partial.

    """

    def __init__(self, data, target):
        self.data = data
        self.target = target
        self.cols = np.arange(data.shape[1])

    def __call__(self, index):
        return np.abs(self.data[index, self.cols].mean() - self.target)


def calc_means(data, seasons=None):
    """Calculate the averages of all n-year segments, and store them in a
    dataframe

    The dataframe columns correspond to the n-year segments (simply
    numbered 0 to ndata-1).

    The dataframe index consists of the variable, season and period
    (control or future), and the individual datasets (runs) as a
    multi-index.

    The values then, of course, are the averages for that specific
    n-year segment, variable, season, period and dataset.

    """

    if seasons is None:
        seasons = ALLSEASONS
    means = {}
    std = {}
    for var in data:
        for season in seasons:
            for period in ['control', 'future']:
                item = data[var][season][period]
                for i, cubelist in enumerate(item):
                    means[(var, season, period, i)] = [cube.data.mean() for cube in cubelist]
                    std[(var, season, period, i)] = [cube.data.std() for cube in cubelist]

    df_mean = pd.DataFrame.from_dict(means, orient='index')
    df_mean.index = pd.MultiIndex.from_tuples(means.keys(),
                                              names=['var', 'season', 'period', 'run'])
    df_mean = df_mean.sort_index()

    df_std = pd.DataFrame.from_dict(std, orient='index')
    df_std.index = pd.MultiIndex.from_tuples(means.keys(),
                                             names=['var', 'season', 'period', 'run'])
    df_std = df_std.sort_index()

    return df_mean, df_std


def prepare_data(dataset, variables, period, control_period, nsections):
    """Prepare the data for a scenario

    - segment the data into nsections.

    - calculate the means for variables and seasons of interest, for
      all individual ensemble runs, for each n-year period.

    """

    data = {}
    ndata = set()
    for var in variables:
        cubes = dataset.loc[dataset['var'] == var, 'cube']
        logger.debug("Segmenting %s data into %d sections", var, nsections)
        data[var] = segment_data(cubes, period, control_period, nsections)
        season = list(data[var].keys())[0]
        segment = data[var][season]['future']
        ndata.add(len(segment))

    assert len(ndata) == 1, "Datasets are not the same length for different variables"
    ndata = ndata.pop()
    indices = create_indices(ndata, nsections)

    logger.debug("start")
    means, _ = calc_means(data)
    logger.debug("stop")

    return data, indices, means, ndata


def calculate_s1(means, indices, target, var='pr', season='djf'):
    """Calculate the subset S1: winter precipitation change equals `target`"""

    # Calculate the procentual change in winter precipitation for all sections individually,
    # With respect to an overall-averaged control period
    controlmean = means.loc[(var, season, 'control'), :]
    controlmean = controlmean.mean().mean()
    data = means.loc[(var, season, 'future'), :].values
    data = 100 * (data - controlmean) / controlmean

    calculation = Calculation(data, target)
    tstart = datetime.now()
    values = list(map(calculation, indices))
    values = np.array(values)

    tstart = datetime.now()
    order = np.argsort(values)

    selind = {'future': indices[order]}

    data = means.loc[(var, season, 'control'), :].values
    data = data - controlmean
    calculation = Calculation(data, 0)
    tstart = datetime.now()
    values = list(map(calculation, indices))
    values = np.array(values)

    tstart = datetime.now()
    order = np.argsort(values)

    selind['control'] = indices[order]

    return selind


def calculate_s2(means, indices, scenarios):
    """Calculate the subset S2: select percentile ranges for average
    precipitation and temperatures"""

    s2_indices = {}
    columns = np.arange(indices['control'].shape[1])
    for period in ['control', 'future']:
        selection = np.ones(indices[period].shape[0], dtype=np.bool)
        s2_indices[period] = indices[period][selection]
        for scenario in scenarios:
            var = scenario['var']
            season = scenario['season']
            logger.debug("Subsetting with percentiles for %s, %s, %s", var, season, period)
            values = means.loc[(var, season, period), :].values[s2_indices[period], columns]
            # Calculate mean along the columns, i.e., one 30-year period
            mean = values.mean(axis=1)
            logger.debug("Min, max, mean, median values: %f  %f  %f  %f",
                         mean.min(), mean.max(), mean.mean(), np.median(mean))
            p_range = scenario[period]
            logger.debug("Percentile range = %s", p_range)
            low, high = np.percentile(mean, p_range)
            logger.debug("Percentile values = %f  --  %f", low, high)
            selection = (low <= mean) & (mean <= high)
            s2_indices[period] = s2_indices[period][selection]
            logger.debug("Subsetting down to %d samples", selection.sum())

    return s2_indices


def calculate_s3(indices_dict, penalties, nstep3=NSTEP3, nsample=NSAMPLE,
                 minimum_penalty=None):

    """Calculate the subset S3: find a subset with the least re-use of
    segments, using random sampling"""

    if minimum_penalty is None:
        minimum_penalty = np.finfo(float).eps

    s3_indices = {}
    rng = np.random.default_rng()   # pylint: disable=no-member
    for period, indices in indices_dict.items():
        n = len(indices)  # pylint: disable=invalid-name
        ncomb = factorial(n) // factorial(nstep3) // factorial(n - nstep3)
        logger.debug("Number of combinations for the %s period = %d", period, ncomb)
        best = np.inf
        for i in range(nsample):
            choice = rng.choice(indices, size=nstep3, replace=False)
            penalty = 0
            for column in choice.T:
                _, counts = np.unique(column, return_counts=True)
                for count in counts:
                    penalty += penalties[count]
            if penalty <= best:
                best = penalty
                s3_indices[period] = choice
            if best < minimum_penalty:
                logger.debug("Minimum penalty reached after %d iterations", i)
                break
        logger.info("Best sample for %s (penalty: %f) out of %d samples: %s",
                    period, best, nsample, s3_indices[period])

    return s3_indices


def find_resamples(indices, means, precip_change, ranges, penalties,
                   nstep1=NSTEP1, nstep3=NSTEP3, nsample=NSAMPLE):
    """Find the (best) resamples

    This does the actual work:

    - run step one (for both future and control), limiting samples by
      required precipitation change

    - run step two, limiting samples by requiring their mean values
      to be in certain percentile change compared to the overall mean
      values of the remaining resamples.

    - run step three, by randomly selecting samples, keeping the
      number of duplicate segments to a minimum.

    Note that the actual data cubes are not required: we have already
    calculated the averages, and those are used in step one and two
    (step three doesn't require any actual data).

    """
    logger.debug("Calculating S1")
    logger.debug("Precipitation change: %.1f", precip_change)

    indices = calculate_s1(means, indices, precip_change)
    indices['control'] = indices['control'][:nstep1]
    indices['future'] = indices['future'][:nstep1]

    indices = calculate_s2(means, indices, ranges)
    logger.debug("The S2 subset has %d & %d indices for the control & future periods, resp.",
                 len(indices['control']), len(indices['future']))

    indices = calculate_s3(indices, penalties, nstep3=nstep3, nsample=nsample)

    return indices


def resample(indices, data, variables, seasons, relative):
    """Perform the actual resampling of data, given the resampled indices"""

    percs = list(map(float, STATS[1:]))

    diffs = {}
    for key, value in indices.items():
        diffs[key] = {}
        tempindices = value['data']
        for var in variables:
            diffs[key][var] = {}
            for season in seasons:
                stats = {}
                for period in ['control', 'future']:
                    cubes = data[key][var][season][period]
                    columns = np.arange(cubes.shape[1])
                    resampled = cubes[tempindices[period], columns]
                    stats[period] = pd.DataFrame({'mean': 0, '5': 0, '10': 0, '25': 0, '50': 0,
                                                  '75': 0, '90': 0, '95': 0},
                                                 index=np.arange(len(resampled)))
                    for i, cubeset in enumerate(resampled):
                        block = np.hstack([cube.data for cube in cubeset])
                        stats[period].loc[i, 'mean'] = block.mean()
                        stats[period].loc[i, STATS[1:]] = np.percentile(block, percs)

                diff = pd.DataFrame()
                for col in STATS:
                    diff[col] = stats['future'][col] - stats['control'][col]
                    if var in relative:
                        diff[col] = 100 * diff[col] / stats['control'][col]
                diffs[key][var][season] = diff
    return diffs


def calc(dataset, steering_table, ranges, penalties,
         nstep1=NSTEP1, nstep3=NSTEP3, nsample=NSAMPLE,
         nsections=NSECTIONS, reference_period=REF_PERIOD,
         relative=RELATIVE):
    """DUMMY DOCSTRING"""


    data = {}
    indices = {}
    means = {}
    variables = dataset['var'].unique()

    for _, row in steering_table.iterrows():
        period = tuple(map(int, row['period']))
        scenario = row['scenario']
        subscenario = row['subscenario']
        epoch = row['epoch']
        mainkey = (str(epoch), scenario, subscenario)
        logger.info("Preparing data for %s_%s - %s %s",
                    scenario, subscenario, epoch, period)
        data[mainkey], indices[mainkey], means[mainkey], _ = prepare_data(
            dataset, variables, period, reference_period, nsections)

    all_indices = {}
    for _, row in steering_table.iterrows():
        period = tuple(map(int, row['period']))
        scenario = row['scenario']
        subscenario = row['subscenario']
        epoch = row['epoch']
        mainkey = (str(epoch), scenario, subscenario)
        precip_change = row['precip_change']

        rrange = ranges[scenario][subscenario][epoch]
        logger.info("Processing %s_%s - %s %s, %.2f pr", scenario, subscenario, epoch,
                    period, precip_change)
        final_indices = find_resamples(indices[mainkey], means[mainkey],
                                       precip_change, rrange, penalties,
                                       nstep1, nstep3, nsample)

        attrs = {
            'scenario': scenario, 'subscenario': subscenario, 'epoch': epoch,
            'period': period, 'reference-period': reference_period,
            'ranges': rrange, 'winter-precip-change': precip_change}
        all_indices[mainkey] = {'data': final_indices, 'meta': attrs}

    diffs = resample(all_indices, data, variables, seasons=['djf', 'mam', 'jja', 'son'],
                     relative=relative)

    return all_indices, diffs


def main(cfg):
    """

    """
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable

    my_files_dict = group_metadata(cfg['input_data'].values(), 'filename')
    #from pprint import pprint
    #pprint(my_files_dict)
    paths = [Path(filename) for filename in my_files_dict.keys()]
    dataset = kcsutils.read_data(paths)
    sel = dataset['model'] == cfg['model']
    dataset = dataset.loc[sel, :]
    temp = []
    for var, group in dataset.groupby('var'):
        group = kcsutils.match(
            group, match_by='ensemble', on_no_match='randomrun',
            historical_key='historical')
        group['shape'] = group['cube'].apply(lambda x: x.shape)
        #with pd.option_context('display.max_rows', None):
        #    print(dataset[['model', 'experiment', 'var']])
        group['var'] = var
        temp.append(kcsutils.concat_cubes(group, historical_key='historical'))
    dataset = pd.concat(temp)
    columns = ['model', 'experiment', 'realization', 'initialization', 'physics', 'prip', 'var', 'shape', 'var']

    filename = (Path(cfg['run_dir']).parent.parent / 'diagnostic1' / 'script2' /
                cfg['steering_table'])
    precip_scenarios = {'L': 4, 'H': 8}
    steering_table = read_steering_table(filename, precip_scenarios)

    path = Path(__file__).parent /  cfg['conditions']
    with open(path) as fh:  # pylint: disable=invalid-name
        conditions = yaml.load(fh, Loader=yaml.Loader)

    penalties = cfg['penalties']
    penalties = {int(key): value for key, value in penalties.items()}
    # Fill up the penalty dict
    for i in range(max(penalties.keys())+1, NSTEP3+1):
        penalties[i] = inf
    indices, diffs = calc(dataset, steering_table, conditions, penalties)
    save_indices_h5(cfg['indices_output'], indices)
    save_resamples(cfg['resamples_output'], diffs)


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
