"""A script to obtain the global mean temperature rise (change).

The global mean temperature rise is one of the two steering variables,
used in the interpolation of the regional climate model.

The minimum data for this script is one dataset with two different experiments:
historical and rcp (Representative Concentration Pathways).
For example:
- {dataset: ACCESS1-0, exp: historical, ensemble: r1i1p1, ...}
- {dataset: ACCESS1-0, exp: rcp45, ensemble: r1i1p1, ...}
"""
from pathlib import Path
from datetime import datetime
import functools
import logging
import warnings
from pprint import pformat
import iris
import pandas as pd
import numpy as np
from esmvaltool.diag_scripts.shared import (group_metadata,
                                            run_diagnostic,
                                            ProvenanceLogger)
import kcsutils
from kcsutils.attributes import get as get_attributes


PERIOD = (1961, 2099)
REFERENCE_PERIOD = (1980, 2009)
MINDATA = {'historical': 20, 'future': 4}

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption': "Global mean temperature rise (change).",
        'domains': ['global'],
        'authors': [
            'kalverla_peter',
            # 'rol_evert',
            'alidoost_sarah',
        ],
        'projects': [
            'esmval',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [],
    }
    return record


def add_extra_coord(cube):
    """Adda some coordinates to cube."""
    if not cube.coords('season'):
        iris.coord_categorisation.add_season(cube, 'time')
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time')
    if not cube.coords('season_year'):
        iris.coord_categorisation.add_season_year(cube, 'time',
                                                name='season_year')


def extract_season(cubes, season):
    """DUMMY DOC-STRING"""
    # Use a class instead of a lambda function, so we can pass the
    # constraint to multiprocessing (which doesn't handle lambdas).
    constraint = iris.Constraint(season=kcsutils.EqualConstraint(season))
    logger.info("Extracting season %s", season)
    for cube in cubes:
        if not cube.coords('season'):
            iris.coord_categorisation.add_season(cube, 'time')
    cubes = list(map(constraint.extract, cubes))
    return cubes


def perform_optional_operations(dataset, yearly=False, season=None):
    """Perform extracting of season, and averaging of years.

    These are optional steps before calculation of reference values.

    """
    if season:
        dataset['cube'] = extract_season(dataset['cube'], season)

    if yearly:
        dataset['cube'] = average_year(dataset['cube'], season=season)
    return dataset


def average_year_cube(cube, season=None):
    """DUMMY DOC-STRING"""
    if season:
        if not cube.coords('season_year'):
            iris.coord_categorisation.add_season_year(cube, 'time')
        mean = cube.aggregated_by('season_year', iris.analysis.MEAN)
    else:
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        mean = cube.aggregated_by('year', iris.analysis.MEAN)
    return mean


def average_year(cubes, season=None):
    """DUMMY DOC-STRING"""
    if season:
        logger.info("Calculating %s averages", season)
    else:
        logger.info("Calculating yearly averages")
    func = functools.partial(average_year_cube, season=season)
    cubes = list(map(func, cubes))
    return cubes


def calculate_reference_values(dataset, yearly=False, season=None,
                          historical_key=None, normby='run'):
    """Calculate reference values.

    model:  *per model*, so that each realization is
    scaled to the reference period following the average *model* reference
    value.

    """
    # Undefined variable default_config
    if not historical_key:
        historical_key = default_config['data']['historical_experiment']

    index = dataset['index_match_run']
    hindex = index[index > -1]
    future_data = dataset[index > -1].copy()
    future_data['match_historical_run'] = dataset.loc[hindex, 'cube'].array

    models = dataset['model'].unique()

    # Obtain reference_values  based on different settings
    if normby == 'model':
        reference_values = reference_value_by_model(future_data, models)
        dataset['reference_value'] = dataset['model'].map(reference_values)
    elif normby == 'experiment':
        reference_values = reference_value_by_experiment(future_data, models)
        dataset['reference_value'] = pd.concat(reference_values)
    else:
        reference_values = reference_value_by_run(future_data, models)
        dataset['reference_value'] = pd.concat(reference_values)
    return dataset


def calc_mean(cubes, mindata, model):
    """DUMMY DOC-STRING"""
    constraint = kcsutils.make_year_constraint_all_calendars(*REFERENCE_PERIOD)
    averages = []
    for cube in cubes:
        calendar = cube.coord('time').units.calendar
        excube = constraint[calendar].extract(cube)
        if excube is None:
            logger.warning("A cube of %s does not support time range: %s",
                            model, cube.coord('time'))
            continue
        ndata = len(excube.coord('time').points)
        # Not enough of data? Ignore
        if ndata < mindata:
            logger.warning("A cube of %s has only %d data points for its time range",
                            model, ndata)
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning,
                message="Collapsing a non-contiguous coordinate. "
                "Metadata may not be fully descriptive for 'year'")
            averages.append((excube.collapsed('time', iris.analysis.MEAN), ndata))
    return averages


def calc_refvalue(cubes, histcubes, model, yearly=False, season=None):
    """DUMMY DOC-STRING"""
    if yearly:
        mindata = MINDATA
    elif season:  # in ['djf', 'mam', 'jja', 'son']:
        # Three months a year
        mindata = {key: 3*value for key, value in MINDATA.items()}
    else:
        # Twelve months a year
        mindata = {key: 12*value for key, value in MINDATA.items()}

    avs = {}
    avs['historical'] = calc_mean(histcubes, mindata['historical'], model)
    avs['future'] = calc_mean(cubes, mindata['future'], model)
    if not avs['historical'] or not avs['future']:  # Too few data to calculate a decent bias
        logger.warning("%s does not have enough data to compute a reference", model)
        return None
    logger.info("Calculating time-weighted reference value")
    ndata = {}
    mean = {}
    for key, values in avs.items():
        n = len(values)
        # Weighted time average for each section
        ndata[key] = sum(value[1] for value in values) / n
        mean[key] = sum(value[0].data for value in values) / n
    logger.debug("Reference data values: %s , with weights %s", pformat(mean), pformat(ndata))
    value = ((mean['historical'] * ndata['historical'] +
                mean['future'] * ndata['future']) /
                (ndata['historical'] + ndata['future']))

    return value


def reference_value_by_model(dataset, models):
    values = {}
    for model in models:
        selected_dataset = dataset[dataset['model'] == model]
        cubes = selected_dataset['cube']
        histcubes = selected_dataset['match_historical_run']
        values[model] = calc_refvalue(cubes, histcubes, model)
    return values


def reference_value_by_experiment(dataset, models):
    values = []
    for model in models:
        value = {}
        selected_dataset = dataset[dataset['model'] == model]
        for exp, group in selected_dataset.groupby('experiment'):
            cubes = group['cube']
            histcubes = group['match_historical_run']
            value[exp] = calc_refvalue(cubes, histcubes, model)
        experiments = dataset.loc[
                (dataset['model'] == model),
                'experiment']
        values.append(experiments.map(value))
    return values


def reference_value_by_run(dataset, models):
    values = []
    for model in models:
        value = {}
        selected_dataset = dataset[dataset['model'] == model]
        for index, row in selected_dataset.iterrows():
            cubes = [row['cube']]
            histcubes = [row['match_historical_run']]
            value[index] = calc_refvalue(cubes, histcubes, model)
        values.append(pd.Series(value))
    return values


def normalize(dataset, relative=False, normby='run'):
    """Normalize cubes to their reference values

    If the value is a relative value, i.e., a percentual change, set
    the 'relative' parameter to `True`.

    This returns a new, modified dataset.

    """

    dataset['matched_exp'] = ''
    if normby != 'model':
        # Add the (double/triple/etc) historical runs for the future experiments,
        # and get rid of the old historical runs
        # We do this by
        # - selecting the indices & reference values for the matching runs
        # - copy the relevant data from the dataset
        #   We need to copy the Iris cubes, otherwise we'll have
        #   multiple references to the same cube
        # - Set the reference values and the matching future indices
        #   for the newly copied historical data
        # - Concatenate the current and new datasets, and drop all
        #   rows that don't have a reference value (original
        #   historical rows)

        sel = dataset['index_match_run'] > -1
        indices = dataset.loc[sel, 'index_match_run']

        reference_values = dataset.loc[sel, 'reference_value']
        hist_data = dataset.loc[indices, :]
        hist_data['cube'] = hist_data['cube'].apply(lambda cube: cube.copy())
        hist_data['reference_value'] = reference_values.array
        hist_data['index_match_run'] = dataset.loc[sel].index.array
        hist_data['matched_exp'] = dataset.loc[sel, 'experiment'].array
        dataset = pd.concat([dataset, hist_data])
        dataset.dropna(axis=0, subset=['reference_value'], inplace=True)

    # Normalization to the reference period
    for cube, value in zip(dataset['cube'], dataset['reference_value']):
        cube.data -= value
        if relative:
            cube.data /= value
            cube.data *= 100
            cube.units = '%'

    return dataset


def calculate_percentiles(dataset, average_experiments=False):
    """Calculate the percentile yearly change distribution for the input data"""
    time = []
    percentile_as_list = []
    for year in list(range(*PERIOD)):
        constraint = iris.Constraint(year=kcsutils.EqualConstraint(year))
        if average_experiments:
            data = []
            for _, group in dataset.groupby(['model', 'experiment', 'matched_exp']):
                cubes = list(filter(None, map(constraint.extract, group['cube'])))
                if cubes:
                    data.append(np.mean([cube.data for cube in cubes]))
        else:
            cubes = list(filter(None, map(constraint.extract, dataset['cube'])))
            data = [cube.data for cube in cubes]
        mean = np.mean(data)
        percentile = np.percentile(data, [5, 10, 25, 50, 75, 90, 95], overwrite_input=True)
        percentile_as_dict = dict(zip(['mean', '5', '10', '25', '50', '75', '90', '95'],
                                    [mean] + percentile.tolist()))
        percentile_as_list.append(percentile_as_dict)
        time.append(datetime(year, 1, 1))

    percentile_as_dataframe = pd.DataFrame(percentile_as_list, index=pd.DatetimeIndex(time))
    return percentile_as_dataframe


def main(cfg):
    """Return percentiles of temperature at each year.

    Four steps:
    1- obtain average values per each year,
    2- obtain the total average that is called reference_value,
    3- normalize temperature using refernce_value,
    4- calculate percentile for 5-95 perentage per each year.

    """
    # Group input files based on filenames
    files = group_metadata(cfg['input_data'].values(), 'filename')

    # Make two lists of paths and cubes for the function get_attrs()
    paths = []
    cubes = []
    for filename in files.keys():
        path = Path(filename)
        cube = iris.load_cube(str(path))
        add_extra_coord(cube)
        cubes.append(cube)
        paths.append(path)

    # Get the attributes, and create a dataframe with cubes & attributes
    dataset = get_attributes(
        cubes, paths, info_from=('attributes', 'filename'),
        attributes=None, filename_pattern=None
    )

    dataset = kcsutils.match(
        dataset, match_by='ensemble', on_no_match='randomrun',
        historical_key='historical')

    # Perform extracting of season, or averaging of years (optional).
    dataset = perform_optional_operations(dataset, yearly=False, season=None)

    # Calculate reference values
    dataset = calculate_reference_values(
        dataset, yearly=False, season=None,
        historical_key="historical", normby='run'
    )

    # Normalize data using reference values
    dataset = normalize(dataset, relative=False, normby='run')

    # Get the percentiles
    percentiles = calculate_percentiles(dataset, average_experiments=False)

    # Write output as a csv file
    percentiles.to_csv(cfg['output'], index_label="date")

    # Make provenance records
    provenance = create_provenance_record()
    provenance['ancestors'] = [path.name for path in paths]
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(cfg['output'], provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
