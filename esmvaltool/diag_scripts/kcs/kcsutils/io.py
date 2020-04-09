"""DUMMY DOCSTRING"""

import pathlib
import logging
import pandas as pd
import iris
import iris.cube
import iris.coord_categorisation
import iris.analysis
from iris.util import unify_time_units
try:
    from iris.util import equalise_attributes
except ImportError:   # Iris 2
    from iris.experimental.equalise_cubes import equalise_attributes
import iris.exceptions
from .constraints import CoordConstraint
from .attributes import get as get_attrs


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


def read_data(paths, info_from=('attributes', 'filename'),
              attributes=None, filename_pattern=None):
    """DUMMY DOC-STRING"""
    cubes = [iris.load_cube(str(path)) for path in paths]

    for cube in cubes:
        if not cube.coords('season'):
            iris.coord_categorisation.add_season(cube, 'time')
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        if not cube.coords('season_year'):
            iris.coord_categorisation.add_season_year(cube, 'time',
                                                      name='season_year')
    # Get the attributes, and create a dataframe with cubes & attributes
    dataset = get_attrs(
        cubes, paths, info_from=info_from,
        attributes=attributes, filename_pattern=filename_pattern)

    return dataset


def load_cube(paths, variable_name=None):
    """Read datasets from paths into Iris cubes.

    Combines cubes if there are more than one dataset in the same file.

    Returns a list of lists. Inner lists corresponds to the areas (in
    order), outer lists corresponds to the paths

    """

    if isinstance(paths, (str, pathlib.Path)):
        if variable_name:
            cubes = iris.load_cubes(str(paths), constraints=variable_name)
        else:
            cubes = iris.load_cubes(str(paths))
    else:
        if variable_name:
            cubes = iris.load([str(path) for path in paths], constraints=variable_name)
        else:
            cubes = iris.load([str(path) for path in paths])
    # Select only the cubes with 3/4D data (time, lat, long, height)
    cubes = iris.cube.CubeList([cube for cube in cubes if len(cube.coords()) >= 3])

    if len(cubes) == 0:
        return None
    equalise_attributes(cubes)
    unify_time_units(cubes)

    try:
        cube = cubes.concatenate_cube()
    except iris.exceptions.ConcatenateError as exc:
        logger.warning("%s for %s", exc, str(paths))
        logger.warning("Using only the first cube of [%s]", cubes)
        cube = cubes[0]  # iris.load always returns a cubelist, so just take the first element
    return cube


def extract_areas(cube, areas=None, targetgrid=None, average_area=True, gridscheme='area'):
    """Regrid, extract and average multiple areas from a cube for a given variable"""

    if areas is None:
        areas = ['global']
    if isinstance(areas, str):
        areas = [areas]

    if targetgrid is not None:
        if gridscheme == 'area':
            scheme = iris.analysis.AreaWeighted()
        elif gridscheme == 'linear':
            scheme = iris.analysis.Linear()
        else:
            scheme = iris.analysis.Linear()
        gridcube = cube.regrid(targetgrid, scheme)
    else:
        gridcube = cube

    results = []
    for area in areas:
        excube = gridcube.copy()
        if area is not None:
            if isinstance(area, iris.Constraint):
                excube = excube.extract(area)
            elif isinstance(area, dict) and 'latitude' in area and 'longitude' in area:
                # To do: this if-statement needs to be rewritten; note the pylint/flake8 comments
                if (isinstance(area['latitude'], (list, tuple)) and len(area['latitude']) == 2 and
                    isinstance(area['longitude'], (list, tuple)) and len(area['longitude']) == 2):  # pylint: disable=bad-continuation,line-too-long  # noqa: E129, E501
                    long_constraint = CoordConstraint(area['longitude'][0], area['longitude'][1])
                    lat_constraint = CoordConstraint(area['latitude'][0], area['latitude'][1])
                    constraint = iris.Constraint(longitude=long_constraint, latitude=lat_constraint)
                    excube = excube.extract(constraint)
                else:
                    coords = [('latitude', area['latitude']), ('longitude', area['longitude'])]
                    excube = excube.interpolate(coords, iris.analysis.Linear())
        if excube is None:
            results.append(None)
            continue

        iris.coord_categorisation.add_season(excube, 'time')
        iris.coord_categorisation.add_season_year(excube, 'time')
        iris.coord_categorisation.add_year(excube, 'time')

        if average_area and (len(excube.coord('latitude').points) > 1 or
                             len(excube.coord('longitude').points) > 1):
            weights = iris.analysis.cartography.area_weights(excube)
            excube_meanarea = excube.collapsed(['latitude', 'longitude'],
                                               iris.analysis.MEAN, weights=weights)
        else:
            excube_meanarea = excube
        results.append(excube_meanarea)

    return results


def concat_cubes(dataset, historical_key=None):
    """Concatenate cubes into a dataset spanning the full time frame"""

    if historical_key is None:
        historical_key = 'history'

    concatenated = pd.DataFrame(columns=dataset.columns)
    for model, group in dataset.groupby('model'):
        future_sel = group['experiment'] != historical_key
        for _, row in group.loc[future_sel, :].iterrows():
            cube = row['cube']
            matchid = row['index_match_run']
            histcube = dataset.loc[matchid, 'cube']
            time = cube.coord('time')
            start = time.units.num2date(time.points)[0]
            time = histcube.coord('time')
            end = time.units.num2date(time.points)[-1]
            if end > start:
                logger.warning("Historical experiment ends past the start of future experiment: "
                               "trimming historical dataset to match %s - %s r%di%dp%d",
                               model, row['experiment'],
                               row['realization'], row['initialization'], row['physics'])
                logger.debug("Historical end: %s. Future start: %s", end, start)
                logger.info("Constraining historical run to end before %s", start)
                # pylint: disable=cell-var-from-loop
                constraint = iris.Constraint(time=lambda cell: cell.point < start)
                histcube = constraint.extract(histcube)
                time = histcube.coord('time')
                end = time.units.num2date(time.points)[-1]
            # Since we may have matched on different realizations, set them to
            # equal, otherwise the concatenation procedure will fail
            try:
                histcube.replace_coord(cube.coord('realization'))
            except iris.exceptions.CoordinateNotFoundError:
                pass
            cubes = iris.cube.CubeList([histcube, cube])
            equalise_attributes(cubes)
            unify_time_units(cubes)
            try:
                row['cube'] = cubes.concatenate_cube()
            except iris.exceptions.ConcatenateError as exc:
                logger.warning("DATA SKIPPED: Error concatenating cubes: %s - %s r%di%dp%d: %s",
                               model, row['experiment'], row['realization'],
                               row['initialization'], row['physics'], exc)
                continue
            logger.info("Concatenated %s - %s r%di%dp%d", model, row['experiment'],
                        row['realization'], row['initialization'], row['physics'])
            # By adding the rows with concatenated cubes,
            # we can construct a new DataFrame with only concatenated cubes,
            # but with all the essential attributes from the input dataset.
            concatenated = concatenated.append(row)
    concatenated.reset_index(inplace=True)
    logger.info("concatenated a total of %d realizations", len(concatenated))

    return concatenated
