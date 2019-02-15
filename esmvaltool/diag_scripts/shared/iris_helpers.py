"""Convenience functions for :mod:`iris` objects."""
import logging

import iris
import numpy as np

from ._base import group_metadata

logger = logging.getLogger(__name__)


def check_coordinate(cubes, coord_name):
    """Compare coordinate of cubes and raise error if not identical.

    Parameters
    ----------
    cubes : list of iris.cube.Cube
        Cubes to be compared.
    coord_name : str
        Name of the coordinate.

    Returns
    -------
    numpy.array
        Points of the coordinate.

    Raises
    ------
    ValueError
        Given coordinate differs for the input cubes.

    """
    coord = None
    for cube in cubes:
        try:
            new_coord = cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            raise ValueError("'{}' is not a coordinate of cube {}".format(
                coord_name, cube))
        if coord is None:
            coord = new_coord
        else:
            if new_coord != coord:
                raise ValueError("Expected cubes with identical coordinates "
                                 "'{}', got {} and {}".format(
                                     coord_name, new_coord, coord))
    logger.debug("Successfully checked coordinate '%s' of cubes %s",
                 coord_name, cubes)
    return coord.points


def convert_to_iris(dict_):
    """Change all appearances of `short_name` to `var_name`.

    Parameters
    ----------
    dict_ : dict
        Dictionary to convert.

    Returns
    -------
    dict
        Converted dictionary.

    """
    dict_ = dict(dict_)
    if 'short_name' in dict_:
        dict_['var_name'] = dict_.pop('short_name')
    return dict_


def iris_project_constraint(projects, cfg, negate=False):
    """Create `iris.Constraint` to select specific projects from data.

    Parameters
    ----------
    projects : list of str
        Projects to be selected.
    cfg : dict
        Diagnostic script configuration.
    negate : bool, optional (default: False)
        Negate constraint (`False`: select all elements that fit `projects`,
        `True`: select all elements that do NOT fit `projects`).

    Returns
    -------
    iris.Constraint
        constraint for coordinate `dataset`.

    """
    datasets = []
    grouped_data = group_metadata(cfg['input_data'].values(), 'project')
    for project in projects:
        for data in grouped_data[project]:
            datasets.append(data['dataset'])

    def project_constraint(cell):
        """Constraint function."""
        if negate:
            return cell not in datasets
        return cell in datasets

    return iris.Constraint(dataset=project_constraint)


def match_dataset_coordinates(cubes):
    """Compare dataset coordinates of cubes and match them if necessary.

    Parameters
    ----------
    cubes : list of iris.cube.Cube
        Cubes to be compared.

    Returns
    -------
    list of iris.cube.Cube
        Transformed cubes.

    """
    common_elements = None
    new_cubes = []
    for cube in cubes:
        if common_elements is None:
            common_elements = set(cube.coord('dataset').points)
        else:
            common_elements = common_elements.intersection(
                set(cube.coord('dataset').points))
    common_elements = list(common_elements)
    for cube in cubes:
        cube = cube.extract(iris.Constraint(dataset=common_elements))
        sorted_idx = np.argsort(cube.coord('dataset').points)
        new_cubes.append(cube[sorted_idx])
    check_coordinate(new_cubes, 'dataset')
    logger.debug(
        "Successfully matched 'dataset' coordinates of cubes %s to "
        "%s", cubes, sorted(common_elements))
    return new_cubes


def unify_1d_cubes(cubes, coord_name):
    """Unify 1D cubes by transforming them to identical coordinate.

    Use longest coordinate as reference and transform other cubes to it by
    adding missing values.

    Parameters
    ----------
    cubes : dict of iris.cube.Cube
        Cubes to be processed.
    coord_name : str
        Name of the coordinate.

    Returns
    -------
    dict of iris.cube.Cube
        Transformed cubes.

    Raises
    ------
    ValueError
        Cubes are not 1D, coordinate name differs or not all cube coordinates
        are subsets of longest coordinate.

    """
    ref_coord = None
    for cube in cubes.values():
        if cube.ndim != 1:
            raise ValueError("Dimension of cube {} is not 1".format(cube))
        try:
            new_coord = cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            raise ValueError("'{}' is not a coordinate of cube {}".format(
                coord_name, cube))
        if ref_coord is None:
            ref_coord = new_coord
        else:
            if ref_coord.shape[0] < new_coord.shape[0]:
                ref_coord = new_coord
    if coord_name == 'time':
        iris.util.unify_time_units(cubes.values())

    # Transform all cubes
    new_cubes = {}
    for (key, cube) in cubes.items():
        coord = cube.coord(coord_name)
        if not np.all(np.isin(coord.points, ref_coord.points)):
            raise ValueError(
                "Coordinate '{}' of cube\n{}\nis not subset of reference "
                "coordinate (longest coordinate in list of cubes)".format(
                    coord_name, cube))
        new_data = np.full(ref_coord.shape, np.nan)
        indices = np.where(np.in1d(ref_coord.points, coord.points))
        new_data[indices] = np.ma.filled(cube.data, np.nan)
        new_cube = iris.cube.Cube(
            np.ma.masked_invalid(new_data),
            aux_coords_and_dims=[(ref_coord, 0)])
        new_cube.metadata = cube.metadata
        new_cubes[key] = new_cube
    check_coordinate(new_cubes.values(), 'year')
    logger.debug("Successfully unified 1D coordinate '%s' of cubes %s to %s",
                 coord_name, cubes.values(), ref_coord)
    return new_cubes
