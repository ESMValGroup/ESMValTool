"""Convenience functions for :mod:`iris` objects."""
import logging

import iris
import numpy as np

from ._base import group_metadata

logger = logging.getLogger(__name__)


def check_dataset_coordinates(cubes):
    """Compare dataset coordinates of cubes and raise error if not identical.

    Parameters
    ----------
    cubes : list of iris.cube.Cube
        Cubes to be compared.

    Returns
    -------
    numpy.array
        All datasets.

    Raises
    ------
    ValueError
        Coordinate `dataset` differs for the input cubes.

    """
    coord = None
    for cube in cubes:
        if coord is None:
            coord = cube.coord('dataset')
        else:
            if cube.coord('dataset') != coord:
                raise ValueError("Expected cubes with identical coordinates "
                                 "'dataset', got {} and {}".format(
                                     cube.coord('dataset'), coord))
    return coord.points


def iris_project_constraint(projects, cfg):
    """Create `iris.Constraint` to select specific projects from data.

    Parameters
    ----------
    projects : list of str
        Projects to be selected.
    cfg : dict
        Diagnostic script configuration.

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

    # Constraint function
    def project_constraint(cell):
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
    logger.debug("Matched 'dataset' coordinates of cubes to %s",
                 sorted(common_elements))
    check_dataset_coordinates(new_cubes)
    return new_cubes
