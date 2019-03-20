"""Convenience functions for :mod:`iris` objects."""
import logging
from pprint import pformat

import iris
import numpy as np

from ._base import group_metadata

logger = logging.getLogger(__name__)


def _transform_coord_to_ref(cubes, ref_coord):
    """Transform coordinates of cubes to reference."""
    ref_coord = iris.coords.DimCoord.from_coord(ref_coord)
    coord_name = ref_coord.name()
    cubes_iterable = enumerate(cubes)
    new_cubes = list(range(len(cubes)))
    if isinstance(cubes, dict):
        cubes_iterable = cubes.items()
        new_cubes = {}
    for (key, cube) in cubes_iterable:
        coord = cube.coord(coord_name)
        if not np.all(np.isin(coord.points, ref_coord.points)):
            raise ValueError(
                "Coordinate '{}' of cube\n{}\nis not subset of reference "
                "coordinate (longest coordinate in iterable of cubes)".format(
                    coord_name, cube))
        new_data = np.full(ref_coord.shape, np.nan)
        indices = np.where(np.in1d(ref_coord.points, coord.points))
        new_data[indices] = np.ma.filled(cube.data, np.nan)
        new_cube = iris.cube.Cube(
            np.ma.masked_invalid(new_data),
            dim_coords_and_dims=[(ref_coord, 0)])
        for aux_coord in cube.coords(dim_coords=False):
            if aux_coord.shape in ((), (1, )):
                new_cube.add_aux_coord(aux_coord, [])
        new_cube.metadata = cube.metadata
        new_cubes[key] = new_cube
    check_coordinate(new_cubes, coord_name)
    logger.debug("Successfully unified 1D coordinates '%s' to %s", coord_name,
                 ref_coord)
    logger.debug("of cubes")
    logger.debug(pformat(cubes))
    if isinstance(cubes, iris.cube.CubeList):
        return iris.cube.CubeList(new_cubes)
    return new_cubes


def check_coordinate(cubes, coord_name):
    """Compare coordinate of cubes and raise error if not identical.

    Parameters
    ----------
    cubes : iterable or dict of iris.cube.Cube
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
    if isinstance(cubes, dict):
        cubes = list(cubes.values())
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
    logger.debug("Successfully checked coordinate '%s' of cubes", coord_name)
    logger.debug(pformat(cubes))
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
        for data in grouped_data.get(project, {}):
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
    cubes : iterable or dict of iris.cube.Cube
        Cubes to be compared.

    Returns
    -------
    iterable or dict of iris.cube.Cube
        Transformed cubes.

    """
    common_elements = None
    cubes_iterable = cubes
    if isinstance(cubes, dict):
        cubes_iterable = cubes.values()

    # Get common elements
    for cube in cubes:
        if common_elements is None:
            common_elements = set(cube.coord('dataset').points)
        else:
            common_elements = common_elements.intersection(
                set(cube.coord('dataset').points))
    common_elements = list(common_elements)

    # Save new cubes
    cubes_iterable = enumerate(cubes)
    new_cubes = list(range(len(cubes)))
    if isinstance(cubes, dict):
        cubes_iterable = cubes.items()
        new_cubes = {}
    for (key, cube) in cubes_iterable:
        cube = cube.extract(iris.Constraint(dataset=common_elements))
        sorted_idx = np.argsort(cube.coord('dataset').points)
        new_cubes[key] = cube[sorted_idx]
    check_coordinate(new_cubes, 'dataset')
    logger.debug("Successfully matched 'dataset' coordinate to %s",
                 sorted(common_elements))
    logger.debug("of cubes")
    logger.debug(pformat(cubes))
    if isinstance(cubes, iris.cube.CubeList):
        return iris.cube.CubeList(new_cubes)
    return new_cubes


def unify_1d_cubes(cubes, coord_name):
    """Unify 1D cubes by transforming them to identical coordinate.

    Use longest coordinate as reference and transform other cubes to it by
    adding missing values.

    Parameters
    ----------
    cubes : iterable or dict of iris.cube.Cube
        Cubes to be processed.
    coord_name : str
        Name of the coordinate.

    Returns
    -------
    iterable or dict of iris.cube.Cube
        Transformed cubes.

    Raises
    ------
    ValueError
        Cubes are not 1D, coordinate name differs or not all cube coordinates
        are subsets of longest coordinate.

    """
    ref_coord = None
    cubes_iterable = cubes
    if isinstance(cubes, dict):
        cubes_iterable = cubes.values()

    # Get reference coordinate
    for cube in cubes_iterable:
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
        iris.util.unify_time_units(cubes_iterable)

    # Transform all cubes
    return _transform_coord_to_ref(cubes, ref_coord)
