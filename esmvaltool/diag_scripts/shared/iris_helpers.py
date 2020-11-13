"""Convenience functions for :mod:`iris` objects."""
import logging
from pprint import pformat

import iris
import numpy as np

from ._base import group_metadata

logger = logging.getLogger(__name__)


def _transform_coord_to_ref(cubes, ref_coord):
    """Transform coordinates of cubes to reference."""
    try:
        # Convert AuxCoord to DimCoord if necessary and possible
        ref_coord = iris.coords.DimCoord.from_coord(ref_coord)
    except ValueError:
        pass
    if not np.array_equal(np.unique(ref_coord.points), np.sort(
            ref_coord.points)):
        raise ValueError(
            f"Expected unique coordinate '{ref_coord.name()}', got "
            f"{ref_coord}")
    coord_name = ref_coord.name()
    new_cubes = iris.cube.CubeList()
    for cube in cubes:
        coord = cube.coord(coord_name)
        if not np.all(np.isin(coord.points, ref_coord.points)):
            raise ValueError(
                f"Coordinate {coord} of cube\n{cube}\nis not subset of "
                f"reference coordinate {ref_coord}")
        new_data = np.full(ref_coord.shape, np.nan)
        indices = np.where(np.in1d(ref_coord.points, coord.points))
        new_data[indices] = np.ma.filled(cube.data, np.nan)
        new_cube = iris.cube.Cube(np.ma.masked_invalid(new_data))
        if isinstance(ref_coord, iris.coords.DimCoord):
            new_cube.add_dim_coord(ref_coord, 0)
        else:
            new_cube.add_aux_coord(ref_coord, 0)
        for aux_coord in cube.coords(dim_coords=False):
            if aux_coord.shape in ((), (1, )):
                new_cube.add_aux_coord(aux_coord, [])
        new_cube.metadata = cube.metadata
        new_cubes.append(new_cube)
    check_coordinate(new_cubes, coord_name)
    logger.debug("Successfully unified coordinate '%s' to %s", coord_name,
                 ref_coord)
    logger.debug("of cubes")
    logger.debug(pformat(cubes))
    return new_cubes


def check_coordinate(cubes, coord_name):
    """Compare coordinate of cubes and raise error if not identical.

    Parameters
    ----------
    cubes : iris.cube.CubeList
        Cubes to be compared.
    coord_name : str
        Name of the coordinate.

    Returns
    -------
    numpy.array
        Points of the coordinate.

    Raises
    ------
    iris.exceptions.CoordinateNotFoundError
        Coordinate ``coord_name`` is not a coordinate of one of the cubes.
    ValueError
        Given coordinate differs for the input cubes.

    """
    coord = None
    for cube in cubes:
        try:
            new_coord = cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            raise iris.exceptions.CoordinateNotFoundError(
                f"'{coord_name}' is not a coordinate of cube\n{cube}")
        if coord is None:
            coord = new_coord
        else:
            if new_coord != coord:
                raise ValueError(
                    f"Expected cubes with identical coordinates "
                    f"'{coord_name}', got {new_coord} and {coord}")
    logger.debug("Successfully checked coordinate '%s' of cubes", coord_name)
    logger.debug(pformat(cubes))
    return coord.points


def convert_to_iris(dict_):
    """Change all appearances of ``short_name`` to ``var_name``.

    Parameters
    ----------
    dict_ : dict
        Dictionary to convert.

    Returns
    -------
    dict
        Converted dictionary.

    Raises
    ------
    KeyError
        :obj:`dict` contains keys``'short_name'`` **and** ``'var_name'``.

    """
    dict_ = dict(dict_)
    if 'short_name' in dict_:
        if 'var_name' in dict_:
            raise KeyError(
                f"Cannot replace 'short_name' by 'var_name', dictionary "
                f"already contains 'var_name' (short_name = "
                f"'{dict_['short_name']}', var_name = '{dict_['var_name']}')")
        dict_['var_name'] = dict_.pop('short_name')
    return dict_


def get_mean_cube(datasets):
    """Get mean cube of a list of datasets.

    Parameters
    ----------
    datasets : list of dict
        List of datasets (given as metadata :obj:`dict`).

    Returns
    -------
    iris.cube.Cube
        Mean cube.

    """
    cubes = iris.cube.CubeList()
    for dataset in datasets:
        path = dataset['filename']
        cube = iris.load_cube(path)
        prepare_cube_for_merging(cube, path)
        cubes.append(cube)
    mean_cube = cubes.merge_cube()
    if len(cubes) > 1:
        mean_cube = mean_cube.collapsed(['cube_label'], iris.analysis.MEAN)
    mean_cube.remove_coord('cube_label')
    return mean_cube


def iris_project_constraint(projects, input_data, negate=False):
    """Create :class:`iris.Constraint` to select specific projects from data.

    Parameters
    ----------
    projects : list of str
        Projects to be selected.
    input_data : list of dict
        List of dataset metadata used to extract all relevant datasets
        belonging to given ``projects``.
    negate : bool, optional (default: False)
        Negate constraint (``False``: select all elements that fit
        ``projects``, `True``: select all elements that do **not** fit
        ``projects``).

    Returns
    -------
    iris.Constraint
        constraint for coordinate ``dataset``.

    """
    datasets = []
    grouped_data = group_metadata(input_data, 'project')
    for project in projects:
        datasets.extend([d['dataset'] for d in grouped_data.get(project, [])])

    def project_constraint(cell):
        """Constraint function."""
        if negate:
            return cell not in datasets
        return cell in datasets

    return iris.Constraint(dataset=project_constraint)


def intersect_dataset_coordinates(cubes):
    """Compare dataset coordinates of cubes and match them if necessary.

    Use intersection of coordinate 'dataset' of all given cubes and remove
    elements which are not given in all cubes.

    Parameters
    ----------
    cubes : iris.cube.CubeList
        Cubes to be compared.

    Returns
    -------
    iris.cube.CubeList
        Transformed cubes.

    Raises
    ------
    iris.exceptions.CoordinateNotFoundError
        Coordinate ``dataset`` is not a coordinate of one of the cubes.
    ValueError
        At least one of the cubes contains a ``dataset`` coordinate with
        duplicate elements or the cubes do not share common elements.

    """
    common_elements = None

    # Get common elements
    for cube in cubes:
        try:
            coord_points = cube.coord('dataset').points
        except iris.exceptions.CoordinateNotFoundError:
            raise iris.exceptions.CoordinateNotFoundError(
                f"'dataset' is not a coordinate of cube\n{cube}")
        if len(set(coord_points)) != len(coord_points):
            raise ValueError(
                f"Coordinate 'dataset' of cube\n{cube}\n contains duplicate "
                f"elements")
        if common_elements is None:
            common_elements = set(coord_points)
        else:
            common_elements = common_elements.intersection(set(coord_points))
    common_elements = list(common_elements)

    # Save new cubes
    new_cubes = iris.cube.CubeList()
    for cube in cubes:
        cube = cube.extract(iris.Constraint(dataset=common_elements))
        if cube is None:
            raise ValueError(f"Cubes {cubes} do not share common elements")
        sorted_idx = np.argsort(cube.coord('dataset').points)
        new_cubes.append(cube[sorted_idx])
    check_coordinate(new_cubes, 'dataset')
    logger.debug("Successfully matched 'dataset' coordinate to %s",
                 sorted(common_elements))
    logger.debug("of cubes")
    logger.debug(pformat(cubes))
    return new_cubes


def prepare_cube_for_merging(cube, cube_label):
    """Prepare single :class:`iris.cube.Cube` in order to merge it later.

    Parameters
    ----------
    cube : iris.cube.Cube
        Cube to be pre-processed.
    cube_label : str
        Label for the new scalar coordinate ``cube_label``.

    """
    cube.attributes = {}
    cube.cell_methods = ()
    for coord in cube.coords(dim_coords=True):
        coord.attributes = {}
    for coord in cube.coords(dim_coords=False):
        cube.remove_coord(coord)
    cube_label_coord = iris.coords.AuxCoord(cube_label,
                                            var_name='cube_label',
                                            long_name='cube_label')
    cube.add_aux_coord(cube_label_coord, [])


def unify_1d_cubes(cubes, coord_name):
    """Unify 1D cubes by transforming them to identical coordinates.

    Use union of all coordinates as reference and transform other cubes to it
    by adding missing values.

    Parameters
    ----------
    cubes : iris.cube.CubeList
        Cubes to be processed.
    coord_name : str
        Name of the coordinate.

    Returns
    -------
    iris.cube.CubeList
        Transformed cubes.

    Raises
    ------
    ValueError
        Cubes are not 1D, coordinate name differs or not all cube coordinates
        are subsets of longest coordinate.

    """
    ref_coord = None

    # Get reference coordinate
    for cube in cubes:
        if cube.ndim != 1:
            raise ValueError(f"Dimension of cube\n{cube}\nis not 1")
        try:
            new_coord = cube.coord(coord_name)
        except iris.exceptions.CoordinateNotFoundError:
            raise iris.exceptions.CoordinateNotFoundError(
                f"'{coord_name}' is not a coordinate of cube\n{cube}")
        if not np.array_equal(np.unique(new_coord.points),
                              np.sort(new_coord.points)):
            raise ValueError(
                f"Coordinate '{coord_name}' of cube\n{cube}\n is not unique, "
                f"unifying not possible")
        if ref_coord is None:
            ref_coord = new_coord
        else:
            new_points = np.union1d(ref_coord.points, new_coord.points)
            ref_coord = ref_coord.copy(new_points)
    if coord_name == 'time':
        iris.util.unify_time_units(cubes)

    # Transform all cubes
    return _transform_coord_to_ref(cubes, ref_coord)


def var_name_constraint(var_name):
    """:class:`iris.Constraint` using ``var_name`` of an :mod:`iris.cube.Cube`.

    Parameters
    ----------
    var_name : str
        Short name (``var_name`` in :mod:`iris`) for the constraint.

    Returns
    -------
    iris.Constraint
        Constraint to select only cubes with correct ``var_name``.

    """
    return iris.Constraint(cube_func=lambda c: c.var_name == var_name)
