"""Common code for radiation budget diagnostic scripts."""

import iris


def load_data(filenames):
    """Return the loaded cubes.

    Parameters
    ----------
    filenames : list of strings
        The filenames to load.

    Returns
    -------
    :class:`iris.cube.Cube`
        The loaded cubes.
    """
    cubes = iris.load(filenames)
    return cubes


def get_filenames(group):
    """Return all the filenames for the group.

    Parameters
    ----------
    group : list(dict)
        The grouped metadata describing preprocessed data.

    Returns
    -------
    list
        All the filenames for the group.
    """
    filenames = [item["filename"] for item in group]
    return filenames
