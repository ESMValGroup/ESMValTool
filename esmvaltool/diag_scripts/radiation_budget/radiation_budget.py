import logging

import iris

from esmvaltool.diag_scripts.shared import run_diagnostic


def load_cubes(filenames):
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
    logger = logging.getLogger(__name__)
    cubes = iris.load(filenames)
    for cube in cubes:
        logger.info(cube)
    return cubes


def main(config):
    """Radiation budget for HadGEM3 vs UKESM1.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = config["input_data"]
    filenames = input_data.keys()
    _ = load_cubes(filenames)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
