import logging

import iris

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic


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
        logger.info(f"{cube.standard_name}, {cube.units}, {cube.data}")
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


def main(config):
    """Radiation budget for HadGEM3 vs UKESM1.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    for dataset, group in datasets.items():
        # 'dataset' is the name of the dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info(f"Processing data for {dataset}")
        filenames = get_filenames(group)
        _ = load_cubes(filenames)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
