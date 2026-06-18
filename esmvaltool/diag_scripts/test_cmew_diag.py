"""Test diagnostic to investigate variable processing."""

import logging
from pathlib import Path

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)


def loop_over_models(cfg):
    """Loop over each section in the config."""
    # Iterate over the data in the cfg
    for section in cfg["input_data"].values():

        # Read relevant facets
        dataset = section["dataset"]
        logger.info("Looking at %s", dataset)

        nc_filepath = section["filename"]
        logger.info("Found filepath %s", nc_filepath)

        # See if there's more than one cube (e.g. for wind, levels)
        cubelist = iris.load(nc_filepath)
        logger.info("Length of cubelist: %d", len(cubelist))

        # Load the cube(s)
        cubes = iris.load_cubes(nc_filepath)

        for cube in cubes:
            logger.info("Cube: %s", cube.name)
            plot_timeseries(cube, cfg)


def plot_timeseries(cube, cfg):
    """Plot timeseries."""
    # Plot the data
    fig, ax = plt.subplots()
    qplt.plot(cube.coord("time").points(), cube.data(), axes=ax)
    # Save the figure (also closes it)
    save_figure(
        cube.name,
        {},
        cfg,
        figure=fig,
        close=True,
    )


def main(cfg):
    """Create two plots and one csv file per diagnostic."""
    # Log config dictionary
    logger.info("cfg:\n%s", cfg)
    loop_over_models(cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
