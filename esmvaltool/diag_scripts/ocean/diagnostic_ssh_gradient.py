"""
Sea Surface Height gradient diagnostic
======================================

Diagnostic to produce maps of the horizontal gradient magnitude of Sea Surface
Height (SSH, ``zos``). The gradient magnitude highlights strong surface
currents and frontal structures and is used to evaluate the Gulf Stream
separation bias against a reference dataset (e.g. ORAS5).

This diagnostic expects the preprocessor to have regridded ``zos`` onto a
regular latitude/longitude grid, **keeping the time axis**. The gradient
magnitude is computed for every time step and then time-averaged, i.e.
``mean(|grad(SSH)|)`` rather than ``|grad(mean(SSH))|`` -- these differ for this
nonlinear quantity, and the former preserves frontal sharpness. An appropriate
preprocessor is::

  preprocessors:
    create_ssh_regridded:
      regrid:
        target_grid: 2x2
        scheme: linear

The gradient magnitude is computed as::

  |grad(SSH)| = sqrt( (d SSH / dx)^2 + (d SSH / dy)^2 )

with the horizontal metric terms ``dx = R cos(lat) dlon`` and ``dy = R dlat``.

Note that this diagnostic may not function on machines with no access to the
internet, as cartopy may try to download coastline shapefiles. See
``diagnostic_maps.py`` for how to provide them offline via
``auxiliary_data_dir``.

Author: Kevin Debeire (DLR).
"""

import logging
import sys
from pathlib import Path

import cartopy
import iris
import iris.analysis.calculus
import iris.analysis.cartography
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_data,
    save_figure,
)

# This part sends debug statements to stdout
logger = logging.getLogger(Path(__file__).name)
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

EARTH_RADIUS = 6371000.0  # m


def compute_ssh_gradient_magnitude(cube):
    """Compute the time-mean horizontal SSH gradient magnitude.

    The gradient magnitude is computed for every time step (if the cube has a
    time dimension) and then averaged over time, i.e. ``mean(|grad(SSH)|)``.

    Parameters
    ----------
    cube: iris.cube.Cube
        SSH field with ``latitude`` and ``longitude`` coordinates on a regular
        grid, optionally with a leading ``time`` dimension.

    Returns
    -------
    iris.cube.Cube
        2D cube of the time-mean gradient magnitude (units ``m m-1``) on the
        interior grid (one point smaller than the input in each horizontal
        direction).
    """
    # Finite-difference gradients (each reduces one horizontal dim by 1). The
    # leading `...` keeps any time dimension so the gradient is per time step.
    # A circular longitude coordinate (as produced by regridding to a global
    # target grid) makes `differentiate` wrap and keep the same length; disable
    # it so longitude reduces by one point, consistent with latitude.
    cube = cube.copy()
    lon_coord = cube.coord("longitude")
    if getattr(lon_coord, "circular", False):
        lon_coord.circular = False
    grad_x = iris.analysis.calculus.differentiate(cube, "longitude")
    grad_y = iris.analysis.calculus.differentiate(cube, "latitude")

    # Align both gradients (and the target cube) onto the common interior grid.
    gx_data = grad_x.data[..., :-1, :]
    gy_data = grad_y.data[..., :, :-1]
    ssh_grad = cube[..., :-1, :-1]

    deg_to_rad = np.pi / 180.0
    dy = EARTH_RADIUS * deg_to_rad
    cos_lat = iris.analysis.cartography.cosine_latitude_weights(ssh_grad)
    dx = EARTH_RADIUS * cos_lat * deg_to_rad

    ssh_grad.data = np.sqrt((gx_data / dx) ** 2 + (gy_data / dy) ** 2)
    ssh_grad.var_name = "ssh_grad_mag"
    ssh_grad.long_name = "SSH Gradient Magnitude"
    ssh_grad.units = "m m-1"

    # Time-average the magnitude, i.e. mean(|grad|), not |grad(mean)|.
    if ssh_grad.coords("time", dim_coords=True):
        ssh_grad = ssh_grad.collapsed("time", iris.analysis.MEAN)
    return ssh_grad


def make_ssh_gradient_plot(cfg, metadata, filename):
    """Compute and plot the SSH gradient magnitude for one dataset.

    Parameters
    ----------
    cfg: dict
        The opened global config dictionary, passed by ESMValTool.
    metadata: dict
        The metadata dictionary for this dataset.
    filename: str
        The preprocessed SSH climatology file.
    """
    cube = iris.load_cube(filename)
    ssh_grad = compute_ssh_gradient_magnitude(cube)

    dataset = metadata["dataset"]
    basename = dataset + "_ssh_gradient"
    provenance_record = {
        "caption": "SSH gradient magnitude map for " + dataset + ".",
        "statistics": ["mean"],
        "domain": ["global"],
        "plot_types": ["map"],
        "ancestors": [filename],
    }

    # Save the plotted data as a NetCDF file in the work directory.
    save_data(basename, provenance_record, cfg, ssh_grad)

    # Plot the map and save the figure (png/pdf as configured).
    plt.figure(figsize=(10, 5))
    qplt.contourf(ssh_grad, 25, linewidth=0, rasterized=True)
    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning("Not able to add coastlines")
    plt.title(dataset + " SSH gradient magnitude")
    save_figure(basename, provenance_record, cfg, dpi=150)


def main(cfg):
    """Load the config file and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        The opened global config dictionary, passed by ESMValTool.
    """
    cartopy.config["data_dir"] = cfg["auxiliary_data_dir"]

    for index, metadata_filename in enumerate(cfg["input_files"]):
        logger.info("metadata filename:\t%s", metadata_filename)
        metadatas = diagtools.get_input_files(cfg, index=index)
        for filename in sorted(metadatas.keys()):
            logger.info("-----------------")
            logger.info("model filenames:\t%s", filename)
            make_ssh_gradient_plot(cfg, metadatas[filename], filename)

    logger.info("Success")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
