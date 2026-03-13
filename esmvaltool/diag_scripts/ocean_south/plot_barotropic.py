"""diagnostic script to plot"""

import iris
import os
import logging
import numpy as np
import iris.plot as iplt
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
    group_metadata,
    select_metadata,
)
from esmvalcore.preprocessor import (
    climate_statistics,
    regrid
)


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def circumpolar_map():
    fig = plt.figure(figsize = (12, 8))
    ax = plt.axes(projection = ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -80, -40], crs = ccrs.PlateCarree())
    ax.set_facecolor('lightgrey')
    # Map the plot boundaries to a circle
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform = ax.transAxes)

    return fig, ax


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""

    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["polar"],
        "plot_types": ["map"],
        "authors": [
            "chun_felicity",
        ],
        "references": [
            "access-nri",
        ],
        "ancestors": ancestor_files,
    }
    return record


def main(cfg):
    """# in script cumsum, units handling, climate mean"""

    input_data = cfg["input_data"].values()

    for dataset in input_data:
        # Load the data
        input_file = dataset['filename']
        name = dataset['dataset']
        cube = iris.load_cube(input_file)

        # units cumulative sum and climate mean
        p0 = 1035.0
        cube = cube.collapsed('depth', iris.analysis.SUM)
        cube.data = cube.data / p0  # divide by density for volume
        cube.data = cube.data / 1e6  # then divide by 10^6 for Sv
        cube = regrid(cube, target_grid='0.5x0.5', scheme="linear")
        cube.data = cube.data.cumsum(1)  #latitude index 1
        logger.info(cube.summary(shorten=True))
        cube.units = 'Sv'
        cube = climate_statistics(cube, period='full', operator='mean')

        # plot
        fig, ax = circumpolar_map()

        p1 = iplt.contourf(cube,
                    levels = np.arange(-50, 150, 10),
                    extend = 'both',
                    cmap = cmocean.cm.speed,
                    )
        plt.colorbar(p1, orientation='vertical', label='$\psi$ (Sv)')

        plt.title(f'Barotropic streamfunction from {name}')
        # Save output
        prov_record = get_provenance_record(
            f'Barotropic streamfunction from {name}.',
            [input_file],
        )
        save_figure(
            f"barotropic_streamfunction_{name}",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
            )

if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
