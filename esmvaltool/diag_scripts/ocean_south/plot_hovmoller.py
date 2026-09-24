"""diagnostic script to plot"""

import iris
import os
import logging
import numpy as np
import iris.plot as iplt
from matplotlib.gridspec import GridSpec
from matplotlib import ticker
import matplotlib.pyplot as plt

import cmocean
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
    group_metadata,
    select_metadata,
)
from esmvalcore.preprocessor import (
    area_statistics,
    anomalies
)


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def plot_hovmoller(fsize = 14):
    """Setup of figures properties."""
    plt.rcParams['font.size'] = fsize
    plt.rcParams['xtick.labelsize'] = fsize-2
    plt.rcParams['ytick.labelsize'] = fsize-2
    
    fig = plt.figure(figsize = (8, 6))
    grid = GridSpec(100, 100)
    
    ax = [fig.add_subplot(grid[:30, :45]), fig.add_subplot(grid[:30, 48:93]),
          fig.add_subplot(grid[32:, :45]), fig.add_subplot(grid[32:, 48:93])]
    
    for i in range(len(ax)):
        # ax[i].xaxis.set_major_formatter(date_format)
        ax[i].tick_params(axis='x', labelrotation=45)
        if i != 0 and i != 2:
            ax[i].set_yticklabels([]) 

    return fig, ax


def hovmoller(ds_cube, startyr):
    """Calculate anomalies summed over area for depth vs time plot."""
    cube = anomalies(ds_cube, period='full', reference={'start_year':startyr, 'start_month':1,'start_day':1,
                                                            'end_year':startyr, 'end_month':12, 'end_day':31})
    # anomalies doesn't retain cell measure required for area statistics
    cube.add_cell_measure(ds_cube.cell_measure('cell_area'), data_dims=[2,3])
    # sum anomalies
    cube = area_statistics(cube, operator='sum')

    #divide by total area cube (anomalies are summed.. weighted by cell_area(multiplied)
    tarea = totalarea_depth(ds_cube)
    return cube / tarea


def totalarea_depth(cube):
    """Calculate total area at each depth level."""
    tmp1 = cube[0] #get first time step
    # divide by itself for 1s mask 3-dimensional
    tmp1.data = tmp1.data / tmp1.data
    # check cell_measure exists?
    if cube.cell_measure('cell_area'):
        return area_statistics(tmp1, operator='sum')


def plot_details(ax, cf_temp, cf_salt):
    for i in range(len(ax)):
        if i < 2:
            ax[i].set_ylim(500, 0)
            ax[i].set_xticklabels([])
        else:
            ax[i].set_xlabel("")
            ax[i].set_ylim(5000, 500)
    ax[0].set_ylabel("Depth [m]")
    ax[1].set_ylabel("")
    ax[2].set_ylabel("Depth [m]")
    ax[3].set_ylabel("")

    # Colorbars
    bar = plt.axes([0.11, 0.97, 0.35, 0.02])
    cbar_1 = plt.colorbar(cf_temp, cax = bar, orientation = 'horizontal', extend='both', format= '%.2f')  
    cbar_1.set_label("Temperature [$\degree$C]")

    bar = plt.axes([0.50, 0.97, 0.35, 0.02])
    cbar_2 = plt.colorbar(cf_salt, cax = bar, orientation = 'horizontal', extend='both', format= '%.2f')  
    cbar_2.set_label("Salinity [psu]")

    for cbar in [cbar_1, cbar_2]:
        tick_locator = ticker.MaxNLocator(nbins=3, prune='both') ## The ticker needs to called within the loop
        cbar.locator = tick_locator
        cbar.update_ticks()

def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""

    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["global"],
        "plot_types": ["vert"],
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
    """Create Hovmoller diagrams."""

    input_data = cfg["input_data"].values()

    # group by dataset
    ds_groups = group_metadata(
        input_data,
        "dataset",
    )

    levels_tempsal = {'thetao':np.arange(-0.3, 0.31, 0.01), 'so':np.arange(-0.03, 0.031, 0.001)}
    cmapcm = {'thetao':cmocean.cm.balance, 'so':cmocean.cm.curl}
    axi = {'thetao':0, 'so':1}
    # for each ds plt ts, sal 
    for grp, var_attr in ds_groups.items():
        fig, ax = plot_hovmoller(fsize = 14)
        files=[]
        for metadata in var_attr:  #[thetao, so]
            name = metadata["dataset"]
            shortname = metadata["short_name"]
            input_file = metadata['filename']
            cube = iris.load_cube(metadata['filename'])
            hov_cube = hovmoller(cube, metadata['start_year'])

            # plot twice
            for i in [axi[shortname], axi[shortname]+2]:
                p1 = iplt.contourf(hov_cube,
                            levels = levels_tempsal[shortname],
                            extend = 'both',
                            cmap = cmapcm[shortname],
                            axes = ax[i], #0,2 temp
                            )
            files.append(input_file)

            if shortname == 'so':
                cf_salt = p1
            else:
                cf_temp = p1
        
        plot_details(ax, cf_temp, cf_salt)
        # Save output
        prov_record = get_provenance_record(
            f'Depth-Time Temperature and Salinity from {name}.',
            files,
        )
        save_figure(
            f"hovmoller_{name}",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
