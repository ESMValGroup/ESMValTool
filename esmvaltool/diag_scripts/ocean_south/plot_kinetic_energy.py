"""diagnostic script to plot"""

import iris
import os
import logging
import numpy as np
import iris.plot as iplt
import matplotlib.pyplot as plt

import cmocean
import cartopy.crs as ccrs
import cartopy.feature as feature
import matplotlib.path as mpath
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
    group_metadata,
    save_data,
)
from esmvalcore.preprocessor import (
    add_supplementary_variables,
    axis_statistics,
    regrid,
)


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def ke_calc(u_vel, v_vel, thkcel):
    KE = 0.5*(u_vel**2 + v_vel**2)
    KE = add_supplementary_variables(KE, [thkcel])
    return KE

def kinetic_energys(u_vel, v_vel, thkcel):
    """Calculate kinetic energy from u and v components."""
    
    u_mean = axis_statistics(u_vel, operator='mean', axis='t')
    v_mean = axis_statistics(v_vel, operator='mean', axis='t')
    mke = ke_calc(u_mean, v_mean, thkcel)
    mke = mke.collapsed('depth', iris.analysis.SUM, weights='cell_thickness')
    mke.long_name = 'Mean Kinetic Energy'
    # eke
    u_transient = u_vel - u_mean
    v_transient = v_vel - v_mean
    eke = ke_calc(u_transient, v_transient, thkcel) 
    eke = eke.collapsed('depth', iris.analysis.SUM, weights='cell_thickness')
    eke = axis_statistics(eke, operator='mean', axis='t')
    eke.long_name = 'Eddy Kinetic Energy'
    # tke
    tke = ke_calc(u_vel, v_vel, thkcel)
    tke = axis_statistics(tke, operator='mean', axis='t')
    tke = tke.collapsed('depth', iris.analysis.SUM, weights='cell_thickness')
    tke.long_name = 'Total Kinetic Energy'
    
    KEall = [regrid(ke, '0.5x0.5', scheme='nearest') for ke in [tke, mke, eke]]  # regrid to map
    return KEall # tke,mke,eke

def circumpolar_map(fig, i): # make 3 subplots
    # fig = plt.figure(figsize = (12, 8)) #move
    # ax = plt.axes(projection = ccrs.SouthPolarStereo())
    ax = fig.add_subplot(i, projection = ccrs.SouthPolarStereo()) # i 221
    ax.set_extent([-180, 180, -80, -50], crs = ccrs.PlateCarree())
    ax.set_facecolor('lightgrey')
    # Map the plot boundaries to a circle
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform = ax.transAxes)

    return ax


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""

    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["shpolar"],
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
    """Create kinetic energy diagrams."""

    input_data = cfg["input_data"].values()

    # group by dataset
    ds_groups = group_metadata(
        input_data,
        "dataset",
    )

    #
    for grp, var_attr in ds_groups.items():
        fig = plt.figure(figsize = (12, 10))
        
        files, data_cubes = [], {}
        for metadata in var_attr:  #{vo,uo, thkcello}
            
            shortname = metadata["short_name"]
            input_file = metadata['filename']
            cube = iris.load_cube(metadata['filename'])
            data_cubes[shortname] = cube
            files.append(input_file)
            thkness = cube.ancillary_variable('cell_thickness')

        tke, mke, eke = kinetic_energys(data_cubes['uo'], data_cubes['vo'], thkness)

        plt_titles = ['Total Kinetic Energy (TKE)', 'Mean Kinetic Energy (MKE)', 'Eddy Kinetic Energy (EKE)']
        # plot TKE, MKE, EKE (3subplts)
        for i, ke_cube in enumerate([tke, mke, eke]):
            axs = circumpolar_map(fig, 221+i)
            cf1 = iplt.contourf(ke_cube, levels=np.arange(0,26,0.5), extend='max', axes=axs, cmap = cmocean.cm.ice)
            axs.set_title(plt_titles[i])

            data_prov = get_provenance_record(
                f'{plt_titles[i]} for {grp}.',
                files,
            ) # save data out
            save_data(f"{plt_titles[i][-4:-1]}_{grp}", data_prov, cfg, ke_cube)
        
        fig.colorbar(cf1, cax=fig.add_axes([0.85, 0.08, 0.02, 0.4]), label='m$^3$ s$^{-2}$',ticks=np.arange(0,26,5), shrink=0.6) # share colorbar
        # Save figure
        prov_record = get_provenance_record(
            f'Kinetic Energy {grp}.',
            files,
        )
        save_figure(
            f"kinetic_energy_{grp}",
            prov_record,
            cfg,
            figure=fig,
            dpi=300,
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
