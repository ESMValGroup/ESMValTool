""" Single model diagnostics
1. Solve the Poisson solver over the ocean
2. Produce and save plots
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from poisson_solver_oht import spherical_ocean_poisson
from scipy.ndimage import binary_fill_holes

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

# Initialise logger
logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(plot_type, ancestor_files):
    # Create a provenance record describing the diagnostic data and plot.

    record = {
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': [plot_type],
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def call_poisson(cube, mask):
    # Remove average of flux field to account for storage term
    data = cube.data.copy()
    grid = iris.analysis.cartography.area_weights(cube, normalize=True)
    data_mean = np.average(data, weights=np.ma.array(grid, mask=mask))
    data -= data_mean
    data *= (1 - mask)

    logger.info("Calling spherical_poisson")
    poisson, mht, dx, dy = spherical_ocean_poisson(logger,
                                                   forcing=data *
                                                   (6371e3**2.0),
                                                   mask=(1 - mask))
    logger.info("Ending spherical_poisson")

    # Poisson data cube
    p_cube = cube.copy()
    p_cube.rename('poisson')
    p_cube.units = 'W'
    p_cube.data = poisson[1:-1, 1:-1]

    # Flux data cube
    cube.rename('flux')
    cube.units = 'W.m-2'
    cube.data = data

    return p_cube, mht, cube, dx, dy


def masks(cfg):
    mask = iris.load_cube(cfg['mask'])

    plt.figure()
    plt.contourf(mask.data)
    plt.colorbar()
    plt.savefig('tst.png')

    # Close any seas not connected to the oceans
    # e.g. Mediterranean Sea at some resolutions
    (y, x) = np.shape(mask)
    wrap_mask = np.zeros([y, x + 2])
    wrap_mask[:, 1:-1] = mask.data
    wrap_mask[:, 0] = mask.data[:, -1]
    wrap_mask[:, -1] = mask.data[:, 0]

    wrap_mask = binary_fill_holes(wrap_mask)
    mask.data = wrap_mask[:, 1:-1]

    # Return wrapped mask
    return mask


def calc_nsf(dataset, mask):
    # Calculate NSF from the input fluxes
    # NSF = (rlds - rlus) + (rsds - rsus) - (hfls + hfss)

    nsf_dict = {}
    for attributes in dataset:
        input_file = attributes['filename']
        variable_name = attributes['short_name']
        logger.info("Loading %s", input_file)
        cube = iris.load_cube(input_file)
        nsf_dict[variable_name] = cube.data

    nsf = cube.copy()
    nsf.rename('Net Surface Flux')
    nsf.data = (nsf_dict['rlds'] - nsf_dict['rlus']) + \
        (nsf_dict['rsds'] - nsf_dict['rsus']) - \
        (nsf_dict['hfss'] + nsf_dict['hfls'])

    # Regrid NSF to the same grid as mask
    nsf = nsf.regrid(mask, iris.analysis.AreaWeighted())
    if 'variable_id' in nsf.attributes:
        del nsf.attributes['variable_id']

    return nsf


def compute_p_field_diagnostic(model_name, model_dict, cfg):
    # Define mask
    mask = masks(cfg)

    # Determine number of experiments
    exp_dataset = group_metadata(model_dict, 'exp', sort='variable_group')

    for exp_nme in exp_dataset:
        # Calculate MHT/solve Poisson equation
        logger.info("Experiment: %s", exp_nme)
        nsf = calc_nsf(exp_dataset[exp_nme], mask)
        p, mht, flux, dx, dy = call_poisson(nsf, mask.data)

        files_dict = {}
        for attributes in exp_dataset[exp_nme]:
            variable_name = attributes['short_name']
            input_file = attributes['filename']
            files_dict[variable_name] = input_file

        ancestor_files = files_dict
        basename = model_name + '_' + exp_nme

        # Plotting
        mht_plot(mht, flux, ancestor_files, cfg, pltname=basename)
        quiver_plot(p,
                    flux,
                    mask.data,
                    dx,
                    dy,
                    ancestor_files,
                    cfg,
                    pltname=basename)


def mht_plot(mht, cube, ancestor_files, cfg, pltname):
    """MHT plot Produces a single plot which compares the total ocean MHT to
    the Atlantic, Indian and Pacific ocean contributions."""

    plt.figure()
    for x in mht:
        plt.plot(cube.coord('latitude').points, mht[x] / 1e15, label=x)
    plt.xlim(-90, 90)
    ymin, ymax = plt.ylim()
    ylim = np.max([np.abs(ymin), np.abs(ymax)])
    plt.ylim(-ylim, ylim)
    plt.xticks(np.arange(-90, 120, 30))
    plt.hlines(0, -90, 90, color='k', linestyles=':')
    plt.vlines(0, -10, 10, color='k', linestyles=':')
    plt.xlabel('Latitude')
    plt.ylabel('MHT (PW)')
    plt.legend()
    plt.title(pltname)
    plt.tight_layout()

    # Save plot
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=ancestor_files)
    save_figure(pltname + '_mht', provenance_record, cfg)


def quiver_plot(poisson, flux, mask, dx, dy, ancestor_files, cfg, pltname):
    """Scalar p-field/NSF maps Produces two plots per dataset.

    1. The scalar p-field and the input flux field as subplots.
    2. Scalar p-field plotted individually with arrows showing
    the direction/magnitude of the heat transport.
    In all maps the land points are masked out and the land
    boundaries are plotted.
    """

    # Determine levels
    plt.figure()
    plt.contourf(
        poisson.coord('longitude').points,
        poisson.coord('latitude').points,
        np.ma.array(poisson.data / 1e15, mask=mask))
    cbar = plt.colorbar()
    cmin, cmax = cbar.mappable.get_clim()
    plt.close()

    levels1 = np.linspace(cmin, cmax, 21)

    # Poisson p-field + flux field
    plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(1, 2)

    ax1 = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
    cb1 = ax1.contourf(poisson.coord('longitude').points,
                       poisson.coord('latitude').points,
                       np.ma.array(poisson.data / 1e15, mask=mask),
                       levels=11,
                       transform=ccrs.PlateCarree(central_longitude=0))
    plt.gca().coastlines()
    ax1.set_xticks(np.arange(-180, 190, 60))
    ax1.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    ax1.set_yticks(np.arange(-90, 100, 30))
    ax1.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
    plt.colorbar(cb1, orientation='horizontal', label='Scalar p-field (PW)')

    ax2 = plt.subplot(gs[0, 1], projection=ccrs.PlateCarree())
    cb2 = ax2.contourf(flux.coord('longitude').points,
                       flux.coord('latitude').points,
                       np.ma.array(flux.data, mask=mask),
                       levels=10,
                       transform=ccrs.PlateCarree(central_longitude=0),
                       cmap='plasma')
    plt.gca().coastlines()
    ax2.set_xticks(np.arange(-180, 190, 60))
    ax2.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    ax2.set_yticks(np.arange(-90, 100, 30))
    ax2.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
    plt.colorbar(cb2, orientation='horizontal', label='Flux (W/m2)')

    plt.suptitle(pltname, y=0.9)
    plt.tight_layout()

    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=ancestor_files)
    save_figure(pltname + '_poisson_flux', provenance_record, cfg)

    # Poisson p-field + quiver
    plt.figure(figsize=(10, 6))
    ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())

    x, y = np.meshgrid(
        poisson.coord('longitude').points,
        poisson.coord('latitude').points)
    v, u = dy / 1e15, dx / 1e15

    cb1 = ax1.contourf(poisson.coord('longitude').points,
                       poisson.coord('latitude').points,
                       np.ma.array(poisson.data / 1e15, mask=mask),
                       levels=levels1,
                       transform=ccrs.PlateCarree(central_longitude=0))
    plt.gca().coastlines()
    ax1.quiver(np.ma.array(x, mask=mask)[9::18, 6::12],
               np.ma.array(y, mask=mask)[9::18, 6::12],
               np.ma.array(u, mask=mask)[9::18, 6::12],
               np.ma.array(v, mask=mask)[9::18, 6::12],
               pivot='mid',
               color='k')

    ax1.set_xticks(np.arange(-180, 190, 60))
    ax1.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    ax1.set_yticks(np.arange(-90, 100, 30))
    ax1.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])

    plt.colorbar(cb1, orientation='horizontal', label='Scalar p-field (PW)')
    plt.title(pltname)

    save_figure(pltname + '_poisson', provenance_record, cfg)


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Group data by dataset
    logger.info("Group input data by model, sort by experiment")
    model_dataset = group_metadata(input_data, 'dataset', sort='exp')

    # Solve Poisson equation and produce plots
    for model_name in model_dataset:
        logger.info("Processing model data: %s", model_name)
        compute_p_field_diagnostic(model_name, model_dataset[model_name], cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
