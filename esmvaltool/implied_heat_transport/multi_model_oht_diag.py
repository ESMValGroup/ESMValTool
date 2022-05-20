""" Multi model diagnostics
1. Solve the Poisson equation for each dataset
2. Calculate/process observational dataset
3. Compare simulations to observational data
4. Produce comparison plots
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib import gridspec
from poisson_solver_oht import spherical_ocean_poisson
from scipy.ndimage import binary_fill_holes

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

# Initialise logger
logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(plot_type, ancestor_files):
    # Create a provenance record describing the diagnostic data and plot

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

    # MHT data cube list
    mht_list = iris.cube.CubeList([])
    for x in mht:
        basin = cube.copy()
        basin = basin.collapsed('longitude', iris.analysis.MEAN)
        basin.rename(x)
        basin.units = 'W'
        basin.data = mht[x]
        mht_list.append(basin)
        del basin

    # Gradient cube list
    dx_cube = cube.copy()
    dx_cube.rename('dx gradient')
    dx_cube.units = 'W'
    dx_cube.data = dx

    dy_cube = cube.copy()
    dy_cube.rename('dy gradient')
    dy_cube.units = 'W'
    dy_cube.data = dy

    grad_list = iris.cube.CubeList([dx_cube, dy_cube])
    del dx_cube, dy_cube

    # Flux data cube
    cube.rename('flux')
    cube.units = 'W.m-2'
    cube.data = data

    return p_cube, mht_list, cube, grad_list


def masks(cfg):
    mask = iris.load_cube(cfg['mask'])

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


def calc_nsf(dataset, mask, cfg):
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

    # Define data dicts
    p_dict = {}
    mht_dict = {}
    flux_dict = {}
    grad_dict = {}
    files_dict = {}

    # Determine number of experiments
    exp_dataset = group_metadata(model_dict, 'exp', sort='variable_group')

    for exp_nme in exp_dataset:
        logger.info("Experiment: %s", exp_nme)
        key = model_name + '_' + exp_nme

        # Calculate NSF/solve Poisson equation
        nsf = calc_nsf(exp_dataset[exp_nme], mask, cfg)
        p, mht, flux, gradient = call_poisson(nsf, mask.data)

        # Save data to dict
        p_dict[key] = p
        mht_dict[key] = mht
        flux_dict[key] = flux
        grad_dict[key] = gradient

        files_dict = {}
        for attributes in exp_dataset[exp_nme]:
            variable_name = attributes['short_name']
            input_file = attributes['filename']
            files_dict[variable_name] = input_file
            filename = attributes['dataset'] + '_' + \
                attributes['mip'] + '_' + attributes['exp'] + '_' + \
                attributes['ensemble'] + '_nsf_' + \
                str(attributes['start_year']) + '-' + \
                str(attributes['end_year'])

        # Save data cubes
        provenance_record = get_provenance_record(plot_type='map',
                                                  ancestor_files=exp_dataset)
        save_data(filename + '_poisson', provenance_record, cfg, p)
        save_data(filename + '_mht', provenance_record, cfg, mht)

    return p_dict, mht_dict, flux_dict, grad_dict, files_dict, mask


def plot_multi_model_diagnostics(poisson, mht, flux, gradient, mask, files,
                                 cfg):
    # Plot comparison of models to observational dataset
    mht_plot(mht, files, cfg)
    quiver_plot(poisson, gradient, mask.data, files, cfg)
    quiver_plot(flux, gradient, mask.data, files, cfg)
    return


def mht_plot(mht, ancestor_files, cfg):
    """MHT plot Produces two plots comparing the estimated MHT of the different
    models to the observational dataset provided.

    1. '_obs_ocean_mht.png'
        Gives the MHT for the given obs dataset, then shows the difference
        between the model and obs for each basin in separate subplots. All the
        difference subplots are presented with the same y-limits to allow for
        easy comparison.
    2. '_ocean_mht.png'
        Four subplots showing the MHT for each model in each basin. The
        observational results are presented as a thick black solid line. Allows
        for easy comparison between the shapes of the MHT curves, but doesn't
        give the same detail as the difference plots.
    For both plots each model is presented in a different colour/linestyle
    configuration.
    """

    # MHT plot: direct comparison to observations
    # Color/linestyle cycler to ensure no models are the same
    cc = (cycler(linestyle=['-', '--', '-.']) * cycler('color', [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]))

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    for i, x in enumerate(['global', 'atlantic', 'pacific', 'indian']):
        ax = plt.subplot(2, 2, i + 1)

        for y in mht:
            cube = mht[y].extract_cube(x)
            if y == cfg['obs_name']:
                # Plot observational data separately (black line)
                ax.plot(cube.coord('latitude').points,
                        cube.data / 1e15,
                        label=y,
                        color='k',
                        lw=2,
                        zorder=100)
                ax.set_prop_cycle(cc)
            else:
                ax.plot(cube.coord('latitude').points,
                        cube.data / 1e15,
                        label=y)
            plt.xlim(-90, 90)
            plt.xticks(np.arange(-90, 120, 30))
            plt.hlines(0, -90, 90, color='k', linestyles=':')
            if i == 0:
                ymin, ymax = plt.ylim()
                ylim = np.max([np.abs(ymin), np.abs(ymax)])
            plt.ylim(-ylim, ylim)
            plt.vlines(0, -ylim, ylim, color='k', linestyles=':')
            plt.title(x)
            plt.xlabel('Latitude')
            plt.ylabel('MHT (PW)')

    # Plot legend to right of figure outside subplots
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
    plt.tight_layout()

    # Save figure
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=ancestor_files)
    save_figure('multi_model_ocean_mht', provenance_record, cfg)

    # MHT plot: difference compared to observations
    # Color/linestyle cycler to ensure no models are the same
    cc = (cycler(linestyle=['-', '--', '-.']) * cycler(
        'color',
        ['#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']))

    if cfg['obs_name'] in mht:
        obs_name = cfg['obs_name']
        vlim = 0

        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        ax1 = plt.subplot(2, 3, 1)  # Obs
        # Subplot for each basin
        for i, x in enumerate(['global', 'atlantic', 'pacific', 'indian']):
            if i > 1:
                ax = plt.subplot(2, 3, i + 3)
            else:
                ax = plt.subplot(2, 3, i + 2)
            ax.set_prop_cycle(cc)
            obs = mht[obs_name].extract_cube(x)
            for y in mht:
                if y == obs_name:
                    # Obs plotted in separate subplot
                    ax1.plot(cube.coord('latitude').points,
                             obs.data / 1e15,
                             label=x)
                    ax1.set_xlim(-90, 90)
                    ax1.set_xticks(np.arange(-90, 120, 30))
                    ax1.axhline(0, -90, 90, color='k', linestyle=':')

                    if i == 3:
                        dmin, dmax = ax1.get_ylim()
                        dlim = np.max([np.abs(dmin), np.abs(dmax)])
                        ax1.set_ylim(-dlim, dlim)
                        ax1.axvline(0, -dlim, dlim, color='k', linestyle=':')
                        del dmin, dmax, dlim

                    ax1.set_xlabel('Latitude')
                    ax1.set_ylabel('MHT (PW)')
                    ax1.set_title(y)
                    ax1.legend()
                    continue

                cube = mht[y].extract_cube(x)
                ax.plot(cube.coord('latitude').points,
                        (cube.data - obs.data) / 1e15,
                        label=y)
                plt.xlim(-90, 90)
                plt.xticks(np.arange(-90, 120, 30))
                plt.hlines(0, -90, 90, color='k', linestyles=':')
                if i == 0:
                    ymin, ymax = plt.ylim()
                    ylim = np.max([np.abs(ymin), np.abs(ymax)])
                    vlim = np.max([vlim, ylim])
                # Get model names for legend
                handles, labels = ax.get_legend_handles_labels()

            # Keep y-limits the same across all difference plots
            plt.ylim(-vlim, vlim)
            plt.vlines(0, -vlim, vlim, color='k', linestyles=':')
            plt.title(x)
            plt.xlabel('Latitude')
            plt.ylabel('MHT (PW)')

        # Plot legend of models under the observational plot
        ax2 = plt.subplot(2, 3, 4)
        plt.axis('off')
        ax2.legend(handles, labels, loc='center left')
        plt.tight_layout()

        # Save plot
        provenance_record = get_provenance_record(
            plot_type='zonal', ancestor_files=ancestor_files)
        save_figure('multi_model_' + obs_name + '_ocean_mht',
                    provenance_record, cfg)


def quiver_plot(data, gradient, mask, ancestor_files, cfg):
    """Maps of scalar p-fields/input flux fields Produces four plots of maps
    comparing the observational data to models.

    Creates the
    same two plots for both the NSF maps and the scalar p-field.
    1. 'multi_model_ocean_flux/poisson.png'
        Plots the full maps for the observational data followed by all of the
        models. All of the subplots share the same colorbar allowing for direct
        comparison of the maps.
    2. 'multi_model_obs_ocean_flux.png'
        Shwows the full map for the observational data with its own colorbar,
        followed by all of the difference between the model maps and the
        observational map. The model difference subplots all share the same
        colorbar allowing for better comparison of the difference between
        models.
    """

    # Plot 1.
    # Determine colorbar levels
    vmin, vmax = 0, 0
    for x in data:
        # Different factor depending on p-field or NSF
        cube = data[x]
        if cube.name() == 'poisson':
            div = 1e15
        else:
            div = 1

        plt.figure()
        plt.contourf(cube.data, levels=10)
        cbar = plt.colorbar()
        cmin, cmax = cbar.mappable.get_clim()
        vmin = np.min([vmin, cmin])
        vmax = np.max([vmax, cmax])
        plt.close()

        del cube

    levels1 = np.linspace(vmin / div, vmax / div, 21)
    del vmin, vmax

    # Size of plot parameters
    # Change these if plots are not displaying properly
    if len(data) <= 4:
        ncols = 2
        nrows = int(np.ceil(len(data) / 2))
    elif len(data) <= 6:
        ncols = 3
        nrows = 2
    else:
        ncols = 4
        nrows = int(np.ceil(len(data) / 4))

    # Create figure axes
    height_ratio = np.ones(nrows + 1)
    height_ratio[:nrows] *= 7

    plt.figure(figsize=(15, 4 * nrows))
    gs = gridspec.GridSpec(nrows + 1, ncols, height_ratios=height_ratio)

    j = 0
    cube = 0
    # Plot each model individually
    for i, model in enumerate(data):
        # Iterate over subplots
        if i >= ncols:
            if i % ncols == 0:
                j += 1
            k = i - (j * ncols)
        else:
            k = i

        del cube
        cube = data[model]

        ax = plt.subplot(gs[j:j + 1, k:k + 1], projection=ccrs.PlateCarree())
        cb = ax.contourf(
            cube.coord('longitude').points,
            cube.coord('latitude').points,
            np.ma.array(cube.data / div, mask=mask),
            levels=levels1,
            transform=ccrs.PlateCarree(central_longitude=0),
            cmap='viridis' if cube.name() == 'poisson' else 'plasma')
        plt.gca().coastlines()

        ax.set_xticks(np.arange(-180, 190, 60))
        ax.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
        ax.set_yticks(np.arange(-90, 100, 30))
        ax.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
        ax.set_title(model)

    # Plot colorbar at bottom spanning length of plot
    ax = plt.subplot(gs[-1, :])
    plt.colorbar(cb,
                 cax=ax,
                 orientation='horizontal',
                 label='Scalar p-field (PW)'
                 if cube.name() == 'poisson' else 'Flux (w/m2)')

    # Adjust and save plot
    plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.1)
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=ancestor_files)
    save_figure('mulit_model_ocean_' + cube.name(), provenance_record, cfg)

    # Plot 2.
    # Difference between models and obs
    if cfg['obs_name'] in data:
        obs_name = cfg['obs_name']

        # Determine colorbar levels
        vmin, vmax = 0, 0
        for x in data:
            cube = data[x]
            obs = data[obs_name].data - np.average(data[obs_name].data)
            # Define different factors for p-field/NSF maps
            if cube.name() == 'poisson':
                div = 1e15
            else:
                div = 1
            cube = cube.data - np.average(cube.data)

            plt.figure()
            plt.contourf(cube - obs, levels=10)
            cbar = plt.colorbar()
            cmin, cmax = cbar.mappable.get_clim()
            vmin = np.min([vmin, cmin])
            vmax = np.max([vmax, cmax])
            plt.close()

            del cube

        vlim = np.max([np.abs(vmin), np.abs(vmax)])
        levels1 = np.linspace(-vlim / div, vlim / div, 21)
        del vmin, vmax

        # Size of plot parameters
        # Change these if plots are not displaying properly
        if len(data) <= 5:
            ncols = 2
            nrows = int(np.ceil((len(data) - 1) / 2))
            xlen = 15
        elif len(data) <= 7:
            ncols = 3
            nrows = 2
            xlen = 15
        else:
            ncols = 4
            nrows = int(np.ceil((len(data) - 1) / 4))
            xlen = 20

        # Create figure axes
        height_ratio = np.ones(nrows + 1)
        height_ratio[:nrows] *= 7

        plt.figure(figsize=(xlen, 4 * nrows))
        gs = gridspec.GridSpec(nrows + 1,
                               ncols + 1,
                               height_ratios=height_ratio)

        j = 0
        for i, model in enumerate(data):
            if model == obs_name:
                # Plot observational map separately
                cube = data[model]
                ax = plt.subplot(gs[nrows - 1:nrows, 0],
                                 projection=ccrs.PlateCarree())
                cb = ax.contourf(
                    cube.coord('longitude').points,
                    cube.coord('latitude').points,
                    np.ma.array(cube.data / div, mask=mask),
                    levels=21,
                    transform=ccrs.PlateCarree(central_longitude=0),
                    cmap='viridis' if cube.name() == 'poisson' else 'plasma')
                plt.gca().coastlines()

                ax.set_xticks(np.arange(-180, 190, 60))
                ax.set_xticklabels(
                    ['180', '120W', '60W', '0', '60E', '120E', '180'])
                ax.set_yticks(np.arange(-90, 100, 30))
                ax.set_yticklabels(
                    ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
                ax.set_title(model)

                ax = plt.subplot(gs[-1, :1])
                plt.colorbar(cb,
                             cax=ax,
                             orientation='horizontal',
                             label='Scalar p-field (PW)'
                             if cube.name() == 'poisson' else 'Flux (w/m2)')
                continue

            if (i - 1) >= ncols:
                if (i - 1) % ncols == 0:
                    j += 1
                k = i - (j * ncols)
            else:
                k = i

            del cube
            # Plot difference between model and obs
            cube = data[model]
            tmp_cube = cube.data - np.average(cube.data)

            ax = plt.subplot(gs[j:j + 1, k:k + 1],
                             projection=ccrs.PlateCarree())
            cb = ax.contourf(cube.coord('longitude').points,
                             cube.coord('latitude').points,
                             np.ma.array((tmp_cube - obs) / div, mask=mask),
                             levels=levels1,
                             transform=ccrs.PlateCarree(central_longitude=0),
                             cmap='coolwarm')
            plt.gca().coastlines()

            ax.set_xticks(np.arange(-180, 190, 60))
            ax.set_xticklabels(
                ['180', '120W', '60W', '0', '60E', '120E', '180'])
            ax.set_yticks(np.arange(-90, 100, 30))
            ax.set_yticklabels(
                ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
            ax.set_title(model + ' - ' + obs_name)

        # Add subplot for colorbar shared between difference plots
        ax = plt.subplot(gs[-1, 1:])
        plt.colorbar(cb,
                     cax=ax,
                     orientation='horizontal',
                     label='Scalar p-field (PW)'
                     if cube.name() == 'poisson' else 'Flux (w/m2)')

        # Save plot
        plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.1)
        provenance_record = get_provenance_record(
            plot_type='map', ancestor_files=ancestor_files)
        save_figure('mulit_model_' + obs_name + '_ocean_' + cube.name(),
                    provenance_record, cfg)


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Dictionaries for data cubes
    model_p = {}
    model_mht = {}
    model_flux = {}
    model_grad = {}
    model_files = {}

    # Observational comparison data
    if cfg['obs_data'] is not False:
        mask = masks(cfg)

        # Regrid to mask shape
        obs_cube = iris.load_cube(cfg['obs_data'])
        obs_cube = obs_cube.regrid(mask, iris.analysis.AreaWeighted())
        obs_name = cfg['obs_name']

        # Solve Poisson equation
        model_p[obs_name], model_mht[obs_name], \
            model_flux[obs_name], model_grad[obs_name] = \
            call_poisson(obs_cube, mask.data)

        # Save data cubes
        provenance_record = get_provenance_record(
            plot_type='map', ancestor_files=cfg['obs_data'])
        save_data(obs_name + '_poisson', provenance_record, cfg,
                  model_p[obs_name])
        save_data(obs_name + '_mht', provenance_record, cfg,
                  model_mht[obs_name])
        save_data(obs_name + '_flux', provenance_record, cfg,
                  model_flux[obs_name])
        del provenance_record

        model_files[obs_name] = cfg['obs_name']

    # Group data by dataset
    logger.info("Group input data by model, sort by experiment")
    model_dataset = group_metadata(input_data, 'dataset', sort='exp')

    for model_name in model_dataset:
        logger.info("Processing model data: %s", model_name)
        # Solve Poisson equation
        p, mht, flux, grad, files, mask = compute_p_field_diagnostic(
            model_name, model_dataset[model_name], cfg)
        model_p.update(p)
        model_mht.update(mht)
        model_flux.update(flux)
        model_grad.update(grad)
        model_files.update(files)

    # Multi-model plots
    plot_multi_model_diagnostics(model_p, model_mht, model_flux, model_grad,
                                 mask, model_files, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
