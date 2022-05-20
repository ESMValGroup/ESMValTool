""" Single model diagnostics
1. Solve the Poisson solver
2. Produce and save plots
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from poisson_solver import spherical_poisson

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


def call_poisson(cube, variable):
    # Remove average of flux field to account for storage term
    data = cube.data.copy()
    grid = iris.analysis.cartography.area_weights(cube, normalize=True)
    data_mean = np.average(data, weights=grid.data)
    data -= data_mean

    logger.info("Calling spherical_poisson")
    poisson, mht = spherical_poisson(logger,
                                     forcing=data * (6371e3**2.0),
                                     tolerance=2.0e-4)
    logger.info("Ending spherical_poisson")

    # Poisson data cube
    p_cube = cube.copy()
    p_cube.rename(variable + '_poisson')
    p_cube.units = 'W'
    p_cube.data = poisson[1:-1, 1:-1]

    # MHT data cube
    mht_cube = cube.copy()
    mht_cube = mht_cube.collapsed('longitude', iris.analysis.MEAN)
    mht_cube.rename(variable + '_mht')
    mht_cube.units = 'W'
    mht_cube.data = mht

    # Flux data cube
    cube.rename(variable + '_flux')
    cube.units = 'W.m-2'
    cube.data = data

    return p_cube, mht_cube, cube


def compute_p_field_diagnostic(model_dict, cfg):
    # Sort experiments by variable
    exp_dataset = group_metadata(model_dict, 'exp', sort='variable_group')

    # Define dicts to hold the data
    p_dict = {}
    mht_dict = {}
    flux_dict = {}
    files_dict = {}

    for exp_nme in exp_dataset:
        logger.info("Experiment: %s", exp_nme)
        p_cube = iris.cube.CubeList([])
        mht_cube = iris.cube.CubeList([])
        flux_cube = iris.cube.CubeList([])

        for attributes in exp_dataset[exp_nme]:
            input_file = attributes['filename']
            variable_name = attributes['short_name']
            logger.info("Loading %s", input_file)
            cube = iris.load_cube(input_file)

            # Model key for dicts
            key = '{0}_{1}'.format(attributes['dataset'], attributes['exp'])
            if key not in files_dict:
                files_dict[key] = {}
            files_dict[key][variable_name] = input_file

            # If flux is upwelling, x -1
            if attributes['short_name'] in [
                    'rlut', 'rlutcs', 'rsut', 'rsutcs'
            ]:
                cube *= -1

            # Call the Poisson solver
            p, mht, flux = call_poisson(cube, variable_name)
            p_cube.append(p)
            mht_cube.append(mht)
            flux_cube.append(cube)

            del p, mht

        # Add iris cubelist to dictionary
        p_dict[key] = p_cube
        mht_dict[key] = mht_cube
        flux_dict[key] = flux_cube
        del key, p_cube, mht_cube, flux_cube

    return p_dict, mht_dict, flux_dict, files_dict


def mht_plot(variables, mht, ancestor_files, cfg, pltname):
    """MHT plot Produces a single plot comparing the estimated MHT due to the
    input variables.

    MHT is presented in PW, plotted against latitude. Up to three
    variables are on each plot.
    """

    # Variables
    var_dict = {
        'rtnt': 'Net TOA',
        'rsnt': 'Net SW',
        'rlut': 'Net LW',
        'rlutcs': '-1 x OLR (clear)',
        'rsut': '-1 x OSR (all)',
        'rsutcs': '-1 x OSR (clear)',
        'netcre': 'NET CRE',
        'swcre': 'SW CRE',
        'lwcre': 'LW CRE'
    }

    # MHT comparison plot
    plt.figure()
    for i, x in enumerate(variables):
        cube = mht.extract_cube(x + '_mht')
        plt.plot(cube.coord('latitude').points,
                 cube.data / 1e15,
                 label=var_dict[x])
        del cube
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
    plt.tight_layout()

    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=ancestor_files)
    save_figure(pltname, provenance_record, cfg)


def quiver_plot(variables, poisson, flux, ancestor_files, cfg, pltname):
    """Maps of scalar p-fields and corresponding flux The maps are plotted such
    that the variables share the same colorbar (one each for the flux maps and
    the scalar p-field maps).

    Scalar p-field is plotted as PW, flux maps are W/m2. The scalar
    p-field maps have arrows imposed over the map which show the
    direction of the implied heat transport. The arrows are also
    normalised across the different variables so the lengths can be
    directly compared (a longer arrow corresponds to a larger heat
    transport).
    """

    # Variables
    var_dict = {
        'rtnt': 'Net TOA',
        'rsnt': 'Net SW',
        'rlut': 'Net LW',
        'rlutcs': '-1 x OLR (clear)',
        'rsut': '-1 x OSR (all)',
        'rsutcs': '-1 x OSR (clear)',
        'netcre': 'NET CRE',
        'swcre': 'SW CRE',
        'lwcre': 'LW CRE'
    }

    # Variables to determine map extent
    # Allows for same colorbar across subplots
    vmin = 0
    vmax = 0
    wmin = 0
    wmax = 0

    for i, p in enumerate(variables):
        # Scalar p-field
        cube = poisson.extract_cube(p + '_poisson')
        tmp = np.average(cube.data)

        plt.figure()
        plt.axes(projection=ccrs.PlateCarree())
        plt.contourf(cube.coord('longitude').points,
                     cube.coord('latitude').points,
                     cube.data - tmp,
                     levels=10,
                     transform=ccrs.PlateCarree(central_longitude=0))
        plt.gca().coastlines()
        cbar = plt.colorbar()
        cmin, cmax = cbar.mappable.get_clim()
        vmin = np.min([vmin, cmin])
        vmax = np.max([vmax, cmax])
        del tmp, cube
        plt.close()

        # Flux
        cube = flux.extract_cube(p + '_flux')
        tmp = np.average(cube.data)

        plt.figure()
        plt.axes(projection=ccrs.PlateCarree())
        plt.contourf(cube.coord('longitude').points,
                     cube.coord('latitude').points,
                     cube.data - tmp,
                     levels=10,
                     transform=ccrs.PlateCarree(central_longitude=0))
        plt.gca().coastlines()
        cbar = plt.colorbar()
        cmin, cmax = cbar.mappable.get_clim()
        wmin = np.min([wmin, cmin])
        wmax = np.max([wmax, cmax])
        del tmp
        plt.close()

    levels1 = np.linspace(vmin / 1e15, vmax / 1e15, 11)
    levels2 = np.linspace(wmin, wmax, 11)

    len_p = len(variables)
    if len_p == 3:
        plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(22, 2)
        gs.update(wspace=0.25, hspace=1.5)
    elif len_p == 2:
        plt.figure(figsize=(10, 6.5))
        gs = gridspec.GridSpec(15, 2)
        gs.update(wspace=0.25, hspace=1.5)

    for i, p in enumerate(variables):
        # Scalar p-field
        cube = poisson.extract_cube(p + '_poisson')
        tmp = np.average(cube.data)

        x, y = np.meshgrid(
            cube.coord('longitude').points,
            cube.coord('latitude').points)
        v, u = np.gradient(cube.data, 1e14, 1e14)
        u = u[1:-1, 1:-1]
        v = v[1:-1, 1:-1]

        ax1 = plt.subplot(gs[i * 7:(i * 7) + 7, 0],
                          projection=ccrs.PlateCarree())
        cb1 = ax1.contourf(cube.coord('longitude').points,
                           cube.coord('latitude').points,
                           (cube.data - tmp) / 1e15,
                           levels=levels1,
                           transform=ccrs.PlateCarree(central_longitude=0))
        plt.gca().coastlines()
        if i == 0:
            R = ax1.quiver(x[15::25, 15::32],
                           y[15::25, 15::32],
                           u[15::25, 15::32],
                           v[15::25, 15::32],
                           pivot='mid',
                           color='w')
            R._init()
        else:
            ax1.quiver(x[15::25, 15::32],
                       y[15::25, 15::32],
                       u[15::25, 15::32],
                       v[15::25, 15::32],
                       scale=R.scale,
                       pivot='mid',
                       color='w')
        ax1.set_xticks(np.arange(-180, 190, 60))
        ax1.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
        ax1.set_yticks(np.arange(-90, 100, 30))
        ax1.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
        ax1.set_title(var_dict[p])
        del cube, tmp, u, v, x, y

        # Flux
        cube = flux.extract_cube(p + '_flux')
        tmp = np.average(cube.data)
        ax2 = plt.subplot(gs[i * 7:(i * 7) + 7, 1],
                          projection=ccrs.PlateCarree())
        cb2 = ax2.contourf(cube.coord('longitude').points,
                           cube.coord('latitude').points,
                           cube.data - tmp,
                           levels=levels2,
                           transform=ccrs.PlateCarree(central_longitude=0),
                           cmap='plasma')
        plt.gca().coastlines()
        ax2.set_xticks(np.arange(-180, 190, 60))
        ax2.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
        ax2.set_yticks(np.arange(-90, 100, 30))
        ax2.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
        ax2.set_title(var_dict[p])
        del tmp

    ax1 = plt.subplot(gs[-1, 0])
    plt.colorbar(cb1,
                 cax=ax1,
                 orientation='horizontal',
                 label='Scalar p-field (PW)')

    ax2 = plt.subplot(gs[-1, 1])
    plt.colorbar(cb2, cax=ax2, orientation='horizontal', label='Flux (W/m2)')

    if len_p == 3:
        plt.subplots_adjust(left=0.1, right=0.94, top=1.0, bottom=0.11)
    elif len_p == 2:
        plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.13)

    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=ancestor_files)
    save_figure(pltname, provenance_record, cfg)


def plot_single_model_diagnostics(poisson, mht, flux, basename, ancestor_files,
                                  cfg):

    # TOA: net, lw, sw
    variables = ['rtnt', 'rsnt', 'rlut']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_toa_mht')
    quiver_plot(variables,
                poisson,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_toa_p_field')

    # OSR: clear vs. all sky
    variables = ['rsut', 'rsutcs']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_osr_mht')
    quiver_plot(variables,
                poisson,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_osr_p_field')

    # CRE: net, lw, sw
    variables = ['netcre', 'swcre', 'lwcre']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_cre_mht')
    quiver_plot(variables,
                poisson,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_cre_p_field')


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""
    input_data = cfg['input_data'].values()

    # Dictionaries for data cubes
    model_p = {}
    model_mht = {}
    model_flux = {}
    model_files = {}

    # Group data by dataset
    logger.info("Group input data by model, sort by experiment")
    model_dataset = group_metadata(input_data, 'dataset', sort='exp')

    # Solve Poisson equation for each dataset
    for model_name in model_dataset:
        if model_name == 'CERES-EBAF':
            # Ignore any observational data
            continue

        logger.info("Processing model data: %s", model_name)
        p, mht, flux, files = compute_p_field_diagnostic(
            model_dataset[model_name], cfg)
        model_p.update(p)
        model_mht.update(mht)
        model_flux.update(flux)
        model_files.update(files)
        del p, mht, flux, files

    # Produce plots
    for model in model_mht:
        output_basename = model
        plot_single_model_diagnostics(model_p[model], model_mht[model],
                                      model_flux[model], output_basename,
                                      model_files[model], cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
