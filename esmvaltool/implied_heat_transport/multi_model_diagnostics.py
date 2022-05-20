""" Multi model diagnostics
1. Solve the Poisson equation for each dataset
2. Compare to CERES-EBAF results
3. Produce comparison plots
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib import gridspec
from poisson_solver import spherical_poisson

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
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


def call_poisson(cube, variable, tol=None):
    logger.info('{0}'.format(cube))

    # Remove average of flux field to account for storage term
    data = cube.data.copy()
    grid = iris.analysis.cartography.area_weights(cube, normalize=True)
    data_mean = np.average(data, weights=grid.data)
    data -= data_mean

    logger.info("Calling spherical_poisson")
    poisson, mht = spherical_poisson(logger,
                                     forcing=data * (6371e3**2.0),
                                     tolerance=tol)
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
    # Sort experiments by variables
    exp_dataset = group_metadata(model_dict, 'exp', sort='variable_group')

    # Define dicts to hold the data
    p_dict = {}
    mht_dict = {}
    flux_dict = {}

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

            # Save data cubes
            provenance_record = get_provenance_record(
                plot_type='map', ancestor_files=exp_dataset)
            filename = attributes['alias'] + '_' + attributes['dataset'] + \
                '_' + attributes['mip'] + '_' + attributes['exp'] + '_' + \
                attributes['ensemble'] + '_' + attributes['variable_group'] \
                + '_' + str(attributes['start_year']) + '-' \
                + str(attributes['end_year'])
            save_data(filename + '_poisson', provenance_record, cfg, p)
            save_data(filename + '_mht', provenance_record, cfg, mht)

            del p, mht, flux, provenance_record, filename

        # Add iris cubelist to dictionary
        p_dict[key] = p_cube
        mht_dict[key] = mht_cube
        flux_dict[key] = flux_cube
        del key, p_cube, mht_cube, flux_cube

    return p_dict, mht_dict, flux_dict


def ceres_p_field_diagnostic(model_dict, cfg):
    # Calculate diagnostics for CERES-EBAF data

    # Load iris cubes: rlut, rlutcs, rsut, rsutcs, rsdt
    for attributes in model_dict:
        input_file = attributes['filename']
        logger.info('{0}'.format(input_file))
        cube = iris.load_cube(input_file)

        if attributes['short_name'] == 'rlut':
            rlut = cube.copy()
            rlut *= -1
            rlut.long_name = attributes['short_name']
        elif attributes['short_name'] == 'rsut':
            rsut = cube.copy()
            rsut *= -1
            rsut.long_name = attributes['short_name']
        elif attributes['short_name'] == 'rlutcs':
            rlutcs = cube.copy()
            rlutcs *= -1
            rlutcs.long_name = attributes['short_name']
        elif attributes['short_name'] == 'rsutcs':
            rsutcs = cube.copy()
            rsutcs *= -1
            rsutcs.long_name = attributes['short_name']
        elif attributes['short_name'] == 'rsdt':
            rsdt = cube.copy()
            rsdt.long_name = attributes['short_name']
        del cube

    # Define additional variables
    rsnt = rsdt + rsut
    rsnt.long_name = 'rsnt'
    rtnt = rsnt + rlut
    rtnt.long_name = 'rtnt'
    swcre = rsut - rsutcs
    swcre.long_name = 'swcre'
    lwcre = rlut - rlutcs
    lwcre.long_name = 'lwcre'
    netcre = swcre + lwcre
    netcre.long_name = 'netcre'

    p_dict = {}
    mht_dict = {}
    flux_dict = {}

    cube_list = [rtnt, rsnt, rlut, rsut, rsutcs, netcre, swcre, lwcre]

    p_cube = iris.cube.CubeList([])
    mht_cube = iris.cube.CubeList([])
    flux_cube = iris.cube.CubeList([])

    for cube in cube_list:
        # Solve Poisson equation for each variable
        p, mht, flux = call_poisson(cube, variable=cube.long_name, tol=2e-4)
        p_cube.append(p)
        mht_cube.append(mht)
        flux_cube.append(cube)
        del p, mht, cube

    # Add data to dictionaries
    key = 'CERES-EBAF'
    p_dict[key] = p_cube
    mht_dict[key] = mht_cube
    flux_dict[key] = flux_cube

    logger.info("Finished processing observational data.")
    return p_dict, mht_dict, flux_dict


def mht_plot(variables, mht, ancestor_files, cfg, pltname):
    """MHT plot Produces two plots comparing the estimated MHT of the different
    models against CERES-EBAF.

    1. '_mht_CERES-EBAF.png'
        The first subplot gives the MHT of CERES-EBAF for all of the given
        variables, then the difference between the model and CERES result in
        each of the other subplots. Each variable has its own subplot for the
        differences. The axes are not forced to the same limits as the
        differences can vary between the variables.
    2. '_mht_comp.png'
        Shows the full MHT estimates for all of the models and CERES-EBAF
        (shown in black). Allows for easy comparison between the shapes of the
        MHT curves, but doesn't give the same detail as the difference plots.
        Each variable has its own subplot.
    """

    # Variables
    var_dict = {
        'rtnt': 'Net TOA',
        'rsnt': 'Net SW',
        'rlut': 'Net LW',
        'rlutcs': '-1 x OLR (clear)',
        'rsut': '-1 x OSR (all)',
        'rsutcs': '-1 x OSR (clear)',
        'rsdt': 'ISR',
        'netcre': 'NET CRE',
        'swcre': 'SW CRE',
        'lwcre': 'LW CRE'
    }
    ''' MHT comparison plot: compares difference with CERES-EBAF '''
    # Different sized plots depending on number of variables
    plt.figure(figsize=(15, 6))
    if len(variables) == 2:
        gs = gridspec.GridSpec(2, 3, height_ratios=[7, 1])
    elif len(variables) == 3:
        gs = gridspec.GridSpec(2, 4, height_ratios=[7, 1])

    # Plot CERES-EBAF results separately
    ax = plt.subplot(gs[:1, :1])
    for i, x in enumerate(variables):
        cube = mht['CERES-EBAF'].extract_cube(x + '_mht')
        ax.plot(cube.coord('latitude').points,
                cube.data / 1e15,
                label=var_dict[x])
        del cube

    ymin, ymax = ax.get_ylim()
    ylim = np.max([np.abs(ymin), np.abs(ymax)])
    ax.set_ylim(-ylim, ylim)
    ax.set_xlim(-90, 90)
    ax.set_xticks(np.arange(-90, 120, 30))
    ax.axhline(0, -90, 90, color='k', ls=':')
    ax.axvline(0, -ylim, ylim, color='k', ls=':')
    ax.set_xlabel('Latitude')
    ax.set_ylabel('MHT (PW)')
    ax.set_title('CERES-EBAF')
    ax.legend()
    vlim = ylim
    del ymin, ymax, ylim

    # Bespoke color cycler to remove those used for the CERES-EBAF
    # Cycles through linestyles for when there are many comparison models
    cc = (cycler(linestyle=['-', '--', '-.']) * cycler('color', [
        '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
        '#17becf'
    ]))

    # Plot difference for all variables/models
    for i, x in enumerate(variables):
        ax = plt.subplot(gs[:1, i + 1:i + 2])
        ax.set_prop_cycle(cc)
        j = 0
        for y in mht:
            if y == 'CERES-EBAF':
                continue
            ceres = mht['CERES-EBAF'].extract_cube(x + '_mht')
            cube = mht[y].extract_cube(x + '_mht')
            ax.plot(cube.coord('latitude').points,
                    (cube.data - ceres.data) / 1e15,
                    label=y)
            del cube
            j += 1

        ymin, ymax = ax.get_ylim()
        ylim = np.max([np.abs(ymin), np.abs(ymax)])
        ax.axhline(0, -90, 90, color='k', ls=':')
        ax.axvline(0, -10, 10, color='k', ls=':')
        ax.set_ylim(-ylim, ylim)
        ax.set_xlim(-90, 90)
        ax.set_xticks(np.arange(-90, 120, 30))
        ax.set_xlabel('Latitude')
        ax.set_title(var_dict[x] + ': (model - CERES-EBAF)')
        del ymin, ymax, ylim

        handles, labels = ax.get_legend_handles_labels()

    # Legend spans the width of the difference subplots
    plt.subplot(gs[1:, 1:])
    plt.axis('off')
    plt.legend(handles, labels, ncol=3, loc='upper center')

    plt.tight_layout()

    # Save plots
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=cfg)
    logger.info('{0}'.format(provenance_record))
    save_figure(pltname + '_CERES-EBAF', provenance_record, cfg)
    plt.close()
    ''' MHT plot: compares full MHT with CERES-EBAF '''
    # Reset color cycler
    cc = (cycler(linestyle=['-', '--', '-.']) * cycler('color', [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]))

    # Different size plots depending on number of variables
    plt.figure(figsize=(15, 6))
    if len(variables) == 2:
        gs = gridspec.GridSpec(2, 2, height_ratios=[7, 1])
    elif len(variables) == 3:
        gs = gridspec.GridSpec(2, 3, height_ratios=[7, 1])

    # Plot MHT for all models including CERES-EBAF (black line)
    for i, x in enumerate(variables):
        ax = plt.subplot(gs[:1, i:i + 1])
        ax.set_prop_cycle(cc)
        for y in mht:
            cube = mht[y].extract_cube(x + '_mht')
            if y == 'CERES-EBAF':
                ax.plot(cube.coord('latitude').points,
                        cube.data / 1e15,
                        'k',
                        ls='-',
                        lw=2,
                        label=y)
            else:
                ax.plot(cube.coord('latitude').points,
                        cube.data / 1e15,
                        label=y)

        ax.set_ylim(-vlim, vlim)
        ax.set_xlim(-90, 90)
        ax.set_xticks(np.arange(-90, 120, 30))
        ax.axhline(0, -90, 90, color='k', ls=':')
        ax.axvline(0, -vlim, vlim, color='k', ls=':')
        ax.set_xlabel('Latitude')
        ax.set_ylabel('MHT (PW)')
        ax.set_title(var_dict[x])

        handles, labels = ax.get_legend_handles_labels()
        del cube

    # Legend spans the width of the full plot
    plt.subplot(gs[1:, :])
    plt.axis('off')
    plt.legend(handles, labels, ncol=3, loc='upper center')

    plt.tight_layout()

    # Save plot
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=cfg)
    logger.info('{0}'.format(provenance_record))
    save_figure(pltname + '_comp', provenance_record, cfg)


def quiver_plot(variables, cubelist, ancestor_files, cfg, pltname):
    """Maps of scalar p-fields or input flux fields Produces a single plot with
    multiple subplots comparing the CERES-EBAF field to the different model
    fields.

    The CERES-EBAF field is shown is full, with all variables sharing
    the same colorbar. The difference between each model and CERES-EBAF
    are shown as separate subplots for each variables. All of the
    difference subplots share the same colorbar. For the scalar p-field
    maps the direction of the implied heat transport is shown by white
    arrows. The arrows are normalised across the different plots and
    variables (separate for the difference plots and CERES-EBAF), so the
    lengths can be directly compared for the magnitude of the heat
    transport.
    """

    # Variables
    var_dict = {
        'rtnt': 'Net TOA',
        'rsnt': 'Net SW',
        'rlut': 'Net LW',
        'rlutcs': '-1 x OLR (clear)',
        'rsut': '-1 x OSR (all)',
        'rsutcs': '-1 x OSR (clear)',
        'rsdt': 'ISR',
        'netcre': 'NET CRE',
        'swcre': 'SW CRE',
        'lwcre': 'LW CRE'
    }

    # Data type
    dtype = pltname.split('_')[-1]
    if dtype == 'flux':
        cbar_label = 'Flux (W/m2)'
        div = 1
    elif dtype == 'poisson':
        cbar_label = 'Scalar p-field (PW)'
        div = 1e15

    # Size of plot parameters
    # Change these if plots are not displaying properly
    if len(cubelist) == 2:
        ncols = 2
        nrows = 1
        xlen = 8
    elif len(cubelist) <= 5:
        ncols = 2
        nrows = int(np.ceil((len(cubelist) - 1) / 2))
        xlen = 10
    elif len(cubelist) <= 7:
        ncols = 3
        nrows = 2
        xlen = 15
    else:
        ncols = 4
        nrows = int(np.ceil((len(cubelist) - 1) / 4))
        xlen = 20

    # Create figure axes
    height_ratio = np.ones((nrows * len(variables)) + 1)
    height_ratio[:(nrows * len(variables))] *= 7

    plt.figure(figsize=(xlen, 5 * nrows))
    gs = gridspec.GridSpec((nrows * len(variables)) + 1,
                           ncols + 1,
                           height_ratios=height_ratio)
    gs.update(hspace=0.5)

    vmin = 0
    vmax = 0
    wmin = 0
    wmax = 0

    # Calculate the plot limits
    # - Same colorbar for CERES plots
    # - Same colorbar for difference plots
    for i, p in enumerate(variables):
        for j, q in enumerate(cubelist):
            if q == 'CERES-EBAF':
                cube = cubelist[q].extract_cube(p + '_' + dtype)
                tmp = np.average(cube.data)

                plt.figure()
                plt.axes(projection=ccrs.PlateCarree())
                plt.contourf(cube.coord('longitude').points,
                             cube.coord('latitude').points,
                             tmp - cube.data,
                             levels=10,
                             transform=ccrs.PlateCarree(central_longitude=0))
                cbar = plt.colorbar()
                cmin, cmax = cbar.mappable.get_clim()
                vmin = np.min([vmin, cmin])
                vmax = np.max([vmax, cmax])
                del tmp, cube
                plt.close()

            else:
                cube = cubelist[q].extract_cube(p + '_' + dtype)
                tmp_cube = cube.data - np.average(cube.data)
                ceres = cubelist['CERES-EBAF'].extract_cube(p + '_' + dtype)
                tmp_ceres = ceres.data - np.average(ceres.data)

                plt.figure()
                plt.axes(projection=ccrs.PlateCarree())
                plt.contourf(cube.coord('longitude').points,
                             cube.coord('latitude').points,
                             tmp_cube - tmp_ceres,
                             levels=10,
                             transform=ccrs.PlateCarree(central_longitude=0))
                cbar = plt.colorbar()
                cmin, cmax = cbar.mappable.get_clim()
                wmin = np.min([wmin, cmin])
                wmax = np.max([wmax, cmax])
                del cube, ceres, tmp_cube, tmp_ceres
                plt.close()

    wlim = np.max([np.abs(wmin), np.abs(wmax)])
    levels1 = np.linspace(vmin / div, vmax / div, 11)  # CERES
    levels2 = np.linspace(-wlim / div, wlim / div, 11)  # difference

    # Plot the maps
    for m, p in enumerate(variables):
        j = 0
        for i, q in enumerate(cubelist):
            # CERES-EBAF data
            if q == 'CERES-EBAF':
                cube = cubelist[q].extract_cube(p + '_' + dtype)
                tmp = np.average(cube.data)
                ax1 = plt.subplot(gs[(m * nrows):(m * nrows) + 1, 0],
                                  projection=ccrs.PlateCarree())
                cb1 = ax1.contourf(
                    cube.coord('longitude').points,
                    cube.coord('latitude').points, (tmp - cube.data) / div,
                    levels=levels1,
                    transform=ccrs.PlateCarree(central_longitude=0))
                plt.gca().coastlines()

                if dtype == 'poisson':
                    # Plot heat transport arrows on p-field
                    x, y = np.meshgrid(
                        cube.coord('longitude').points,
                        cube.coord('latitude').points)
                    v, u = np.gradient(cube.data, 1e14, 1e14)
                    if m == 0:
                        Q = ax1.quiver(x[15::25, 15::32],
                                       y[15::25, 15::32],
                                       u[15::25, 15::32],
                                       v[15::25, 15::32],
                                       pivot='mid',
                                       color='w')
                        Q._init()
                    else:
                        ax1.quiver(x[15::25, 15::32],
                                   y[15::25, 15::32],
                                   u[15::25, 15::32],
                                   v[15::25, 15::32],
                                   scale=Q.scale,
                                   pivot='mid',
                                   color='w')
                    del x, y, u, v

                ax1.set_xticks(np.arange(-180, 190, 60))
                ax1.set_xticklabels(
                    ['180', '120W', '60W', '0', '60E', '120E', '180'])
                ax1.set_yticks(np.arange(-90, 100, 30))
                ax1.set_yticklabels(
                    ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
                ax1.set_title('{0}: CERES-EBAF'.format(var_dict[p]))
                del tmp

            # Difference between models and CERES-EBAF
            else:
                # Subplot details
                if (i + 1) >= ncols:
                    if i % ncols == 0:
                        j += 1
                    k = (i + 1) - (j * ncols)
                else:
                    k = i + 1

                cube = cubelist[q].extract_cube(p + '_' + dtype)
                tmp_cube = cube.data - np.average(cube.data)
                ceres = cubelist['CERES-EBAF'].extract_cube(p + '_' + dtype)
                tmp_ceres = ceres.data - np.average(ceres.data)

                ax2 = plt.subplot(gs[(m * nrows) + j:(m * nrows) + j + 1,
                                     k:k + 1],
                                  projection=ccrs.PlateCarree())
                cb2 = ax2.contourf(
                    cube.coord('longitude').points,
                    cube.coord('latitude').points,
                    (tmp_cube - tmp_ceres) / div,
                    levels=levels2,
                    transform=ccrs.PlateCarree(central_longitude=0),
                    cmap='coolwarm')
                plt.gca().coastlines(color='k')

                if dtype == 'poisson':
                    # Heat transport arrows
                    x, y = np.meshgrid(
                        cube.coord('longitude').points,
                        cube.coord('latitude').points)
                    v, u = np.gradient(tmp_cube - tmp_ceres, 1e14, 1e14)
                    if m == 0:
                        R = ax2.quiver(x[15::25, 15::32],
                                       y[15::25, 15::32],
                                       u[15::25, 15::32],
                                       v[15::25, 15::32],
                                       pivot='mid',
                                       color='w')
                        R._init()
                    else:
                        ax2.quiver(x[15::25, 15::32],
                                   y[15::25, 15::32],
                                   u[15::25, 15::32],
                                   v[15::25, 15::32],
                                   scale=R.scale,
                                   pivot='mid',
                                   color='w')
                    del x, y, u, v

                ax2.set_xticks(np.arange(-180, 190, 60))
                ax2.set_xticklabels(
                    ['180', '120W', '60W', '0', '60E', '120E', '180'])
                ax2.set_yticks(np.arange(-90, 100, 30))
                ax2.set_yticklabels(
                    ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
                if len(cubelist) > 3:
                    ax2.set_title('{0}'.format(q))
                else:
                    ax2.set_title('{0} - CERES-EBAF'.format(q))
                del cube, tmp_cube, ceres, tmp_ceres

    # Plot colorbars
    ax1 = plt.subplot(gs[-1, 0])
    plt.colorbar(cb1, cax=ax1, orientation='horizontal', label=cbar_label)

    ax2 = plt.subplot(gs[-1, 1:])
    plt.colorbar(cb2,
                 cax=ax2,
                 orientation='horizontal',
                 label='Anomaly ' + cbar_label.split(' ')[-1])

    # Plot adjustments
    if len(variables) == 2:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.97, bottom=0.1)
    elif len(variables) < 5:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.97, bottom=0.07)
    elif len(variables) == 5:
        plt.subplots_adjust(left=0.06, right=0.94, top=0.97, bottom=0.07)

    # Save plot
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=cfg)
    save_figure(pltname, provenance_record, cfg)


def plot_multi_model_diagnostics(poisson, mht, flux, ancestor_files, cfg):
    basename = 'multi_model_implied'

    # TOA: net, sw, lw
    variables = ['rtnt', 'rsnt', 'rlut']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_toa_mht')
    quiver_plot(variables,
                poisson,
                ancestor_files,
                cfg,
                pltname=basename + '_toa_poisson')
    quiver_plot(variables,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_toa_flux')

    # OSR: clear vs. all sky
    variables = ['rsut', 'rsutcs']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_osr_mht')
    quiver_plot(variables,
                poisson,
                ancestor_files,
                cfg,
                pltname=basename + '_osr_poisson')
    quiver_plot(variables,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_osr_flux')

    # CRE: net, sw, lw
    variables = ['netcre', 'swcre', 'lwcre']
    mht_plot(variables,
             mht,
             ancestor_files,
             cfg,
             pltname=basename + '_cre_mht')
    quiver_plot(variables,
                poisson,
                ancestor_files,
                cfg,
                pltname=basename + '_cre_poisson')
    quiver_plot(variables,
                flux,
                ancestor_files,
                cfg,
                pltname=basename + '_cre_flux')


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""
    input_data = cfg['input_data'].values()

    # Dictionaries for data cubes
    model_p = {}
    model_mht = {}
    model_flux = {}
    ceres = False

    # Group data by dataset
    logger.info("Group input data by model, sort by experiment")
    model_dataset = group_metadata(input_data, 'dataset', sort='exp')

    for model_name in model_dataset:
        # Treat CERES-EBAF data separately
        if model_name == 'CERES-EBAF':
            ceres = True
            continue

        logger.info("Processing model data: %s", model_name)
        p, mht, flux = compute_p_field_diagnostic(model_dataset[model_name],
                                                  cfg)
        model_p.update(p)
        model_mht.update(mht)
        model_flux.update(flux)
        del p, mht, flux

    if ceres is True:
        model_name = 'CERES-EBAF'
        logger.info("Processing observational data: CERES-EBAF")
        # Need to calculate some fluxes for CERES-EBAF
        p, mht, flux = ceres_p_field_diagnostic(model_dataset[model_name], cfg)
        model_p.update(p)
        model_mht.update(mht)
        model_flux.update(flux)
        del p, mht, flux

    # Plotting multi-model diagnostics
    plot_multi_model_diagnostics(model_p, model_mht, model_flux, input_data,
                                 cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
