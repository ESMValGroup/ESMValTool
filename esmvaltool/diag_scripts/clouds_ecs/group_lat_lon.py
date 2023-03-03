"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
from iris.analysis.stats import pearsonr
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from scipy.stats import bootstrap

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    extract_variables,
    group_metadata,
    run_diagnostic,
    get_diagnostic_filename,
    get_plot_filename,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
    io,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)

VAR_NAMES = {
    'clt': 'total_cloud_fraction',
    'lwp': 'liquid_water_path',
    'clivi': 'ice_water_path',
    'netcre': 'net_cre',
    'swcre': 'sw_cre',
    'lwcre': 'lw_cre',
}
PANEL = {
    'ECS_high': 222, # (0, 1),
    'ECS_med':  223, # (1, 0),
    'ECS_low':  224, # (1, 1),
    'OBS':      221  # (0, 0)
}
PANEL_LABELS = {
    'ECS_high': 'b)', # (0, 1),
    'ECS_med':  'c)', # (1, 0),
    'ECS_low':  'd)', # (1, 1),
    'OBS':      'a)'  # (0, 0)
}
PANDAS_PRINT_OPTIONS = ['display.max_rows', None, 'display.max_colwidth', -1]


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Climatology of {standard_name}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
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


def _get_cube_list(input_files):
    """Get :class:`iris.cube.CubeList` of input files."""
    cubes = iris.cube.CubeList()

    # Input files
    for filename in input_files:
        logger.info("Loading '%s'", filename)
        cube = _load_cube_with_dataset_coord(filename)
        cube.attributes['filename'] = filename
        cubes.append(cube)

    # Check metadata of cubes
    for cube in cubes:
        check_metadata(cube.attributes)

    return cubes


def area_weighted_mean(cube):
    logger.debug("Computing field mean")
    grid_areas = iris.analysis.cartography.area_weights(cube)
    mean = cube.collapsed(['longitude', 'latitude'],
                               iris.analysis.MEAN,
                                           weights=grid_areas)
    return mean


def calculate_bias(model_cube, obs_cube):
    logger.debug("Computing bias")
    diff = model_cube - obs_cube
    bias = area_weighted_mean(diff)
    bias.attributes = model_cube.attributes
    return bias


def calculate_rmsd(model_cube, obs_cube):
    logger.debug("Computing RMSD")
    diff = model_cube - obs_cube
    rmsd = area_weighted_mean(diff**2)**0.5
    rmsd.attributes = model_cube.attributes
    return rmsd


def calculate_corr(model_cube, obs_cube):
    logger.debug("Computing Correlation")
    grid_areas = iris.analysis.cartography.area_weights(model_cube)
    corr = pearsonr(model_cube, obs_cube, weights=grid_areas)
    #corr = pearsonr(model_cube, obs_cube)
    return corr


def compute_diagnostic(filename):
    
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    cube = iris.util.squeeze(cube)
    return cube


def plot_model(cube, attributes, cfg):
    # Plot each model.

    levels = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    if attributes['short_name'] == 'clt':
        levels = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        cmap = 'viridis'
    elif attributes['short_name'] == 'clivi':
        levels = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        cmap = 'viridis'
    elif attributes['short_name'] == 'lwp':
        levels = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        cmap = 'viridis'
    elif attributes['short_name'] == 'netcre':
        levels = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
        cmap = 'bwr' 
    elif attributes['short_name'] == 'lwcre':
        #levels = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
        levels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
        cmap = 'Reds' 
    elif attributes['short_name'] == 'swcre':
        levels = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0]
        cmap = 'Blues_r' 
    plt.axes(projection=ccrs.Robinson())
    im = iplt.contourf(cube, levels=levels, cmap=cmap, extend ='both')
    #plt.clim(0., 100.)
    plt.gca().coastlines()
    colorbar = plt.colorbar(orientation='horizontal')
    colorbar.set_label( cube.var_name + '/' + cube.units.origin)
    if attributes['short_name'] == 'clt':
        ticks = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    elif attributes['short_name'] == 'clivi':
        ticks = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    elif attributes['short_name'] == 'lwp':
        ticks = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    elif attributes['short_name'] == 'netcre':
        ticks = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
    elif attributes['short_name'] == 'lwcre':
        ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    elif attributes['short_name'] == 'swcre':
        ticks = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0]
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([str(tick) for tick in ticks])

    # Appearance
    dataset_name = attributes['dataset']
    title = f'{VAR_NAMES.get(cube.var_name, cube.var_name)} for {dataset_name}'
    filename = ('{}_{}_{}'.format(VAR_NAMES.get(cube.var_name, cube.var_name),
                                  attributes['exp'], dataset_name))

    plt.title(title)
    plot_path = get_plot_filename(filename, cfg)
    plt.savefig(plot_path,
                bbox_inches='tight',
                orientation='landscape')
    logger.info("Wrote %s", plot_path)
    plt.close()


def read_data(groups, cfg):
    logger.debug("Read data")
    cubes = iris.cube.CubeList()
    cubes_out = iris.cube.CubeList()

    for group_name in groups:
        logger.info("Processing variable %s", group_name)

        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = compute_diagnostic(input_file)
            cube.attributes['variable_group'] = group_name
            cube.attributes['dataset'] = attributes['dataset']

            cubes.append(cube)

            if attributes['dataset'] == 'MultiModelMean' or group_name == 'OBS':
                cubes_out.append(cube)
            else:
                if cfg['plot_each_model']:
                    plot_model(cube, attributes, cfg)

    return cubes, cubes_out


def plot_diagnostic(cubes, attributes, cfg):
    """Create diagnostic data and plot it."""

    fig = plt.figure(constrained_layout=True)
    fig.set_figheight(10)
    fig.set_figwidth(14)
    plt.subplots_adjust(left=0.05, bottom=0.21, right=0.95, top=0.94, wspace=0.02, hspace=0.02)

    cmap = 'bwr' 
    if attributes['short_name'] == 'clt':
        levels = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        cmap = 'viridis'
    elif attributes['short_name'] == 'clivi':
        levels = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        cmap = 'viridis'
    elif attributes['short_name'] == 'lwp':
        levels = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        cmap = 'viridis'
    elif attributes['short_name'] == 'netcre':
        levels = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
        cmap = 'bwr' 
    elif attributes['short_name'] == 'lwcre':
        levels = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
        cmap = 'Reds' 
    elif attributes['short_name'] == 'swcre':
        levels = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0]
        cmap = 'Blues_r' 
    elif attributes['short_name'] == 'clt_diff':
        levels = list(np.arange(-30, 31, 2.5))
    elif attributes['short_name'] == 'clivi_diff':
        levels = list(np.arange(-0.1, 0.105, 0.01))
    elif attributes['short_name'] == 'lwp_diff':
        levels = list(np.arange(-0.1, 0.105, 0.01))
    elif attributes['short_name'] in ['netcre_diff', 'lwcre_diff', 'swcre_diff']:
        levels = list(np.arange(-30, 31, 2.5))

    for cube in cubes:
        logger.info("Plotting %s %s of group %s", cube.attributes['dataset'], attributes['short_name'], cube.attributes['variable_group'])
        mean = area_weighted_mean(cube)

        legend = cube.attributes['variable_group']

        ipanel = PANEL.get(legend, None)
        plt.subplot(ipanel, projection=ccrs.Robinson())

        im = iplt.contourf(cube, levels=levels, cmap=cmap, extend ='both')

        plt.gca().coastlines()
        plt.title(legend, fontsize=18)
        if attributes['short_name'] in ['clt', 'netcre']:
            plt.title('mean = {:.1f}      '.format(mean.data), fontsize = 14, loc='right')
        elif attributes['short_name'] in ['clivi', 'lwp']:
            plt.title('mean = {:.3f}      '.format(mean.data), fontsize = 14, loc='right')
        elif attributes['short_name'] in ['clivi_diff', 'lwp_diff']:
            plt.title('bias = {:.3f}      '.format(mean.data), fontsize = 14, loc='right')
        elif attributes['short_name'] in ['clt_diff', 'netcre_diff']:
            plt.title('bias = {:.1f}      '.format(mean.data), fontsize = 14, loc='right')
        else:
            plt.title('{:.1f}      '.format(mean.data), fontsize = 14, loc='right')
        ipanel_label = PANEL_LABELS.get(legend, None)
        plt.title(ipanel_label, fontsize = 22, loc='left')


    title = attributes['long_name']
    fig.suptitle(title, fontsize = 22)
    cbar_ax = fig.add_axes([0.2, 0.18, 0.6, 0.03])
    colorbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    colorbar.set_label( cube.var_name + '/' + cube.units.origin, fontsize = 16)
    if attributes['short_name'] == 'clt':
        ticks = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    elif attributes['short_name'] == 'clivi':
        ticks = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    elif attributes['short_name'] == 'lwp':
        ticks = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    elif attributes['short_name'] == 'netcre':
        ticks = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
    elif attributes['short_name'] == 'lwcre':
        ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
    elif attributes['short_name'] == 'swcre':
        ticks = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0]

    elif attributes['short_name'] == 'clt_diff':
        ticks = list(np.arange(-30,31,5))
    elif attributes['short_name'] == 'clivi_diff':
        ticks = [-0.1, -0.08, -0.06, -0.04, -0.02, 0., 0.02, 0.04, 0.06, 0.08, 0.1]
    elif attributes['short_name'] == 'lwp_diff':
        ticks = [-0.1, -0.08, -0.06, -0.04, -0.02, 0., 0.02, 0.04, 0.06, 0.08, 0.1]
        #ticks = list(np.arange(-0.1, 0.12, 0.02))
    elif attributes['short_name'] in ['netcre_diff', 'lwcre_diff', 'swcre_diff']:
        #ticks = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50]
        ticks = list(np.arange(-30,31,5))

    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([str(tick) for tick in ticks], fontsize = 16)

    # Save the data and the plot
    provenance_record = get_provenance_record(
            attributes, ancestor_files=cfg['input_files'])
    basename = 'map_' + attributes['short_name']

    save_data(basename, provenance_record, cfg, cubes)
    save_figure(basename, provenance_record, cfg)


def get_dataframe(cubes, cube_obs):

    df = pd.DataFrame(columns=['Dataset', 'Group', 'Statistic', 'Value'])
    idf = 0
    
    print(cubes)

    for cube in cubes:
        dataset = cube.attributes['dataset']
        group = cube.attributes['variable_group']
        logger.info("Computing statistics of dataset %s", dataset)

        mean = area_weighted_mean(cube)
        bias = calculate_bias(cube, cube_obs)
        rmsd = calculate_rmsd(cube, cube_obs)
        corr = calculate_corr(cube, cube_obs)

        df.loc[idf] = [dataset, group, 'Mean', mean.data]
        idf = idf + 1
        df.loc[idf] = [dataset, group, 'Bias', bias.data]
        idf = idf + 1
        df.loc[idf] = [dataset, group, 'RMSD', rmsd.data]
        idf = idf + 1
        df.loc[idf] = [dataset, group, 'Corr', corr.data]
        idf = idf + 1

    return df


def write_statistics(df, attributes, cfg):
    df['Value'] = df['Value'].astype(str).astype(float)
     
    basename = "statistic_all_" + attributes['short_name']
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    df.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Data:\n%s", df)

    stat = df.groupby(['Statistic', 'Group'])['Value'].describe()
    basename = "statistic_" + attributes['short_name']
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    stat.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Data:\n%s", df)


def bootstrapping(cubes, cube_obs, all_groups, attributes, cfg):
    logger.info("Bootstrapping")

    for group in all_groups:
        if group != 'OBS':
            logger.info("Processing group %s", group)
            cubes_part = {}
            datasets = []
            for cube in cubes:
                if cube.attributes['variable_group'] == group:
                    dataset = cube.attributes['dataset']
                    cubes_part[dataset] = cube
                    datasets.append(dataset)

            nsample = 1000
            sample_stat = pd.DataFrame(columns=['Mean', 'Bias', 'RMSD', 'Corr'])

            ncubes = len(cubes_part)
            array = list(np.arange(0, ncubes))
            for iboot in range(0, nsample):
                ires = random.choices(array, k=len(array))
                for i, icube in enumerate(ires):
                    if i == 0:
                        cube = cubes_part[datasets[icube]].copy()
                    else:
                        cube += cubes_part[datasets[icube]]
                cube.data = cube.data / ncubes
                mean = area_weighted_mean(cube).data
                bias = calculate_bias(cube, cube_obs).data
                rmsd = calculate_rmsd(cube, cube_obs).data
                corr = calculate_corr(cube, cube_obs).data
                sample_stat.loc[iboot] = [mean, bias, rmsd, corr]
                sample_stat[:] = sample_stat[:].astype(str).astype(float)

            stat = sample_stat.describe()
            basename = "bootstrapping_" + attributes['short_name'] + "_" + group
            csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
            stat.to_csv(csv_path)
            logger.info("Wrote %s", csv_path)


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    cfg.setdefault('plot_each_model', False)
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    input_data = list(cfg['input_data'].values())

    groups = group_metadata(input_data, 'variable_group', sort='dataset')
    attributes = next(iter(extract_variables(cfg).values()))
    all_groups = list(group_metadata(input_data, 'variable_group'))
    ngroups = len(all_groups)
    print('all_groups')
    print(all_groups)

    # Read data
    cubes, cubes_out = read_data(groups, cfg)

    # Plotting climatologies
    plot_diagnostic(cubes_out, attributes, cfg)

    # Compute bias plots
    cube_obs = cubes_out.extract_cube(iris.Constraint
               (cube_func=lambda cube: cube.attributes['variable_group']=='OBS'))

    # Bootstrapping
    bootstrapping(cubes, cube_obs, all_groups, attributes, cfg)

    # Compute statistics
    df = get_dataframe(cubes, cube_obs)

    # write statistics
    write_statistics(df, attributes, cfg)

    # compute bias
    cubes_diff = iris.cube.CubeList()
    attributes['short_name'] = attributes['short_name'] + "_diff"

    for cube in cubes_out:
        if cube.attributes['variable_group'] != 'OBS' or cube.attributes['dataset'] != 'MultiModelMean':
            logger.info("Processing %s of group %s", cube.attributes['dataset'], cube.attributes['variable_group'])
            mean = area_weighted_mean(cube)
            bias = calculate_bias(cube, cube_obs)
            rmsd = calculate_rmsd(cube, cube_obs)
            corr = calculate_corr(cube, cube_obs)
            cube_diff = cube - cube_obs
            cube_diff.attributes = cube.attributes
            cube_diff.var_name = cube.var_name
            cube_diff.attributes['short_name'] = attributes['short_name']
            cubes_diff.append(cube_diff)
            print('{0} : bias = {1}, rmsd = {2}, corr = {3}'
                  .format(cube.attributes['variable_group'], bias.data, rmsd.data, corr.data))

    # Plotting biases
    plot_diagnostic(cubes_diff, attributes, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
