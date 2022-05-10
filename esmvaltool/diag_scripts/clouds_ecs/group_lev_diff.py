"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
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
    'cl': 'cloud_fraction',
    'cli': 'ice_water_content',
    'clw': 'liquid_water_content',
}
LINE_LEGEND = {
    'ECS_high_hist': 'ECS_high',
    'ECS_med_hist': 'ECS_med',
    'ECS_low_hist': 'ECS_low',
}
LINE_COLOR = {
    'ECS_high_hist': 'royalblue',
    'ECS_high_scen': 'royalblue',
    'ECS_med_hist': 'green',
    'ECS_med_scen': 'green',
    'ECS_low_hist': 'orange',
    'ECS_low_scen': 'orange',
    'CMIP6': 'firebrick',
    'CMIP5': 'royalblue',
    'CMIP3': 'darkcyan',
    'OBS': 'black'
}
LINE_DASH = {
    'ECS_high_hist': 'solid',
    'ECS_high_scen': 'dashed',
    'ECS_med_hist': 'solid',
    'ECS_med_scen': 'dashed',
    'ECS_low_hist': 'solid',
    'ECS_low_scen': 'dashed',
    'CMIP6': 'solid',
    'CMIP5': 'solid',
    'CMIP3': 'solid',
    'OBS': 'solid'
}
FIGURE_NUMBER = {
    'ECS_high_hist': 231,
    'ECS_med_hist': 232,
    'ECS_low_hist': 233,
    'ECS_high_scen': 425,
    'ECS_med_scen': 426,
    'ECS_low_scen': 427,
    'OBS': 424
}
FIGURE_NUMBER_DIFF = {
    'ECS_high_hist': 234,
    'ECS_med_hist': 235,
    'ECS_low_hist': 236,
}

def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    #print(attributes)
    caption = ("Average {short_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

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

def _get_multi_model_mean(cubes, var):
    """Compute multi-model mean."""

    logger.debug("Calculating multi-model mean")
    datasets = []
    mmm = []
    for (dataset_name, cube) in cubes.items():
        datasets.append(dataset_name)
        mmm.append(cube.data)
    mmm = np.ma.array(mmm)
    dataset_0 = list(cubes.keys())[0]
    mmm_cube = cubes[dataset_0].copy(data=np.ma.mean(mmm, axis=0))
    attributes = {
        'dataset': 'MultiModelMean',
        'short_name': var,
        'datasets': '|'.join(datasets),
    }
    mmm_cube.attributes = attributes
    print(mmm_cube)
    return  mmm_cube


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    if cube.var_name == 'cli':
        cube.convert_units('g/kg')
    elif cube.var_name == 'clw':
        cube.convert_units('g/kg')

    logger.debug("Reading data")
    cube = iris.util.squeeze(cube)
    return cube


def compute_diff(filename1, filename2):
    """Compute difference between two cubes."""
    logger.debug("Loading %s", filename1)
    cube1 = iris.load_cube(filename1)
    cube2 = iris.load_cube(filename2)

    if cube1.var_name == 'cli':
        cube1.convert_units('g/kg')
        cube2.convert_units('g/kg')
    elif cube1.var_name == 'clw':
        cube1.convert_units('g/kg')
        cube2.convert_units('g/kg')

    cube = cube2 - cube1
    cube.metadata = cube1.metadata
    return cube

def compute_diff_temp(input_data, group, dataset):
    """Compute relative change per temperture change."""

    dataset_name = dataset['dataset']
    var = dataset['short_name']

    input_file_1 = dataset['filename']

    var_data_2 = select_metadata(input_data,
                                 short_name=var,
                                 dataset=dataset_name,
                                 variable_group=group[1]) 
    if not var_data_2:
        raise ValueError(
            f"No '{var}' data for '{dataset_name}' in '{group[1]}' available")

    input_file_2 = var_data_2[0]['filename']

    ta_data_1 = select_metadata(input_data,
                              short_name='ta',
                              dataset=dataset_name,
                              variable_group='ta_'+group[0]) 
    ta_data_2 = select_metadata(input_data,
                              short_name='ta',
                              dataset=dataset_name,
                              variable_group='ta_'+group[1]) 
    if not ta_data_1:
        raise ValueError(
            f"No 'ta' data for '{dataset_name}' in '{group[0]}' available")
    if not ta_data_2:
        raise ValueError(
            f"No 'ta' data for '{dataset_name}' in '{group[1]}' available")
    input_file_ta_1 = ta_data_1[0]['filename']
    input_file_ta_2 = ta_data_2[0]['filename']

    cube = compute_diagnostic(input_file_1)
    cube.data[cube.data < 0.001] = 0.0

    cube_diff = compute_diff(input_file_1, input_file_2)
    cube_ta_diff = compute_diff(input_file_ta_1, input_file_ta_2)

    #cube_diff = cube_diff / cube_ta_diff
    cube_diff = 100. * (cube_diff / cube) / cube_ta_diff

    #cube_diff.units = 'g/kg/K'
    cube_diff.units = '%/K'
    
    return cube_diff


def plot_model(cube, attributes, plot_type, cfg):
    """Plot each single model."""

    plt.ylim(1000.,100.)
    plt.yscale('log')
    plt.yticks([1000., 800., 600., 400., 300., 200., 100.], [1000, 800, 600, 400, 300, 200, 100])
    cube.coord('air_pressure').convert_units('hPa')
    if cube.var_name == 'cl':
        levels = np.linspace(0., 50., 11)
    elif cube.var_name == 'cli':
        levels = np.linspace(0., 0.02, 11)
        #cube.convert_units('g/kg')
    elif cube.var_name == 'clw':
        levels = np.linspace(0., 0.05, 11)
        #cube.convert_units('g/kg')
    qplt.contourf(cube, levels=levels, extend='max')

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


def plot_diagnostic(cube, legend, plot_type, cfg):
    """Create diagnostic data and plot it."""

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
    else:
        nplot = FIGURE_NUMBER.get(legend, legend)

        plt.subplot(nplot)

        plt.ylim(1000.,100.)
        plt.yscale('log')
        plt.yticks([1000., 800., 600., 400., 300., 200., 100.], [1000, 800, 600, 400, 300, 200, 100])
        cube.coord('air_pressure').convert_units('hPa')
        if cube.var_name == 'cl':
            levels = np.linspace(0., 50., 11)
        elif cube.var_name == 'cli':
            levels = np.linspace(0., 0.02, 11)
            cube.convert_units('g/kg')
        elif cube.var_name == 'clw':
            levels = np.linspace(0., 0.05, 11)
            cube.convert_units('g/kg')
        qplt.contourf(cube, levels=levels, extend='max')

        logger.info("Plotting %s", legend)


def plot_diagnostic_diff(cube, legend, plot_type, cfg):
    """Create diagnostic data and plot it."""

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
    else:
        nplot = FIGURE_NUMBER_DIFF.get(legend, legend)

        plt.subplot(nplot)

        plt.ylim(1000.,100.)
        plt.yscale('log')
        plt.yticks([1000., 800., 600., 400., 300., 200., 100.], [1000, 800, 600, 400, 300, 200, 100])
        cube.coord('air_pressure').convert_units('hPa')
        levels = np.linspace(-20., 20., 21)
        if cube.var_name == 'cl':
            levels = np.linspace(-6., 6., 13)
        elif cube.var_name == 'cli':
            levels = np.linspace(-0.01, 0.01, 11)
            cube.convert_units('g/kg')
        elif cube.var_name == 'clw':
            levels = np.linspace(-0.01, 0.01, 11)
            cube.convert_units('g/kg')
        qplt.contourf(cube, levels=levels, extend='both', cmap='coolwarm')
        #qplt.contourf(cube, extend='both', cmap='coolwarm')

        #plt.colorbar(label = '%/K')
        #cbar = plt.colorbar()
        #cbar.set_label('%/K')

        logger.info("Plotting %s", legend)


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    cfg.setdefault('filename_attach', 'base')
    cfg.setdefault('plot_each_model', False)
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    plot_type = cfg['plot_type']

    input_data = list(cfg['input_data'].values())
    all_vars = list(group_metadata(input_data, 'short_name'))

    groups = group_metadata(input_data, 'variable_group', sort='dataset')

    cubes = iris.cube.CubeList()

    plt.figure(constrained_layout=True, figsize=(12, 8))

    ngroups = len(groups)

    print("Number of variable  groups: ", ngroups)

    for group_name in groups:
        if 'hist' in group_name:
            if 'ta_' not in group_name:
                logger.info("Processing variable %s", group_name)

                for attributes in groups[group_name]:
                    logger.info("Loop dataset %s", attributes['dataset'])

                    input_file = attributes['filename']
                    cube = compute_diagnostic(input_file)

                    if cfg['plot_each_model']:
                        plot_model(cube, attributes, plot_type, cfg)

                    if attributes['dataset'] == 'MultiModelMean' or group_name == 'OBS':
                      logger.info("Processing dataset %s", attributes['dataset'])
                      cubes.append(cube)

                      plot_diagnostic(cube, group_name, plot_type, cfg)

                      if plot_type == 'height':
                        title = group_name
                      else:
                        title = attributes['long_name']

                      plt.title(title, fontsize=10)


    for group_name in cfg['group_by']:

        logger.info("Processing group %s", group_name[0])

        dataset_names = []
        cubes_diff = {}

        for dataset in groups[group_name[0]]:
            dataset_name = dataset['dataset']
            var = dataset['short_name']

            if dataset_name != 'MultiModelMean':
                logger.info("Loop dataset %s", dataset_name)
                dataset_names.append(dataset_name)

                cube_diff = compute_diff_temp(input_data, group_name, dataset) 
                
                cubes_diff[dataset_name] = cube_diff

        cube_mmm = _get_multi_model_mean(cubes_diff, var)

        plot_diagnostic_diff(cube_mmm, group_name[0], plot_type, cfg)


        if plot_type == 'height':
          title = group_name[1] + " - " + group_name[0]
        else:
          title = cube_mmm['short_name']

        plt.title(title, fontsize=9)
        #plt.ylabel('Air pressure (hPa)', ncol=1)
        #plt.legend(ncol=1)

    plt.suptitle(var)

    provenance_record = get_provenance_record(
        attributes, ancestor_files=cfg['input_files'])

    basename = 'diff_' + var + '_' + cfg['filename_attach']
    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cubes)

    # And save the plot
    save_figure(basename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
