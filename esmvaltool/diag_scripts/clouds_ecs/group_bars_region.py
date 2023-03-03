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
import pandas as pd
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
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
    'cl': 'cloud_fraction',
    'cli': 'ice_water_content',
    'clw': 'liquid_water_content',
}
PALETTE = {
    'high ECS': 'royalblue',
    'med ECS': 'green',
    'low ECS': 'orange',
}
PANDAS_PRINT_OPTIONS = ['display.max_rows', None, 'display.max_colwidth', -1]


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
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
    dataset_names = []
    #mmm = []
    mmm = {}
    for (dataset, cube) in cubes.items():
        dataset_names.append(dataset)
        #mmm.append(cube.data)
        mmm['dataset'] = cube.data
    #mmm = np.ma.array(mmm)
    mmm = np.ma.masked_invalid(list(mmm.values()))
    mmm_cube = cube.copy(data=np.ma.mean(mmm, axis=0))
    #dataset_0 = list(cubes.keys())[0]
    #print("dataset_0")
    #print(dataset_0)
    #mmm_cube = cubes[dataset_0].copy(data=np.ma.mean(mmm, axis=0))
    #print(mmm_cube)
    attributes = {
        'dataset': 'MultiModelMean',
        'short_name': var,
        'datasets': '|'.join(dataset_names),
    }
    mmm_cube.attributes = attributes
    return  mmm_cube


def read_data(filename):
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

def compute_diff_temp(input_data, group, var, dataset):
    """Compute relative change per temperture change."""

    dataset_name = dataset['dataset']
    var = dataset['short_name']

    input_file_1 = dataset['filename']

    var_data_2 = select_metadata(input_data,
                                 short_name=var,
                                 dataset=dataset_name,
                                 variable_group=var+"_"+group[1]) 
    if not var_data_2:
        raise ValueError(
            f"No '{var}' data for '{dataset_name}' in '{group[1]}' available")

    input_file_2 = var_data_2[0]['filename']

    tas_data_1 = select_metadata(input_data,
                              short_name='tas',
                              dataset=dataset_name,
                              variable_group='tas_'+group[0]) 
    tas_data_2 = select_metadata(input_data,
                              short_name='tas',
                              dataset=dataset_name,
                              variable_group='tas_'+group[1]) 
    if not tas_data_1:
        raise ValueError(
            f"No 'tas' data for '{dataset_name}' in '{group[0]}' available")
    if not tas_data_2:
        raise ValueError(
            f"No 'tas' data for '{dataset_name}' in '{group[1]}' available")
    input_file_tas_1 = tas_data_1[0]['filename']
    input_file_tas_2 = tas_data_2[0]['filename']

    cube = read_data(input_file_1)

    cube_diff = compute_diff(input_file_1, input_file_2)
    cube_tas_diff = compute_diff(input_file_tas_1, input_file_tas_2)

    #cube_diff = 100. * (cube_diff / cube)
    cube_diff = 100. * (cube_diff / iris.analysis.maths.abs(cube)) / cube_tas_diff

    return cube_diff


def check_cube(cube, filename):
    """Check properties of cube."""
    if cube.ndim != 1:
        raise ValueError(
            f"Expected 1D data in file '{filename}', got {cube.ndim:d} cube")
    try:
        cube.coord('dataset')
    except iris.exceptions.CoordinateNotFoundError:
        raise iris.exceptions.CoordinateNotFoundError(
            f"Necessary coordinate 'dataset' not found in file '{filename}'")


def create_data_frame(input_data, cfg):
    """Create data frame."""
    data_frame = pd.DataFrame(columns=['Variable', 'Group', 'Dataset', 'Data'])

    ifile = 0

    all_vars = group_metadata(input_data, 'short_name')
    groups = group_metadata(input_data, 'variable_group', sort='dataset')

    for var in all_vars:
        if var != 'tas':
            logger.info("Processing variable %s", var)

            for group_names in cfg['group_by']:
                logger.info("Processing group %s of variable %s", group_names[0], var)

                for dataset in groups[var + "_" + group_names[0]]:
                    dataset_name = dataset['dataset']
                    #logger.info("Loop dataset %s", dataset_name)

                    if dataset_name not in cfg['exclude_datasets']:
                        cube_diff = compute_diff_temp(input_data, group_names, var, dataset)

                        #if var in ['clivi', 'lwp']:
                        #    cube_diff = 1000 * cube_diff

                        group_name = group_names[0].split('_')[1] + " ECS"

                        data_frame.loc[ifile] = [var, group_name, dataset_name, cube_diff.data]
                        #data_frame.loc[dataset_name] = [var, group_name, cube_diff.data]
                        ifile = ifile + 1


    data_frame['Data'] = data_frame['Data'].astype(str).astype(float)

    # Sort and add numerical index for labels
    #data_frame.index.names = ['variable', 'group', 'dataset']
    #data_frame = data_frame.sort_index()
    #data_frame = _add_numerical_index(data_frame, exclude_datasets)

    return data_frame


def plot_boxplot(data_frame, cfg):

    sns.set_style('darkgrid')
    sns.set(font_scale=2)
    sns.boxplot(data=data_frame, x='Variable', y='Data', hue='Group', palette=PALETTE)
    plt.ylabel('Relative change (%/K)')
    #plt.ylabel('Relative change (%)')
    #plt.ylabel('Absolute change')
    #if data_frame['Variable'] == 'clt':
    #    plt.ylabel('Absolute change (%)')
    #if data_frame['Variable'] in ['clivi', 'lwp']:
    #    plt.ylabel('Absolute change (kg m-3)')
    #if data_frame['Variable'] == 'netcre':
    #    plt.ylabel('Absolute change (W m-2)')
    plt.ylim([-15., 40.])
    plt.title(cfg['title'])

    # Save plot
    plot_path = get_plot_filename('boxplot'+ '_' + cfg['filename_attach'], cfg)
    plt.savefig(plot_path)
    logger.info("Wrote %s", plot_path)
    plt.close()


def main(cfg):
    """Run diagnostic."""
    cfg.setdefault('exclude_datasets', ['MultiModelMean', 'MultiModelP5', 'MultiModelP95'])
    cfg.setdefault('title', 'Test')

    plt.figure(constrained_layout=True, figsize=(12, 8))

    # Get input data
    input_data = list(cfg['input_data'].values())
    groups = group_metadata(input_data, 'variable_group', sort='dataset')

    # Create data frame
    data_frame = create_data_frame(input_data, cfg)

    # Create plot
    plot_boxplot(data_frame, cfg)

    # Save file
    #basename = '-'.join(data_frame.index.levels[0]) + '_'
    #basename += '-'.join(data_frame.columns)
    basename = "boxplot_region_" + cfg['filename_attach']
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    data_frame.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Data:\n%s", data_frame)

#    # Provenance
#    write_provenance(cfg, csv_path, data_frame, input_files)

#    cfg = deepcopy(cfg)
#    cfg.setdefault('title_key', 'dataset')
#    cfg.setdefault('filename_attach', 'base')
#    logger.info("Using key '%s' to create titles for datasets",
#                cfg['title_key'])
#
#    input_data = list(cfg['input_data'].values())
#    all_vars = group_metadata(input_data, 'short_name')
#
#    plt.figure(constrained_layout=True, figsize=(12, 8))
#
#    for var in all_vars:
#
#        if var != 'tas':
#
#            logger.info("Processing variable %s", var)
#
#            groups = group_metadata(all_vars[var], 'variable_group', sort='dataset')
#
#            for group_name in groups:
#                if 'hist' in group_name:
#                    if 'tas_' not in group_name:
#                        logger.info("Processing group %s", group_name)
#
#                        dataset_names = []
#                        cubes = {}
#
#                        for dataset in groups[group_name]:
#
#                            dataset_name = dataset['dataset']
#                            var = dataset['short_name']
#
#                            if dataset_name != 'MultiModelMean':
#
#                                logger.info("Processing dataset %s", dataset_name)
#
#                                input_file = dataset['filename']
#                                cube = read_data(input_file)
#
#                                #cube = cube.collapsed('longitude', iris.analysis.MEAN)
#
#                                if cfg['plot_each_model']:
#                                    plot_model(cube, dataset, plot_type, cfg)
#
#                                cubes[dataset_name] = cube
#
#                        cube_mmm = _get_multi_model_mean(cubes, var)
#
#                        plot_diagnostic(cube, group_name, plot_type, cfg)
#
#                        title = group_name.split('_',1)[1]
#                        plt.title(title, fontsize=10)
#
#                        provenance_record = get_provenance_record(
#                            dataset, ancestor_files=cfg['input_files'])
#                        basename = 'mean_' + group_name + '_' + cfg['filename_attach']
#                        # Save the data used for the plot
#                        save_data(basename, provenance_record, cfg, cube_mmm)


#            for group_name in cfg['group_by']:
#
#                logger.info("Processing group %s of variable %s", group_name[0], var)
#
#                dataset_names = []
#                cubes = {}
#
#                for dataset in groups[var + "_" + group_name[0]]:
#                    dataset_name = dataset['dataset']
#                    var = dataset['short_name']
#
#                    if dataset_name not in ['MultiModelMean', 'MultiModelP5', 'MultiModelP95']:
#                        logger.info("Loop dataset %s", dataset_name)
#                        dataset_names.append(dataset_name)
#
#                        var_data = select_metadata(input_data,
#                                                   short_name=var,
#                                                   dataset=dataset_name,
#                                                   variable_group=var+"_"+group_name[0]) 
#                        if not var_data:
#                            raise ValueError(
#                                f"No '{var}' data for '{dataset_name}' in '{group_name[0]}' available")
#
#                        input_file = var_data[0]['filename']
#
#                        cube = read_data(input_file)
#                        
#                        cubes[dataset_name] = cube

#                cube_mmm = _get_multi_model_mean(cubes_diff, var)
#
#                plot_diagnostic_diff(cube_mmm, group_name[0], plot_type, cfg)
#
#                title = group_name[1] + " - " + group_name[0]
#                plt.title(title, fontsize=9)
#
#                provenance_record = get_provenance_record(
#                   dataset, ancestor_files=cfg['input_files'])
#                #basename = 'diff_' + var + '_' + group_name[0] + '_' + cfg['filename_attach']
#                basename = 'diff_' + var + '_' + group_name[0] + '_' + group_name[1] + '_' + cfg['filename_attach']
#                # Save the data used for the plot
#                save_data(basename, provenance_record, cfg, cube_mmm)
#
#            plt.suptitle(var)
#
#            provenance_record = get_provenance_record(
#                dataset, ancestor_files=cfg['input_files'])
#
#            basename = 'diff_' + var + '_' + cfg['filename_attach']
#            # Save the data used for the plot
#            save_data(basename, provenance_record, cfg, cube_mmm)
#
#            # And save the plot
#            save_figure(basename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
