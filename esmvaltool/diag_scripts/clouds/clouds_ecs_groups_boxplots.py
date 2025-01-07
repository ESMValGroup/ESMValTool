"""Python diagnostic for plotting boxplots."""
import logging
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

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


def get_provenance_record(ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Relative change per degree of warming averaged over the"
               "chosen region.")

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'bock_lisa',
        ],
        'references': [
            'bock24acp',
        ],
        'ancestors': ancestor_files,
    }
    return record


def read_data(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    if cube.var_name == 'cli':
        cube.convert_units('g/kg')
    elif cube.var_name == 'clw':
        cube.convert_units('g/kg')

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
                                 variable_group=var + "_" + group[1])
    if not var_data_2:
        raise ValueError(
            f"No '{var}' data for '{dataset_name}' in '{group[1]}' available")

    input_file_2 = var_data_2[0]['filename']

    tas_data_1 = select_metadata(input_data,
                                 short_name='tas',
                                 dataset=dataset_name,
                                 variable_group='tas_' + group[0])
    tas_data_2 = select_metadata(input_data,
                                 short_name='tas',
                                 dataset=dataset_name,
                                 variable_group='tas_' + group[1])
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

    cube_diff = (100. * (cube_diff / iris.analysis.maths.abs(cube)) /
                 cube_tas_diff)

    return cube_diff


def create_data_frame(input_data, cfg):
    """Create data frame."""
    data_frame = pd.DataFrame(columns=['Variable', 'Group', 'Dataset', 'Data'])

    ifile = 0

    all_vars = group_metadata(input_data, 'short_name')
    groups = group_metadata(input_data, 'variable_group', sort='dataset')

    for var in all_vars:
        if var != 'tas':
            logger.info("Processing variable %s", var)

            if var == 'clivi':
                varname = 'iwp'
            else:
                varname = var

            for group_names in cfg['group_by']:
                logger.info("Processing group %s of variable %s",
                            group_names[0], var)

                for dataset in groups[var + "_" + group_names[0]]:
                    dataset_name = dataset['dataset']

                    if dataset_name not in cfg['exclude_datasets']:
                        cube_diff = compute_diff_temp(input_data, group_names,
                                                      var, dataset)

                        group_name = group_names[0].split('_')[1] + " ECS"

                        data_frame.loc[ifile] = [
                            varname, group_name, dataset_name, cube_diff.data
                        ]
                        ifile = ifile + 1

    data_frame['Data'] = data_frame['Data'].astype(str).astype(float)

    return data_frame


def plot_boxplot(data_frame, input_data, cfg):
    """Create boxplot."""
    sns.set_style('darkgrid')
    sns.set(font_scale=2)
    sns.boxplot(data=data_frame,
                x='Variable',
                y='Data',
                hue='Group',
                palette=PALETTE)
    plt.ylabel('Relative change (%/K)')
    if 'y_range' in cfg:
        plt.ylim(cfg.get('y_range'))
    plt.title(cfg['title'])

    provenance_record = get_provenance_record(
        ancestor_files=[d['filename'] for d in input_data])

    # Save plot
    plot_path = get_plot_filename('boxplot' + '_' + cfg['filename_attach'],
                                  cfg)
    plt.savefig(plot_path)
    logger.info("Wrote %s", plot_path)
    plt.close()

    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)


def main(cfg):
    """Run diagnostic."""
    cfg.setdefault('exclude_datasets',
                   ['MultiModelMean', 'MultiModelP5', 'MultiModelP95'])
    cfg.setdefault('title', 'Test')

    plt.figure(constrained_layout=True, figsize=(12, 8))

    # Get input data
    input_data = list(cfg['input_data'].values())

    # Create data frame
    data_frame = create_data_frame(input_data, cfg)

    # Create plot
    plot_boxplot(data_frame, input_data, cfg)

    # Save file
    basename = "boxplot_region_" + cfg['filename_attach']
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    data_frame.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
