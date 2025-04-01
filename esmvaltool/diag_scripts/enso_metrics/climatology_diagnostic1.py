"""diagnostic script to plot basic climatology metrics

"""

import logging
import os
from pprint import pformat

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    convert_units,
    meridional_statistics,
    zonal_statistics,
)

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(input_data, cfg):
    """Create plots for pair of input data."""
    plt.clf()
    filename = []
    obs_data, model_data = None
    figure = plt.figure(figsize=(10, 6), dpi=300)
    var_units = {'tos': 'degC', 'pr': 'mm/day', 'tauu': '1e-3 N/m2'}

    for dataset in input_data:
        # Load the data
        filep, sname = dataset['filename'], dataset['short_name']
        dtname, proj = dataset['dataset'], dataset['project']

        logger.info("dataset: %s - %s", dtname, dataset['long_name'])

        cube = iris.load_cube(filep)
        # convert units for different variables
        cube = convert_units(cube, units=var_units[sname])

        title = f"Mean {dataset['long_name']}"
        if len(cube.coords('month_number')) == 1:
            cube = sea_cycle_month_stdev(cube, dataset['preprocessor'])
            # plot title
            title = f"{dataset['long_name']} seasonal cycle"

        if proj == 'CMIP6':  # group by models/ for each model with obs
            qplt.plot(cube, label=dtname)
            model_data = cube.data
            filename = [dtname, dataset['variable_group']]
        else:
            qplt.plot(cube, label=f'ref: {dtname}', color='black')
            obs_data = cube.data

    rmse = np.sqrt(np.mean((obs_data - model_data) ** 2))
    metricfile = get_diagnostic_filename('matrix', cfg, extension='csv')
    with open(metricfile, 'a+',  encoding='utf-8') as fileo:
        fileo.write(f"{filename[0]},{filename[1]},{rmse}\n")

    plt.title(title)
    plt.legend()
    plt.grid(linestyle='--')
    plt.ylabel(f"{sname.upper()} ({cube.units})")

    plt.text(0.5, 0.95, f"RMSE: {rmse:.2f} {cube.units}", fontsize=12,
             ha='center', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    if dataset['preprocessor'].startswith('ITCZ'):
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lat))
    else:
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lon))

    return figure, filename


def sea_cycle_month_stdev(cube, preproc):
    """Process standard deviation for seasonal cycle."""
    cube.coord('month_number').guess_bounds()
    cube = cube.collapsed('month_number', iris.analysis.STD_DEV)
    # ITCZ or zonal
    if preproc.startswith('ITCZ'):
        cube = zonal_statistics(cube, 'mean')
    else:
        cube = meridional_statistics(cube, 'mean')

    return cube


def format_lat(x_val, pos):
    """Format latitude for plot axis."""
    if x_val < 0:
        return f'{int(abs(x_val))}°S'
    if x_val > 0:
        return f'{int(x_val)}°N'
    else:
        return '0°'


def format_lon(val, pos):
    """Format longitude for plot axis."""
    if val > 180:
        return f'{int(360 - val)}°W'
    if val == 180:
        return f'{int(val)}°'
    else:
        return f'{int(val)}°E'


def main(cfg):
    """Compute sea ice area for each input dataset."""
    provenance_record = {
        'caption': "ENSO metrics",
        'authors': [
            'chun_felicity',
        ],
        'references': [''],
        'ancestors': list(cfg['input_data'].keys()),
    }
    input_data = cfg['input_data'].values()

    # group by variables
    variable_groups = group_metadata(input_data, 'variable_group',
                                     sort='project')
    # for each select obs and iterate others, obs last
    for grp, var_attr in variable_groups.items():

        logger.info("%s : %d, %s", grp, len(var_attr), pformat(var_attr))
        obs_data = var_attr[-1]

        for metadata in var_attr:
            logger.info('iterate though datasets\n %s', pformat(metadata))
            pairs = [obs_data]
            if metadata['project'] == 'CMIP6':
                pairs.append(metadata)
                fig, filename = plot_level1(pairs, cfg)

                save_figure('_'.join(filename), provenance_record,
                            cfg, figure=fig, dpi=300)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
