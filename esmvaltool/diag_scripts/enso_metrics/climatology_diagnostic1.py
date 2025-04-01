"""diagnostic script to plot ENSO metrics

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
    # input data is 2 - model and obs
    plt.clf()
    figure = plt.figure(figsize=(10, 6), dpi=300)
    var_units = {'tos': 'degC', 'pr': 'mm/day', 'tauu': '1e-3 N/m2'}

    for dataset in input_data:
        # Load the data
        fp, sn, dt, proj = (dataset['filename'], dataset['short_name'],
                            dataset['dataset'], dataset['project'])

        logger.info(f"dataset: {dt} - {dataset['long_name']}")

        cube = iris.load_cube(fp)
        # convert units for different variables
        cube = convert_units(cube, units=var_units[sn])

        title = f"Mean {dataset['long_name']}"
        if len(cube.coords('month_number')) == 1:
            cube = sea_cycle_month_stdev(cube, dataset['preprocessor'])
            # plot title
            title = f"{dataset['long_name']} seasonal cycle"

        if proj == 'CMIP6':  # group by models/ for each model with obs
            qplt.plot(cube, label=dt)
            model_data = cube.data
            filename = [dt, dataset['variable_group']]
        else:
            qplt.plot(cube, label=f'ref: {dt}', color='black')
            obs_data = cube.data

    rmse = np.sqrt(np.mean((obs_data - model_data) ** 2))
    metricfile = get_diagnostic_filename('matrix', cfg, extension='csv')
    with open(metricfile, 'a+') as f:
        f.write(f"{filename[0]},{filename[1]},{rmse}\n")

    plt.title(title)
    plt.legend()
    plt.grid(linestyle='--')
    plt.ylabel(f"{sn.upper()} ({cube.units})")

    plt.text(0.5, 0.95, f"RMSE: {rmse:.2f} {cube.units}", fontsize=12,
             ha='center', transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    if dataset['preprocessor'].startswith('ITCZ'):
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lat))
    else:
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lon))

    return figure, filename


def sea_cycle_month_stdev(cube, preproc):

    cube.coord('month_number').guess_bounds()
    cube = cube.collapsed('month_number', iris.analysis.STD_DEV)
    # ITCZ or zonal
    if preproc.startswith('ITCZ'):
        cube = zonal_statistics(cube, 'mean')
    else:
        cube = meridional_statistics(cube, 'mean')

    return cube


def format_lat(x, pos):
    if x < 0:
        return f'{int(abs(x))}°S'
    elif x > 0:
        return f'{int(x)}°N'
    else:
        return '0°'


def format_lon(x, pos):
    if x > 180:
        return f'{int(360 - x)}°W'
    elif x == 180:
        return f'{int(x)}°'
    else:
        return f'{int(x)}°E'


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
    for grp in variable_groups:
        msg = "{} : {}, {}".format(grp, len(variable_groups[grp]),
                                   pformat(variable_groups[grp]))
        logger.info(msg)
        obs_data = variable_groups[grp][-1]

        for metadata in variable_groups[grp]:
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
