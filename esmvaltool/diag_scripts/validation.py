"""
Validation Diagnostic

This diagnostic uses two datasets (control and experiment),
applies operations on their data, and plots one against the other.
It can optionally use a number of OBS, OBS4MIPS datasets.
"""

import os
import logging

import numpy as np

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata)
from esmvaltool.preprocessor._area_pp import area_slice

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt      # noqa
import iris                          # noqa
import iris.analysis.maths as imath  # noqa
import iris.quickplot as qplt        # noqa


logger = logging.getLogger(os.path.basename(__file__))


def plot_contour(cube, plt_title, file_name):
    """Plot a contour with iris.quickplot (qplot)"""
    qplt.contourf(cube, cmap='RdYlBu_r', bbox_inches='tight')
    plt.title(plt_title)
    plt.gca().coastlines()
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


def plot_latlon_cubes(cube_1, cube_2, cfg, data_names, obs_name=None):
    """
    Plot different things

    cube_1: first cube (dataset: dat1)
    cube_2: second cube (dataset: dat2)
    cfg: configuration dictionary
    data_names: var + '_' + dat1 + '_vs_' + dat2
    """
    plot_name = cfg['analysis_type'] + '_' + data_names + '.png'
    plot_title = cfg['analysis_type'] + ': ' + data_names
    cubes = [cube_1, cube_2]

    # plot difference: cube_1 - cube_2; use iris subtract=idiff
    plot_contour(imath.abs(cube_1) - imath.abs(cube_2),
                 'Difference ' + plot_title,
                 os.path.join(cfg['plot_dir'], 'Difference_' + plot_name))

    # plot each cube
    var = data_names.split('_')[0]
    if not obs_name:
        cube_names = [data_names.split('_')[1], data_names.split('_')[3]]
        for cube, cube_name in zip(cubes, cube_names):
            plot_contour(cube,
                         cube_name + ' ' + cfg['analysis_type'] + ' ' + var,
                         os.path.join(cfg['plot_dir'],
                                      cube_name + '_' + var + '.png'))
    else:
        # obs is always cube_2
        plot_contour(cube_2,
                     obs_name + ' ' + cfg['analysis_type'] + ' ' + var,
                     os.path.join(cfg['plot_dir'],
                                  obs_name + '_' + var + '.png'))


def plot_lon_cubes(cube_1, cube_2, cfg, data_names):
    """Plot cubes vs latitude"""
    var = data_names.split('_')[0]
    cube_names = [data_names.split('_')[1], data_names.split('_')[3]]
    lat_points = cube_1.coord('latitude').points
    plt.plot(lat_points, cube_1.data, label=cube_names[0])
    plt.plot(lat_points, cube_2.data, label=cube_names[1])
    plt.title('Zonal mean on longitude: ' + data_names)
    plt.xlabel('Latitude (deg)')
    plt.ylabel(var)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(cfg['plot_dir'],
                             'Zonal_Means_LON_' + data_names + '.png'))
    plt.close()


def plot_lat_cubes(cube_1, cube_2, cfg, data_names):
    """Plot cubes vs longitude"""
    var = data_names.split('_')[0]
    cube_names = [data_names.split('_')[1], data_names.split('_')[3]]
    lat_points = cube_1.coord('longitude').points
    plt.plot(lat_points, cube_1.data, label=cube_names[0])
    plt.plot(lat_points, cube_2.data, label=cube_names[1])
    plt.title('Zonal mean on latitude: ' + data_names)
    plt.xlabel('Longitude (deg)')
    plt.ylabel(var)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(cfg['plot_dir'],
                             'Zonal_Means_LAT_' + data_names + '.png'))
    plt.close()


def apply_supermeans(ctrl, exper, obs_list):
    """Apply supermeans on data components"""
    ctrl_file = ctrl['filename']
    exper_file = exper['filename']
    logger.info("Loading %s", ctrl_file)
    logger.info("Loading %s", exper_file)
    ctrl_cube = iris.load_cube(ctrl_file)
    exper_cube = iris.load_cube(exper_file)
    ctrl_cube = ctrl_cube.collapsed('time', iris.analysis.MEAN)
    logger.debug("Time-averaged control %s", ctrl_cube)
    exper_cube = exper_cube.collapsed('time', iris.analysis.MEAN)
    logger.debug("Time-averaged experiment %s", exper_cube)
    if obs_list:
        obs_cube_list = []
        for obs in obs_list:
            obs_file = obs['filename']
            logger.info("Loading %s", obs_file)
            obs_cube = iris.load_cube(obs_file)
            obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)
            logger.debug("Time-averaged obs %s", obs_cube)
            obs_cube_list.append(obs_cube)
    else:
        obs_cube_list = None

    return ctrl_cube, exper_cube, obs_cube_list


def apply_analysis(data_set, cfg):
    """Configure stats for any dataset"""
    # see what analysis needs performing
    analysis_type = cfg['analysis_type']

    # if subset on LAT-LON
    if 'lat_lon_slice' in cfg:
        start_longitude = cfg['lat_lon_slice']['start_longitude']
        end_longitude = cfg['lat_lon_slice']['end_longitude']
        start_latitude = cfg['lat_lon_slice']['start_latitude']
        end_latitude = cfg['lat_lon_slice']['end_latitude']
        data_set = area_slice(data_set, start_longitude, end_longitude,
                              start_latitude, end_latitude)

    # if apply mask
    if '2d_mask' in cfg:
        mask_file = os.path.join(cfg['2d_mask'])
        mask_cube = iris.load_cube(mask_file)
        if 'mask_threshold' in cfg:
            thr = cfg['mask_threshold']
            data_set.data = np.ma.masked_array(data_set.data,
                                               mask=(mask_cube.data > thr))
        else:
            logger.warning('Could not find masking threshold')
            logger.warning('Please specify it if needed')
            logger.warning('Masking on 0-values = True (masked value)')
            data_set.data = np.ma.masked_array(data_set.data,
                                               mask=(mask_cube.data == 0))

    # if zonal mean on LON
    if analysis_type == 'zonal_mean_longitude':
        data_set = data_set.collapsed('longitude', iris.analysis.MEAN)

    # if zonal mean on LAT
    if analysis_type == 'zonal_mean_latitude':
        data_set = data_set.collapsed('latitude', iris.analysis.MEAN)

    # if vertical mean
    elif analysis_type == 'vertical_mean':
        data_set = data_set.collapsed('pressure', iris.analysis.MEAN)

    return data_set


def do_preamble(cfg):
    """Execute some preamble functionality"""
    # prepare output dirs
    if not os.path.exists(cfg['plot_dir']):
        os.makedirs(cfg['plot_dir'])
    if not os.path.exists(cfg['work_dir']):
        os.makedirs(cfg['work_dir'])

    # get data
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(
        input_data, 'short_name', sort='dataset')

    return input_data, grouped_input_data


def get_all_datasets(short_name, input_data, cfg):
    """Get control, exper and obs datasets"""
    dataset_selection = select_metadata(
        input_data, short_name=short_name, project='CMIP5')

    # get the obs datasets
    if 'observational_datasets' in cfg.keys():
        obs_selection = [
            select_metadata(input_data, short_name=short_name,
                            dataset=obs_dataset)[0]
            for obs_dataset in cfg['observational_datasets']
        ]
    else:
        obs_selection = []

    # determine CONTROL and EXPERIMENT datasets
    for model in dataset_selection:
        if model['dataset'] == cfg['control_model']:
            logger.info("Control dataset %s", model['dataset'])
            control = model
        elif model['dataset'] == cfg['exper_model']:
            logger.info("Experiment dataset %s", model['dataset'])
            experiment = model

    if obs_selection:
        logger.info("Observations dataset(s) %s",
                    [obs['dataset'] for obs in obs_selection])

    return control, experiment, obs_selection


def main(cfg):
    """Execute the radiation rms diag"""
    logger.setLevel(cfg['log_level'].upper())
    input_data, grouped_input_data = do_preamble(cfg)

    # select variables and their corresponding
    # obs files
    for short_name in grouped_input_data:
        logger.info("Processing variable %s", short_name)

        # control, experiment and obs's and the names
        ctrl, exper, obs = get_all_datasets(short_name, input_data, cfg)
        ctrl_name = ctrl['dataset']
        exper_name = exper['dataset']

        # apply the supermeans
        ctrl, exper, obs_list = apply_supermeans(ctrl, exper, obs)

        # apply the desired analysis
        ctrl = apply_analysis(ctrl, cfg)
        exper = apply_analysis(exper, cfg)

        # plot
        plot_key = short_name + '_' + ctrl_name + '_vs_' + exper_name
        if cfg['analysis_type'] == 'lat_lon':
            plot_latlon_cubes(ctrl, exper, cfg, plot_key)
        elif cfg['analysis_type'] == 'zonal_mean_longitude':
            plot_lon_cubes(ctrl, exper, cfg, plot_key)
        elif cfg['analysis_type'] == 'zonal_mean_latitude':
            plot_lat_cubes(ctrl, exper, cfg, plot_key)

        # apply desired analysis on obs's
        if obs_list:
            for obs_i, obsfile in zip(obs_list, obs):
                obs_analyzed = apply_analysis(obs_i, cfg)
                obs_name = obsfile['dataset']
                plot_key = short_name + '_' + ctrl_name + '_vs_' + obs_name
                if cfg['analysis_type'] == 'lat_lon':
                    plot_latlon_cubes(ctrl, obs_analyzed, cfg,
                                      plot_key, obs_name=obs_name)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
