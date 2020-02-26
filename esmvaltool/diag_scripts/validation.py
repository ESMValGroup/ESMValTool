"""
Validation Diagnostic

This diagnostic uses two datasets (control and experiment),
applies operations on their data, and plots one against the other.
It can optionally use a number of OBS, OBS4MIPS datasets.
"""

import logging
import os

import iris
import iris.analysis.maths as imath
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (apply_supermeans,
                                            get_control_exper_obs,
                                            group_metadata, run_diagnostic)
from esmvalcore.preprocessor import extract_region, extract_season

logger = logging.getLogger(os.path.basename(__file__))


def plot_contour(cube, plt_title, file_name):
    """Plot a contour with iris.quickplot (qplot)"""
    if len(cube.shape) == 2:
        qplt.contourf(cube, cmap='RdYlBu_r', bbox_inches='tight')
    else:
        qplt.contourf(cube[0], cmap='RdYlBu_r', bbox_inches='tight')
    plt.title(plt_title)
    plt.gca().coastlines()
    plt.tight_layout()
    plt.savefig(file_name)
    plt.close()


def plot_latlon_cubes(cube_1, cube_2, cfg, data_names, obs_name=None):
    """
    Plot lat-lon vars for control, experiment, and obs

    Also plot Difference plots (control-exper, control-obs)
    cube_1: first cube (dataset: dat1)
    cube_2: second cube (dataset: dat2)
    cfg: configuration dictionary
    data_names: var + '_' + dat1 + '_vs_' + dat2
    """
    plot_name = cfg['analysis_type'] + '_' + data_names + '.png'
    plot_title = cfg['analysis_type'] + ': ' + data_names
    cubes = [cube_1, cube_2]

    # plot difference: cube_1 - cube_2; use numpy.ma.abs()
    diffed_cube = imath.subtract(cube_1, cube_2)
    plot_contour(diffed_cube, 'Difference ' + plot_title,
                 os.path.join(cfg['plot_dir'], 'Difference_' + plot_name))

    # plot each cube
    var = data_names.split('_')[0]
    if not obs_name:
        cube_names = [data_names.split('_')[1], data_names.split('_')[3]]
        for cube, cube_name in zip(cubes, cube_names):
            plot_contour(
                cube, cube_name + ' ' + cfg['analysis_type'] + ' ' + var,
                os.path.join(cfg['plot_dir'], cube_name + '_' + var + '.png'))
    else:
        # obs is always cube_2
        plot_contour(
            cube_2, obs_name + ' ' + cfg['analysis_type'] + ' ' + var,
            os.path.join(cfg['plot_dir'], obs_name + '_' + var + '.png'))


def plot_zonal_cubes(cube_1, cube_2, cfg, plot_data):
    """Plot cubes data vs latitude or longitude when zonal meaning"""
    # xcoordinate: latotude or longitude (str)
    data_names, xcoordinate, period = plot_data
    var = data_names.split('_')[0]
    cube_names = [data_names.split('_')[1], data_names.split('_')[3]]
    lat_points = cube_1.coord(xcoordinate).points
    plt.plot(lat_points, cube_1.data, label=cube_names[0])
    plt.plot(lat_points, cube_2.data, label=cube_names[1])
    if xcoordinate == 'latitude':
        plt.title(period + ' Zonal Mean for ' + var + ' ' + data_names)
    elif xcoordinate == 'longitude':
        plt.title(period + ' Meridional Mean for ' + var + ' ' + data_names)
    plt.xlabel(xcoordinate + ' (deg)')
    plt.ylabel(var)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    if xcoordinate == 'latitude':
        png_name = 'Zonal_Mean_' + xcoordinate + '_' + data_names + '.png'
    elif xcoordinate == 'longitude':
        png_name = 'Merid_Mean_' + xcoordinate + '_' + data_names + '.png'
    plt.savefig(os.path.join(cfg['plot_dir'], period, png_name))
    plt.close()


def apply_seasons(data_set_dict):
    """Extract seaons and apply a time mean per season"""
    data_file = data_set_dict['filename']
    logger.info("Loading %s for seasonal extraction", data_file)
    data_cube = iris.load_cube(data_file)
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    season_cubes = [extract_season(data_cube, season) for season in seasons]
    season_meaned_cubes = [
        season_cube.collapsed('time', iris.analysis.MEAN)
        for season_cube in season_cubes
    ]

    return season_meaned_cubes


def coordinate_collapse(data_set, cfg):
    """Perform coordinate-specific collapse and (if) area slicing and mask"""
    # see what analysis needs performing
    analysis_type = cfg['analysis_type']

    # if subset on LAT-LON
    if 'lat_lon_slice' in cfg:
        start_longitude = cfg['lat_lon_slice']['start_longitude']
        end_longitude = cfg['lat_lon_slice']['end_longitude']
        start_latitude = cfg['lat_lon_slice']['start_latitude']
        end_latitude = cfg['lat_lon_slice']['end_latitude']
        data_set = extract_region(data_set, start_longitude, end_longitude,
                                  start_latitude, end_latitude)

    # if apply mask
    if '2d_mask' in cfg:
        mask_file = os.path.join(cfg['2d_mask'])
        mask_cube = iris.load_cube(mask_file)
        if 'mask_threshold' in cfg:
            thr = cfg['mask_threshold']
            data_set.data = np.ma.masked_array(
                data_set.data, mask=(mask_cube.data > thr))
        else:
            logger.warning('Could not find masking threshold')
            logger.warning('Please specify it if needed')
            logger.warning('Masking on 0-values = True (masked value)')
            data_set.data = np.ma.masked_array(
                data_set.data, mask=(mask_cube.data == 0))

    # if zonal mean on LON
    if analysis_type == 'zonal_mean':
        data_set = data_set.collapsed('longitude', iris.analysis.MEAN)

    # if zonal mean on LAT
    if analysis_type == 'meridional_mean':
        data_set = data_set.collapsed('latitude', iris.analysis.MEAN)

    # if vertical mean
    elif analysis_type == 'vertical_mean':
        data_set = data_set.collapsed('pressure', iris.analysis.MEAN)

    return data_set


def do_preamble(cfg):
    """Execute some preamble functionality"""
    # prepare output dirs
    time_chunks = ['alltime', 'DJF', 'MAM', 'JJA', 'SON']
    time_plot_dirs = [
        os.path.join(cfg['plot_dir'], t_dir) for t_dir in time_chunks
    ]
    for time_plot_dir in time_plot_dirs:
        if not os.path.exists(time_plot_dir):
            os.makedirs(time_plot_dir)

    # get data
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(
        input_data, 'short_name', sort='dataset')

    return input_data, grouped_input_data


def plot_ctrl_exper(ctrl, exper, cfg, plot_key):
    """Call plotting functions and make plots depending on case"""
    if cfg['analysis_type'] == 'lat_lon':
        plot_latlon_cubes(ctrl, exper, cfg, plot_key)
    elif cfg['analysis_type'] == 'zonal_mean':
        plot_info = [plot_key, 'latitude', 'alltime']
        plot_zonal_cubes(ctrl, exper, cfg, plot_info)
    elif cfg['analysis_type'] == 'meridional_mean':
        plot_info = [plot_key, 'longitude', 'alltime']
        plot_zonal_cubes(ctrl, exper, cfg, plot_info)


def plot_ctrl_exper_seasons(ctrl_seasons, exper_seasons, cfg, plot_key):
    """Call plotting functions and make plots with seasons"""
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    if cfg['analysis_type'] == 'zonal_mean':
        for c_i, e_i, s_n in zip(ctrl_seasons, exper_seasons, seasons):
            plot_info = [plot_key, 'latitude', s_n]
            plot_zonal_cubes(c_i, e_i, cfg, plot_info)
    elif cfg['analysis_type'] == 'meridional_mean':
        for c_i, e_i, s_n in zip(ctrl_seasons, exper_seasons, seasons):
            plot_info = [plot_key, 'longitude', s_n]
            plot_zonal_cubes(c_i, e_i, cfg, plot_info)


def main(cfg):
    """Execute validation analysis and plotting"""
    logger.setLevel(cfg['log_level'].upper())
    input_data, grouped_input_data = do_preamble(cfg)

    # select variables and their corresponding obs files
    for short_name in grouped_input_data:
        logger.info("Processing variable %s", short_name)

        cmip_era = cfg["cmip_era"]
        # get the control, experiment and obs dicts
        ctrl, exper, obs = get_control_exper_obs(short_name, input_data,
                                                 cfg, cmip_era)
        # set a plot key holding info on var and data set names
        plot_key = short_name + '_' + ctrl['dataset'] \
            + '_vs_' + exper['dataset']

        # get seasons if needed then apply analysis
        if cfg['seasonal_analysis']:
            ctrl_seasons = apply_seasons(ctrl)
            exper_seasons = apply_seasons(exper)
            ctrl_seasons = [
                coordinate_collapse(cts, cfg) for cts in ctrl_seasons
            ]
            exper_seasons = [
                coordinate_collapse(exps, cfg) for exps in exper_seasons
            ]
            plot_ctrl_exper_seasons(ctrl_seasons, exper_seasons, cfg, plot_key)

        # apply the supermeans (MEAN on time), collapse a coord and plot
        ctrl, exper, obs_list = apply_supermeans(ctrl, exper, obs)
        ctrl = coordinate_collapse(ctrl, cfg)
        exper = coordinate_collapse(exper, cfg)
        plot_ctrl_exper(ctrl, exper, cfg, plot_key)

        # apply desired analysis on obs's
        if obs_list:
            for obs_i, obsfile in zip(obs_list, obs):
                obs_analyzed = coordinate_collapse(obs_i, cfg)
                obs_name = obsfile['dataset']
                plot_key = short_name + '_CONTROL_vs_' + obs_name
                if cfg['analysis_type'] == 'lat_lon':
                    plot_latlon_cubes(
                        ctrl, obs_analyzed, cfg, plot_key, obs_name=obs_name)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
