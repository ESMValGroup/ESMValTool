"""
Port to Version 2 with implementation of v2-specific changes

Valeriu Predoi, UREAD, July 2018

Porting replicates the functionality to minimum errors.

Original Description from Version 1 Diagnostic:
;;###########################################################################
;; AutoAssess_radiation_rms.py
;; Author: Yoko Tsushima (Met Office, UK)
;; CMUG project
;;###########################################################################
;; Description
;;    This script is the RMS error metric script of
;;    AutoAssess radiation
;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;; [general]
;; plot_clouds:           Switch to plot cloud diags
;; plot_fluxes:           Switch to plot fluxes diags
;; plot_radiation:        Switch to plot radiation diags
;; plot_scatter:          Switch to plot scatter diags
;; plot_background_grid:  Switch to plot background grid
;;
;; plot_total_cover:     Sub map switch for plotting
;; plot_liquid_path:     Sub map switch for plotting
;; plot_ice_path:        Sub map switch for plotting
;; plot_optical_depth:   Sub map switch for plotting
;; plot_flux_maps:       Sub map switch for plotting
;; plot_radiation_maps:  Sub map switch for plotting
;;
;; plot_lat_averages:        General switch valid for all plots
;; plot_lon_averages:        General switch valid for all plots
;; plot_monthly_averages:    General switch valid for all plots
;; plot_sub_areas:           General switch valid for all plots
;;
;; mask_unwanted_values: Switch for masking values
;; mask_limit_low:       Mask values below this value
;; mask_limit_high:      Mask values above this value
;;
;; [SouthernHemisphere]
;; areas:           Control what plots (regions/seasons) are generated
;; sub_areas:       Control what plots (regions/seasons) are generated
;; scatter_areas:   Control what plots (regions/seasons) are generated
;; seasons:         Control what plots (regions/seasons) are generated
;;
;; [SouthernHemisphere_default]
;; # Latitudes [-90, 90 degrees]: 10S = -10, 10N = 10;
;; longtitudes [0, 360 degrees]
;; lat_min: Default area
;; lat_max: Default area
;; lon_min: Default area
;; lon_max: Default area
;;
;; stride:      Color difference interval
;; maxshades:   Color difference interval
;;
;; contour_limits_clt:   (min, max, diff, [dev_min, dev_max]
;; contour_limits_clivi: (min, max, diff, [dev_min, dev_max]
;; contour_limits_clwvi: (min, max, diff, [dev_min, dev_max]
;; contour_limits_hfls:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_hfss:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rlut:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rsut:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rlds:  (min, max, diff, [dev_min, dev_max]
;; contour_limits_rsds:  (min, max, diff, [dev_min, dev_max]
;;
;; colourmap_clouds: Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_model:  Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_diff:   Python matplotlib colormaps, '_r' reversed colormap
;; colourmap_dev:    Python matplotlib colormaps, '_r' reversed colormap
;;
;;
;; # Define area specifications for sub_areas
;; [SouthernHemisphere_northern]
;; lat_min: Define areas for northern parts
;; lat_max: Define areas for northern parts
;; lon_min: Define areas for northern parts
;; lon_max: Define areas for northern parts
;;
;; [SouthernHemisphere_southern]
;; lat_min: Define areas for southern parts
;; lat_max: Define areas for southern parts
;; lon_min: Define areas for southern parts
;; lon_max: Define areas for southern parts
;;
;; # Months to use for each season - 1 is January and so forth.
;; [SouthernHemisphere_season_DJF]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_MAM]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_JJA]
;; season_months: Months to use, 1 i January and so forth
;;
;; [SouthernHemisphere_season_SON]
;; season_months: Months to use, 1 i January and so forth
;;
;; # Define configuration for cloud vs radiation scatter plots
;; [SouthernHemisphere_scatter_default]
;; lat_min: Configuration for cloud vs radiation scatter plots
;; lat_max: Configuration for cloud vs radiation scatter plots
;; lon_min: Configuration for cloud vs radiation scatter plots
;; lon_max: Configuration for cloud vs radiation scatter plots
;; points:  Configuration for cloud vs radiation scatter plots
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;    long_name: Name displayed in plot
;;    units:     Units
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;   landsea.nc from NCLA for land/sea mask and resolution.
;;   This does not have coord system info.
;;
;; Modification history
;;    20180712- autoassess_radiation_rms: porting to v2
;;    20170323-_AutoAssess_radiation_rms: Test finished.
;;    20160819-_test_AutoAssess_radiation_rms: written based on calc_rms code.
;;
;; ###########################################################################
"""

import os
import logging
import iris
from esmvaltool.diag_scripts.autoassess.autoassess_source import rms
from esmvaltool.diag_scripts.autoassess.autoassess_source \
    import valmod_radiation as vm
from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def map_var_to_obs(cfg):
    """Get the preferred OBS files for each var"""
    if cfg['diagnostic_class'] == 'Amon_all':
        var_dict = {
            'pr': ['GPCP-SG'],  # underived
            'rlut': ['CERES-EBAF'],  # underived
            'rsut': ['CERES-EBAF'],  # underived
            'rsutcs': ['CERES-EBAF'],  # underived
            'rlutcs': ['JRA-55'],  # underived
            'prw': ['SSMI-MERIS'],  # underived
            'rldscs': ['JRA-55'],  # underived
            'rtnt': ['CERES-EBAF'],  # derived
            'rsnt': ['CERES-EBAF'],  # derived
            'rsns': ['CERES-EBAF'],  # derived
            'rlns': ['CERES-EBAF'],  # derived
            'lwcre': ['CERES-EBAF'],  # derived
            'swcre': ['CERES-EBAF']  # derived
        }

    return var_dict


def apply_supermeans(ctrl, exper, obs_list):
    """Grab data and apply supermeans"""
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


def apply_rms(data_1, data_2, cfg, component_dict, var_name):
    """Compute RMS for any data1-2 combination"""
    data_names = [model['dataset'] for model in component_dict.values()]
    plot_title = var_name + ': ' + data_names[0] + ' vs ' + data_names[1]
    rms_list = rms.start(data_names[0], data_names[1])
    filename = var_name + '_' + data_names[0] + '_vs_' + data_names[1]
    analysis_type = cfg['analysis_type']
    landsea_mask_file = os.path.join(
        os.path.dirname(__file__), 'autoassess_source', cfg['landsea_mask'])
    landsea_mask_cube = iris.load_cube(landsea_mask_file)
    data1_vs_data2 = vm.perform_equation(data_1, data_2, analysis_type)

    # call to calc_all to compute rms
    rms.calc_all(rms_list, data1_vs_data2, landsea_mask_cube, plot_title,
                 filename)
    rms.end(rms_list, cfg['work_dir'])


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
    var_dict = map_var_to_obs(cfg)
    dataset_selection = select_metadata(
        input_data, short_name=short_name, project='CMIP5')
    obs_selection = [
        select_metadata(input_data, short_name=short_name,
                        dataset=var_to_obs)[0]
        for var_to_obs in var_dict[short_name]
    ]

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

        # control, experiment and obs's
        ctrl, exper, obslist = get_all_datasets(short_name, input_data, cfg)

        # apply the supermeans
        ctrl_sm, exper_sm, obs_sm_list = apply_supermeans(ctrl, exper, obslist)

        # assemble a dict that contains various params depending
        # on the data combinations for RMS computations
        # control-experiment
        data_component_dict = {'ct-ex': {'ctrl': ctrl, 'exper': exper}}
        logger.info("Computing CONTROL-EXPERIMENT RMS...")
        apply_rms(ctrl_sm, exper_sm, cfg, data_component_dict['ct-ex'],
                  short_name)
        if obs_sm_list:
            for obs, obsfile in zip(obs_sm_list, obslist):
                data_component_dict = {
                    'ct-obs': {
                        'ctrl': ctrl,
                        'obs': obsfile
                    },
                    'ex-obs': {
                        'exper': exper,
                        'obs': obsfile
                    }
                }

                # ctrl-obs
                logger.info("Computing CONTROL-OBS RMS...")
                apply_rms(ctrl_sm, obs, cfg, data_component_dict['ct-obs'],
                          short_name)
                # exper-obs
                logger.info("Computing EXPERIMENT-OBS RMS...")
                apply_rms(exper_sm, obs, cfg, data_component_dict['ex-obs'],
                          short_name)
        else:
            # only ctrl-exper
            data_component_dict = {'ct-ex': {'ctrl': ctrl, 'exper': exper}}
            logger.info("Computing CONTROL-EXPERIMENT RMS...")
            apply_rms(ctrl_sm, exper_sm, cfg, data_component_dict['ct-ex'],
                      short_name)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
