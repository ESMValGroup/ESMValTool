"""
Port to Version 2 with implementation of v2-specific changes

Uses: ESMValTool v2, Python3.x
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
import autoassess_source.rms as rms
import autoassess_source.valmod_radiation as vm
from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


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


def apply_rms(data_1, data_2, cfg, component_dict, var_name):
    """Compute RMS for any data1-2 combination"""
    data_names = [model['dataset'] for model in component_dict.values()]
    plot_title = var_name + ': ' + data_names[0] + ' vs ' + data_names[1]
    rms_list = rms.start(data_names[0], data_names[1])
    analysis_type = cfg['analysis_type']
    landsea_mask_file = os.path.join(
        os.path.dirname(__file__), 'autoassess_source', cfg['landsea_mask'])
    landsea_mask_cube = iris.load_cube(landsea_mask_file)
    data1_vs_data2 = vm.perform_equation(data_1, data_2, analysis_type)

    # call to rms.calc_all() to compute rms; rms.end() to write results
    rms.calc_all(rms_list, data1_vs_data2, landsea_mask_cube, plot_title)
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
