"""
Port to Version 2 with implementation of v2-specific changes.

Uses: ESMValTool v2, Python3.x
Valeriu Predoi, UREAD, July 2018

Porting replicates the functionality to minimum errors.

Original Description from Version 1 Diagnostic:
;;###########################################################################
;; AutoAssess_radiation_rms.py
;;###########################################################################
;; Description
;;    This script is the RMS error metric script of
;;    AutoAssess radiation
;; ###########################################################################

This diagnostic uses CMIP5 data; to switch to CMIP6 change _CMIP_TYPE.
"""

import os
import logging
import iris
from esmvaltool.diag_scripts.autoassess._rms_radiation import (start, end,
                                                               calc_all)
from esmvaltool.diag_scripts.autoassess._valmod_radiation import (
    perform_equation)
from esmvaltool.diag_scripts.shared import (
    group_metadata, run_diagnostic, get_control_exper_obs, apply_supermeans)

logger = logging.getLogger(os.path.basename(__file__))

_CMIP_TYPE = 'CMIP5'


def apply_rms(data_1, data_2, cfg, component_dict, var_name):
    """Compute RMS for any data1-2 combination."""
    data_names = [model['dataset'] for model in component_dict.values()]
    plot_title = var_name + ': ' + data_names[0] + ' vs ' + data_names[1]
    rms_list = start(data_names[0], data_names[1])
    analysis_type = cfg['analysis_type']
    landsea_mask_file = os.path.join(
        os.path.dirname(__file__), 'autoassess_source', cfg['landsea_mask'])
    landsea_mask_cube = iris.load_cube(landsea_mask_file)
    data1_vs_data2 = perform_equation(data_1, data_2, analysis_type)

    # call to rms.calc_all() to compute rms; rms.end() to write results
    calc_all(rms_list, data1_vs_data2, landsea_mask_cube, plot_title)
    end(rms_list, cfg['work_dir'])


def do_preamble(cfg):
    """Execute some preamble functionality."""
    # get data
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(
        input_data, 'short_name', sort='dataset')

    return input_data, grouped_input_data


def main(cfg):
    """Execute the radiation rms diag."""
    logger.setLevel(cfg['log_level'].upper())
    input_data, grouped_input_data = do_preamble(cfg)

    # select variables and their corresponding
    # obs files
    for short_name in grouped_input_data:
        logger.info("Processing variable %s", short_name)

        # control, experiment and obs's
        ctrl, exper, obslist = get_control_exper_obs(short_name, input_data,
                                                     cfg, _CMIP_TYPE)

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
