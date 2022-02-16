"""Load functions needed by diags with CONTROL and EXPERIMENT."""
import logging
import os
import sys

import iris
from esmvalcore.preprocessor import climate_statistics

from esmvaltool.diag_scripts.shared import select_metadata

logger = logging.getLogger(os.path.basename(__file__))


def _disentagle_iden_datasets(dat_selection):
    """
    Disentangle identical dataset names for CONTROL and EXPERIMENT.

    This func takes a list of exactly two dictionaries
    that have the same value for `dataset` key and returns a composite
    dataset name assembled from the actual dataset name + parameter that
    is different between the two dicts (e.g. mip or exp or grid etc.)
    """
    orig_ctrl, orig_exper = dat_selection
    dif_pars = [
        key for key in orig_ctrl.keys() & orig_exper
        if orig_ctrl[key] != orig_exper[key]
    ]

    # do not populate dataset name with these (long) params
    unwanted_pars = ['filename', 'alias', 'recipe_dataset_index']
    dif_pars = [
        dif_par for dif_par in dif_pars if dif_par not in unwanted_pars
    ]
    if not dif_pars:
        logger.error("Your CONTROL and EXPERIMENT "
                     "datasets are completely identical, your analysis "
                     "will output garbage, exiting.")
        sys.exit(1)

    # assemble new dataset names
    new_ctrl = [str(orig_ctrl[k]) for k in dif_pars]
    new_ctrl = orig_ctrl['dataset'] + "-" + "-".join(new_ctrl)
    new_exper = [str(orig_exper[k]) for k in dif_pars]
    new_exper = orig_exper['dataset'] + "-" + "-".join(new_exper)

    # recast the new names in the old dicts
    orig_ctrl['dataset'] = new_ctrl
    orig_exper['dataset'] = new_exper

    return orig_ctrl, orig_exper


def get_control_exper_obs(short_name, input_data, cfg, cmip_type):
    """
    Get control, exper and obs datasets.

    This function is used when running recipes that need
    a clear distinction between a control dataset, an experiment
    dataset and have optional obs (OBS, obs4MIPs etc) datasets;
    such recipes include recipe_validation, and all the autoassess
    ones;
    short_name: variable short name
    input_data: dict containing the input data info
    cfg: config file as used in this module
    """
    # select data per short name and CMIP type
    dataset_selection = select_metadata(
        input_data, short_name=short_name, project=cmip_type)

    # get the obs datasets if specified in recipe
    if 'observational_datasets' in cfg:
        obs_selection = [
            select_metadata(
                input_data, short_name=short_name, dataset=obs_dataset)[0]
            for obs_dataset in cfg['observational_datasets']
        ]
    else:
        obs_selection = []

    # print out OBS's
    if obs_selection:
        logger.info("Observations dataset(s) %s",
                    [obs['dataset'] for obs in obs_selection])

    # determine CONTROL and EXPERIMENT datasets

    # corner case: they could be the same dataset name
    if cfg['control_model'] == cfg['exper_model']:
        logger.info("Identical Control/Experiment dataset names: %s",
                    dataset_selection[0]['dataset'])
        control, experiment = _disentagle_iden_datasets(dataset_selection)
        return control, experiment, obs_selection

    # if they're not the same dataset, fire away
    for model in dataset_selection:
        if model['dataset'] == cfg['control_model']:
            logger.info("Control dataset %s", model['dataset'])
            control = model
        elif model['dataset'] == cfg['exper_model']:
            logger.info("Experiment dataset %s", model['dataset'])
            experiment = model

    return control, experiment, obs_selection


# apply supermeans: handy function that loads CONTROL, EXPERIMENT
# and OBS (if any) files and applies climate_statistics() to mean the cubes
def apply_supermeans(ctrl, exper, obs_list):
    """
    Apply supermeans on data components ie MEAN on time.

    This function is an extension of climate_statistics() meant to ease the
    time-meaning procedure when dealing with CONTROL, EXPERIMENT and OBS
    (if any) datasets.
    ctrl: dictionary of CONTROL dataset
    exper: dictionary of EXPERIMENT dataset
    obs_lis: list of dicts for OBS datasets (0, 1 or many)

    Returns: control and experiment cubes and list of obs cubes
    """
    ctrl_file = ctrl['filename']
    exper_file = exper['filename']
    ctrl_cube = iris.load_cube(ctrl_file)
    exper_cube = iris.load_cube(exper_file)
    ctrl_cube = climate_statistics(ctrl_cube)
    exper_cube = climate_statistics(exper_cube)
    if obs_list:
        obs_cube_list = []
        for obs in obs_list:
            obs_file = obs['filename']
            obs_cube = iris.load_cube(obs_file)
            obs_cube = climate_statistics(obs_cube)
            obs_cube_list.append(obs_cube)
    else:
        obs_cube_list = None

    return ctrl_cube, exper_cube, obs_cube_list
