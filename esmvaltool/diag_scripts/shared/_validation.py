"""Load functions needed by diags with CONTROL and EXPERIMENT."""

import logging
import os

import iris
from esmvalcore.preprocessor import climate_statistics

from esmvaltool.diag_scripts.shared import select_metadata

logger = logging.getLogger(os.path.basename(__file__))


def get_control_exper_obs(short_name, input_data, cfg, cmip_type=None):
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
    cmip_type: optional, CMIP project type (CMIP5 or CMIP6)
    """
    # select data per short name and optional CMIP type
    if not cmip_type:
        dataset_selection = select_metadata(input_data, short_name=short_name)
    else:
        dataset_selection = select_metadata(
            input_data, short_name=short_name, project=cmip_type
        )

    # get the obs datasets if specified in recipe
    if "observational_datasets" in cfg:
        obs_selection = [
            select_metadata(
                input_data, short_name=short_name, dataset=obs_dataset
            )[0]
            for obs_dataset in cfg["observational_datasets"]
        ]
    else:
        obs_selection = []

    # print out OBS's
    if obs_selection:
        logger.info(
            "Observations dataset(s) %s",
            [obs["dataset"] for obs in obs_selection],
        )

    # make sure the chosen datasets for control and exper are available
    alias_selection = []
    for model in dataset_selection:
        try:
            dataset_name = model["alias"].split("_")[1]
        except IndexError:
            dataset_name = model["alias"]
        alias_selection.append(dataset_name)

    if cfg["control_model"] not in alias_selection:
        raise ValueError(
            f"Control dataset {cfg['control_model']} not in datasets"
        )

    if cfg["exper_model"] not in alias_selection:
        raise ValueError(
            f"Experiment dataset {cfg['exper_model']} not in datasets"
        )

    # pick control and experiment dataset
    for model in dataset_selection:
        if cfg["control_model"] in model["alias"].split("_"):
            logger.info("Control dataset %s", model["alias"])
            control = model
        elif cfg["exper_model"] in model["alias"].split("_"):
            logger.info("Experiment dataset %s", model["alias"])
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
    ctrl_file = ctrl["filename"]
    exper_file = exper["filename"]
    ctrl_cube = iris.load_cube(ctrl_file)
    exper_cube = iris.load_cube(exper_file)
    ctrl_cube = climate_statistics(ctrl_cube)
    exper_cube = climate_statistics(exper_cube)
    if obs_list:
        obs_cube_list = []
        for obs in obs_list:
            obs_file = obs["filename"]
            obs_cube = iris.load_cube(obs_file)
            obs_cube = climate_statistics(obs_cube)
            obs_cube_list.append(obs_cube)
    else:
        obs_cube_list = None

    return ctrl_cube, exper_cube, obs_cube_list
