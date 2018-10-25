"""Run the CMORization module."""
import logging
import os
import sys
import subprocess
import shutil

from ._task import write_ncl_settings


logger = logging.getLogger(__name__)


def _assemble_datasets(raw_obs, obs_list):
    """Get my datasets as dictionary keyed on Tier."""
    # check for desired datasets only (if any)
    # if not, walk all over rawobs dir
    # assume a RAWOBS/TierX/DATASET input structure
    tiers = []
    datasets = {}

    # get all available tiers in source dir
    for _, tier, _ in os.walk(raw_obs, followlinks=True):
        tiers.append(tier)

    # if user specified obs list
    if obs_list:
        for tier in tiers[0]:
            datasets[tier] = []
            for dataset_name in obs_list.split(','):
                if os.path.isdir(os.path.join(raw_obs, tier, dataset_name)):
                    datasets[tier].append(dataset_name)

    # otherwise go through the whole raw_obs dir
    else:
        for tier in tiers[0]:
            datasets[tier] = []
            for _, dats, _ in os.walk(os.path.join(raw_obs, tier),
                                      followlinks=True):
                datasets[tier].append(dats)
            datasets[tier] = datasets[tier][0]

    return datasets


def _write_ncl_infofile(project_info, dataset, output_dir):
    """Write the information needed by the ncl reformat script."""
    info = {
        'input_dir_path': project_info[dataset]['indir'],
        'output_dir_path': project_info[dataset]['outdir'],
    }

    filename = os.path.join(output_dir, dataset + '_info.ncl')
    if not os.path.isfile(filename):
        write_ncl_settings(info, filename)
    return filename


def _run_ncl_script(in_dir,
                    out_dir,
                    dataset,
                    reformat_script):
    """Run the NCL cmorization mechanism."""
    project = {}
    project[dataset] = {}
    project[dataset]['indir'] = in_dir
    project[dataset]['outdir'] = out_dir
    _write_ncl_infofile(project, dataset, out_dir)
    ncl_call = ['ncl', os.path.basename(reformat_script)]
    logger.info("Executing cmd: %s", ' '.join(ncl_call))
    process = subprocess.Popen(ncl_call, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    output, _ = process.communicate()
    for oline in str(output).split('\\n'):
        logger.info('[NCL] ' + oline)


def _run_pyt_script(in_dir, out_dir):
    """Run the Python cmorization mechanism."""
    import py_cmor
    py_cmor.cmorization(in_dir, out_dir)


def cmor_reformat(config, obs_list):
    """Run the oldskool v1 cmorization scripts."""
    logger.info("Running the CMORization scripts.")

    # master directory
    raw_obs = config["rootpath"]["RAWOBS"][0]

    # set the reformat scripts dir
    reformat_scripts = os.path.join(os.path.dirname(__file__),
                                    'cmor/cmorize_obs')

    # datsets dictionary of Tier keys
    datasets = _assemble_datasets(raw_obs, obs_list)
    logger.info("Processing datasets %s", datasets)

    # assemble i/o information
    project_info = {}

    # loop through tier/datasets to be cmorized
    for tier, _ in datasets.items():
        for dataset in datasets[tier]:
            project_info[dataset] = {}
            reformat_script_root = os.path.join(reformat_scripts,
                                                'cmorize_obs_' + dataset)

            # in-data dir; build out-dir tree
            in_data_dir = os.path.join(raw_obs, tier, dataset)
            out_data_dir = os.path.join(config['output_dir'], tier, dataset)
            if not os.path.isdir(out_data_dir):
                os.makedirs(out_data_dir)

            # all operations are done in the working dir now
            os.chdir(out_data_dir)

            # figure out what language the script is in
            if os.path.isfile(reformat_script_root + '.ncl'):
                reformat_script = reformat_script_root + '.ncl'
                logger.info("CMORizing dataset %s using NCL script %s",
                            dataset, reformat_script)
                # copy over the reformat script
                shutil.copyfile(reformat_script,
                                os.path.join(out_data_dir, reformat_script))
                # call the ncl script
                _run_ncl_script(in_data_dir,
                                out_data_dir,
                                dataset,
                                reformat_script)
            elif os.path.isfile(reformat_script_root + '.py'):
                py_reformat_script = reformat_script_root + '.py'
                logger.info("CMORizing dataset %s using Python script %s",
                            dataset, py_reformat_script)
                # copy over the reformat script
                shutil.copyfile(py_reformat_script,
                                os.path.join(out_data_dir, 'py_cmor.py'))
                sys.path.append(out_data_dir)
                _run_pyt_script(in_data_dir, out_data_dir)
            else:
                logger.info("No need to CMORize, could not find CMOR script.")
