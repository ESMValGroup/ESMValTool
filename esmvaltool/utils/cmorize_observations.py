"""
Run the CMORization module as a utility executable.

This utility allows the user to call and execute CMOR reformatting
scripts (support for NCL and Python at the moment), that will use
two I/O variables passed by this utility: an input directory as
specified in config-user.yml by the RAWOBS key, and an output dir
created in the form of output_dir/CMOR_DATE_TIME/TierTIER/DATASET.
The user can specify a list of DATASETS that the CMOR reformatting
can by run on by using -o (--obs-list-cmorize) command line argument.
The CMOR reformatting scripts are to be found in:
esmvaltool/cmor/cmorize_obs

    Usage
    ------
        cmorize_observations --help
        cmorize_observations -c config-user.yml (for CMORization of
            all datasets in RAWOBS)
        cmorize_observations -c config-user.yml -o DATASET1,DATASET2...
            (for CMORization of select datasets)
        cmorize_observations -c config-user.yml -o DATA1,DATA2 -l LOGLEVEL
            (to set the log level: debug, info, warning, error)

"""
import argparse
import logging
import os
import sys
import datetime
import subprocess
import shutil

from .._task import write_ncl_settings
from .._config import read_config_user_file

logger = logging.getLogger(__name__)

HEADER = r"""
______________________________________________________________________
          _____ ____  __  ____     __    _ _____           _
         | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
         |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
         | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
         |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
______________________________________________________________________

""" + __doc__


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


def _write_ncl_settings(project_info, dataset, run_dir, reformat_script):
    """Write the information needed by the ncl reformat script."""
    settings = {
        'cmorization_script': reformat_script,
        'input_dir_path': project_info[dataset]['indir'],
        'output_dir_path': project_info[dataset]['outdir'],
    }
    settings_filename = os.path.join(run_dir, dataset, 'settings.ncl')
    if not os.path.isdir(os.path.join(run_dir, dataset)):
        os.makedirs(os.path.join(run_dir, dataset))
    # write the settings file
    write_ncl_settings(settings, settings_filename)
    return settings_filename


def _run_ncl_script(in_dir,
                    out_dir,
                    run_dir,
                    dataset,
                    reformat_script):
    """Run the NCL cmorization mechanism."""
    project = {}
    project[dataset] = {}
    project[dataset]['indir'] = in_dir
    project[dataset]['outdir'] = out_dir
    settings_file = _write_ncl_settings(project, dataset, run_dir,
                                        reformat_script)
    # put settings in environment
    env = dict(os.environ)
    env['settings'] = settings_file

    # call NCL
    ncl_call = ['ncl', os.path.basename(reformat_script)]
    logger.info("Executing cmd: %s", ' '.join(ncl_call))
    process = subprocess.Popen(ncl_call, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, env=env)
    output, _ = process.communicate()
    for oline in str(output).split('\\n'):
        logger.info('[NCL] %s', oline)


def _run_pyt_script(in_dir, out_dir):
    """Run the Python cmorization mechanism."""
    import py_cmor
    py_cmor.cmorization(in_dir, out_dir)


def execute_cmorize():
    """Run it as executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-o',
        '--obs-list-cmorize',
        type=str,
        help='List of obs datasets to cmorize')
    parser.add_argument(
        '-c',
        '--config-file',
        default=os.path.join(os.path.dirname(__file__), 'config-user.yml'),
        help='Config file')
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        choices=['debug', 'info', 'warning', 'error'])
    args = parser.parse_args()

    # get and read config file
    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    # Read user config file
    if not os.path.exists(config_file):
        logger.error("Config file %s does not exist", config_file)

    # read the file in
    config_user = read_config_user_file(config_file, 'CMOR')

    # set the run dir to hold the settings and log files
    run_dir = os.path.join(config_user['output_dir'], 'run')
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)

    # set logging for screen and file output
    root_logger = logging.getLogger()
    out_fmt = "%(asctime)s %(levelname)-8s %(name)s,%(lineno)s\t%(message)s"
    logging.basicConfig(
        filename=os.path.join(run_dir, 'cmorization_log.txt'),
        filemode='a',
        format=out_fmt,
        datefmt='%H:%M:%S',
        level=logging.DEBUG)
    root_logger.setLevel(args.log_level.upper())
    logfmt = logging.Formatter(out_fmt)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logfmt)
    root_logger.addHandler(console_handler)

    # print header
    logger.info(HEADER)

    # run
    timestamp1 = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d %H:%M:%S"

    logger.info(
        "Starting the CMORization Tool at time: %s UTC",
        timestamp1.strftime(timestamp_format))

    logger.info(70 * "-")
    logger.info("INPUTDIR  = %s", config_user["rootpath"]["RAWOBS"][0])
    logger.info("OUTPUTDIR = %s", config_user["output_dir"])
    logger.info(70 * "-")

    # call the reformat function
    if args.obs_list_cmorize:
        obs_list = args.obs_list_cmorize
    else:
        obs_list = []
    _cmor_reformat(config_user, obs_list)

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info(
        "Ending the CMORization Tool at time: %s UTC",
        timestamp2.strftime(timestamp_format))
    logger.info(
        "Time for running the CMORization scripts was: %s",
        timestamp2 - timestamp1)


def _cmor_reformat(config, obs_list):
    """Run the cmorization routine."""
    logger.info("Running the CMORization scripts.")

    # master directory
    raw_obs = config["rootpath"]["RAWOBS"][0]

    # set the reformat scripts dir
    reformat_scripts = os.path.join(os.path.dirname(__file__),
                                    '../cmor/cmorize_obs')
    run_dir = os.path.join(config['output_dir'], 'run')
    # datsets dictionary of Tier keys
    datasets = _assemble_datasets(raw_obs, obs_list)
    logger.info("Processing datasets %s", datasets)

    # loop through tier/datasets to be cmorized
    for tier, _ in datasets.items():
        for dataset in datasets[tier]:
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
                shutil.copy2(reformat_script, out_data_dir)
                # call the ncl script
                _run_ncl_script(in_data_dir,
                                out_data_dir,
                                run_dir,
                                dataset,
                                reformat_script)
            elif os.path.isfile(reformat_script_root + '.py'):
                py_reformat_script = reformat_script_root + '.py'
                logger.info("CMORizing dataset %s using Python script %s",
                            dataset, py_reformat_script)
                # copy over the reformat script
                shutil.copy2(py_reformat_script,
                             os.path.join(out_data_dir, 'py_cmor.py'))
                sys.path.append(out_data_dir)
                _run_pyt_script(in_data_dir, out_data_dir)
            else:
                logger.info("No need to CMORize, could not find CMOR script.")


if __name__ == '__main__':
    execute_cmorize()
