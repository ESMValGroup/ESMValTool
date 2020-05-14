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
esmvalcore.cmor/cmorizers/obs
"""
import argparse
import datetime
import importlib
import logging
import os
import subprocess
from pathlib import Path

import esmvalcore
from esmvalcore._config import configure_logging, read_config_user_file
from esmvalcore._task import write_ncl_settings

from esmvaltool.cmorizers.obs.utilities import read_cmor_config

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

    # get all available tiers in source dir
    tiers = ['Tier{}'.format(i) for i in range(2, 4)]
    tiers = [
        tier for tier in tiers if os.path.exists(os.path.join(raw_obs, tier))
    ]
    datasets = {tier: [] for tier in tiers}

    # if user specified obs list
    if obs_list:
        for dataset_name in obs_list:
            for tier in tiers:
                if os.path.isdir(os.path.join(raw_obs, tier, dataset_name)):
                    datasets[tier].append(dataset_name)
                    break
            else:
                logger.warning("Could not find raw data %s in %s/%s",
                               dataset_name, raw_obs, tier)

    # otherwise go through the whole raw_obs dir
    else:
        for tier in tiers:
            for dats in os.listdir(os.path.join(raw_obs, tier)):
                datasets[tier].append(dats)

    return datasets


def _write_ncl_settings(project_info, dataset, run_dir, reformat_script,
                        log_level):
    """Write the information needed by the ncl reformat script."""
    settings = {
        'cmorization_script': reformat_script,
        'input_dir_path': project_info[dataset]['indir'],
        'output_dir_path': project_info[dataset]['outdir'],
        'config_user_info': {
            'log_level': log_level
        },
    }
    settings_filename = os.path.join(run_dir, dataset, 'settings.ncl')
    if not os.path.isdir(os.path.join(run_dir, dataset)):
        os.makedirs(os.path.join(run_dir, dataset))
    # write the settings file
    write_ncl_settings(settings, settings_filename)
    return settings_filename


def _run_ncl_script(in_dir, out_dir, run_dir, dataset, reformat_script,
                    log_level):
    """Run the NCL cmorization mechanism."""
    logger.info("CMORizing dataset %s using NCL script %s",
                dataset, reformat_script)
    project = {}
    project[dataset] = {}
    project[dataset]['indir'] = in_dir
    project[dataset]['outdir'] = out_dir
    settings_file = _write_ncl_settings(project, dataset, run_dir,
                                        reformat_script, log_level)
    esmvaltool_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(reformat_script))))

    # put settings in environment
    env = dict(os.environ)
    env['settings'] = settings_file
    env['esmvaltool_root'] = esmvaltool_root
    env['cmor_tables'] = str(Path(esmvalcore.cmor.__file__).parent / 'tables')
    logger.info("Using CMOR tables at %s", env['cmor_tables'])
    # call NCL
    ncl_call = ['ncl', reformat_script]
    logger.info("Executing cmd: %s", ' '.join(ncl_call))
    process = subprocess.Popen(ncl_call,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               env=env)
    output, err = process.communicate()
    for oline in str(output.decode('utf-8')).split('\n'):
        logger.info('[NCL] %s', oline)
    if err:
        logger.info('[NCL][subprocess.Popen ERROR] %s', err)


def _run_pyt_script(in_dir, out_dir, dataset, user_cfg):
    """Run the Python cmorization mechanism."""
    module_name = 'esmvaltool.cmorizers.obs.formatters.datasets.{}'.format(
        dataset.lower().replace("-", "_"))
    module = importlib.import_module(module_name)
    logger.info("CMORizing dataset %s using Python script %s",
                dataset, module.__file__)
    cmor_cfg = read_cmor_config(dataset)
    module.cmorization(in_dir, out_dir, cmor_cfg, user_cfg)


def main():
    """Run it as executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o',
                        '--obs-list-cmorize',
                        type=str,
                        help='List of obs datasets to cmorize. \
              If no list provided: CMORization of \
              all datasets in RAWOBS; \
              -o DATASET1,DATASET2... : \
              for CMORization of select datasets.')
    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')
    parser.add_argument('-d', '--download', action='store_true')
    parser.add_argument('--startdate',
                        '-s',
                        help='starting date as YYYYMMDD')
    parser.add_argument('--enddate',
                        '-e',
                        help='end date as YYYYMMDD')
    args = parser.parse_args()

    # get and read config file
    config_file = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    # Read user config file
    if not os.path.exists(config_file):
        logger.error("Config file %s does not exist", config_file)

    # read the file in
    config_user = read_config_user_file(config_file, 'cmorize_obs', options={})

    # set the run dir to hold the settings and log files
    run_dir = os.path.join(config_user['output_dir'], 'run')
    if not os.path.isdir(run_dir):
        os.makedirs(run_dir)

    # configure logging
    log_files = configure_logging(
        output_dir=run_dir, console_log_level=config_user['log_level'])
    logger.info("Writing program log files to:\n%s", "\n".join(log_files))

    # print header
    logger.info(HEADER)

    # run
    timestamp1 = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d %H:%M:%S"

    logger.info("Starting the CMORization Tool at time: %s UTC",
                timestamp1.strftime(timestamp_format))

    logger.info(70 * "-")
    logger.info("input_dir  = %s", config_user["rootpath"]["RAWOBS"][0])
    # check if the inputdir actually exists
    if not args.download and \
       not os.path.isdir(config_user["rootpath"]["RAWOBS"][0]):
        logger.error("Directory %s does not exist",
                     config_user["rootpath"]["RAWOBS"][0])
        raise ValueError
    logger.info("output_dir = %s", config_user["output_dir"])
    logger.info(70 * "-")

    # call the reformat function
    if args.obs_list_cmorize:
        obs_list = args.obs_list_cmorize.split(',')
    else:
        obs_list = []
    if args.download:
        if not obs_list:
            logger.error(
                "In order to download automatically, you must provide "
                "the desired datasets"
            )
            raise ValueError
        start_date = datetime.datetime.strptime(args.startdate, '%Y%m%d')
        end_date = datetime.datetime.strptime(args.enddate, '%Y%m%d')
        _download(config_user, obs_list, start_date, end_date)

    _format(config_user, obs_list)

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info("Ending the CMORization Tool at time: %s UTC",
                timestamp2.strftime(timestamp_format))
    logger.info("Time for running the CMORization scripts was: %s",
                timestamp2 - timestamp1)


def _download(config, obs_list, start_date, end_date):
    """Automatic download"""
    logger.info("Downloading RAW data")
    # master directory
    failed_datasets = []
    for dataset in obs_list:
        dataset_module = dataset.lower().replace('-', '_')
        logger.info("Download module: %s", dataset_module)
        try:
            downloader = importlib.import_module(
                f'.{dataset_module}',
                package='esmvaltool.cmorizers.obs.downloaders.datasets'
            )
        except ImportError:
            logger.exception('Could not find cmorizer for %s', dataset)
            failed_datasets.append(dataset)
        else:
            downloader.download_dataset(config, dataset, start_date, end_date)

    if failed_datasets:
        raise Exception(
            'Could not find downloader for %s datasets ' %
            ' '.join(failed_datasets)
        )


def _format(config, obs_list):
    """Run the cmorization routine."""
    logger.info("Running the CMORization scripts.")

    # master directory
    raw_obs = config["rootpath"]["RAWOBS"][0]

    # set the reformat scripts dir
    reformat_scripts = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'formatters', 'datasets'
    )
    logger.info("Using cmorizer scripts repository: %s", reformat_scripts)
    run_dir = os.path.join(config['output_dir'], 'run')
    # datsets dictionary of Tier keys
    datasets = _assemble_datasets(raw_obs, obs_list)
    if not datasets:
        logger.warning("Check input: could not find required %s in %s",
                       obs_list, raw_obs)
    logger.info("Processing datasets %s", datasets)

    # loop through tier/datasets to be cmorized
    failed_datasets = []
    for tier in datasets:
        for dataset in datasets[tier]:
            reformat_script_root = os.path.join(
                reformat_scripts,
                dataset.lower().replace('-', '_'),
            )
            # in-data dir; build out-dir tree
            in_data_dir = os.path.join(raw_obs, tier, dataset)
            logger.info("Input data from: %s", in_data_dir)
            out_data_dir = os.path.join(config['output_dir'], tier, dataset)
            logger.info("Output will be written to: %s", out_data_dir)
            if not os.path.isdir(out_data_dir):
                os.makedirs(out_data_dir)

            # all operations are done in the working dir now
            os.chdir(out_data_dir)

            # figure out what language the script is in
            logger.info("Reformat script: %s", reformat_script_root)
            if os.path.isfile(reformat_script_root + '.ncl'):
                reformat_script = reformat_script_root + '.ncl'
                _run_ncl_script(
                    in_data_dir,
                    out_data_dir,
                    run_dir,
                    dataset,
                    reformat_script,
                    config['log_level'],
                )
            elif os.path.isfile(reformat_script_root + '.py'):
                _run_pyt_script(in_data_dir, out_data_dir, dataset, config)
            else:
                logger.error('Could not find cmorizer for %s', dataset)
                failed_datasets.append(dataset)
                raise Exception(
                    'Could not find cmorizers for %s datasets ' %
                    ' '.join(failed_datasets)
                )


class DataCommand():

    def _start(self, config_file):
        # get and read config file
        config_file = os.path.abspath(
            os.path.expandvars(os.path.expanduser(config_file)))

        # Read user config file
        if not os.path.exists(config_file):
            logger.error("Config file %s does not exist", config_file)

        # read the file in
        config_user = read_config_user_file(config_file, 'cmorize_obs')

        # set the run dir to hold the settings and log files
        run_dir = os.path.join(config_user['output_dir'], 'run')
        if not os.path.isdir(run_dir):
            os.makedirs(run_dir)

        # configure logging
        log_files = configure_logging(
            output=run_dir, console_log_level=config_user['log_level'])
        logger.info("Writing program log files to:\n%s", "\n".join(log_files))

        # print header
        logger.info(HEADER)

        # run
        timestamp1 = datetime.datetime.utcnow()
        timestamp_format = "%Y-%m-%d %H:%M:%S"

        logger.info("Starting the CMORization Tool at time: %s UTC",
                    timestamp1.strftime(timestamp_format))

        logger.info(70 * "-")
        logger.info("input_dir  = %s", config_user["rootpath"]["RAWOBS"][0])
        # check if the inputdir actually exists
        if not os.path.isdir(config_user["rootpath"]["RAWOBS"][0]):
            logger.error("Directory %s does not exist",
                         config_user["rootpath"]["RAWOBS"][0])
            raise ValueError
        logger.info("output_dir = %s", config_user["output_dir"])
        logger.info(70 * "-")
        return config_user

    def _parse_datasets(self, datasets):
        if datasets:
            return datasets.split(',')
        return []

    def download(self, datasets, config_file, start_date=None, end_date=None,
                 overwrite=False):
        config_user = self._start(config_file)
        datasets = self._parse_datasets(datasets)
        if not datasets:
            logger.error(
                "In order to download automatically, you must provide "
                "the desired datasets"
            )
            raise ValueError
        start_date = datetime.datetime.strptime(str(start_date), '%Y%m%d')
        end_date = datetime.datetime.strptime(str(end_date), '%Y%m%d')
        _download(config_user, datasets, start_date, end_date)

    def format(self, dataset=None, config_file=None):
        config_user = self._start(config_file)
        datasets = self._parse_datasets(datasets)
        if not datasets:
            logger.error(
                "In order to download automatically, you must provide "
                "the desired datasets"
            )
            raise ValueError
        start_date = datetime.datetime.strptime(str(start_date), '%Y%m%d')
        end_date = datetime.datetime.strptime(str(end_date), '%Y%m%d')
        _format(config_user, datasets)

    def prepare(self, datasets, config_file, start_date=None, end_date=None,
                overwrite=False):
        config_user = self._start(config_file)
        datasets = self._parse_datasets(datasets)
        if not datasets:
            logger.error(
                "In order to download automatically, you must provide "
                "the desired datasets"
            )
            raise ValueError
        start_date = datetime.datetime.strptime(str(start_date), '%Y%m%d')
        end_date = datetime.datetime.strptime(str(end_date), '%Y%m%d')
        _download(config_user, datasets, start_date, end_date)
        _format(config_user, datasets)


if __name__ == '__main__':
    main()
