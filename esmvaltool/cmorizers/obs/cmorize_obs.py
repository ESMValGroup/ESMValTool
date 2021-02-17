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
import numbers
import os
import subprocess
import time
from pathlib import Path

import yaml

import esmvalcore.cmor

from .utilities import read_cmor_config

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


def _normalize_path(path):
    """Normalize paths.

    Expand ~ character and environment variables and convert path to absolute.

    Parameters
    ----------
    path: str
        Original path

    Returns
    -------
    str:
        Normalized path
    """
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def read_config_user_file(config_file):
    """Read config user file and store settings in a dictionary."""
    config_file = _normalize_path(config_file)

    # Check user config file
    if not os.path.exists(config_file):
        print(f"Config file {config_file} does not exist")

    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    for key in cfg['rootpath']:
        root = cfg['rootpath'][key]
        if isinstance(root, str):
            cfg['rootpath'][key] = [_normalize_path(root)]
        else:
            cfg['rootpath'][key] = [_normalize_path(path) for path in root]

    # insert a directory date_time_recipe_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join(('cmorize_obs', now))
    output_dir = _normalize_path(cfg.get('output_dir', 'esmvaltool_output'))

    # set up create subdirectories
    cfg['output_dir'] = os.path.join(output_dir, new_subdir)
    cfg['work_dir'] = os.path.join(cfg['output_dir'], 'work')
    cfg['run_dir'] = os.path.join(cfg['output_dir'], 'run')

    return cfg


def configure_logging(output_dir, log_level='INFO'):
    """Configure logging.

    Parameters
    ----------
    output_dir : str, optional
        Output directory for the log files.
    log_level : str, optional
        If `None`, use the default (INFO).

    """
    fmt = '%(asctime)s UTC [%(process)d] %(levelname)-7s %(message)s'
    debugfmt = ('%(asctime)s UTC [%(process)d] %(levelname)-7s '
                '%(name)s:%(lineno)s %(message)s')

    logging.Formatter.converter = time.gmtime
    logging.captureWarnings(True)

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Configure stdout
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(fmt))
    console_handler.setLevel(log_level.upper())
    root_logger.addHandler(console_handler)

    # Create directory for log files
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Configure main_log.txt
    log_file = Path(output_dir) / 'main_log.txt'
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(logging.Formatter(fmt))
    log_handler.setLevel(log_level.upper())
    root_logger.addHandler(log_handler)

    # Configure main_log_debug.txt
    debug_file = Path(output_dir) / 'main_log_debug.txt'
    debug_handler = logging.FileHandler(debug_file)
    debug_handler.setFormatter(logging.Formatter(debugfmt))
    debug_handler.setLevel(logging.DEBUG)
    root_logger.addHandler(debug_handler)

    # print header and info
    logger.info(HEADER)
    logger.info("Writing program log to:\n%s", log_file)
    logger.info("Writing debug program log to:\n%s", debug_file)


def _py2ncl(value, var_name=''):
    """Format a structure of Python list/dict/etc items as NCL."""
    txt = var_name + ' = ' if var_name else ''
    if value is None:
        txt += '_Missing'
    elif isinstance(value, str):
        txt += '"{}"'.format(value)
    elif isinstance(value, (list, tuple)):
        if not value:
            txt += '_Missing'
        else:
            if isinstance(value[0], numbers.Real):
                type_ = numbers.Real
            else:
                type_ = type(value[0])
            if any(not isinstance(v, type_) for v in value):
                raise ValueError(
                    "NCL array cannot be mixed type: {}".format(value))
            txt += '(/{}/)'.format(', '.join(_py2ncl(v) for v in value))
    elif isinstance(value, dict):
        if not var_name:
            raise ValueError(
                "NCL does not support nested dicts: {}".format(value))
        txt += 'True\n'
        for key in value:
            txt += '{}@{} = {}\n'.format(var_name, key, _py2ncl(value[key]))
    else:
        txt += str(value)
    return txt


def write_ncl_settings(settings, filename, mode='wt'):
    """Write a dictionary with generic settings to NCL file."""
    logger.debug("Writing NCL configuration file %s", filename)

    def _ncl_type(value):
        """Convert some Python types to NCL types."""
        typemap = {
            bool: 'logical',
            str: 'string',
            float: 'double',
            int: 'int64',
            dict: 'logical',
        }
        for type_ in typemap:
            if isinstance(value, type_):
                return typemap[type_]
        raise ValueError("Unable to map {} to an NCL type".format(type(value)))

    lines = []
    for var_name, value in sorted(settings.items()):
        if isinstance(value, (list, tuple)):
            # Create an NCL list that can span multiple files
            lines.append('if (.not. isdefined("{var_name}")) then\n'
                         '  {var_name} = NewList("fifo")\n'
                         'end if\n'.format(var_name=var_name))
            for item in value:
                lines.append('ListAppend({var_name}, new(1, {type}))\n'
                             'i = ListCount({var_name}) - 1'.format(
                                 var_name=var_name, type=_ncl_type(item)))
                lines.append(_py2ncl(item, var_name + '[i]'))
        else:
            # Create an NCL variable that overwrites previous variables
            lines.append('if (isvar("{var_name}")) then\n'
                         '  delete({var_name})\n'
                         'end if\n'.format(var_name=var_name))
            lines.append(_py2ncl(value, var_name))

    with open(filename, mode) as file:
        file.write('\n'.join(lines))
        file.write('\n')


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
        for dataset_name in obs_list.split(','):
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
    logger.info("CMORizing dataset %s using NCL script %s", dataset,
                reformat_script)
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
    module_name = 'esmvaltool.cmorizers.obs.cmorize_obs_{}'.format(
        dataset.lower().replace("-", "_"))
    module = importlib.import_module(module_name)
    logger.info("CMORizing dataset %s using Python script %s", dataset,
                module.__file__)
    cmor_cfg = read_cmor_config(dataset)
    module.cmorization(in_dir, out_dir, cmor_cfg, user_cfg)


def main():
    """Run it as executable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-o',
        '--obs-list-cmorize',
        type=str,
        help=('List of obs datasets to cmorize. If no list provided: '
              'CMORization of all datasets in RAWOBS; '
              '-o DATASET1,DATASET2... : for CMORization of select datasets.'),
    )
    parser.add_argument(
        '-c',
        '--config-file',
        default=str(Path.home() / '.esmvaltool' / 'config-user.yml'),
        help=('ESMValTool configuration file '
              '(default: ~/.esmvaltool/config-user.yml).'),
    )
    args = parser.parse_args()

    # read the file in
    config_user = read_config_user_file(args.config_file)

    # configure logging
    configure_logging(config_user['run_dir'], config_user['log_level'])

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

    # call the reformat function
    if args.obs_list_cmorize:
        obs_list = args.obs_list_cmorize
    else:
        obs_list = []
    _cmor_reformat(config_user, obs_list)

    # End time timing
    timestamp2 = datetime.datetime.utcnow()
    logger.info("Ending the CMORization Tool at time: %s UTC",
                timestamp2.strftime(timestamp_format))
    logger.info("Time for running the CMORization scripts was: %s",
                timestamp2 - timestamp1)


def _cmor_reformat(config, obs_list):
    """Run the cmorization routine."""
    logger.info("Running the CMORization scripts.")

    # master directory
    raw_obs = config["rootpath"]["RAWOBS"][0]

    # set the reformat scripts dir
    reformat_scripts = os.path.dirname(os.path.abspath(__file__))
    logger.info("Using cmorizer scripts repository: %s", reformat_scripts)
    # datasets dictionary of Tier keys
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
                'cmorize_obs_' + dataset.lower().replace('-', '_'),
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
                    config['run_dir'],
                    dataset,
                    reformat_script,
                    config['log_level'],
                )
            elif os.path.isfile(reformat_script_root + '.py'):
                _run_pyt_script(in_data_dir, out_data_dir, dataset, config)
            else:
                logger.error('Could not find cmorizer for %s', dataset)
                failed_datasets.append(dataset)
                raise Exception('Could not find cmorizers for %s datasets ' %
                                ' '.join(failed_datasets))


if __name__ == '__main__':
    main()
