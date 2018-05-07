"""ESMValTool configuration."""
import datetime
import logging
import logging.config
import os

import yaml

logger = logging.getLogger(__name__)


def read_config_user_file(config_file, namelist_name):
    """Read config user file and store settings in a dictionary."""
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    # set defaults
    defaults = {
        'write_plots': True,
        'write_netcdf': True,
        'exit_on_warning': False,
        'max_data_filesize': 100,
        'output_file_type': 'ps',
        'output_dir': './output_dir',
        'save_intermediary_cubes': False,
        'max_parallel_tasks': 1,
        'run_diagnostic': True,
        'drs': {},
    }

    for key in defaults:
        if key not in cfg:
            logger.warning("No %s specification in config file, "
                           "defaulting to %s", key, defaults[key])
            cfg[key] = defaults[key]

    # expand ~ to /home/username in directory names and normalize paths
    cfg['output_dir'] = os.path.abspath(os.path.expanduser(cfg['output_dir']))

    for key in cfg['rootpath']:
        cfg['rootpath'][key] = os.path.abspath(
            os.path.expanduser(cfg['rootpath'][key]))

    # insert a directory date_time_namelist_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((namelist_name, now))
    cfg['output_dir'] = os.path.join(cfg['output_dir'], new_subdir)

    # create subdirectories
    cfg['preproc_dir'] = os.path.join(cfg['output_dir'], 'preproc')
    cfg['work_dir'] = os.path.join(cfg['output_dir'], 'work')
    cfg['plot_dir'] = os.path.join(cfg['output_dir'], 'plots')
    cfg['run_dir'] = os.path.join(cfg['output_dir'], 'run')

    return cfg


def read_config_developer_file(cfg_file=None):
    """Read the developer's configuration file."""
    if cfg_file is None:
        cfg_file = os.path.join(
            os.path.dirname(__file__),
            'config-developer.yml',
        )

    with open(cfg_file, 'r') as file:
        cfg = yaml.safe_load(file)

    return cfg


def configure_logging(cfg_file=None, output=None, console_log_level=None):
    """Set up logging"""
    if cfg_file is None:
        cfg_file = os.path.join(
            os.path.dirname(__file__), 'config-logging.yml')

    if output is None:
        output = os.getcwd()

    cfg_file = os.path.abspath(cfg_file)
    with open(cfg_file) as file_handler:
        cfg = yaml.safe_load(file_handler)

    log_files = []
    for handler in cfg['handlers'].values():
        if 'filename' in handler:
            if not os.path.isabs(handler['filename']):
                handler['filename'] = os.path.join(output, handler['filename'])
            log_files.append(handler['filename'])
        if console_log_level is not None and 'stream' in handler:
            if handler['stream'] in ('ext://sys.stdout', 'ext://sys.stderr'):
                handler['level'] = console_log_level.upper()

    logging.config.dictConfig(cfg)

    return log_files


CFG = read_config_developer_file()


def get_project_config(project):
    """Get developer-configuration for project."""
    logger.debug("Retrieving %s configuration", project)
    return CFG[project]


def cmip5_model2inst(model):
    """Return the institute given the model name in CMIP5."""
    logger.debug("Retrieving institute for CMIP5 model %s", model)
    return CFG['CMIP5']['institute'][model]


def cmip5_mip2realm_freq(mip):
    """Return realm and frequency given the mip in CMIP5."""
    logger.debug("Retrieving realm and frequency for CMIP5 mip %s", mip)
    return CFG['CMIP5']['realm_frequency'][mip]
