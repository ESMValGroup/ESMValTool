"""ESMValTool configuration."""
import datetime
import logging
import logging.config
import os
import time
from distutils.version import LooseVersion

import iris
import six
import yaml

from .cmor.table import read_cmor_tables

logger = logging.getLogger(__name__)

CFG = {}
CFG_USER = {}


def use_legacy_iris():
    """Return True if legacy iris is used."""
    return LooseVersion(iris.__version__) < LooseVersion("2.0.0")


def read_config_user_file(config_file, recipe_name):
    """Read config user file and store settings in a dictionary."""
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    # set defaults
    defaults = {
        'write_plots': True,
        'write_netcdf': True,
        'compress_netcdf': False,
        'exit_on_warning': False,
        'max_data_filesize': 100,
        'output_file_type': 'ps',
        'output_dir': './output_dir',
        'auxiliary_data_dir': './auxiliary_data',
        'save_intermediary_cubes': False,
        'remove_preproc_dir': False,
        'max_parallel_tasks': 1,
        'run_diagnostic': True,
        'profile_diagnostic': False,
        'config_developer_file': None,
        'drs': {},
    }

    for key in defaults:
        if key not in cfg:
            logger.info(
                "No %s specification in config file, "
                "defaulting to %s", key, defaults[key])
            cfg[key] = defaults[key]

    cfg['output_dir'] = _normalize_path(cfg['output_dir'])
    cfg['auxiliary_data_dir'] = _normalize_path(cfg['auxiliary_data_dir'])

    cfg['config_developer_file'] = _normalize_path(
        cfg['config_developer_file'])

    for key in cfg['rootpath']:
        root = cfg['rootpath'][key]
        if isinstance(root, six.string_types):
            cfg['rootpath'][key] = [_normalize_path(root)]
        else:
            cfg['rootpath'][key] = [_normalize_path(path) for path in root]

    # insert a directory date_time_recipe_usertag in the output paths
    now = datetime.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    new_subdir = '_'.join((recipe_name, now))
    cfg['output_dir'] = os.path.join(cfg['output_dir'], new_subdir)

    # create subdirectories
    cfg['preproc_dir'] = os.path.join(cfg['output_dir'], 'preproc')
    cfg['work_dir'] = os.path.join(cfg['output_dir'], 'work')
    cfg['plot_dir'] = os.path.join(cfg['output_dir'], 'plots')
    cfg['run_dir'] = os.path.join(cfg['output_dir'], 'run')

    # Save user configuration in global variable
    for key, value in six.iteritems(cfg):
        CFG_USER[key] = value

    # Read developer configuration file
    cfg_developer = read_config_developer_file(cfg['config_developer_file'])
    for key, value in six.iteritems(cfg_developer):
        CFG[key] = value
    read_cmor_tables(CFG)

    return cfg


def get_config_user_file():
    """Return user configuration dictionary."""
    return CFG_USER


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
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


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
    """Set up logging."""
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
    logging.Formatter.converter = time.gmtime
    logging.captureWarnings(True)

    return log_files


def get_project_config(project):
    """Get developer-configuration for project."""
    logger.debug("Retrieving %s configuration", project)
    return CFG[project]


def get_institutes(variable):
    """Return the institutes given the dataset name in CMIP5."""
    dataset = variable['dataset']
    project = variable['project']
    logger.debug("Retrieving institutes for dataset %s", dataset)
    return CFG.get(project, {}).get('institutes', {}).get(dataset, [])


def replace_mip_fx(fx_file):
    """Replace MIP so to retrieve correct fx files."""
    default_mip = 'Amon'
    if fx_file not in CFG['CMIP5']['fx_mip_change']:
        logger.warning(
            'mip for fx variable %s is not specified in '
            'config_developer.yml, using default (%s)', fx_file, default_mip)
    new_mip = CFG['CMIP5']['fx_mip_change'].get(fx_file, default_mip)
    logger.debug("Switching mip for fx file finding to %s", new_mip)
    return new_mip


TAGS_CONFIG_FILE = os.path.join(
    os.path.dirname(__file__), 'config-references.yml')


def _load_tags(filename=TAGS_CONFIG_FILE):
    """Load the refence tags used for provenance recording."""
    logger.debug("Loading tags from %s", filename)
    with open(filename) as file:
        return yaml.safe_load(file)


TAGS = _load_tags()


def get_tag_value(section, tag):
    """Retrieve the value of a tag."""
    if section not in TAGS:
        raise ValueError("Section '{}' does not exist in {}".format(
            section, TAGS_CONFIG_FILE))
    if tag not in TAGS[section]:
        raise ValueError(
            "Tag '{}' does not exist in section '{}' of {}".format(
                tag, section, TAGS_CONFIG_FILE))
    return TAGS[section][tag]


def replace_tags(section, tags):
    """Replace a list of tags with their values."""
    return tuple(get_tag_value(section, tag) for tag in tags)
