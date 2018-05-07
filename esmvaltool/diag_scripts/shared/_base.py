"""Convenience functions for running a diagnostic script."""
import contextlib
import logging
import os
import sys

import yaml

logger = logging.getLogger(__name__)


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg


def _get_input_data_files(cfg):
    """Get a dictionary containing all data input files."""
    input_files = {}
    for filename in cfg['input_files']:
        with open(filename) as file:
            metadata = yaml.safe_load(file)
            input_files.update(metadata)

    return input_files


@contextlib.contextmanager
def run_diagnostic():
    """Run a diagnostic."""
    # Implemented as context manager so we can support clean up actions later
    cfg = get_cfg()

    # Set up logging
    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    logging.getLogger().setLevel(cfg['log_level'].upper())

    # Read input metadata
    cfg['input_data'] = _get_input_data_files(cfg)

    # Create output directories
    if cfg['write_netcdf']:
        os.makedirs(cfg['work_dir'])
    if cfg['write_plots']:
        os.makedirs(cfg['plot_dir'])

    yield cfg
