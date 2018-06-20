"""Convenience functions for running a diagnostic script."""
import argparse
import contextlib
import logging
import os
import shutil
import sys
import time

import yaml

logger = logging.getLogger(__name__)


def get_cfg(filename=None):
    """Read diagnostic script configuration from settings.yml."""
    if filename is None:
        filename = sys.argv[1]
    with open(filename) as file:
        cfg = yaml.safe_load(file)
    return cfg


def _get_input_data_files(cfg):
    """Get a dictionary containing all data input files."""
    input_files = {}
    for filename in cfg['input_files']:
        if os.path.basename(filename) == 'metadata.yml':
            with open(filename) as file:
                metadata = yaml.safe_load(file)
                input_files.update(metadata)

    return input_files


@contextlib.contextmanager
def run_diagnostic():
    """Run a diagnostic."""
    # Implemented as context manager so we can support clean up actions later
    parser = argparse.ArgumentParser(description="Diagnostic script")
    parser.add_argument('filename', help="Path to settings.yml")
    parser.add_argument(
        '-f',
        '--force',
        help=("Force emptying the output directories"
              "(useful when re-running the script)"),
        action='store_true',
    )
    parser.add_argument(
        '-l',
        '--log-level',
        help=("Set the log-level"),
        choices=['debug', 'info', 'warning', 'error'],
    )
    args = parser.parse_args()

    cfg = get_cfg(args.filename)

    # Set up logging
    if args.log_level:
        cfg['log_level'] = args.log_level

    logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                        "%(name)s,%(lineno)s\t%(message)s")
    logging.Formatter.converter = time.gmtime
    logging.getLogger().setLevel(cfg['log_level'].upper())

    # Read input metadata
    cfg['input_data'] = _get_input_data_files(cfg)

    logger.info("Starting diagnostic script %s with configuration:\n%s",
                cfg['script'], yaml.safe_dump(cfg))

    # Create output directories
    output_directories = []
    if cfg['write_netcdf']:
        output_directories.append(cfg['work_dir'])
    if cfg['write_plots']:
        output_directories.append(cfg['plot_dir'])

    existing = [p for p in output_directories if os.path.exists(p)]

    if existing:
        if args.force:
            for output_directory in existing:
                logger.info("Removing %s", output_directory)
                shutil.rmtree(output_directory)
        else:
            logger.error(
                "Script will abort to prevent accidentally overwriting your "
                "data in these directories:\n%s\n"
                "Use -f or --force to force emptying the output directories.",
                '\n'.join(existing))

    for output_directory in output_directories:
        logger.info("Creating %s", output_directory)
        os.makedirs(output_directory)

    yield cfg

    logger.info("End of diagnostic script run.")
