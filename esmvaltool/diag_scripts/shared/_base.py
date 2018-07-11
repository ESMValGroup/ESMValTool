"""Convenience functions for running a diagnostic script."""
import argparse
import contextlib
import glob
import logging
import os
import shutil
import sys
import time
from collections import OrderedDict

import yaml

logger = logging.getLogger(__name__)


def select_metadata(metadata, **attributes):
    """Select specific metadata describing preprocessed data.

    Parameters
    ----------
    metadata : :obj:`list` of :obj:`dict`
        A list of metadata describing preprocessed data.
    **attributes :
        Keyword arguments specifying the required variable attributes and
        their values.
        Use the value '*' to select any variable that has the attribute.

    Returns
    -------
    :obj:`list` of :obj:`dict`
        A list of matching metadata.

    """
    selection = []
    for attribs in metadata:
        if all(
                a in attribs and (
                    attribs[a] == attributes[a] or attributes[a] == '*')
                for a in attributes):
            selection.append(attribs)
    return selection


def group_metadata(metadata, attribute, sort=None):
    """Group metadata describing preprocessed data by attribute.

    Parameters
    ----------
    metadata : :obj:`list` of :obj:`dict`
        A list of metadata describing preprocessed data.
    attribute : str
        The attribute name that the metadata should be grouped by.
    sort :
        See `sorted_group_metadata`.

    Returns
    -------
    :obj:`dict` of :obj:`list` of :obj:`dict`
        A dictionary containing the requested groups. If sorting is requested,
        an `OrderedDict` will be returned.

    """
    groups = {}
    for attributes in metadata:
        key = attributes.get(attribute)
        if key not in groups:
            groups[key] = []
        groups[key].append(attributes)

    if sort:
        groups = sorted_group_metadata(groups, sort)

    return groups


def sorted_metadata(metadata, sort):
    """Sort a list of metadata describing preprocessed data.

    Sorting is done on strings and is not case sensitive.

    Parameters
    ----------
    metadata : :obj:`list` of :obj:`dict`
        A list of metadata describing preprocessed data.
    sort : :obj:`str` or :obj:`list` of :obj:`str`
        One or more attributes to sort by.

    Returns
    -------
    :obj:`list` of :obj:`dict`
        The sorted list of variable metadata.

    """
    if isinstance(sort, str):
        sort = [sort]

    def normalized_variable_key(attributes):
        """Define a key to sort the list of attributes by."""
        return tuple(str(attributes.get(k, '')).lower() for k in sort)

    return sorted(metadata, key=normalized_variable_key)


def sorted_group_metadata(metadata_groups, sort):
    """Sort grouped metadata.

    Sorting is done on strings and is not case sensitive.

    Parameters
    ----------
    metadata_groups : :obj:`dict` of :obj:`list` of :obj:`dict`
        Dictionary containing the groups of metadata.
    sort : :obj:`bool` or :obj:`str` or :obj:`list` of :obj:`str`
        One or more attributes to sort by or True to just sort the groups but
        not the lists.

    Returns
    -------
    :obj:`OrderedDict` of :obj:`list` of :obj:`dict`
        A dictionary containing the requested groups.

    """
    if sort is True:
        sort = []

    def normalized_group_key(key):
        """Define a key to sort the OrderedDict by."""
        return '' if key is None else str(key).lower()

    groups = OrderedDict()
    for key in sorted(metadata_groups, key=normalized_group_key):
        groups[key] = sorted_metadata(metadata_groups[key], sort)

    return groups


def get_cfg(filename=None):
    """Read diagnostic script configuration from settings.yml."""
    if filename is None:
        filename = sys.argv[1]
    with open(filename) as file:
        cfg = yaml.safe_load(file)
    return cfg


def _get_input_data_files(cfg):
    """Get a dictionary containing all data input files."""
    metadata_files = []
    for filename in cfg['input_files']:
        if os.path.isdir(filename):
            metadata_files.extend(
                glob.glob(os.path.join(filename, '*metadata.yml')))
        elif os.path.basename(filename) == 'metadata.yml':
            metadata_files.append(filename)

    input_files = {}
    for filename in metadata_files:
        with open(filename) as file:
            metadata = yaml.safe_load(file)
            input_files.update(metadata)

    return input_files


@contextlib.contextmanager
def run_diagnostic():
    """Run a diagnostic.

    Example
    -------
    See esmvaltool/diag_scripts/examples/diagnostic.py for an example of how to
    start your diagnostic.

    """
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
        '-i',
        '--ignore-existing',
        help=("Force running the script, even if output files exists."
              "(useful when re-running the script, use at your own risk)"),
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
        elif not args.ignore_existing:
            logger.error(
                "Script will abort to prevent accidentally overwriting your "
                "data in these directories:\n%s\n"
                "Use -f or --force to force emptying the output directories "
                "or use -i or --ignore-existing to ignore existing output "
                "directories.", '\n'.join(existing))

    for output_directory in output_directories:
        logger.info("Creating %s", output_directory)
        if args.ignore_existing and os.path.exists(output_directory):
            continue
        os.makedirs(output_directory)

    yield cfg

    logger.info("End of diagnostic script run.")
