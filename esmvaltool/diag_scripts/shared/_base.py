"""Convenience functions for running a diagnostic script."""
import argparse
import contextlib
import glob
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import yaml

logger = logging.getLogger(__name__)


def get_plot_filename(basename, cfg):
    """Get a valid path for saving a diagnostic plot.

    Parameters
    ----------
    basename: str
        The basename of the file.
    cfg: dict
        Dictionary with diagnostic configuration.

    Returns
    -------
    str:
        A valid path for saving a diagnostic plot.
    """
    return os.path.join(
        cfg['plot_dir'],
        f"{basename}.{cfg['output_file_type']}",
    )


def get_diagnostic_filename(basename, cfg, extension='nc'):
    """Get a valid path for saving a diagnostic data file.

    Parameters
    ----------
    basename: str
        The basename of the file.
    cfg: dict
        Dictionary with diagnostic configuration.
    extension: str
        File name extension.

    Returns
    -------
    str:
        A valid path for saving a diagnostic data file.
    """
    return os.path.join(
        cfg['work_dir'],
        f"{basename}.{extension}",
    )


def save_figure(basename, provenance, cfg, figure=None, close=True, **kwargs):
    """Save a figure to file.

    Parameters
    ----------
    basename: str
        The basename of the file.
    provenance: dict
        The provenance record for the figure.
    cfg: dict
        Dictionary with diagnostic configuration.
    figure: matplotlib.figure.Figure
        Figure to save.
    close: bool
        Close the figure after saving.
    **kwargs:
        Keyword arguments to pass to :obj:`matplotlib.figure.Figure.savefig`.

    See Also
    --------
    ProvenanceLogger: For an example provenance record that can be used
        with this function.
    """
    if cfg.get('output_file_type') is None:
        extensions = ('png', 'pdf')
    elif isinstance(cfg['output_file_type'], str):
        extensions = (cfg['output_file_type'], )
    else:
        extensions = cfg['output_file_type']

    for ext in extensions:
        filename = Path(cfg['plot_dir']) / ext / f"{basename}.{ext}"
        filename.parent.mkdir(exist_ok=True)
        logger.info("Plotting analysis results to %s", filename)
        fig = plt if figure is None else figure
        fig.savefig(filename, **kwargs)
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(filename, provenance)

    if close:
        plt.close(figure)


def save_data(basename, provenance, cfg, cube, **kwargs):
    """Save the data used to create a plot to file.

    Parameters
    ----------
    basename: str
        The basename of the file.
    provenance: dict
        The provenance record for the data.
    cfg: dict
        Dictionary with diagnostic configuration.
    cube: iris.cube.Cube
        Data cube to save.
    **kwargs:
        Extra keyword arguments to pass to :obj:`iris.save`.

    See Also
    --------
    ProvenanceLogger: For an example provenance record that can be used
        with this function.
    """
    if 'target' in kwargs:
        raise ValueError(
            "Please use the `basename` argument to specify the output file")

    filename = get_diagnostic_filename(basename, cfg)
    logger.info("Saving analysis results to %s", filename)
    iris.save(cube, target=filename, **kwargs)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


class ProvenanceLogger:
    """Open the provenance logger.

    Parameters
    ----------
    cfg: dict
        Dictionary with diagnostic configuration.

    Example
    -------
        Use as a context manager::

            record = {
                'caption': "This is a nice plot.",
                'statistics': ['mean'],
                'domain': ['global'],
                'plot_type': ['zonal'],
                'authors': [
                    'first_author',
                    'second_author',
                ],
                'references': [
                    'author20journal',
                ],
                'ancestors': [
                    '/path/to/input_file_1.nc',
                    '/path/to/input_file_2.nc',
                ],
            }
            output_file = '/path/to/result.nc'

            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, record)
    """
    def __init__(self, cfg):
        """Create a provenance logger."""
        self._log_file = os.path.join(cfg['run_dir'],
                                      'diagnostic_provenance.yml')

        if not os.path.exists(self._log_file):
            self.table = {}
        else:
            with open(self._log_file, 'r') as file:
                self.table = yaml.safe_load(file)

    def log(self, filename, record):
        """Record provenance.

        Parameters
        ----------
        filename: str
            Name of the file containing the diagnostic data.
        record: dict
            Dictionary with the provenance information to be logged.

            Typical keys are:
                - ancestors
                - authors
                - caption
                - domain
                - plot_type
                - references
                - statistics

        Note
        ----
            See the provenance `documentation`_ for more information.

        .. _documentation: https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance
        """  # noqa
        if isinstance(filename, Path):
            filename = str(filename)
        if filename in self.table:
            raise KeyError(
                "Provenance record for {} already exists.".format(filename))

        self.table[filename] = record

    def _save(self):
        """Save the provenance log to file."""
        dirname = os.path.dirname(self._log_file)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(self._log_file, 'w') as file:
            yaml.safe_dump(self.table, file)

    def __enter__(self):
        """Enter context."""
        return self

    def __exit__(self, *_):
        """Save the provenance log before exiting context."""
        self._save()


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
        if all(a in attribs and (
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
        A dictionary containing the requested groups.
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
    :obj:`dict` of :obj:`list` of :obj:`dict`
        A dictionary containing the requested groups.
    """
    if sort is True:
        sort = []

    def normalized_group_key(key):
        """Define a key to sort by."""
        return '' if key is None else str(key).lower()

    groups = {}
    for key in sorted(metadata_groups, key=normalized_group_key):
        groups[key] = sorted_metadata(metadata_groups[key], sort)

    return groups


def extract_variables(cfg, as_iris=False):
    """Extract basic variable information from configuration dictionary.

    Returns `short_name`, `standard_name`, `long_name` and `units` keys for
    each variable.

    Parameters
    ----------
    cfg : dict
        Diagnostic script configuration.
    as_iris : bool, optional
        Replace `short_name` by `var_name`, this can be used directly in
        :mod:`iris` classes.

    Returns
    -------
    dict
        Variable information in :obj:`dict`s (values) for each `short_name`
        (key).
    """
    keys_to_extract = [
        'short_name',
        'standard_name',
        'long_name',
        'units',
    ]

    # Extract variables
    input_data = cfg['input_data'].values()
    variable_data = group_metadata(input_data, 'short_name')
    variables = {}
    for (short_name, data) in variable_data.items():
        data = data[0]
        variables[short_name] = {}
        info = variables[short_name]
        for key in keys_to_extract:
            if key in data:
                info[key] = data[key]

        # Replace short_name by var_name if desired
        if as_iris:
            info['var_name'] = info.pop('short_name')
            if info['standard_name'] == '':
                info['standard_name'] = None

    return variables


def variables_available(cfg, short_names):
    """Check if data from certain variables is available.

    Parameters
    ----------
    cfg : dict
        Diagnostic script configuration.
    short_names : list of str
        Variable `short_names` which should be checked.

    Returns
    -------
    bool
        `True` if all variables available, `False` if not.
    """
    input_data = cfg['input_data'].values()
    available_short_names = list(group_metadata(input_data, 'short_name'))
    for var in short_names:
        if var not in available_short_names:
            return False
    return True


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
    """Run a Python diagnostic.

    This context manager is the main entry point for most Python diagnostics.

    Example
    -------
    See esmvaltool/diag_scripts/examples/diagnostic.py for an extensive
    example of how to start your diagnostic.

    Basic usage is as follows, add these lines at the bottom of your script::

        def main(cfg):
            # Your diagnostic code goes here.
            print(cfg)

        if __name__ == '__main__':
            with run_diagnostic() as cfg:
                main(cfg)

    The `cfg` dict passed to `main` contains the script configuration that
    can be used with the other functions in this module.
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
    logging.captureWarnings(True)
    logging.getLogger().setLevel(cfg['log_level'].upper())

    # Read input metadata
    cfg['input_data'] = _get_input_data_files(cfg)

    logger.info("Starting diagnostic script %s with configuration:\n%s",
                cfg['script'], yaml.safe_dump(cfg))

    # Clean run_dir and output directories from previous runs
    default_files = {
        'diagnostic_provenance.yml',
        'log.txt',
        'profile.bin',
        'resource_usage.txt',
        'settings.yml',
    }

    output_directories = (cfg['work_dir'], cfg['plot_dir'])
    old_content = [
        p for p in output_directories
        if Path(p).exists() and any(Path(p).iterdir())
    ]
    old_content.extend(p for p in glob.glob(f"{cfg['run_dir']}{os.sep}*")
                       if not os.path.basename(p) in default_files)

    if old_content:
        if args.force:
            for content in old_content:
                logger.info("Removing %s", content)
                if os.path.isfile(content):
                    os.remove(content)
                else:
                    shutil.rmtree(content)
        elif not args.ignore_existing:
            raise FileExistsError(
                "Script will abort to prevent accidentally overwriting "
                "your data in the following output files or directories:"
                "\n%s\n Use -f or --force to force emptying the output "
                "directories or use -i or --ignore-existing to ignore "
                "existing output directories." % '\n'.join(old_content))

    # Create output directories
    for output_directory in output_directories:
        if not os.path.isdir(output_directory):
            logger.info("Creating %s", output_directory)
            os.makedirs(output_directory)

    provenance_file = os.path.join(cfg['run_dir'], 'diagnostic_provenance.yml')
    if os.path.exists(provenance_file):
        logger.info("Removing %s from previous run.", provenance_file)
        os.remove(provenance_file)

    yield cfg

    logger.info("End of diagnostic script run.")
