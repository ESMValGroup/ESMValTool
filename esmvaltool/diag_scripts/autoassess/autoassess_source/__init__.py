"""
Initialisation code for AutoAssess.
"""
from __future__ import print_function
import argparse
import os.path
import sys


__version__ = '2018.05.24'
__author__ = 'Met Office ported to ESMValTool'
__license__ = 'BSD 3-clause license'
__copyright__ = 'Crown Copyright 2017, the Met Office'

#######################################################

_ROOT = os.path.abspath(os.path.dirname(__file__))  # path of this file


def query_template_dir():
    """Return path to MOOSE query templates."""
    path = os.path.join(_ROOT, 'moose_queries')
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def _area_dir(area):
    """Return path to the directory of an assessment area."""
    path = os.path.join(_ROOT, 'assessment_areas', area)
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def metrics_from_obs_dir(area):
    """Return path to directory with metrics calculated from observations."""
    area = area.lower()
    path = os.path.join(_area_dir(area), 'assessment_plot_obs')
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def metric_selections_dir(area):
    """Return path to directory with metrics selections."""
    area = area.lower()
    path = os.path.join(_area_dir(area), 'assessment_plot_selections')
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def area_page_template_dir(area):
    """Return path to directory with area page template."""
    area = area.lower()
    path = os.path.join(_area_dir(area), 'assessment_page_template')
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def overview_page_template_dir():
    """Return path to directory with overview page template."""
    path = os.path.join(_ROOT, 'overview_page_template')
    if not os.path.isdir(path):
        raise Exception('Directory {0} not found.'.format(path))
    return path


def create_dir(path, safe=False):
    """
    Make sure directory exists. If safe=True directory must not exist already.
    ``create_dir`` is equivalent to os.makedirs(path, exists_ok=True) in
    Python3.

    :param str path: Full path to be created.
    :param bool safe: If True raise OSError if path is existing directory.
    :raises: OSError if cannot make directory at path, and path is not an
        existing directory.
             OSError if safe=True and directory exists at path.
    """
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            if safe:
                msg = (
    '######################################################################\n'
    'In order to re-run the task, the already existing directory\n'
    'has to be deleted.\n'
    '######################################################################\n'
                    )
                print(msg, file=sys.stderr)
                # VPREDOI
                pass
            else:
                pass
        else:
            raise


def argparser(description, pos_args, opt_args, short=False):
    '''
    Generate an argument parser based on given configuration.

    The pos_args and opt_args arguments are lists with each entry a tuple
    containing the name of each argument, and a dictionary defining the
    keywords to be used in ``argparse.add_argument()`` for the argument.

    :param str desc: Description for the argument parser, usually __doc__
    :param list pos_args: Postional argument configuration
    :param list opt_args: Optional argument configuration
    :param bool short: Include shortened single letter optional arguments
    :returns: Parser object
    :rtype: argparse.ArgumentParser
    '''

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Add positional arguments to the parser
    for (arg, kwargs) in pos_args:
        parser.add_argument(arg, **kwargs)

    # Add optional arguments to the parser
    # Option to have shortened single letter optional arguments
    for (arg, kwargs) in opt_args:
        if short:
            parser.add_argument('-'+arg[0], '--'+arg, **kwargs)
        else:
            parser.add_argument('--'+arg, **kwargs)

    # Returned parsed options
    return parser
