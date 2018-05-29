"""
Initialisation code for AutoAssess.
"""
from __future__ import print_function
import argparse
import os.path
import sys


def create_dir(path, safe=False):
    """
    Make a dir
    """
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
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
