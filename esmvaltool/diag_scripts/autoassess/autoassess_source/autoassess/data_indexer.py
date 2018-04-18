#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Command line tool to recursively load all PP files in the given directory with
Iris, and save the pickled Iris cube list in the given directory. Also important
meta data is added to the cube attributes.

By providing `Iris cubes` data is made available in a file format independent
way. As the loading process can take a lot of time for files without a header
describing the entire file (for example PP files), the result is save as
Pickled string in `cubes.pickle` inside the top directory.

Currently limited to `*.pp` files.
"""

import argparse
import fnmatch
import os
import os.path
import sys

import iris

try:
    import cPickle as pickle
except ImportError:
    import pickle


class ExistingAttributeValue(Exception):
    pass


def parse_args(args):
    """
    Parse arguments in a function to facilitate testing. Contains all command
    line options.

    :param list args: Command line arguments from sys.argv.
    :returns: Checked command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', required=True, help='Directory tree with PP files.')
    parser.add_argument('--suite-id', required=True, help='Model run ID to be added to cube attributes.')
    args = parser.parse_args(args)

    return args


def iris_load_dir(data_dir):
    """
    Load recursively all PP files in `data_dir` with Iris.
    Load all files in a directory into a single cube, thereby avoid the merging
    of different meaning periods into a single cube.

    :param str data_dir: Path of directory with climate data.
    :returns: Iris cube list of all fields in `data_dir`.
    :rtype: Iris cube list.
    """
    cube_list = iris.cube.CubeList()
    for dirpath, dirnames, filenames in os.walk(data_dir):  # pylint: disable=unused-variable
        uris = [os.path.join(dirpath, f) for f in fnmatch.filter(filenames, '*.pp') +
                                                  fnmatch.filter(filenames, '*.PP')]
        cube_list.extend(iris.load(uris))
    return cube_list


def pickle_to_dir(data_dir, cube_list):
    """
    Pickle Iris cube list to `cubes.pickle` in `data_dir`.

    :param str data_dir: Directory path.
    :param obj cube_list: Python object to be pickled.
    """
    with open(os.path.join(data_dir, 'cubes.pickle'), 'w') as fh:
        pickle.dump(cube_list, fh)


def unpickle_from_dir(data_dir):
    """
    Un-Pickle `cubes.pickle` in `data_dir`.

    :param str data_dir: Directory path.
    :returns: Iris cube list.
    """
    with open(os.path.join(data_dir, 'cubes.pickle'), 'r') as fh:
        cube_list = pickle.load(fh)
    return cube_list


def add_suite_id_to_cube(cubes, suite_id):
    """Add given key, value pairs to the attributes of each cube in cube list.

    :param cubes: Iris cube list
    :param suite_id: Model run suite ID
    :returns: Iris cube list with new MODEL_RUN_ID attribute, and ammended
        history for each cube.
    """
    new_attrib = 'MODEL_RUN_ID'
    for cube in cubes:
        if new_attrib in cube.attributes and cube.attributes[new_attrib] != suite_id:
            raise ExistingAttributeValue('Attribute: {}, Existing: {}, New: {}'.format(
                new_attrib,
                cube.attributes[new_attrib],
                suite_id))
        else:
            cube.attributes[new_attrib] = suite_id
            # update history
            if not 'history' in cube.attributes:
                cube.attributes['history'] = ''

            cube.attributes['history'] += 'Attribute MODEL_RUN_ID added by AutoAssess.\n'
    return cubes


def data_indexer():
    args = parse_args(sys.argv[1:])
    cube_list = iris_load_dir(args.data_dir)
    cube_list = add_suite_id_to_cube(cube_list, args.suite_id)
    pickle_to_dir(args.data_dir, cube_list)

data_indexer.__doc__ = __doc__  # use module docstring


if __name__ == '__main__':
    data_indexer()
