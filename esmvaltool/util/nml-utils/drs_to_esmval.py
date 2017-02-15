#! /usr/bin/env python

# Recursively traverse from the given root dir and rewrite
# any matched CMIP5-style files to ESMValTool <model>-tags

import os
import pdb
import re
import sys
from argparse import ArgumentParser

pattern_separator = "_"


def fname2modeltag(fname, fdir):
    # Fix possible (and erroneous) underscore instead of dash in years
    fname = re.sub("([0-9]{4}[0-9]+)_([0-9]{4}[0-9]+).nc", r"\1-\2.nc", fname)
    regex_yrs = re.compile("([0-9]{4})[0-9]+-([0-9]{4})[0-9]+.nc")
    cmip5_drs = ["var", "mip", "name", "exp", "ens", "yrs"]
    fname_dict = dict(zip(cmip5_drs, fname.split("_")))
    fname_dict['start_year'] = regex_yrs.search(fname_dict["yrs"]).group(1)
    fname_dict['end_year'] = regex_yrs.search(fname_dict["yrs"]).group(2)

    entries = ["name", "mip", "exp", "ens", "start_year", "end_year"]
    rearranged_entries = [fname_dict[key] for key in entries]
    rearranged_entries.append(fdir + "  </model>")
    rearranged_entries.insert(0, "<model> CMIP5")
    return "  ".join(rearranged_entries)

# Check command arguments.
description = """drs-to-esmval
Recursively searches the given <data_root_path> for matching entries
and rewrites these to ESMValTool <model>-tags."""

parser = ArgumentParser(description=description)
parser.add_argument("-d", "--data-root-folder",
                    dest="root_folders",
                    type=str,
                    nargs='+',
                    required=True,
                    help="root folder for CMIP input files. Set to 'stdin' to read from stdin")
parser.add_argument("-v", "--variable",
                    default=".*?",
                    help="variable to match for")
parser.add_argument("-m", "--mip",
                    default=".*?",
                    help="mip to match for")
parser.add_argument("-n", "--name",
                    default=".*?",
                    help="model name to match for")
parser.add_argument("-x", "--experiment",
                    default=".*?",
                    help="experiment to match for")
parser.add_argument("-e", "--ensemble",
                    default=".*?",
                    help="ensemble to match for")
args = parser.parse_args()

keys = ["variable", "mip", "name", "experiment", "ensemble"]
entries = dict([(key, args.__dict__[key]) for key in keys])

# Replace comma separated entires, '<entry1>,<entry2>', with
# '(<entry1>|<entry2>)' to fit with re-expressions
for key in entries.keys():
    entries[key] = '(' + '|'.join(entries[key].split(',')) + ')'

pattern_to_match = pattern_separator.join([entries[key] for key in keys]) + '.*.nc'

# Read list of full paths to CMIP5-files from stdin
if args.root_folders == "stdin":
    regex_stdin = re.compile("(.*)/(.*nc)")
    for stdin in sys.stdin.read().split('\n'):
        if len(stdin) == 0:
            continue
        dirName = regex_stdin.search(stdin).group(1)
        fname = regex_stdin.search(stdin).group(2)
        if re.search(pattern_to_match, fname) is not None:
            modeltag = fname2modeltag(fname, dirName)
            print modeltag

# Read list of full paths to CMIP5-files by os.walk on 'root_dir'
else:
    for root_folder in args.root_folders:
        for dirName, subdirList, fileList in os.walk(root_folder, topdown=False):
            for fname in fileList:
                if re.search(pattern_to_match, fname) is not None:
                    modeltag = fname2modeltag(fname, dirName)
                    print modeltag
