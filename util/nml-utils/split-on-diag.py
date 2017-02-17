#! /usr/bin/env python

# Parse namelist and split it on diags (=independent chunks)
#

import glob
import os
import pdb
import re
import xml.sax
import library_nml

from argparse import ArgumentParser

import sys


# Check command arguments.
description = """split-on-diag
Rewrites each 'diag'-section in a namelist to separate namelists
This is under the assumption that they are independent and can be
parallelized.
"""

parser = ArgumentParser(description=description)
parser.add_argument("-n", "--nml",
                    dest="namelist",
                    type=str,
                    required=True,
                    help="Namelist to split")
parser.add_argument("-d", "--target-directory",
                    dest="target_directory",
                    type=str,
                    default="tmp/",
                    help="Director to write output namelists")
args = parser.parse_args()


# Parse input namelist into project_info-dictionary.
namelist = library_nml.splitDiags()
parser = xml.sax.make_parser()
parser.setContentHandler(namelist)
parser.parse(args.namelist)
nml_sections = namelist.namelist_sections

os.makedirs(args.target_directory)
nml_name = os.path.basename(re.search("(.*).xml", args.namelist).group(1))
for num, diag in enumerate(nml_sections['diags']):
    nml_out = os.path.join(args.target_directory,
                           nml_name + "-" + str(num).zfill(2) + ".xml")
    fnml = open(nml_out, "w")
    fnml.write(nml_sections['header'].encode('utf8'))
    fnml.write(diag.encode('utf8'))
    fnml.write(nml_sections['footer'].encode('utf8'))
