#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:55:26 2017

@author: bmueller
"""

from ESMValMD import nclFileMD
from optparse import OptionParser
import sys

usage = "%prog Metadata_file.txt"
description = """Wrapper to write METAdata with ncl scripts!"""


parser = OptionParser(usage=usage, description=description)
options, args = parser.parse_args()
if len(args) != 1:
    parser.print_help()
    sys.exit(0)

nclFileMD(args[0])

# The txt-file looks as follows:

# ./mobile_home_flying.jpg            | = name of file to be changed
# both                                | = type of METAdata (both, xml, meta)
# mobile,home,flying                  | = tags, comma seperated
# This is a flying mobile home.       | = caption in plain text
# #IDtest                             | = Block ID, for report structuring
