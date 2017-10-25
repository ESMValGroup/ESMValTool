#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:55:26 2017

@author: bmueller
"""

from ESMValMD import ESMValMD
from optparse import OptionParser
import sys

usage = "%prog filename [Options] "
description = """Wrapper to write METAdata with ncl scripts!"""

parser = OptionParser(usage=usage, description=description)

parser.add_option("-d", "--dtype", dest="datatype", default="both",
                  help="type of METAdata (both, xml, meta)")
parser.add_option("-t", "--tags", dest="tags", default="",
                  help="tags for METAdata")
parser.add_option("-c", "--cap", dest="caption", default="NO CAPTION.",
                  help="caption for METAdata")
parser.add_option("-b", "--block", dest="block", default="#ID0",
                  help="block ID for report structure")

options, args = parser.parse_args()

if len(args) != 1:
    parser.print_help()
    sys.exit(0)

ESMValMD(options.datatype,
         args[0],
         options.tags.split(","),
         options.caption,
         options.block)
