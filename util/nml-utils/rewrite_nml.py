#! /usr/bin/env python

# Parse namelist and split it on diags (=independent chunks)
#

import glob
import os
import pdb
import re
import shutil
import xml.sax
import library_nml

from argparse import ArgumentParser

import sys


# Check command arguments.
description = """rewrite-nml.py
Wrapper to xml.sax to rewrites certain entries in the namelist.
"""

parser = ArgumentParser(description=description)
parser.add_argument("-n", "--nml",
                    dest="namelist",
                    type=str,
                    required=True,
                    help="Input namelist to rewrite")
parser.add_argument("-o", "--output-namelist",
                    dest="output_namelist",
                    type=str,
                    help="Use this as output instead of overwriting input")
parser.add_argument("-m", "--add-model-entry",
                    dest="model_entry",
                    type=str,
                    default=None,
                    help="Add this line (str) as a new <model>-tag in <MODELs>")
parser.add_argument("-v", "--variable",
                    dest="variable",
                    type=str,
                    default=None,
                    help="When adding models, select only diags matching this variable")
parser.add_argument("-r", "--replacements-in-global",
                    dest="repl_dict",
                    type=str,
                    default=None,
                    help="Comma separated assigments to update in GLOBALS section")
args = parser.parse_args()


nml_name = re.search("(.*).xml", args.namelist).group(1)
nml_tmp = nml_name + "-tmp.xml"
nml_target = nml_name + ".xml"

if args.output_namelist is not None:
    nml_target = args.output_namelist

# Add models to <MODELS>-tags
if args.model_entry is not None and args.variable is None:
    models_start = '<MODELS>' + '\n'
    models_end = '</MODELS>' + '\n'

    namelist = library_nml.addModels()
    parser = xml.sax.make_parser()
    parser.setContentHandler(namelist)
    parser.parse(args.namelist)
    nml_sections = namelist.namelist_sections

    shutil.copy(nml_name + '.xml', nml_name + '-backup.xml')
    ftmp = open(nml_tmp, "w")
    ftmp.write(nml_sections['header'].encode('utf8'))
    ftmp.write(models_start.encode('utf8'))

    ftmp.write(nml_sections['MODELS'] + '\n'
               + '<model> '
               + args.model_entry
               + '</model>'
               + '\n'.encode('utf8'))

    ftmp.write(models_end.encode('utf8'))
    ftmp.write(nml_sections['footer'].encode('utf8'))
    ftmp.close()
    shutil.move(nml_tmp, nml_target)

# Add models to specific diag-entries
if args.model_entry is not None and args.variable is not None:
    namelist = library_nml.addModel(args.variable)
    parser = xml.sax.make_parser()
    parser.setContentHandler(namelist)
    parser.parse(args.namelist)
    nml_sections = namelist.namelist_sections

    shutil.copy(nml_name + '.xml', nml_name + '-backup.xml')
    ftmp = open(nml_tmp, "w")
    if len(nml_sections['header']) > 0:
        ftmp.write(nml_sections['header'].encode('utf8'))
    
        ftmp.write('\n'
                   + '    <model> '
                   + args.model_entry
                   + '   </model>'.encode('utf8'))
    ftmp.write(nml_sections['footer'].encode('utf8'))
    ftmp.close()
    shutil.move(nml_tmp, nml_target)

# set global settings
if args.repl_dict is not None:
    replacement_dict = dict(item.split("=") for item in
                            [e.strip() for e in args.repl_dict.split(",")])
    namelist = library_nml.setGlobal(replacement_dict)
    parser = xml.sax.make_parser()
    parser.setContentHandler(namelist)
    parser.parse(args.namelist)

    shutil.copy(nml_name + '.xml', nml_name + '-backup.xml')
    ftmp = open(nml_tmp, "w")
    ftmp.write(namelist.str.encode('utf8'))
    ftmp.close()
    shutil.move(nml_tmp, nml_target)
