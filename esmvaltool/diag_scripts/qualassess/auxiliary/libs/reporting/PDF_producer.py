#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 14:49:51 2018

@author: bmueller
"""

import subprocess
import os
import shutil

dir2sources = "." + os.sep

source_list = [x[0] for x in os.walk(dir2sources) if "source" in x[0]]

print("PDF production started.")

for source in source_list:

    print "processing " + source

    build = source.replace("source", "build")

    os.environ['SOURCEDIR'] = source
    os.environ['BUILDDIR'] = build

    try:
        subprocess.call("make latexpdf", shell=True)

        # move pdf to the parent directory and rename to report_xxx.pdf
        parent_dir = os.path.dirname(os.path.normpath(source))
        report_title = "report_" + ("_").join(source.split("_")[-3:-1])
        pdfname = parent_dir + os.sep + report_title + ".pdf"
        shutil.move(
            build +
            os.sep +
            "latex" +
            os.sep +
            "ESMValToolC3S_511Report.pdf",
            pdfname.strip(os.sep))
    except BaseException:
        print("processing " + source + " not possible. Please investigate!")

print("PDF production finished.")
