#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
report_cleanup

* cleans directories that were produced by a former reporting job
* is run before a new report is started
* run "python report_cleanup.py" for manual cleanup

Created on Fri May 12 10:13:39 2017

@author: bmueller
"""

from path import path


def cleanup_report(template_dir, image_dir, do_print):
    def remove_files(l_dir, pattern):
        """
        removes files in dir with pattern
        """
        files = l_dir.walkfiles(pattern)
        for l_file in files:
            l_file.remove()
            if do_print:
                print "Removed {} file".format(l_file)

    # cleanup html files
    d = path(template_dir)
    remove_files(d, "*.html")

    # cleanup image directories
    d = path(image_dir)
    for i in d.walk():
        if i.isdir():
            remove_files(i, "*")
            i.rmdir()
            if do_print:
                print "Removed {} directory".format(i)


if __name__ == '__main__':
    do_print = True
    template_dir = '../templates/'
    image_dir = '../static/images/'
    cleanup_report(template_dir, image_dir, do_print)
