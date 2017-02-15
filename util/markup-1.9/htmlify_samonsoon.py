#! /usr/bin/env python

import markup
import htmlify
import os
import pdb
import re
import sys
import shutil
from operator import attrgetter
from optparse import OptionParser
from markup import oneliner as e

project_type = re.compile("(CMIP5|OBS)")
svn_url = "https://svn.dlr.de/ESM-Diagnostic/source/branches/EMBRACE/WAMonsoon/"
diagnostic_folder = "diag_scripts/"

# Check command arguments.
usage = "%prog <figure_root_path>"
description = """htmlify
Arranges a set of figures in root-path in a
nice html structure for overview browsing. 
Usually run from the top figure folder"""


parser = OptionParser(usage=usage, description=description)
parser.add_option("-t", "--target-folder",
                  default="normal",
                  dest="target_folder",
                  help="target folder for html files, use any valid folder name")
parser.add_option("-l", "--link-folder",
                  dest="link_folder",
                  default=[],
                  action="append",
                  help="folders to link from target html index-file."\
                  + "These folders should be created as separate targets.")
options, args = parser.parse_args()

htmlify.mkdir_p(options.target_folder)
html_file = os.path.join(options.target_folder, "index.html")

if len(args) == 0:
    parser.print_help()
    sys.exit(0)

fig_path = args[0]
run_path = os.path.dirname(os.path.realpath(sys.argv[0]))
diagnostics = [dir for dir in os.listdir(fig_path)
               if os.path.isdir(os.path.join(fig_path, dir))]
diagnostics = [dir for dir in diagnostics if dir not in ['normal', 'debug', 'netcdf'] ]

plots = []
# Read and (very) simple parsing of figure filenames
for diag in diagnostics:
    root_path = os.path.join(os.path.abspath(fig_path), diag)
    panels = htmlify.Panels(options.target_folder)
    singles = htmlify.Singles(options.target_folder)
    for path, dirs, files in os.walk(root_path):
        if re.search('netcdf', path):
            continue
        for file in files:
            if project_type.search(file) is None:
                panels.append(htmlify.Panel(diag, file))
            else:
                singles.append(htmlify.Single(diag, file))

        plots.append(htmlify.Diags(root_path, diag, panels, singles))

page = markup.page(mode='strict_html')
shutil.copy(os.path.join(run_path, 'htmlify.css'), options.target_folder)
page.init(css=os.path.join("..", 'htmlify.css'))

# Sort plots
plots = sorted(plots, key=attrgetter('diag'))
for plot in plots:
    plot.panels.sort('descr')

# Write html with figures
page.h1("Diagnostics", id="the_top")
for diag in plots:
    page.ul()
    page.li(e.a(diag.diag, href='#' + diag.diag))
    page.ul.close()

for diag in plots:
    # Work with panel plots
    for panel in diag.panels.files:
        page.div(class_='diag', id=diag.diag)
        page.img(src=os.path.join("..", fig_path, diag.diag, panel.fig),
                 id=panel.fig,
                 class_="figure")

        # Create single plots
        single_page = markup.page(mode='strict_html')
        single_page.init(css=os.path.join(fig_path, 'htmlify.css'))
        remove_err_vars = [s for s in diag.singles.files if s.var == panel.var]
        htmlify.write_single(options.target_folder, 
                             options.link_folder, 
                             diag.diag, 
                             remove_err_vars,
                             panel, 
                             single_page, 
                             run_path, 
                             fig_path, 
                             panel.fig)

        htmlify.link_to_single(page, remove_err_vars, panel)

        page.a("Link to diag script",
               href=svn_url + diagnostic_folder + diag.diag + ".ncl",
               class_='link')
        page.add(" (DLR password required)")
        page.br()
        for link in options.link_folder:
            rewrite_linked_filename = getattr(vars()['htmlify'], options.target_folder + "_" + link)()
            link_file = rewrite_linked_filename.switch(options.target_folder,
                                                       link,
                                                       panel.fig)
            page.a("Link to " + link + " version",
            href=os.path.join("..", link, "index.html") + "#" + link_file,
            class_='link')
            page.br()
        page.a("Back to the top",
               href='#the_top',
               class_='link',
               style='display:block;width:auto;')
        page.div.close()

for diag in plots:
    if len(diag.panels.files) == 0:
        for single in diag.singles.files:
            page.div(class_='diag_single', id=diag.diag)
            page.img(src=os.path.join("..", fig_path, diag.diag, single.fig),
                     id=single.fig,
                     class_="figure")
            page.a("Link to diag script",
                   href=svn_url + diagnostic_folder + diag.diag + ".ncl",
                   class_='link')
            page.add(" (DLR password required)")
            page.br()
            page.a("Back to the top",
                   href='#the_top',
                   class_='link',
                   style='display:block;width:auto;')
            page.div.close()

f = open(html_file, "w")
for line in page.content:
    f.write("%s\n" % line)
f.close()
