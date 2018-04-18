#! /usr/bin/env python2.7

import copy
import datetime
import glob
import gzip
import os
import pdb
import re
import shutil

import markup
import globalvar
import valmod as vm

import ftp_valnote_mod
import make_means
import speed
import summary


def valnote_top(css, title, jscript, links):
    page = markup.page()
    css = os.path.basename(css)
    # Concern: the order of the scripts is important (showplot.js must come
    # before plots.js), but dictionary has no guarenteed order.  Luckily,
    # current result is correct, but could change.  If it does, solution might
    # be to hack markup.py, to put plots.js script inside html page, or to
    # make scripts order-independent (if possible).

    page.init(title=title, css=css,
              script={'plots.js': 'javascript', jscript: 'javascript'},
              bodyattrs={'onLoad': 'updateFields()'},
              footer=datestamp())

    page.h1(title)

    page.table(class_='mivnchead')
    page.tr()
    for link, text in links:
        page.td(markup.oneliner.a(text, href=link))
    page.tr.close()
    page.table.close()

    return page


def datestamp():
    today = datetime.date.today()
    datestamp = today.strftime("%A %d %B %Y")
    stamp = 'Validation Note created by valnote.py on <b>%s</b>' % datestamp
    return stamp


def valnote_table(page_top, valorder_file, outfile, strict=False):
    '''
    Make a plot table

    page_top = (markup obj) Information on the top of the html page
    outfile = (str) the output filename for this page
    valnote_dir = (str) the valnote directory (useful for finding control files)
    strict = (bool) only show the plots that are in the valorder file
    '''

    # Initialise variables
    ncols = 2   # Number of columns
    i = -1      # Current column

    # Initialise the page
    page = copy.deepcopy(page_top)
    page.table(class_="mivnctable")

    # Generate a list of thumbnail plots to show
    thumb_list = glob.glob('thumbs/*.png')

    # Load the valorder file (needed for plot titles)
    valorder_dict = vm.read_info_file(valorder_file)

    for thumb in sorted(thumb_list):

        foundplot = True

        # Deduce the filename of the large image
        thumb_split = thumb.split('/')
        filename = thumb_split[-1][:-4]
        fileend = '.png'
        if globalvar.pub:
            fileend = '.eps'
        plot = os.path.join('images', filename+fileend)

        # Get the plot title from the valorder file
        filename_split = filename.split('_')
        plot_number_str = filename_split[0]
        try:
            plot_title = valorder_dict[plot_number_str]
        except KeyError:

            if strict:
                foundplot=False
            else:
                # If the plot title was not in the valorder file simply use the file name without underscores or .png
                plot_title = filename.replace('.png','').replace('_',' ')

        # If we have found the plot in the valorder file or strict is False then add plot to the page
        if foundplot:

            # Increment the column and decide if we need to start or end a row
            i += 1
            startrow = i % ncols == 0
            endrow = (i+1) % ncols == 0
            if startrow:
                page.tr()

            # Add the plot to the page
            page.td()
            page.p(plot_title)
            img = markup.oneliner.img(src=thumb)
            page.a(img, href=plot)
            page.td.close()

            # End the row if needed
            if endrow:
                page.tr.close()

    # Close the table
    page.table.close()

    # Output to file
    with open(outfile, 'w') as ofile:
        print >>ofile, page


def valnote_browser(page_top, outfile):

    page = copy.deepcopy(page_top)

    page.form(target="", id="selecta")

    page.addcontent("Field:")
    page.select(name="field", id="field", onClick="updateSeasons(this)",
                style="width:150")
    page.select.close()

    page.addcontent("Season:")
    page.select(name="season", id="season", onClick="updateObs(this)",
                style="width:70")
    page.select.close()

    page.addcontent("Obs:")
    page.select(name="obs", id="obs", onClick="updateProc(this)",
                style="width:90")
    page.select.close()

    page.select(name="proc", id="proc",
                style="width:70")
    page.select.close()

    page.input(type="button", value="Display plot", onclick="showplot(this)")

    page.form.close()

    page.div(id="plot_box")
    page.div.close()

    with open(outfile, 'w') as ofile:
        print >>ofile, page


def write_jsplots():
    '''
    This routine looks at the list of files in the images directory and assigns
    them to particular variables, seasons and observations. This information is
    then put into the plots.js file.
    '''

    # Initialise a list that contains all the plot info
    jsplots = []

    # Get a list of files
    plot_list = glob.glob('images/*.png')

    # Loop through the files
    for (i, plot) in enumerate(sorted(plot_list)):

        # Split the string into two halves. The first contains the number,
        # variable name and season. The second contains the observations.
        if plot.count('_v_') == 1:
            sections_list = plot.split('_v_')
            first_part = sections_list[0]
            second_part = sections_list[1]
        else:
            first_part = plot.split('.png')[0]
            second_part = "none.png"

        # Get the variable name after the first underscore and before the last
        # underscore in this section.
        where_first_under = first_part.find('_')
        where_last_under = first_part.rfind('_')
        variable_name = first_part[where_first_under+1:where_last_under]
        variable_name = variable_name.replace('_', ' ')

        # Get the season which is after the last underscore in this section.
        season_name = first_part[where_last_under+1:]

        # Get the observations which is in the second section
        where_dot = second_part.find('.')
        obs_name = second_part[:where_dot]
        obs_name = obs_name.replace('_', ' ')

        # Append all the relevant information to the jsplots list
        jsplots.append('plots[%i] = new Plot("%s", "%s", "%s", "%s", "%s")\n' %
                       (i, plot, variable_name, season_name, obs_name, "mean"))

    # Add all the information on the plots to the plots.js file
    nplot = len(jsplots)
    with open("plots.js", 'w') as jfile:
        jfile.write("var nplot = %i\n" % (nplot,))
        jfile.write("var plots = new Array(nplot)\n")
        jfile.writelines(jsplots)


def create_valnote(outdir, images_dir, thumbs_dir, valnote_dir, title,
                   docurl=None, exper=None, control=None, csv_dir=None):
    '''
    docurl is the URL containing the documentation for the validation note
    title is the title of the validation note

    Optimized for Scientific Desktops with 1280x1024 pixel screens
    '''

    # Change the current working directory to outdir
    os.chdir(outdir)

    # Setup file names
    table = "table.html"
    browser = "browser.html"
    meansfile = "means_global.html"
    css = "valnote.css"
    jscript = "showplot.js"
    if docurl is None:
        docurl = "file.html"

    # Setup path of supplimental html files
    html_dir = os.path.join(valnote_dir, 'html')
    html_dir = os.path.abspath(html_dir)

    # Setup array of text files
    html_aux = [os.path.join(html_dir, css), os.path.join(html_dir, jscript)]

    # Copy text files to outdir
    for hfile in html_aux:
        try:
            shutil.copy(hfile, outdir)
        except IOError:
            print "Warning: failed to copy %s to %s" % (hfile, outdir)

    print "  Creating validation note in", outdir

    # Convert the summary tables from csv to html
    # Return the html filename of the first summary page
    summary_page = "summary_global.html"

    # Choose the links to appear in the header table
    links = [(browser, "Plot browser"),
             (table, "All plots"),
             (summary_page, "Summary table"),
             (docurl, "Documentation"),
             (meansfile, "Global means"),
             ('speed.html', "Speed test")]

    valorderfile_short = globalvar.valorder_file.replace('.dat', '_short.dat')
    if os.path.isfile(valorderfile_short):
        table_short = table.replace('.html', '_short.html')
        links.insert(2, (table_short, "Core plots"))
    else:
        table_short = None
        
    # Make a link to the land surface list
    valorderfile_land = globalvar.valorder_file.replace('.dat','_land.dat')
    if os.path.isfile(valorderfile_land):
        table_land = table.replace('.html','_land.html')
        links.insert(3, (table_land, "Land surface") )
    else:
        table_land = None

    # Create the main validation note pages
    page_top = valnote_top(css, title, jscript, links)
    valnote_browser(page_top, browser)
    write_jsplots()
    valnote_table(page_top, globalvar.valorder_file, table)
    if os.path.isfile(valorderfile_short):
        valnote_table(page_top, valorderfile_short, table_short, strict=True)
    else:
        print "Error: Could not find file "+valorderfile_short
    if os.path.isfile(valorderfile_land):
        valnote_table(page_top, valorderfile_land, table_land, strict=True)
    else:
        print "Error: Could not find file "+valorderfile_land
    try:
        speed.create_page(page_top, exper, control)
    except:
        print "Error creating speed page"

    # Create the summary pages
    summary.csv2html(page_top, valorderfile=globalvar.valorder_file,
                     csv_dir=csv_dir, outdir=outdir)

    # Create the means table page
    make_means.csv2html(page_top, csv_dir=csv_dir)

    # Set up URL strings
    bfile = os.path.abspath(browser)
    main = '/spice/project/gmed/valnote'
    main_net = '/net/spice/project/gmed/valnote'
    main_url = 'http://www-hc/~hadvg/valnote'

    # Copy the validation note to the collaboration twiki
    extweb = 'x'
    if bfile.startswith(main) or bfile.startswith(main_net):
        basename = os.path.basename(outdir)
        extweb = ftp_valnote_mod.check_netrc(basename)
    else:
        print 'Not mirroring to collaboration wiki as file path not recognised'
        print 'file = ', bfile

    # Work out the URL of the valnote and print it to screen
    if bfile.startswith(main):
        url = bfile.replace(main, main_url)
    elif bfile.startswith(main_net):
        url = bfile.replace(main_net, main_url)
    else:
        url = 'file://' + bfile
    print
    print "Internal validation note URL is:"
    print url

    if extweb != 'x':
        print
        print "External validation note URL is:"
        print "http://collab.metoffice.gov.uk/twiki/bin/viewfile/Static/"+extweb+"/valid/"+basename+"/browser.html"

# main script for testing
if __name__ == '__main__':
    exper = 'aliur'
    control = 'akkvi'
    outdir = "/project/gmed/valnote/%s_v_%s" % (exper, control)
    images_dir = os.path.abspath(os.path.join(outdir, 'images'))
    thumbs_dir = os.path.abspath(os.path.join(outdir, 'thumbs'))
    valnote_dir = '/home/h03/hadco/cma/valnote'
    title = 'Test title'
    docurl = 'file.html'
    globalvar.valorder_file = os.path.join(valnote_dir,
                                           'control',
                                           'valorder.dat')

    create_valnote(outdir, images_dir, thumbs_dir, valnote_dir, title,
                   docurl=docurl, exper=exper, control=control)
