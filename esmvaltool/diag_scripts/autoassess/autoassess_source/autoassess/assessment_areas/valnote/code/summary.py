'''
Module to generate summary RMSE files
'''

import copy
import csv
import glob
import os
import pdb
import re

import markup
from markup import oneliner as line
import globalvar

# compile regular expression for isccp format for use by two routines
# Todo: could make this more restrictive by including e.g. Thin|Medium|Thick
#        but need to check syntax
isccp_style = re.compile('-top [MT][a-z]+ Cloud')


def isccp_reshape(isccp_file):
    # input file has format description, value
    # only one description per line
    # - expect each description to be repeated 3 times

    # output file is to have description, value, value, value
    # with no descriptions repeated
    # want to preserve order of descriptions

    # open the file and read all lines.  This leaves the position at the
    # end of the file (ready for appending)

    with open(isccp_file, 'r') as infile:
        cvs_reader = csv.reader(infile)

        descriptions = []  # descriptions is a list of strings
        lines = []         # lines is a list of lists
        start = False
        for desc, value in cvs_reader:

            # Allow for any number of header lines
            # Last header line begins with "Description"
            if not start:
                if desc == "Description":
                    start = True
                continue

            try:
                i = descriptions.index(desc)
            except ValueError:
                # desc not found before: append to lists
                # perhaps need to check is not header row
                # (here or before 'try')
                descriptions.append(desc)
                lines.append([desc, value])
            else:
                # desc is in list: append value to correct list
                lines[i].append(value)

    return lines


def isccp_append():
    # find all isccp csv files, find matching main file,
    # for each file, re-shape and append the isccp values

    isccp_files = glob.glob('summary_*_isccp.csv')

    for ifile in isccp_files:
        mainfile = ifile.replace('_isccp', '')
        assert os.path.isfile(mainfile), 'No main csv file corresponding to '+ifile

        # TODO: Can this be done with context manager form?
        outfile = open(mainfile, 'r+')
        lines = outfile.readlines()

        skip = False
        for line in lines:
            if isccp_style.search(line):
                print 'ISSCP summary already present in main file: '+mainfile
                print 'skipping...'
                skip = True
                outfile.close()
                break

        if not skip:
            lines_append = isccp_reshape(ifile)

            # append a newline in case main file doesn't end in one
            # (only happens if main step fails uncleanly - otherwise this will
            # add a blank line before isccp, but that's OK,
            # csv.reader ignores it)
            outfile.write('\n')
            csv_writer = csv.writer(outfile)
            csv_writer.writerows(lines_append)
            outfile.close()


def html_init(descrip, region, fileinfo, page_top):

    # Initialise the page
    page = copy.deepcopy(page_top)

    rms1 = descrip['rms1']

    index = rms1.find('expt:') + 5
    exp = rms1[index:index+7]

    index = rms1.find('control:') + 8
    ctl = rms1[index:index+7]

    title = "Validation note summary for %s v %s - %s" % (exp, ctl, region)
    page.init(title=title)
    page.region = region

    page.table()
    page.tr()
    for dfile in fileinfo:
        page.td(line.a(dfile['linkname'], href=dfile['html']))
    page.tr.close()
    page.table.close()

    page.p()
    page.h2(title)
    page.p()

    page.table(border=1)
    page.tr()

    page.td('Variable')
    page.td('PNG')
    page.td('RMS1 (expt v control)')
    page.td('RMS2 (control v obs)')
    page.td('RMS3 (expt v obs)')
    page.td('RMS3-RMS2 (expt-control)')
    page.td('(RMS3-RMS2) /RMS2(%)')
    page.td('better(+) /worse(-)')

    page.tr.close()

    return page


def do_sums(d, valcodes, page):
    # do the sums to get the extra columns and the --/++ stuff
    desc = d['description']
    rms1 = d['rms1']
    rms2 = d['rms2']
    rms3 = d['rms3']

    # work out the derived fields from the rms numbers
    try:
        # Convert the strings to floats with 'eval'
        rms1 = eval(rms1)
        rms2 = eval(rms2)
        rms3 = eval(rms3)
    except SyntaxError:
        print "Warning: problem with values in CSV file for:", page.region
        print "Description:", desc
        print "RMS1:", rms1
        print "RMS2:", rms2
        print "RMS3:", rms3
        if globalvar.debug:
            pdb.set_trace()
        else:
            print "Skipping this field in RMS table"
            return
    except (TypeError, NameError):
        # Probably due to blank/NaN in some columns (e.g. where no obs)
        # Just write what we've got and leave the derived values blank
        betworse = None
        diff = None
        diff_pct = None
    else:
        # do the sums
        diff = rms3 - rms2

        if rms2 != 0:
            diff_pct = int(100.0 * diff / rms2)
        else:
            diff_pct = 0.0

        # todo: put green / red in for the ++ and --
        # do via style sheet
        if diff_pct <= -10:
            betworse = '++'
        elif diff_pct < 0:
            betworse = '+'
        elif diff_pct == 0:
            betworse = '&nbsp;'
        elif diff_pct < 10:
            betworse = '-'
        elif diff_pct >= 10.0:
            betworse = '--'
        else:
            raise Exception('Problem with diff_pct number: not valid float?')

        diff_pct = int(diff_pct)

    # Get the png file from the dictionary
    png_file = d['png_file']
    png_file_path = os.path.join('images', png_file)

    # write the html file using the markup module
    # (could optionally put the same info into a csv file - hopefully no need)

    page.tr()
    page.td(desc)
    page.td(line.a(png_file, href=png_file_path))

    page.td(rms1)
    page.td(rms2)
    page.td(rms3)

    page.td(diff)
    page.td(diff_pct)
    page.td(betworse)

    page.tr.close()


def csv2html(page_top, valorderfile=None, csv_dir=None, outdir=None):

    # read the valorder file into a dictionary
    # then invert the dict because we want to map descriptions to codes
    with open(valorderfile, 'r') as f:
        valorder = dict([line.strip().split(':', 1)
                         for line in f.readlines()
                         if ':' in line and line[0] != '#'])
    valcodes = dict([[v, k] for k, v in valorder.items()])

    # append the ISCCP files onto the main files
    isccp_append()

    csvfiles = glob.glob(os.path.join(csv_dir, 'summary*.csv'))
    # want certain regions to come first, followed by others in alphabetical
    # order so - sort alphabetically, and then bring desired regions to front
    # of list
    csvfiles.sort()
    order = ['global', 'north', 'south', 'tropical_ocean', 'tropical_land']
    order.reverse()  # so that 1st element is brought to front last, and so on
    for region in order:
        cfile = 'summary_'+region+'.csv'
        if cfile in csvfiles:
            csvfiles.remove(cfile)
            csvfiles.insert(0, cfile)

    # create 'fileinfo': a list of dictionaries
    fileinfo = list()
    for cfile in csvfiles:
        if not cfile.endswith('_isccp.csv'):
            csv_file_end = os.path.split(cfile)[1]
            region = csv_file_end.replace('.csv', '')
            region = region.replace('summary', '')
            region = region.replace('_', ' ')
            filedict = {'region': region,
                        'linkname': 'Summary:'+region,
                        'csv': cfile,
                        'html': csv_file_end.replace('.csv', '.html')}
            fileinfo.append(filedict)

    # loop over csvfiles, converting them all.
    for filedict in fileinfo:

        print "Creating "+filedict['linkname']+" html page"

        # loop over the plots, doing sums and writing to html file
        # read the csv file
        # get description, rms1, rms2, rms3
        columns = ['description', 'png_file', 'rms1', 'rms2', 'rms3']
        with open(filedict['csv'], 'r') as ifile:
            dictreader = csv.DictReader(ifile, columns)

            start = False
            for d in dictreader:

                # Allow for any number of header lines
                # Last header line begins with "Description"
                if not start:
                    if d['description'] == "Description":
                        page = html_init(d, filedict['region'], fileinfo,
                                         page_top)
                        start = True
                    continue

                # Calculate derived stats and write to html page
                do_sums(d, valcodes, page)

            page.table.close()

        with open(filedict['html'], 'w') as ofile:
            print >>ofile, page
    
    # Return the URL to the first page (if it exists)    
    summary_page_url = 'unknown.html'
    if fileinfo:
        summary_page_url = fileinfo[0]['html']
    return summary_page_url


# main script for testing
if __name__ == '__main__':
    isccp_append()
    linkfile = csv2html(valorderfile='./valorder.dat')
    print linkfile
