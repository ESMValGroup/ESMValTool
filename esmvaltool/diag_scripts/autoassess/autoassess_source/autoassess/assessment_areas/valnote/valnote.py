#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
'''
Validation Notes
'''

import argparse
import os
import pdb
import sys
import traceback

from code import globalvar
from code import make_means
from code import make_plots
from code import valmod as vm
from code import vnc

# Links to
DOCLINKS = {
    'GJ': 'https://code.metoffice.gov.uk/trac/GA/wiki/GAJobs',
    'GM': 'https://code.metoffice.gov.uk/trac/gmed/ticket',
    'GA': 'https://code.metoffice.gov.uk/trac/GA/ticket',
    'GL': 'https://code.metoffice.gov.uk/trac/GL/ticket',
    'GO': 'https://code.metoffice.gov.uk/trac/GO/ticket',
    'TW': 'http://collab.metoffice.gov.uk/twiki/bin/view/Project/CAPTIVATE',
}

DOCMSG = '''Use XXtttt where XX is one of:
GJ = links to GA trac wiki page under GAJobs where tttt is the wiki page;
GM = links to GMED ticket number where tttt is the ticket number;
GA = links to GA ticket number where tttt is the ticket number;
GL = links to GL ticket number where tttt is the ticket number;
GO = links to GO ticket number where tttt is the ticket number;
TW = links to collaboration twiki page where tttt is the twiki page;'''


def run(options):

    # Define global variables
    globalvar.debug = options.debug
    globalvar.pub = options.pub
    globalvar.pkl = options.pkl
    globalvar.store_regrids = options.store_regrids

    # Deduce input directories
    if options.control_dir:
        control_dir = os.path.abspath(options.control_dir)
    else:
        control_dir = os.path.join(options.path, 'control')
    print 'Using control files from: ' + control_dir

    # Set the valorder file
    if options.valorder:
        globalvar.valorder_file = os.path.abspath(options.valorder)
    else:
        globalvar.valorder_file = os.path.join(control_dir, 'valorder.dat')
    print 'Using valorder file: ' + globalvar.valorder_file

    # Set the source file
    if options.source_file:
        globalvar.source_file = os.path.abspath(options.source_file)
    else:
        globalvar.source_file = 'source_file.dat'
    print 'Using source file: ' + globalvar.source_file

    # Load the source_file
    make_plots.source_dict = vm.read_control_file(globalvar.source_file)

    # Get the names of the experiment and control keys for future use.
    # This is a bit of a fudge and needs sorting out.
    exper_key = make_plots.source_dict.keys()[0]
    control_key = make_plots.source_dict.keys()[1]
    for key in make_plots.source_dict.keys():
        if key[0] == 'e':
            exper_key = key
        if key[0] == 'c':
            control_key = key

    # Discover what the experiment and control jobids are from the
    #  source_file.dat
    exper = make_plots.source_dict[exper_key]['jobid'][0]
    control = make_plots.source_dict[control_key]['jobid'][0]

    # Decide what steps we are doing
    steps = options.steps.upper()
    do_main = 'M' in steps
    do_means = 'R' in steps
    do_vnc = 'V' in steps

    # Set up the title
    if options.title is None:
        options.title = '{0} vs {1}'.format(exper, control)
        print 'Warning: you have not specified a title (see --help)'
        print '  Using automatic title: ' + options.title

    # Generate a documentation link from the --link input
    if options.link:
        link_pref = options.link[0:2]
        link_suff = options.link[2:]
        if link_pref in DOCLINKS:
            options.docurl = os.path.join(DOCLINKS[link_pref], link_suff)
        else:
            print 'You have not specified --link using the correct format.'
            print DOCMSG

    # Display warning if no documentation link supplied or generated.
    if options.docurl is None:
        print 'Warning: you have not supplied a documentation URL (see --help)'
        print "  You have documented your run, haven't you?"
        options.docurl = 'file.html'

    # If this is for a publication append _pub to the end of the output directory
    if globalvar.pub:
        options.out_dir += '_pub'

    # Define output sub-directories
    outdir = os.path.abspath(options.out_dir)
    images_dir = os.path.join(outdir, 'images')
    thumbs_dir = os.path.join(outdir, 'thumbs')
    csv_dir = os.path.join(outdir, 'csv')

    # Make output (sub)directories
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists(images_dir):
        os.makedirs(images_dir)
    if not os.path.exists(thumbs_dir):
        os.makedirs(thumbs_dir)
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)
    if not os.path.exists(options.tmp_dir):
        os.makedirs(options.tmp_dir)

    # Load the data into one big dictionary
    if do_main or do_means:
        print 'Loading model and obs data'
        make_plots.obs_dict = vm.read_control_file('obs_file.dat',
                                                   ctldir=control_dir)
        # Prefix obs files with main obs directory
        for key in make_plots.obs_dict.keys():
            for skey in make_plots.obs_dict[key].keys():
                if skey not in ['name', 'constraint_limiter']:
                    make_plots.obs_dict[key][skey] = [os.path.join(options.obs_dir, fname) for fname in make_plots.obs_dict[key][skey]]
        make_plots.load_all_data(do_main, options.tmp_dir)

        print 'Loading extra data'
        extras_dict = vm.read_info_file(os.path.join(control_dir,
                                                     'extras_file.dat'))
        # Prefix extras files with main extras directory
        for key in extras_dict.keys():
            extras_dict[key] = os.path.join(options.extras_dir, extras_dict[key])
        vm.load_extra_data(extras_dict)

    # Run the relevant bits of code
    if do_main:
        try:
            make_plots.main(images_dir, control_dir, thumbs_dir, csv_dir,
                            options.tmp_dir, num_to_save=options.num_to_save)
        except:
            debug_code()

    if do_means:
        try:
            make_means.main(control_dir, csv_dir)
        except:
            debug_code()

    if do_vnc:
        # Generate the main validation note web page
        print 'Creating validation note webpage...'
        try:
            vnc.create_valnote(outdir, images_dir, thumbs_dir, options.path,
                               options.title, docurl=options.docurl,
                               exper=exper, control=control, csv_dir=csv_dir)
        except:
            debug_code()


def debug_code():
    '''
    Code to start an interactive session at the point of the
    crash if debug is turned on.
    '''

    import code.globalvar as globalvar
    type, value, tb = sys.exc_info()
    traceback.print_exc()

    if globalvar.debug:
        print 'Debugging is turned on'
        print 'Starting interactive session at the point of the crash.'
        print 'Commands are:'
        print '? = help'
        print 'w = where am I'
        print 'q = quit debugger'
        pdb.post_mortem(tb)
    else:
        exit(1)


def parse_args(cli_args):
    """
    Parse arguments in a function to facilitate testing. Contains all command
    line options.

    :param list cli_args: Command line arguments from sys.argv.
    :returns: Checked command line arguments.
    :rtype: argparse.Namespace
    """

    if 'SCRATCH' in os.environ:
        scratch = os.environ['SCRATCH']
    else:
        scratch = None

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--out-dir', required=True, help='Parent directory for output validations.')
    parser.add_argument('--title', required=True, help='Title for validation note web page')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--link', help='Generates a link to documentation. ' + DOCMSG)
    group.add_argument('--url', dest='docurl', help='Full URL of run documentation page. Alternative to --link LINK.')

    parser.add_argument('--source', dest='source_file', help='Name of source file (defaults to ./source_file.dat)')
    parser.add_argument('--valorder', default=None, help='Name of valorder file if different from standard plot list')
    parser.add_argument('--steps', default='MRV', help='Perform only some steps.  STEPS is a series of characters identifying individual steps: M = main plots, R = means table (previously called radiation table), V = Create validation note webpage')

    parser.add_argument('--tmp-dir', default=scratch, help='Directory in which to locate the temporary working directory')
    parser.add_argument('--control-dir', default=None, help='Directory that contains all the control files. If omitted it will use the files in valnote/control alongside this executable. Use this if you want to customise your plots (e.g. change contour levels and colours).')
    parser.add_argument('--obs-dir', default='/project/cma/clim', help='Directory that contains all the observation files.')
    parser.add_argument('--ancil-dir', default='/project/cma/ancil', help='Directory that contains all the ancillary files.')

    parser.add_argument('--num-to-save', default=0, type=int, help='Plot number (as listed in valorder) to save in pkl file for use in test_colours.py (for use in code development)')
    parser.add_argument('--debug', default=False, action='store_true', help='Switch on debugging. This stops the code at breakpoints when there are problems (using the pdb debugger).')
    parser.add_argument('--pkl', default=False, action='store_true', help='Generate and use pkl files in temporary directory to speed up reruns.')
    parser.add_argument('--pub', default=False, action='store_true', help='Validation note is for a publication. Do not output the UMUI job ids, Rose suite ids, observation dates or field names in the final plots. Plot titles only state model name as specified in the source_file.dat and the observation name. Title font size is increased to 20. Output in eps (as well as png). Also appends _pub to the final validation note folder (so as to not overwrite the original). Please use alongside a short valorder file (pointed to using --valorder) to limit the validation note to just the plots you want in the publication.')
    parser.add_argument('--store-regrids', default=False, action='store_true', help='Store regridded data for future use. Speeds up the validation note but uses more memory. Default=True')

    # Return parsed args
    return parser.parse_args(cli_args)


# **************************************************
# Code below is only called if this file is run as a
#  script from the command line
# **************************************************

def main():
    args = parse_args(sys.argv[1:])

    # Add new variables to args that are useful
    args.path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    args.tmp_dir = os.path.join(args.tmp_dir, 'valnote')
    args.extras_dir = os.path.join(args.ancil_dir, 'masks')

    run(args)


if __name__ == '__main__':
    main()
