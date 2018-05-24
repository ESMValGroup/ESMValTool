#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Create normalised assessment criteria plots (NAC plots).

The produced metrics (metrics.csv) of both runs are compared to metrics produced
from observations. Additionally, NAC plots are produced for subset of related
metrics.
"""


from __future__ import division, print_function
import argparse
import fnmatch
import glob
import os.path
import subprocess as subp
import sys
from autoassess import argparser, metrics_from_obs_dir, metric_selections_dir

# Create script argument parser
POS_ARGS = []
OPT_ARGS = [
    ('area', dict(required=True, help='Name of the assessment area.')),
    ('area-title', dict(required=True, help='Display name of the assessment area.')),
    ('ref-suite-id', dict(required=True, help='Rose suite ID of reference model run.')),
    ('exp-suite-id', dict(required=True, help='Rose suite ID of experiment model run.')),
    ('ref-desc', dict(required=True, help='Description of reference model run.')),
    ('exp-desc', dict(required=True, help='Description of experiment model run.')),
    ('out-dir', dict(required=True, help='Directory with AutoAssess output.')),
    ('plot-dir', dict(required=True, help='Output directory for NAC plots.')),
]
ARG_PARSE = argparser(__doc__, POS_ARGS, OPT_ARGS)


def main():
    """ """
    # Deal with script arguments
    args = ARG_PARSE.parse_args(sys.argv[1:])
    args.out_dir = os.path.abspath(args.out_dir)
    args.plot_dir = os.path.abspath(args.plot_dir)

    # Set managable variables for control and reference suite ids
    ref = args.ref_suite_id
    exp = args.exp_suite_id

    # Set plot legends and title
    title = '"{0}"'.format(args.area_title)
    leg_template = '"{0} ({1})"'
    leg_ref = leg_template.format(args.ref_desc, ref)
    leg_exp = leg_template.format(args.exp_desc, exp)

    # get metrics from observations and metric subset files from AutoAssess
    obs_dir = metrics_from_obs_dir(args.area)
    selections_dir = metric_selections_dir(args.area)

    # Set output dir and metric file locations
    aa_out_dir = os.path.join(args.out_dir, exp + '_vs_' + ref)
    file_ref = os.path.join(aa_out_dir, args.area, ref, 'metrics.csv')
    file_exp = os.path.join(aa_out_dir, args.area, exp, 'metrics.csv')

    # Identify observation file - should only be one
    obs_files = [file_ for file_ in os.listdir(obs_dir)
                 if fnmatch.fnmatch(file_, '*.csv')]
    if len(obs_files) == 0:
        obs_file = ''
    elif len(obs_files) > 1:
        raise MultipleObsFilesError(
            'More than one file csv file found in {0}'.format(obs_dir)
        )
    else:
        obs_file = obs_files[0]
        obs_file = os.path.join(obs_dir, obs_file)

    # Loop over sets of metrics to produce plot
    for selection_path in glob.glob(os.path.join(selections_dir, '*.csv')):
        selection_filename = os.path.basename(selection_path)
        plot_path = os.path.join(args.plot_dir, selection_filename.split('.')[0])

        # Create shell command
        command = ' '.join([
            'aa_nac_plot',
            '--plot', plot_path,
            '--title', title,
            '--ref', leg_ref,
            '--exp', leg_exp,
            '--file-ref', file_ref,
            '--file-exp', file_exp,
            '--file-ord', selection_path,
            '--file-obs', obs_file
        ])

        # Run shell command
        try:
            subp.check_output(command, stderr=subp.STDOUT, shell=True)
        except subp.CalledProcessError as e:
            print('### Shell output ###')
            print(e.output, file=sys.stderr)
            print('####################')
            raise


class MultipleObsFilesError(Exception):
    pass


if __name__ == '__main__':
    main()
