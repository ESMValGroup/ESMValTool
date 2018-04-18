#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Create HTML page from an AutoAssess template. Placeholders for assessment title,
area name, suite IDs, and dates.

A path to another template can also be specified.
"""
from __future__ import division, print_function
import argparse
import os.path
import sys
from jinja2 import Environment, FileSystemLoader
from autoassess import area_page_template_dir
from autoassess import overview_page_template_dir


def parse_args(cli_args):
    """
    Parse arguments in a function to facilitate testing. Contains all command
    line options.

    :param list cli_args: Command line arguments from sys.argv.
    :returns: Checked command line arguments.
    :rtype: argparse.Namespace
    """

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--title', required=True,
                        help='Description of the model comparison.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--overview', action='store_true',
                       help='Create overview page.')
    group.add_argument('--area',
                       help='Create page for assessment area.')

    parser.add_argument('--template',
                        help='Path to jinja2 template.')

    parser.add_argument('--ref-suite-id', required=True,
                        help='Rose suite ID of reference model run.')
    parser.add_argument('--exp-suite-id', required=True,
                        help='Rose suite ID of experiment model run.')
    parser.add_argument('--start-date', required=True,
                        help='Start of assessment period (YYYY/MM/DD).')
    parser.add_argument('--end-date', required=True,
                        help='End of assessment period (YYYY/MM/DD).')

    # Return parsed args
    return parser.parse_args(cli_args)


def load_template_from_path(template_path):
    """
    Load and return jinja2 template from the given path.

    :param str template_path:
    :returns: Jinja2 template
    :rtype: jinja2.environment.Template
    """
    template_dir = os.path.dirname(template_path)
    template_filename = os.path.basename(template_path)
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template(template_filename)
    return template


def main():
    """
    Load template and fill in placeholders with arguments given at the
    command line.
    """
    TEMPLATE = 'jinja2_template.html'
    args = parse_args(sys.argv[1:])
    if not args.area:
        args.area = ''

    if args.template:
        template = load_template_from_path(args.template)
    else:
        if args.overview:
            template_path = os.path.join(overview_page_template_dir(), TEMPLATE)
            template = load_template_from_path(template_path)
        elif args.area:
            template_path = os.path.join(area_page_template_dir(args.area), TEMPLATE)
            template = load_template_from_path(template_path)
        else:
            raise Exception()

    print(template.render(title=args.title,
                          area=args.area.title(),  # capitalise each word
                          ref_suite_id=args.ref_suite_id,
                          exp_suite_id=args.exp_suite_id,
                          start_date=args.start_date,
                          end_date=args.end_date))

if __name__ == '__main__':
    main()
