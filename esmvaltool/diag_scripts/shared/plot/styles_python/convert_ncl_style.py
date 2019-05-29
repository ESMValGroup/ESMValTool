#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to convert nvl style files to python style files.

Example
-------
The script is simply called by::

    $ python convert_ncl_style.py

Global attributes
-----------------
INPUT_FILE : str
    Name of the ncl style file.
OUTPUT_FILE : str
    Name of the new python style file (yml format).

"""

import os

import yaml

# Global variables
INPUT_FILE = 'cmip6.style'
OUTPUT_FILE = 'cmip6.yml'
BASE_DIR = os.path.dirname(os.path.realpath(__file__))

HEADER_FILE = 'style_header'
DATASET = 'dataset'
COLOR = 'color'
DASH = 'dash'
THICKNESS = 'thick'
MARK = 'mark'
AVG_STD = 'avgstd'
FILLING = 'facecolor'

INFORMATION = [DATASET, COLOR, DASH, THICKNESS, MARK, AVG_STD]


def read_line(line):
    """Read line of the ncl style file."""
    # Read information
    info_dict = {}
    for (idx, line_elem) in enumerate(line):
        info = line_elem.strip()
        option = INFORMATION[idx]

        # Convert color to hex string
        if option == COLOR:
            color = info.split(',')
            info = '#'
            for col in color:
                col = "{:02x}".format(int(col))
                info += col

        # Convert mark index to matplotlib marker
        elif option == MARK:
            # Filling
            info_dict[FILLING] = info_dict[COLOR] if info == '16' else 'none'

            # Shape
            shape = {
                '0': 'x',
                '1': '.',
                '2': '+',
                '3': 'x',
                '4': 'o',
                '5': 'x',
                '6': 's',
                '7': '^',
                '8': 'v',
                '9': 'D',
                '10': '<',
                '11': '>',
                '12': '*',
                '13': 'h',
                '14': '.',
                '15': 'x',
                '16': 'o',
            }
            info = shape.get(info, 'o')

        # Convert dash index to matplotlib dash marker
        elif option == DASH:
            dash = {
                '0': '-',
                '1': '--',
                '2': ':',
                '3': '-.',
                '4': '-.',
                '5': '--',
                '6': '--',
                '7': '-.',
                '8': '-.',
                '9': '-.',
                '10': '-.',
                '11': '--',
                '12': '--',
                '13': '--',
                '14': '--',
                '15': '--',
                '16': '--',
            }
            info = dash.get(info, '-')

        # Convert str to int
        elif option in (AVG_STD, THICKNESS):
            info = int(info)

        # Add information
        info_dict[option] = info
    return info_dict


def read_ncl_style(file_name):
    """Read ncl style file."""
    output = []
    with open(file_name, 'r') as file_:
        for line in file_:
            line = line.strip()

            # Ignore commentary lines
            if line.startswith('#'):
                continue

            # Get lines with valid information (seperated by '|')
            line = line.split('|')
            if len(line) != len(INFORMATION):
                continue

            # Read line
            output.append(read_line(line))
    print("Read '{}'".format(file_name))

    # Convert list to dictionary
    output_dict = {}
    for info in output:
        dataset = info[DATASET]
        del info[DATASET]
        output_dict[dataset] = info

    return output_dict


def write_yml_file(dataset_info, file_name):
    """Write configuration file."""
    header_path = os.path.join(BASE_DIR, HEADER_FILE)
    with open(file_name, 'w') as outfile:
        with open(header_path, 'r') as header_file:
            header = header_file.read()
        outfile.write(
            header.format(
                output_file=OUTPUT_FILE, script=os.path.basename(__file__)))
        yaml.safe_dump(dataset_info, outfile, default_flow_style=False)
    print("Wrote '{}'".format(file_name))


# Execute script if called directly
if __name__ == '__main__':
    INPUT_PATH = os.path.normpath(
        os.path.join(BASE_DIR, '..', 'styles', INPUT_FILE))
    OUTPUT_PATH = os.path.join(BASE_DIR, OUTPUT_FILE)
    STYLES = read_ncl_style(INPUT_PATH)
    write_yml_file(STYLES, OUTPUT_PATH)
