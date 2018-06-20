#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script to convert nvl style files to python style files.

Example
-------
The script is simply called by::

    $ python import_ncl_style.py

Global attributes
-----------------
INPUT_FILE : str
    Path to the ncl style file.
OUTPUT_FILE : str
    Path to the new python style file.

"""


import numpy as np


# Global variables
INPUT_FILE = '../styles/cmip5.style'
OUTPUT_FILE = 'cmip5.style'
HEADER_FILE = 'style_header'
MODEL = 'model'
COLOR = 'color'
DASH = 'dash'
THICKNESS = 'thick'
MARK = 'mark'
AVG_STD = 'avgstd'
FILLING = 'facecolor'
INFORMATION = [MODEL, COLOR, DASH, THICKNESS, MARK, AVG_STD]


def convert_mark(info):
    """Convert numerical mark info to matplotlib mark."""
    if info == '0':
        return 'x'
    elif info == '1':
        return '.'
    elif info == '2':
        return '+'
    elif info == '3':
        return 'x'
    elif info == '4':
        return 'o'
    elif info == '5':
        return 'x'
    elif info == '6':
        return 's'
    elif info == '7':
        return '^'
    elif info == '8':
        return 'v'
    elif info == '9':
        return 'D'
    elif info == '10':
        return '<'
    elif info == '11':
        return '>'
    elif info == '12':
        return '*'
    elif info == '13':
        return 'h'
    elif info == '14':
        return '.'
    elif info == '15':
        return 'x'
    elif info == '16':
        return 'o'
    else:
        return 'o'


def read_ncl_style(file_name):
    """Read ncl style file."""
    output = []
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()

            # Ignore commentary lines
            if line.startswith('#'):
                continue

            # Get lines with valid information (seperated by '|')
            line = line.split('|')
            if len(line) != len(INFORMATION):
                continue

            # Read information
            infos = {}
            for (idx, line_elem) in enumerate(line):
                info = line_elem.strip()
                option = INFORMATION[idx]

                # Convert color to hex string
                if option == COLOR:
                    color = info.split(',')
                    info = '#'
                    for col in color:
                        col = np.int(col.strip())
                        info += format(col, '02X')

                # Convert mark index to matplotlib marker
                if option == MARK:
                    # Filling
                    if info == '16':
                        infos.update({FILLING: infos[COLOR]})
                    else:
                        infos.update({FILLING: 'none'})

                    # Shape
                    info = convert_mark(info)

                # Add information
                infos.update({INFORMATION[idx]: info})
            output.append(infos)

    return output


def write_config_file(model_infos, file_name):
    """Write configuration file."""
    with open(file_name, 'w') as file:
        with open(HEADER_FILE, 'r') as header_file:
            header = header_file.read()
        file.write(header)

        # Write model styles
        for model_info in model_infos:
            file.write("\n")
            for info in INFORMATION:
                if info == MODEL:
                    file.write("[{0}]\n".format(model_info[info]))
                else:
                    file.write("{0} = {1}\n".format(info, model_info[info]))
            file.write("{0} = {1}\n".format(FILLING, model_info[FILLING]))


# Execute script if called directly
if __name__ == '__main__':
    STYLES = read_ncl_style(INPUT_FILE)
    write_config_file(STYLES, OUTPUT_FILE)
