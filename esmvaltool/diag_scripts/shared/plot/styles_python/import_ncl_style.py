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


def read_style(file_name):
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
            for i in range(len(line)):
                info = line[i].strip()
                option = INFORMATION[i]

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
                    if info == '0':
                        info = 'x'
                    elif info == '1':
                        info = '.'
                    elif info == '2':
                        info = '+'
                    elif info == '3':
                        info = 'x'
                    elif info == '4':
                        info = 'o'
                    elif info == '5':
                        info = 'x'
                    elif info == '6':
                        info = 's'
                    elif info == '7':
                        info = '^'
                    elif info == '8':
                        info = 'v'
                    elif info == '9':
                        info = 'D'
                    elif info == '10':
                        info = '<'
                    elif info == '11':
                        info = '>'
                    elif info == '12':
                        info = '*'
                    elif info == '13':
                        info = 'h'
                    elif info == '14':
                        info = '.'
                    elif info == '15':
                        info = 'x'
                    elif info == '16':
                        info = 'o'
                    else:
                        info = 'o'

                # Add information
                infos.update({INFORMATION[i]: info})
            output.append(infos)

    return output


def write_config_file(model_infos, file_name):
    """Write configuration file."""
    with open(file_name, 'w') as file:
        with open(HEADER_FILE, 'r') as header_file:
            header = header_file.read()
        file.write(header)

        # Write model styles
        for model_infos in model_infos:
            file.write("\n")
            for info in INFORMATION:
                if info == MODEL:
                    file.write("[{0}]\n".format(model_infos[info]))
                else:
                    file.write("{0} = {1}\n".format(info, model_infos[info]))
            file.write("{0} = {1}\n".format(FILLING, model_infos[FILLING]))


# Execute script if called directly
if __name__ == '__main__':
    styles = read_style(INPUT_FILE)
    write_config_file(styles, OUTPUT_FILE)
