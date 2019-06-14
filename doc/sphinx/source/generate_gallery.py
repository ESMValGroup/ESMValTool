#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create gallery with all available recipes."""

import os

RECIPE_DIR = 'recipes'
OUT_PATH = os.path.abspath('gallery.rst')
HEADER = ('#######\nGallery\n#######\n\n'
          'This section shows example plots produced by ESMValTool. For more '
          'information, click on the footnote below the image.\n\n')
WIDTH = ':width: 90%'
FIGURE_STR = '.. figure::'
IMAGE_STR = ' image:: '
TABLE_SEP = ('+---------------------------------------------------'
             '+---------------------------------------------------+\n')
EMPTY_TABLE = ('|                                                   '
               '|                                                   |\n')
CELL_WIDTH = 50


def _get_figure_index(file_content):
    """Get index of figure in text."""
    if FIGURE_STR in file_content:
        return file_content.index(FIGURE_STR) + len(FIGURE_STR)
    if IMAGE_STR in file_content:
        return file_content.index(IMAGE_STR) + len(IMAGE_STR)
    raise ValueError("File does not contain image")


def _get_next_row(filenames, file_contents):
    """Get next row."""
    figure_idx = [_get_figure_index(content) for content in file_contents]
    figure_paths = [
        file_contents[idx][fig_idx:].split('\n')[0].strip()
        for (idx, fig_idx) in enumerate(figure_idx)
    ]
    subst = [
        '|{}|'.format(os.path.splitext(filename)[0]) for filename in filenames
    ]
    link = [file_contents[0].split()[1][1:-1]]
    if figure_paths[1] == '':
        subst[1] = ''
        link.append('')
    else:
        link.append(file_contents[1].split()[1][1:-1])

    # Build table
    row = ''
    refs = ''
    row += TABLE_SEP
    row += '| {}| {}|\n'.format(subst[0].ljust(CELL_WIDTH),
                                subst[1].ljust(CELL_WIDTH))
    row += EMPTY_TABLE
    left_col = '[#]_'.ljust(CELL_WIDTH)
    if figure_paths[1] == '':
        right_col = ''.ljust(CELL_WIDTH)
    else:
        right_col = '[#]_'.ljust(CELL_WIDTH)
    row += '| {}| {}|\n'.format(left_col, right_col)

    # Build refs
    for (idx, path) in enumerate(figure_paths):
        if path == '':
            continue
        refs += f'.. {subst[idx]} image:: {path}\n'
        refs += f'   {WIDTH}\n'
        refs += '\n'
        refs += f'.. [#] :ref:`{link[idx]}`\n'
        refs += '\n'

    return (row, refs)


def main():
    """Generate gallery for recipe plots."""
    print(f"Generating gallery at {OUT_PATH}")
    left_col = True
    table = ''
    refs = ''
    filenames = []
    file_contents = []
    for filename in sorted(os.listdir(RECIPE_DIR)):
        if not filename.startswith('recipe_'):
            continue
        if not filename.endswith('.rst'):
            continue
        with open(os.path.join(RECIPE_DIR, filename), 'r') as in_file:
            recipe_file = in_file.read()
        if (FIGURE_STR not in recipe_file and IMAGE_STR not in recipe_file):
            print(f"INFO: {filename} does not contain an image, skipping")
            continue
        if not recipe_file.startswith('..'):
            print(f"INFO: {filename} does not contain reference at top, "
                  "skipping")
            continue

        # Get next row
        if left_col:
            left_col = False
            filenames = [filename]
            file_contents = [recipe_file]
            continue
        else:
            left_col = True
            filenames.append(filename)
            file_contents.append(recipe_file)
            new_row = _get_next_row(filenames, file_contents)
            table += new_row[0]
            refs += new_row[1]

    # Last row
    if len(filenames) == 1:
        filenames.append('')
        file_contents.append(f'{FIGURE_STR}\n')
        new_row = _get_next_row(filenames, file_contents)
        table += new_row[0]
        refs += new_row[1]
    table += TABLE_SEP
    table += '\n'

    # Write file
    whole_file = HEADER + table + refs
    with open(OUT_PATH, 'w') as out_file:
        print(whole_file, file=out_file)
    print(f"Wrote {OUT_PATH}")


if __name__ == '__main__':
    main()
