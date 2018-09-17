"""Assemble energy budget results in one table figure."""
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import numpy as np


def render_mpl_table(table, header_rows, header_columns,
                     col_width, highlight_cells, **kwargs):
    """
    Create table using matplotlib.

    Inspired by: http://stackoverflow.com/a/39358752
    """
    # some default values
    row_height = 0.6
    font_size = 14
    header_color = '#40466e'
    highlight_color = '#c0848e'
    row_colors = ['#f1f1f2', '#ffffff']
    edge_color = '#ffffff'
    bbox = [0, 0, 1, 1]

    table = np.array(table)
    # create empty plot; no axes
    size = (np.array(table.shape[::-1]) + np.array([0, 1])) \
        * np.array([col_width, row_height])
    fig, ax = plt.subplots(figsize=size)
    ax.axis('off')

    # plain table with data
    if header_rows:
        mpl_table = ax.table(cellText=table[1:], bbox=bbox,
                             colLabels=table[0], **kwargs)
    else:
        mpl_table = ax.table(cellText=table[:], bbox=bbox, **kwargs)

    # fontsize
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    # format table
    # -> http://matplotlib.org/1.3.1/users/text_props.html
    # -> https://github.com/matplotlib/matplotlib/blob/master/
    # lib/matplotlib/table.py
    # VPREDOI
    # python3 uses items() instead of iteritems()
    for k, cell in mpl_table._cells.items():
        # line color
        cell.set_edgecolor(edge_color)

        if k[0] < header_rows:               # header rows
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        elif k[1] < header_columns:           # header columns
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        elif k in highlight_cells:
            # highlighted cells (headers have priority)
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(highlight_color)
        else:                                 # format data cells
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
    return ax
