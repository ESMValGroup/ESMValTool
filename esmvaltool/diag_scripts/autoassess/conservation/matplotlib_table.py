"""Assemble energy budget results in one table figure."""
import numpy as np
import matplotlib.pyplot as plt


# need to clean this function up, it is UGLY
def render_mpl_table(table, header_rows=0, header_columns=0,
                     highlight_cells=None,
                     col_width=3.0, row_height=0.6, font_size=14,
                     header_color='#40466e', highlight_color='#c0848e',
                     row_colors=['#f1f1f2', '#ffffff'], edge_color='#ffffff',
                     bbox=[0, 0, 1, 1], **kwargs):
    """
    Create table using matplotlib.

    Inspired by: http://stackoverflow.com/a/39358752

    Example:

    table = [
        ['',     '',     '',     ''],
        ['TRIP', 'bla',  '2200', ''],
        ['',     'blub', '2100', ''],
        ['',     'foo',  '1500', ''],
        ['',     'foo',  '1500', ''],
        ['',     'foo',  '1500', ''],
        ['',     'foo',  '1500', ''],
        ['',     'DFs',  '8888', '']
    ]

    render_mpl_table(array, header_rows=2, header_columns=2,
                     highlight_cells=[(6,2), (7,2)],
                     col_width=2.4, font_size=14, row_height=0.6)

    plt.savefig('test.png', format='png')
    plt.show()
    """
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
