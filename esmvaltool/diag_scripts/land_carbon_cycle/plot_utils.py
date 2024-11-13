"""Provide additional plotting utilities."""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


def ax_clr():
    """Remove all the axis lines."""
    rem_ax_line(rem_list=['top', 'right', 'left', 'bottom'])
    rem_ticks(which_ax='both')


def ax_clr_x(axfs=7, axlw=0.3, nticks=3):
    """Remove all the axis lines except the left one."""
    rem_ax_line(rem_list=['top', 'right', 'bottom'])
    rem_ticks(which_ax='x')
    plt.gca().tick_params(axis='y', labelsize=axfs)
    put_ticks(which_ax='y', axlw=axlw, nticks=nticks)


def ax_clr_y(axfs=7, axlw=0.3, nticks=3):
    """Remove all the axis lines except the bottom one."""
    rem_ax_line(rem_list=['top', 'right', 'left'])
    rem_ticks(which_ax='y')
    put_ticks(which_ax='x', axlw=axlw, nticks=nticks)
    plt.gca().tick_params(axis='x', labelsize=axfs)


def ax_clr_xy(axfs=7, axlw=0.3, nticks=3):
    """Remove the top and right axis."""
    rem_ax_line(rem_list=['top', 'right'])
    put_ticks(which_ax='y', axlw=axlw, nticks=nticks)
    plt.gca().tick_params(axis='both', labelsize=axfs)


def ax_orig(axfs=7, axlw=0.3, nticks=3):
    """
    Remove the top and right axis line.

    Set the axis linewidth to axlw.
    """
    rem_ax_line(rem_list=('top', 'right'), axlw=axlw)
    put_ticks(which_ax='both', axlw=axlw, nticks=nticks)
    plt.gca().tick_params(axis='both', labelsize=axfs)


def draw_line_legend(ax_fs=8):
    """
    Draws a legend for line plots.

    It puts it outside the plot area in
    x-direction
    """
    leg = plt.legend(loc=(1.00974, .06),
                     fontsize=ax_fs,
                     ncol=1,
                     columnspacing=0.05,
                     fancybox=True,
                     handlelength=1.5)
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('#eeeeee')
    leg.legendPatch.set_alpha(0.45)
    texts = leg.get_texts()
    plt.setp(texts, fontsize=ax_fs * 0.9)
    return leg


def get_colomap(cmap_nm, bounds__, lowp=0.05, hip=0.95):
    """
    Get the list of colors from any official colormaps in matplotlib.

    It returns the number of colors based on the number of items in the bounds.
    Bounds is a list of boundary for each color.
    """
    cmap__ = mpl.cm.get_cmap(cmap_nm)
    clist_v = np.linspace(lowp, hip, len(bounds__) - 1)
    rgba_ = [cmap__(_cv) for _cv in clist_v]
    return rgba_


def mk_colo_tau(axcol_,
                bounds__,
                cm2,
                cblw=0.1,
                cbrt=0,
                cbfs=9,
                cbtitle='',
                tick_locs=(),
                ex_tend='both',
                cb_or='horizontal',
                spacing='uniform'):
    """
    Plot the colorbar to the axis given by axcol_.

    Uses arrows on two sides.
    """
    axco1 = plt.axes(axcol_)
    col_bar = mpl.colorbar.ColorbarBase(axco1,
                                        cmap=cm2,
                                        norm=mpl.colors.BoundaryNorm(
                                            bounds__, cm2.N),
                                        boundaries=bounds__,
                                        orientation=cb_or,
                                        drawedges=False,
                                        extend=ex_tend,
                                        ticks=tick_locs,
                                        spacing=spacing)
    col_bar.ax.tick_params(labelsize=cbfs, size=2, width=0.3)
    # hack the lines of the colorbar to make them white, the same color of
    # background so that the colorbar looks broken.
    col_bar.outline.set_alpha(0.)
    col_bar.outline.set_color('white')
    col_bar.outline.set_linewidth(0 * cblw)
    for ti_ck in col_bar.ax.get_yticklabels():
        ti_ck.set_fontsize(cbfs)
        ti_ck.set_rotation(cbrt)

    for ti_ck in col_bar.ax.get_xticklabels():
        ti_ck.set_fontsize(cbfs)
        ti_ck.set_rotation(cbrt)
        ti_ck.set_y(-0.02)
    if cbtitle != '':
        col_bar.ax.set_title(cbtitle, fontsize=1.3 * cbfs)
    col_bar.update_ticks()
    return col_bar


def mk_colo_cont(axcol_,
                 bounds__,
                 cm2,
                 cblw=0.1,
                 cbrt=0,
                 cbfs=9,
                 nticks=10,
                 cbtitle='',
                 col_scale='linear',
                 tick_locs=(),
                 ex_tend='both',
                 cb_or='horizontal',
                 spacing='uniform'):
    """
    Plot the colorbar to the axis given by axcol_.

    Uses arrows on two sides.
    """
    axco1 = plt.axes(axcol_)
    if col_scale == 'linear':
        col_bar = mpl.colorbar.ColorbarBase(axco1,
                                            cmap=cm2,
                                            norm=mpl.colors.BoundaryNorm(
                                                bounds__, cm2.N),
                                            boundaries=bounds__,
                                            orientation=cb_or,
                                            drawedges=False,
                                            extend=ex_tend,
                                            ticks=bounds__[1:-1],
                                            spacing=spacing)
        if not tick_locs:
            tick_locator = mpl.ticker.MaxNLocator(nbins=nticks,
                                                  min_n_ticks=nticks)
        else:
            tick_locator = mpl.ticker.FixedLocator(tick_locs)
        col_bar.locator = tick_locator
        col_bar.update_ticks()
    if col_scale == 'log':
        col_bar = mpl.colorbar.ColorbarBase(axco1,
                                            cmap=cm2,
                                            norm=mpl.colors.BoundaryNorm(
                                                bounds__, cm2.N),
                                            boundaries=bounds__,
                                            orientation=cb_or,
                                            extend=ex_tend,
                                            drawedges=False,
                                            ticks=bounds__[1:-1],
                                            spacing=spacing)
        if cb_or == 'horizontal':
            tick_locs_ori = col_bar.ax.xaxis.get_ticklocs()
        else:
            tick_locs_ori = col_bar.ax.yaxis.get_ticklocs()
        tick_locs_bn = []
        for _tl in tick_locs:
            tl_ind = np.argmin(np.abs(bounds__[1:-1] - _tl))
            tick_locs_bn = np.append(tick_locs_bn, tick_locs_ori[tl_ind])
        if cb_or == 'horizontal':
            col_bar.ax.xaxis.set_ticks(tick_locs_bn)
            col_bar.ax.xaxis.set_ticklabels(tick_locs)
        else:
            col_bar.ax.yaxis.set_ticks(tick_locs_bn)
            col_bar.ax.yaxis.set_ticklabels(tick_locs)
    col_bar.ax.tick_params(labelsize=cbfs, size=2, width=0.3)
    col_bar.outline.set_alpha(0.)
    col_bar.outline.set_color('white')
    col_bar.outline.set_linewidth(0 * cblw)
    for ti_ck in col_bar.ax.get_yticklabels():
        ti_ck.set_fontsize(cbfs)
        ti_ck.set_rotation(cbrt)

    for ti_ck in col_bar.ax.get_xticklabels():
        ti_ck.set_fontsize(cbfs)
        ti_ck.set_rotation(cbrt)
        ti_ck.set_y(-0.02)
    if cbtitle != '':
        col_bar.ax.set_title(cbtitle, fontsize=1.3 * cbfs)
    return col_bar


def put_ticks(nticks=5, which_ax='both', axlw=0.3):
    """Put the ticks on given locations and sets the width of axis lines."""
    if which_ax == 'x':
        plt.gca().xaxis.set_ticks_position('bottom')
        lines = plt.gca().get_xticklines()
        labels = plt.gca().get_xticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKDOWN)
        for label in labels:
            label.set_y(-0.02)
        plt.gca().xaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))

    if which_ax == 'y':
        plt.gca().yaxis.set_ticks_position('left')
        lines = plt.gca().get_yticklines()
        labels = plt.gca().get_yticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKLEFT)
            line.set_linewidth(axlw)
        plt.gca().yaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))
    if which_ax == 'both':
        plt.gca().yaxis.set_ticks_position('left')
        lines = plt.gca().get_yticklines()
        labels = plt.gca().get_yticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKLEFT)
            line.set_linewidth(axlw)
        plt.gca().yaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))
        plt.gca().xaxis.set_ticks_position('bottom')
        lines = plt.gca().get_xticklines()
        for line in lines:
            line.set_marker(mpl.lines.TICKDOWN)
        plt.gca().xaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))


def rem_ticks(which_ax='both'):
    """Remove ticks from either x or y axis and preserves the lines."""
    if which_ax in ('x', 'both'):
        plt.gca().set_xticklabels([])
        plt.gca().xaxis.set_ticks_position("none")
    if which_ax in ('y', 'both'):
        plt.gca().set_yticklabels([])
        plt.gca().yaxis.set_ticks_position("none")


def rem_ax_line(rem_list=('top', 'right'), axlw=0.4):
    """
    Remove the axis lines.

    It uses the list of which lines to remove
    rem_list can be 'left', 'right', 'top', 'bottom'
    """
    for loc, spine in plt.gca().spines.items():
        if loc in rem_list:
            spine.set_position(('outward', 0))
            spine.set_linewidth(0.)
        else:
            spine.set_linewidth(axlw)


def rotate_labels(which_ax='both', rot=0, axfs=6):
    """
    Rotate the ticks labels to rot.

    Also sets it fontsize to axfs
    """
    if which_ax in ('x', 'both'):
        _, labels = plt.xticks()
        plt.setp(labels, rotation=rot, fontsize=axfs)
    if which_ax in ('y', 'both'):
        _, labels = plt.yticks()
        plt.setp(labels, rotation=rot, fontsize=axfs)


def save_figure(plot_path, _extr_art=None):
    """Write the figure to a file."""
    plt.savefig(plot_path,
                bbox_inches='tight',
                bbox_extra_artists=_extr_art,
                dpi=450)
