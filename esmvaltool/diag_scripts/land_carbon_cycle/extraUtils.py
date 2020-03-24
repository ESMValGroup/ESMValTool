import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl


def remove_invalid(tmp, fill_value=-9999.):
    '''
    removes the invalid non-numeric values from the input array and fills it with fill_value. Also removes all large and small values with magnitude beyond 1e15
    '''
    tmp = np.ma.masked_outside(tmp, -1e15, 1e15).filled(fill_value)
    whereisNan = np.isnan(tmp)
    tmp[whereisNan] = fill_value
    whereisNan = np.isinf(tmp)
    tmp[whereisNan] = fill_value
    return tmp


# axis lines, ticks, and fontsizes
def put_ticks(nticks=5, which_ax='both', axlw=0.3):
    '''
    puts the ticks on given locations and sets the width of axis lines
    '''
    ax = plt.gca()
    if which_ax == 'x':
        ax.xaxis.set_ticks_position('bottom')
        lines = ax.get_xticklines()
        labels = ax.get_xticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKDOWN)
        for label in labels:
            label.set_y(-0.02)
        ax.xaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))

    if which_ax == 'y':
        ax.yaxis.set_ticks_position('left')
        lines = ax.get_yticklines()
        labels = ax.get_yticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKLEFT)
            line.set_linewidth(axlw)
        ax.yaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))
    if which_ax == 'both':
        ax.yaxis.set_ticks_position('left')
        lines = ax.get_yticklines()
        labels = ax.get_yticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKLEFT)
            line.set_linewidth(axlw)
        ax.yaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))
        ax.xaxis.set_ticks_position('bottom')
        lines = ax.get_xticklines()
        labels = ax.get_xticklabels()
        for line in lines:
            line.set_marker(mpl.lines.TICKDOWN)
        ax.xaxis.set_major_locator(
            plt.MaxNLocator(nbins=nticks, min_n_ticks=nticks))
    return


def rem_ticks(which_ax='both'):
    '''
    removes ticks from either x or y axis and preserves the lines
    '''
    ax = plt.gca()
    if which_ax == 'x' or which_ax == 'both':
        ax.set_xticklabels([])
        ax.xaxis.set_ticks_position("none")
    if which_ax == 'y' or which_ax == 'both':
        ax.set_yticklabels([])
        ax.yaxis.set_ticks_position("none")
    return


def rem_axLine(rem_list=['top', 'right'], axlw=0.4):
    '''
    removes the axis lines from the list of which lines to remove
    rem_list can be 'left', 'right', 'top', 'bottom'
    '''
    ax = plt.gca()
    for loc, spine in ax.spines.items():
        if loc in rem_list:
            spine.set_position(('outward', 0))
            spine.set_linewidth(0.)
        else:
            spine.set_linewidth(axlw)
    return


def ax_clr(axfs=7):
    '''
    removes all the axis lines
    '''
    rem_axLine(rem_list=['top', 'right', 'left', 'bottom'])
    rem_ticks(which_ax='both')
    return


def ax_clrX(axfs=7, axlw=0.3, nticks=3):
    '''
    removes all the axis lines except the left one
    '''
    rem_axLine(rem_list=['top', 'right', 'bottom'])
    rem_ticks(which_ax='x')
    plt.gca().tick_params(axis='y', labelsize=axfs)
    put_ticks(which_ax='y', axlw=0.3, nticks=nticks)
    return


def ax_clrY(axfs=7, axlw=0.3, nticks=3):
    '''
    removes all the axis lines except the bottom one
    '''
    rem_axLine(rem_list=['top', 'right', 'left'])
    rem_ticks(which_ax='y')
    put_ticks(which_ax='x', axlw=0.3, nticks=nticks)
    plt.gca().tick_params(axis='x', labelsize=axfs)
    return


def ax_clrXY(axfs=7, axlw=0.3, nticks=3):
    '''
    removes the top and right axis
    '''
    rem_axLine(rem_list=['top', 'right'])
    put_ticks(which_ax='y', axlw=0.3, nticks=nticks)
    plt.gca().tick_params(axis='both', labelsize=axfs)
    return


def ax_orig(axfs=7, axlw=0.3, nticks=3):
    '''
    removes the top and right axis and sets the axis linewidth to a thin one
    '''
    rem_axLine(rem_list=['top', 'right'], axlw=axlw)
    put_ticks(which_ax='both', axlw=0.3, nticks=nticks)
    plt.gca().tick_params(axis='both', labelsize=axfs)
    return


def rotate_labels(which_ax='both', rot=0, axfs=6):
    '''
    rotates the ticks labels to rot and sets it fontsize to axfs
    '''
    if which_ax == 'x' or which_ax == 'both':
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=rot, fontsize=axfs)
    if which_ax == 'y' or which_ax == 'both':
        locs, labels = plt.yticks()
        plt.setp(labels, rotation=rot, fontsize=axfs)
    return


def draw_line_legend(ax_fs=8):
    '''
    draws a legend for line plots and puts it outside the plot area in x-direction
    '''
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


# MAPs
def get_colomap(cmap_nm, bounds__, lowp=0.05, hip=0.95):
    '''
    Get the list of colors from any official colormaps in matplotlib. 
    It returns the number of colors based on the number of items in the bounds. 
    Bounds is a list of boundary for each color.
    '''
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
                nticks=10,
                cbtitle='',
                tick_locs=[],
                ex_tend='both',
                cb_or='horizontal',
                spacing='uniform'):
    '''
    Plots the colorbar to the axis given by axcol_. Uses arrows on two sides.
    '''
    axco1 = plt.axes(axcol_)
    cb = mpl.colorbar.ColorbarBase(axco1,
                                   cmap=cm2,
                                   norm=mpl.colors.BoundaryNorm(
                                       bounds__, cm2.N),
                                   boundaries=bounds__,
                                   orientation=cb_or,
                                   drawedges=False,
                                   extend=ex_tend,
                                   ticks=tick_locs,
                                   spacing=spacing)
    cb.ax.tick_params(labelsize=cbfs, size=2, width=0.3)
    # hack the lines of the colorbar to make them white, the same color of background so that the colorbar looks broken.
    cb.outline.set_alpha(0.)
    cb.outline.set_color('white')
    cb.outline.set_linewidth(0 * cblw)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(cbfs)
        t.set_rotation(cbrt)

    for t in cb.ax.get_xticklabels():
        t.set_fontsize(cbfs)
        t.set_rotation(cbrt)
        t.set_y(-0.02)
    if cbtitle != '':
        cb.ax.set_title(cbtitle, fontsize=1.3 * cbfs)
    cb.update_ticks()
    return cb


def mk_colo_cont(axcol_,
                 bounds__,
                 cm2,
                 cblw=0.1,
                 cbrt=0,
                 cbfs=9,
                 nticks=10,
                 cbtitle='',
                 col_scale='linear',
                 tick_locs=[],
                 ex_tend='both',
                 cb_or='horizontal',
                 spacing='uniform'):
    '''
    Plots the colorbar to the axis given by axcol_. Uses arrows on two sides.
    '''

    axco1 = plt.axes(axcol_)
    if col_scale == 'linear':
        cb = mpl.colorbar.ColorbarBase(axco1,
                                       cmap=cm2,
                                       norm=mpl.colors.BoundaryNorm(
                                           bounds__, cm2.N),
                                       boundaries=bounds__,
                                       orientation=cb_or,
                                       drawedges=False,
                                       extend=ex_tend,
                                       ticks=bounds__[1:-1],
                                       spacing=spacing)
        if len(tick_locs) == 0:
            tick_locator = mpl.ticker.MaxNLocator(nbins=nticks,
                                                  min_n_ticks=nticks)
        else:
            tick_locator = mpl.ticker.FixedLocator(tick_locs)
        cb.locator = tick_locator
        cb.update_ticks()
    if col_scale == 'log':
        cb = mpl.colorbar.ColorbarBase(axco1,
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
            tick_locs_ori = cb.ax.xaxis.get_ticklocs()
        else:
            tick_locs_ori = cb.ax.yaxis.get_ticklocs()
        tick_locs_bn = []
        for _tl in tick_locs:
            tlInd = np.argmin(np.abs(bounds__[1:-1] - _tl))
            tick_locs_bn = np.append(tick_locs_bn, tick_locs_ori[tlInd])
        if cb_or == 'horizontal':
            cb.ax.xaxis.set_ticks(tick_locs_bn)
            cb.ax.xaxis.set_ticklabels(tick_locs)
        else:
            cb.ax.yaxis.set_ticks(tick_locs_bn)
            cb.ax.yaxis.set_ticklabels(tick_locs)
    cb.ax.tick_params(labelsize=cbfs, size=2, width=0.3)
    cb.outline.set_alpha(0.)
    cb.outline.set_color('white')
    cb.outline.set_linewidth(0 * cblw)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(cbfs)
        t.set_rotation(cbrt)

    for t in cb.ax.get_xticklabels():
        t.set_fontsize(cbfs)
        t.set_rotation(cbrt)
        t.set_y(-0.02)
    if cbtitle != '':
        cb.ax.set_title(cbtitle, fontsize=1.3 * cbfs)
    return cb
